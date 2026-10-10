# project-ostrea-chil

Genome-guided transcriptome analysis of the Chilean flat oyster (*Ostrea chilensis*, "Och"), using Oxford Nanopore long-read RNA-seq from oysters collected at three locations in southern Chile, aligned against a newly assembled and annotated diploid (two-haplotype) reference genome.

The repository holds the analysis code (RStudio + R Markdown notebooks and shell workflows), a small amount of committed data, and the genome annotation report. Large inputs and all pipeline outputs live outside git — see [Data availability](#data-availability).

---

## Biological question

Three questions drive the work in this repo:

1. **Reference choice and haplotype structure** — how do the two haplotype assemblies (`Och_HapA`, `Och_HapB`) and the merged assembly relate to one another (synteny), and which is the right alignment target?
2. **Population/site differences in gene expression** — do oysters from the three sampling locations differ in transcript abundance, and in which genes?
3. **Annotation completeness** — how much of the assembled transcriptome is captured by the upstream gene set, how much is novel, and can the unmatched fraction be annotated (Swiss-Prot / nr / InterPro / KEGG, plus candidate lncRNAs)?

## Samples

Twelve barcoded ONT RNA-seq libraries, four per site:

| Site prefix | Source directory | Barcodes |
|---|---|---|
| `Pum` | Pumalin | barcode02, barcode12_1, barcode13, barcode13_1 |
| `Qui` | Quilhua | barcode01, barcode01_1, barcode02_1, barcode12 |
| `Rio` | Rio_Pudeto | barcode03, barcode03_1, barcode14, barcode14_1 |

Site prefix is the grouping variable throughout: differential expression is run as the three pairwise contrasts `Pum_vs_Qui`, `Pum_vs_Rio`, `Qui_vs_Rio`.

One public Illumina paired-end run (`SRR30335149`) is also aligned, in [10-illumina.Rmd](code/10-illumina.Rmd), as an independent check on the annotation-aware HISAT2 index.

## Reference genome

Assembled and annotated outside this repository. The accompanying report ([CGBIH240109_Genome_Annotation-Report.pdf](reports/CGBIH240109_Genome_Annotation-Report.pdf), C. Gallardo-Escárate, 2024-02-06) gives these headline figures:

- Two haplotype assemblies, `Och_HapA_assembly` and `Och_HapB_assembly`, 10 sequences each; ~790 Mb total per haplotype
- **69.4%** repeat content (largest classes: unknown 28.7%, satellite/tandem 14.5%, LINE 10.1%)
- **21,444** predicted genes / mRNAs; mean mRNA length 7,532 bp; 13,708 multi-exon vs 7,736 single-exon
- 2,923 non-coding RNAs (miRNA 1,329; snoRNA 913; tRNA 341; …)
- Functional annotation coverage: nr 80.6%, Uniprot 75.8%, KEGG 75.1%, InterPro 61.6%, GO 37.2%

Reference files referenced by the code (see [Data availability](#data-availability)):

| File | Role |
|---|---|
| `Och_HapA_assembly.fa`, `Och_HapB_assembly.fa` | the two haplotype assemblies |
| `merged_out.fasta` (+ `.fai`) | merged assembly — alignment target for the current pipeline |
| `merged_out.fasta.TBtools.fa` (+ `.fai`) | renamed/reformatted merged assembly used by the `04`/`08` series |
| `GN.gene.gff3`, `GN.gene.gtf` | upstream structural annotation (GTF produced via `gffread`) |
| `uniprot_sprot_r2024_03.fasta` | Swiss-Prot release used for all BLAST/DIAMOND annotation |
| `SPannot-GO.tsv` | Swiss-Prot → GO/protein-name table for joining BLAST hits |

## Repository layout

```
code/       numbered R Markdown notebooks + numbered shell workflows
data/       inputs (most gitignored; large reference files not committed)
output/     per-notebook results, one subdirectory per notebook (gitignored)
reports/    genome annotation report
bu.sh       rsync backup of the repo to gannet
```

Convention: notebook `NN-name.Rmd` writes to `output/NN-name/`, and all paths inside notebooks are relative to `code/` (i.e. `../data/...`, `../output/...`). Open [project-ostrea-chil.Rproj](project-ostrea-chil.Rproj) so the working directory resolves correctly.

## Analysis workflow

| Notebook | What it does |
|---|---|
| [00-qc.Rmd](code/00-qc.Rmd) | ONT read QC and trimming: FastQC + MultiQC, NanoFilt (`-q 10 -l 500`), Porechop adapter removal, NanoPlot; compiles per-sample NanoStats into one CSV |
| [01-repro-annot.Rmd](code/01-repro-annot.Rmd) | Builds the Swiss-Prot BLAST database; `blastp` of an *Ostrea* protein set (`PO2457_Ostrea_lurida.protein.fasta`) against Swiss-Prot at `1e-20`, aimed at characterizing reproduction-associated proteins. Output: [Olur-uniprot_blastp.tab](output/01-repro-annot/Olur-uniprot_blastp.tab) |
| [02-RNAseq.Rmd](code/02-RNAseq.Rmd) | Concatenates per-barcode FASTQs into 12 sample files; `minimap2 -ax splice -k14 --MD` to `Och_HapA`; StringTie (`-L`) assemble → merge (guided by `GN.gene.gtf`) → quantify (`-e -B`) → `prepDE.py` count matrices; extracts transcript FASTA and runs `blastx` vs Swiss-Prot; exploratory `bambu` run |
| [03-annotation.Rmd](code/03-annotation.Rmd) | Splits the pipe-delimited Swiss-Prot accession field out of the blastx table and left-joins it to `SPannot-GO.tsv` to attach protein names and GO terms |
| [04-align-v1.Rmd](code/04-align-v1.Rmd) | Same alignment/assembly against `merged_out.fasta.TBtools.fa`; BAM → bedGraph → bigWig track generation; DESeq2 size-factor normalization, limma pairwise contrasts, variable-gene heatmaps, per-gene expression plots |
| [04.5-align-v1-tr.Rmd](code/04.5-align-v1-tr.Rmd) | Transcript-level (rather than gene-level) version of the `04` count matrix, plus three PCA variants (top-500 by variance, top-500 by per-gene ANOVA, all genes) |
| [05-synteny.Rmd](code/05-synteny.Rmd) | `minimap2 -x asm5` whole-genome alignments producing PAFs for HapA vs HapB, and each haplotype vs the merged assembly |
| [06-Annot-novel.Rmd](code/06-Annot-novel.Rmd) | `megablast` of assembled transcripts against NCBI `nt` to characterize transcripts with no Swiss-Prot hit |
| [07-lncRNA-discovery.Rmd](code/07-lncRNA-discovery.Rmd) | lncRNA discovery pipeline: StringTie merge → `gffcompare` → keep class codes `u`/`x`/`o`/`i` and length > 199 bp → `bedtools getfasta` → CPC2 coding-potential filter → lncRNA FASTA/BED/GTF + summary stats. **Not yet ported to this project** — see [Known issues](#known-issues) |
| [08-align-v1-contrain.Rmd](code/08-align-v1-contrain.Rmd) | Annotation-*constrained* quantification (`stringtie -e -G GN.gene.gtf`, no novel transcripts) as a contrast to the de novo assembly in `04`; DE analysis, DIAMOND and NCBI `blastx` annotation of significant genes, UpSet plot of contrast overlap |
| [09-align-v1-clean.Rmd](code/09-align-v1-clean.Rmd) | **Current main pipeline.** Porechop-cleaned reads → `merged_out.fasta` → StringTie assemble/merge/quantify → count matrices → DESeq2 normalization → limma pairwise DE (`adj.P.Val < 0.05`, `abs(logFC) > 1`) → UpSet plot; bigWig tracks; exports a GTF subset containing only DE genes |
| [10-illumina.Rmd](code/10-illumina.Rmd) | HISAT2 index built with extracted exons and splice sites from `GN.gene.gtf`; aligns public Illumina run `SRR30335149`; BAM → sorted bedGraph → bigWig |
| [11-extract_sequences.sh](code/11-extract_sequences.sh) | `gffread` extraction of transcript FASTA from the `09` merged GTF |
| [11-extract_and_blast_workflow.sh](code/11-extract_and_blast_workflow.sh) | Extraction plus DIAMOND `blastx` against Swiss-Prot |
| [11-identify_unmatched.sh](code/11-identify_unmatched.sh) | Set-difference of assembled transcripts vs Swiss-Prot hits; writes the unmatched FASTA |
| [11-simple_annotation.sh](code/11-simple_annotation.sh) | Length statistics and ORF screening of the unmatched transcripts |
| [11-annotate_unmatched_workflow.sh](code/11-annotate_unmatched_workflow.sh) | Broader annotation attempt on the unmatched set: DIAMOND vs nr, plus InterPro / HMMER / KEGG result directories |
| [12-compare-genomes.py](code/12-compare-genomes.py) | Downloads Assembly 1.0, `Och_HapA`, and `Och_HapB`; reports indexed assembly statistics; runs all three pairwise `minimap2 -x asm5 --secondary=no` alignments; summarizes identity and reciprocal coverage; and draws SVG dot plots |
| [13-wgbs-submit.sh](code/13-wgbs-submit.sh) → [13.1](code/13.1-wgbs-genome-prep.sh) / [13.2](code/13.2-wgbs-trim-align.sh) / [13.3](code/13.3-wgbs-qc-summary.sh) | **WGBS DNA methylation (SLURM, klone).** 9 independent oysters (1 Pum, 4 Qui, 4 Rio; [data/13-wgbs-samples.tsv](data/13-wgbs-samples.tsv)) out of 15 delivered libraries — duplicated FASTQs excluded, see [Known issues](#known-issues), directional 2×151 PE. Downloads `merged_out.fasta` + `GN.gene.gff3` and builds the Bismark index; per-sample array: fastp → Bismark (`--score_min L,0,-0.6`) → dedup → methylation extraction → strand-merged CpG coverage; Bismark/MultiQC reports + `alignment_summary.tsv`. Shared settings in [13.0-wgbs-config.sh](code/13.0-wgbs-config.sh) |
| [13.4-wgbs-methylation.Rmd](code/13.4-wgbs-methylation.Rmd) | methylKit: coverage/global methylation QC, PCA, feature annotation (exon/intron/promoter/intergenic), pairwise DMLs and 1 kb DMRs; gene-body methylation tested with limma on M-values using the same three contrasts as the RNA-seq DE |
| [13.5-wgbs-gene-annotation.sh](code/13.5-wgbs-gene-annotation.sh) | Per-gene functional annotation of `GN.gene.gff3` for 13.6: `gffread` proteins → DIAMOND blastp vs current Swiss-Prot (E ≤ 1e-8) → best hit + GO (BP/MF/CC) per GN gene. Needs a `wgbs-annot` conda env (`diamond`, `gffread`) |
| [13.6-wgbs-methylation-landscape.Rmd](code/13.6-wgbs-methylation-landscape.Rmd) | Genome-wide methylome, pooling all oysters in the sample table: global and per-CpG methylation; exon/intron/promoter/downstream/intergenic composition, CpG O/E, methylation and enrichment of methylated CpGs; 100 kb chromosome profiles vs gene density; gene-body methylation (bimodality, CpG O/E, gene structure, metagene TSS→TES profile, exon/intron rank, among-oyster variability); methylation vs annotation status, GO BP term tests, and broad functional groups |
| [13.7-wgbs-bs-snp-check.sh](code/13.7-wgbs-bs-snp-check.sh) + [.py](code/13.7-wgbs-bs-snp-check.py), [13.7-…-batch.sh](code/13.7-wgbs-bs-snp-check-batch.sh) | Bisulfite-aware SNP check: is a methylation difference epigenetic, or a SNP / repeat artefact? Splits each oyster's BAM by Bismark strand and counts each base only where bisulfite can't change it; reports variant sites differentiated between the compared populations (pooled reads, ≥ 2 carrier oysters), CpG-destroying changes (C>T at C, G>A at G), per-oyster CpG methylation, depth, and share of MAPQ < 10 reads. One gene: `bash code/13.7-wgbs-bs-snp-check.sh GN020540 2000 Qui,Rio`; a DML/DMR table: `sbatch code/13.7-wgbs-bs-snp-check-batch.sh output/13-wgbs/13.4-methylation/DMR_Qui_vs_Rio.tsv`, which classes each region. Of the 305 Qui vs Rio DMRs: 208 no genetic signal, 41 CpG-destroying SNP, 36 repeat-like, 20 other differentiated variant; GN020540 is clean |
| [13.7-wgbs-snps.sh](code/13.7-wgbs-snps.sh) | SNPs from the bisulfite BAMs, to flag C/T variation that mimics methylation loss: `biscuit bsstrand` (adds strand tags to the Bismark BAMs) → joint `biscuit pileup` → confident SNPs (QUAL ≥ 20) → SNPs at reference CpGs → overlap with the 13.4 DMLs/DMRs. Needs a `wgbs-snp` conda env (`biscuit`, `samtools`, `bcftools`, `bedtools`) |

Note that the shell scripts in the `11-` and `13.x` series use paths relative to the **repository root** (`output/...`, `data/...`), unlike the notebooks. Run them from the repo root. Start the WGBS chain with `bash code/13-wgbs-submit.sh --test` (one sample) to check M-bias and mapping before the full run.

### Reading order

For the current state of the analysis, read `00` → `09` → `11-*`. The `02` and `04`/`04.5`/`08` series are earlier passes kept for provenance; they differ mainly in which reference they align to and whether StringTie is allowed to assemble novel transcripts.

## Compute environment

The notebooks are written for the Roberts Lab `raven` server and call tools by absolute path, e.g.:

```
/home/shared/hisat2-2.2.1/
/home/shared/samtools-1.12/
/home/shared/stringtie-2.2.1.Linux_x86_64/
/home/shared/ncbi-blast-2.11.0+/bin/   (also 2.15.0+ for nt searches)
/home/shared/diamond-2.1.8
/home/shared/bedtools2/bin/
/home/shared/gffcompare-0.12.6.Linux_x86_64/
/home/shared/CPC2_standalone-1.0.1/
/opt/anaconda/anaconda3/bin/           (conda, minimap2, NanoFilt, porechop, bedGraphToBigWig)
```

Thread counts are hardcoded per chunk (20–48). R analyses use `tidyverse`, `DESeq2`, `limma`, `pheatmap`, `UpSetR`, `rtracklayer`, `DT`, `Biostrings`, and (exploratory) `bambu`. Running elsewhere means editing tool paths and thread counts — `07-lncRNA-discovery.Rmd` is the only notebook that centralizes them as variables, and it carries a commented-out block for the `klone` HPC as a model for how to do this.

## Data availability

Committed to the repo:

| Path | Notes |
|---|---|
| [data/Och_contigs_proteins.fa](data/Och_contigs_proteins.fa) | 36,163 contigs from a trimmed-read assembly (55 MB). Despite the filename these are **nucleotide** sequences |
| [data/blast_contig_Och.xlsx](data/blast_contig_Och.xlsx) | BLAST results for those contigs |
| [output/01-repro-annot/Olur-uniprot_blastp.tab](output/01-repro-annot/Olur-uniprot_blastp.tab) | 27-line blastp excerpt (the only committed pipeline output) |
| [reports/CGBIH240109_Genome_Annotation-Report.pdf](reports/CGBIH240109_Genome_Annotation-Report.pdf) | genome annotation report |

**Not** in the repo, and required to run anything: the genome assemblies, the upstream annotation GFF/GTF, the raw ONT FASTQs (read in `02-RNAseq.Rmd` from `/home/shared/16TB_HDD_01/valentina/data/`), the Illumina FASTQs, the Swiss-Prot release, the BLAST/DIAMOND databases, and every `output/` subdirectory — `.gitignore` excludes `output`, `blastdb`, `data/*assembly*`, `data/P*`, and `data/uniprot_sprot_r2024_03.fast*`.

`bu.sh` rsyncs the working directory (excluding SAM files, bisulfite intermediates, and dotfiles) to `gannet:/volume1/v1_web/owlshell/bu-github/`, which is where the large outputs and rendered notebooks are served from.

## Reproducing

```bash
git clone https://github.com/RobertsLab/project-ostrea-chil.git
```

Then, on a machine with the tool paths above:

1. Stage the reference files listed under [Reference genome](#reference-genome) into `data/`.
2. Open `project-ostrea-chil.Rproj` in RStudio and run [00-qc.Rmd](code/00-qc.Rmd) chunk by chunk to produce `output/00-qc/trimmed/`.
3. Run [09-align-v1-clean.Rmd](code/09-align-v1-clean.Rmd) for alignment through differential expression.
4. Run the `11-*` shell scripts from the repository root for annotation of the resulting transcript set.

The three reference versions can be downloaded and compared independently with
Python 3, `curl`, and `minimap2`:

```bash
python3 code/12-compare-genomes.py --threads 16
```

This stores the downloaded FASTA/FAI files in `data/` and all generated tables,
compressed PAF alignments, logs, and SVG dot plots in
`output/12-compare-genomes/`. To inspect assembly sizes without downloading the
multi-gigabyte FASTA files or running alignments, use:

```bash
python3 code/12-compare-genomes.py --stats-only
```

Chunks are meant to be run interactively, not knit end to end: several are long-running, some are alternatives to one another, and `01-repro-annot.Rmd` sets `eval = FALSE` globally by design.

## Known issues

Open items a new reader should know about before trusting or extending the code:

- **WGBS: 6 of the 15 delivered FASTQ pairs duplicate other samples and are excluded.** All 15 come from one lane (`LH00129:489:23JFWNLT4:8`). Sampling 2,000 R1 read IDs per file and searching every other file ([output/13-wgbs/raw_read_overlap.tsv](output/13-wgbs/raw_read_overlap.tsv)) shows: `RP5` is `RP3` reordered (same 33,380,360 reads); all sampled `Qui5` reads are in `Qui1`; the five Pumalin files are overlapping subsets of one library (all `PUM5` reads are in `PUM1`, `PUM3` and `PUM4`; all `PUM2` reads in `PUM1`), and `PUM4` also carries 36-42% of the `Qui1`/`Qui5` reads. `Qui2-4`, `RP1`, `RP2` and `RP4` share no reads with any other file. The sample table keeps one file per distinct library: Pum1 (superset of PUM2/3/5), Qui1 (superset of Qui5), Qui2-4, Rio1-4. 13.4 therefore tests Qui vs Rio only (4 vs 4, `min_per_group = 3`); Pumalin (n = 1) waits for re-delivered libraries. 13.2 outputs for the dropped samples under `output/13-wgbs/` are stale. Ask the sequencing facility to re-demultiplex the lane.

- **The reference target is inconsistent across notebooks.** `02` aligns to `Och_HapA_assembly.fa`, `04`/`04.5`/`08` to `merged_out.fasta.TBtools.fa`, and `09`/`10`/`11-*` to `merged_out.fasta`. Results are not directly comparable across those series. `09` is the current one.
- **`07-lncRNA-discovery.Rmd` is an unadapted template.** Its variables still point at an *Acropora pulchra* project (`../data/31-Apul-lncRNA-discovery`, `Apulcra-genome.fa`, gannet URLs for that project), and one chunk references an unrelated `../output/17-Ptuh-lncRNA/` path. The pipeline logic is sound but the paths must be repointed at the *O. chilensis* genome and the `09` BAMs before it will run.
- **Sample→group assignment is positional.** `sample_info` in `04.5`, `08`, and `09` hardcodes `c("A","A","A","A","B",...)` against whatever order `prepDE.py` emitted; if column order changes, groups are silently wrong. The limma `group` factor is derived from the `Pum`/`Qui`/`Rio` name prefixes and is safe — the two should be reconciled to a single explicit metadata table.
- **Existing DE results predate the limma fix.** `04`, `04.5`, `08`, and `09` previously fit `lmFit` to linear-scale DESeq2-normalized counts; they now use `log2(norm_counts + 1)` with `eBayes(trend = TRUE)` (limma-trend). Both p-values and fold changes shift, and the `abs(logFC) > 1` filter only means 2-fold under the new version — anything already written to `output/` needs regenerating before it can be compared or reported.
- **No environment capture.** There is no `renv.lock`, conda YAML, or recorded tool-version list; versions must be inferred from the hardcoded paths.
- **Repo hygiene.** `.DS_Store` is tracked and unignored; `project-ostrea-chil.Rproj` is listed in `.gitignore` yet tracked; `output` is gitignored yet one output file is tracked; and 65 MB of `data/` binaries are committed.
- **`data/Och_contigs_proteins.fa` is misnamed** — the contents are nucleotide contigs, not proteins.
- **`01-repro-annot.Rmd` operates on a file named for *Ostrea lurida***, inside an *O. chilensis* project, and writes `Olur-`-prefixed output. Worth confirming whether that protein set is genuinely *O. lurida* (a comparative outgroup) or a mislabeled *O. chilensis* annotation.

## Repository

<https://github.com/RobertsLab/project-ostrea-chil> — Roberts Lab, University of Washington (SAFS).
