#!/usr/bin/env bash
#SBATCH --job-name=13.7-wgbs-snps
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# SNPs from the bisulfite BAMs, to separate genetic C/T variation from methylation.
# A C->T SNP at a CpG (G->A on the other strand) reads as an unmethylated CpG, so
# population-specific SNPs can show up as DMLs/DMRs in 13.4.
#   per sample:  Bismark dedup BAM -> biscuit bsstrand (adds the YD strand tag
#                biscuit needs)
#   all samples: biscuit pileup (bisulfite-aware joint genotyping) -> VCF
#   then:        confident SNPs (QUAL >= MIN_QUAL), the subset at CpGs, and how
#                many 13.4 DMLs/DMRs they touch.
# Samples are the rows of data/13-wgbs-samples.tsv (duplicated libraries excluded).
#
# Needs biscuit, samtools, bcftools and bedtools, which are not in `myflow`. One-time setup:
#   conda create -n wgbs-snp -c conda-forge -c bioconda biscuit samtools bcftools htslib bedtools
#
# Idempotent: skips any step whose output already exists.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate wgbs-snp

SNP="${OUT}/13.7-snps"
TAGGED="${SNP}/yd-bams"   # strand-tagged BAM copies; removed once the VCF exists
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MIN_MAPQ=20               # Bismark: 42 = unique; biscuit's default (40) suits bwa-meth
MIN_QUAL=20               # site QUAL for the confident set (biscuit's PASS is ~QUAL >= 5)
mkdir -p "${SNP}" "${LOGS}"

mapfile -t samples < <(awk -F'\t' 'NR > 1 {print $1}' "${SAMPLES}")
echo "[$(date)] ${#samples[@]} samples: ${samples[*]}"

vcf="${SNP}/wgbs_snps.all.vcf.gz"

if [[ ! -s "${vcf}.tbi" ]]; then
  # ---- 1. strand tags ------------------------------------------------------------
  # Bismark records strand in XG/XR; biscuit pileup reads YD, which bsstrand -c
  # infers from C>T / G>A counts. Bismark BAMs carry no AS tag, so biscuit's
  # alignment-score filter (-a) never applies. One single-threaded job per sample.
  mkdir -p "${TAGGED}"
  for s in "${samples[@]}"; do
    out="${TAGGED}/${s}.yd.bam"
    [[ -s "${out}.bai" ]] && continue
    (
      biscuit bsstrand -c "${GENOME_FA}" "${OUT}/dedup/${s}.sorted.bam" "${out}" \
        > "${TAGGED}/${s}.bsstrand.log" 2>&1
      samtools index "${out}"
    ) &
  done
  wait
  for s in "${samples[@]}"; do
    [[ -s "${TAGGED}/${s}.yd.bam.bai" ]] || { echo "ERROR: bsstrand failed for ${s}" >&2; exit 1; }
  done
  echo "[$(date)] strand tags done"

  # ---- 2. joint pileup / genotyping ------------------------------------------------
  # VCF sample columns follow the BAM order here (= sample table order).
  bams=(); for s in "${samples[@]}"; do bams+=("${TAGGED}/${s}.yd.bam"); done
  biscuit pileup -@ "${THREADS}" -m "${MIN_MAPQ}" \
    -o "${SNP}/wgbs_snps.all.vcf" "${GENOME_FA}" "${bams[@]}"
  bgzip -f -@ "${THREADS}" "${SNP}/wgbs_snps.all.vcf"
  tabix -f -p vcf "${vcf}"
  rm -rf "${TAGGED}"   # ~3 GB per sample; re-creatable from the dedup BAMs
  echo "[$(date)] pileup done"
fi

# ---- 3. confident SNPs ----------------------------------------------------------------
# Biscuit writes every covered site; keep variant sites only. ALT "N" (ambiguous alt,
# e.g. T vs unconverted C) is kept on purpose: those sites can mimic methylation loss.
snps="${SNP}/wgbs_snps.q${MIN_QUAL}.vcf.gz"
if [[ ! -s "${snps}.tbi" ]]; then
  bcftools view -f PASS -i "QUAL >= ${MIN_QUAL} && ALT != \".\"" \
    --threads "${THREADS}" -Oz -o "${snps}" "${vcf}"
  tabix -f -p vcf "${snps}"
fi

# SNP sites as BED, and the subset at a CpG (the SNP base is the C or the G of a
# reference CpG), using the reference base on each side of the site.
snp_bed="${SNP}/wgbs_snps.q${MIN_QUAL}.bed"
cpg_bed="${SNP}/wgbs_snps.q${MIN_QUAL}.CpG.bed"
if [[ ! -s "${cpg_bed}" ]]; then
  bcftools query -f '%CHROM\t%POS0\t%END\t%REF\t%ALT\t%QUAL\n' "${snps}" > "${snp_bed}"
  # Window = [SNP-1, SNP+2) clipped at the contig start; name records where the SNP is.
  awk -F'\t' -v OFS='\t' '{ s = ($2 > 0) ? $2 - 1 : 0; print $1, s, $2 + 2, $1 ":" $2 ":" ($2 - s + 1) }' "${snp_bed}" \
    | bedtools getfasta -fi "${GENOME_FA}" -bed - -name -tab \
    | awk -F'\t' -v OFS='\t' '{
        split($1, k, "::"); n = split(k[1], p, ":")   # chrom may not contain ":" here
        i = p[n]; pos = p[n - 1]; seq = toupper($2)
        b = substr(seq, i, 1); nxt = substr(seq, i + 1, 1); prv = (i > 1) ? substr(seq, i - 1, 1) : ""
        if ((b == "C" && nxt == "G") || (b == "G" && prv == "C")) print p[1], pos, pos + 1
      }' \
    | sort -k1,1 -k2,2n -u > "${cpg_bed}"
fi
echo "[$(date)] $(wc -l < "${snp_bed}") SNPs (QUAL >= ${MIN_QUAL}); $(wc -l < "${cpg_bed}") at CpGs"

# ---- 4. overlap with 13.4 calls ---------------------------------------------------------
# A DML is affected if a SNP hits either base of its CpG; a DMR (1 kb tile) if any CpG
# inside it carries a SNP. The affected calls are written out for filtering in 13.4.
M="${OUT}/13.4-methylation"
summary="${SNP}/dm_snp_overlap.tsv"
printf 'file\tcalls\twith_CpG_SNP\tpct\n' > "${summary}"
for bed in "${M}"/DML_*.bed "${M}"/DMR_*.bed; do
  [[ -s "${bed}" ]] || continue
  name=$(basename "${bed}" .bed)
  if [[ "${name}" == DML_* ]]; then
    awk -F'\t' -v OFS='\t' '{print $1, $2, $3 + 1}' "${bed}"   # both CpG bases
  else
    cut -f1-3 "${bed}"
  fi > "${SNP}/${name}.query.bed"
  bedtools intersect -u -a "${SNP}/${name}.query.bed" -b "${cpg_bed}" > "${SNP}/${name}.with_CpG_SNP.bed"
  n=$(wc -l < "${SNP}/${name}.query.bed"); hit=$(wc -l < "${SNP}/${name}.with_CpG_SNP.bed")
  awk -v f="${name}" -v n="${n}" -v h="${hit}" 'BEGIN { printf "%s\t%d\t%d\t%.1f\n", f, n, h, (n ? 100 * h / n : 0) }' >> "${summary}"
  rm -f "${SNP}/${name}.query.bed"
done
column -t "${summary}"
echo "[$(date)] done"
