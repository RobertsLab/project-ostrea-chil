#!/usr/bin/env bash
#SBATCH --job-name=13.5-gene-annot
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# Per-gene functional annotation of GN.gene.gff3, for 13.6 (methylation vs function).
# The upstream genome report annotated these genes, but no per-gene table ships with
# it, so this rebuilds one: GN proteins -> DIAMOND blastp vs Swiss-Prot -> GO terms.
# Independent of the WGBS alignments; can run any time after 13.1 has staged the GFF.
#
# Needs gffread + DIAMOND, which are not in `myflow`. One-time setup:
#   conda create -n wgbs-annot -c conda-forge -c bioconda diamond gffread
#
# Idempotent: skips any step whose output already exists.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate wgbs-annot

ANNOT="${OUT}/13.5-gene-annotation"
SPROT="${ANNOT}/swissprot"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
EVALUE="1e-8"                     # same threshold as the genome annotation report
mkdir -p "${SPROT}" "${LOGS}"

# ---- 1. GN protein sequences --------------------------------------------------------
# One mRNA per gene in GN.gene.gff3; -S writes stops as '*', which DIAMOND accepts.
prot="${ANNOT}/GN.proteins.faa"
[[ -s "${prot}" ]] || gffread "${GFF}" -g "${GENOME_FA}" -S -y "${prot}"
echo "GN proteins: $(grep -c '^>' "${prot}")"

# ---- 2. Swiss-Prot sequences + GO/protein-name table (current release) -----------------
if [[ ! -s "${SPROT}/uniprot_sprot.dmnd" ]]; then
  curl -fL --retry 3 -o "${SPROT}/uniprot_sprot.fasta.gz" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
  diamond makedb --in "${SPROT}/uniprot_sprot.fasta.gz" -d "${SPROT}/uniprot_sprot" -p "${THREADS}"
  curl -fL --retry 3 -o "${SPROT}/reldate.txt" \
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt"
fi

go_tab="${SPROT}/uniprot_sprot_GO.tsv.gz"
if [[ ! -s "${go_tab}" ]]; then
  # go_p/go_f/go_c are "term name [GO:nnnnnnn]; ..." lists, so 13.6 gets names without GO.db.
  curl -fL --retry 3 -o "${go_tab}" \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=%28reviewed%3Atrue%29&fields=accession%2Cid%2Cprotein_name%2Corganism_name%2Cgo_p%2Cgo_f%2Cgo_c%2Cgo_id"
fi
[[ -s "${go_tab}" ]] || { echo "ERROR: Swiss-Prot GO table download failed" >&2; exit 1; }

# ---- 3. DIAMOND blastp, best hit per protein -------------------------------------------
hits="${ANNOT}/GN-sprot_blastp.tsv"
if [[ ! -s "${hits}" ]]; then
  diamond blastp \
    --query "${prot}" --db "${SPROT}/uniprot_sprot" \
    --more-sensitive --evalue "${EVALUE}" --max-target-seqs 1 \
    --threads "${THREADS}" \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --out "${hits}"
fi

# ---- 4. one row per gene: GN gene ID, best hit accession, GO ---------------------------
# sseqid is sp|ACCESSION|NAME; mRNA IDs are GN000001.1 -> gene GN000001.
out="${ANNOT}/GN_gene_annotation.tsv"
{
  printf 'gene_id\tmrna_id\taccession\tpident\taln_len\tqlen\tslen\tevalue\tbitscore\tentry_name\tprotein_name\torganism\tgo_bp\tgo_mf\tgo_cc\tgo_ids\n'
  awk -F'\t' -v OFS='\t' '
    NR == FNR { if (FNR > 1) ann[$1] = $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8; next }
    !seen[$1]++ {
      split($2, s, "|"); acc = s[2]
      gene = $1; sub(/\.[0-9]+$/, "", gene)
      a = (acc in ann) ? ann[acc] : "\t\t\t\t\t\t"
      print gene, $1, acc, $3, $4, $5, $6, $7, $8, a
    }' <(gzip -dc "${go_tab}") "${hits}"
} > "${out}"

n_genes=$(grep -c '^>' "${prot}")
n_hit=$(( $(wc -l < "${out}") - 1 ))
n_go=$(awk -F'\t' 'NR > 1 && $16 != ""' "${out}" | wc -l)
printf 'genes\t%s\nswissprot_hit\t%s\nwith_GO\t%s\n' "${n_genes}" "${n_hit}" "${n_go}" \
  | tee "${ANNOT}/annotation_summary.tsv"
echo "Swiss-Prot release: $(head -1 "${SPROT}/reldate.txt" 2>/dev/null)"
