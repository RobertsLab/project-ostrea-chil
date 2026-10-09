#!/usr/bin/env bash
#SBATCH --job-name=14.1-orthologs
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=output/14-meth-compare/logs/%x_%j.out

# One-to-one orthologs between O. chilensis and O. lurida, for 14.2.
#   per species: GFF + genome -> proteins (gffread) -> longest isoform per gene,
#                headers renamed to gene IDs
#   both:        DIAMOND blastp in each direction -> reciprocal best hits
# Species, genomes and GFFs come from data/14-species.tsv (first two rows), the
# same table 14.2 reads, so gene IDs agree between the two steps.
#
# Run from the repository root:
#   mkdir -p output/14-meth-compare/logs && sbatch code/14.1-meth-compare-orthologs.sh
# Needs gffread, diamond and samtools: uses the `wgbs-annot` env shared with 13.5
# (see that script for the one-time conda create). Set CONDA_ENV= (empty) to use
# whatever is already on PATH instead.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

SPECIES_TSV="data/14-species.tsv"
OUT="output/14-meth-compare/orthologs"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
EVALUE="1e-10"
CONDA_ROOT="/mmfs1/gscratch/srlab/sr320/miniforge3"
CONDA_ENV="${CONDA_ENV-wgbs-annot}"

[[ -f "${SPECIES_TSV}" ]] || { echo "ERROR: run from the repository root" >&2; exit 1; }
if [[ -n "${CONDA_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_ROOT}/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
fi
for tool in gffread diamond samtools; do
  command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not on PATH" >&2; exit 1; }
done
mkdir -p "${OUT}"

# ---- 1. proteins, one per gene ---------------------------------------------------
# gffread names proteins by transcript ID; map transcript -> gene from the mRNA
# lines (ID=, Parent=) and keep the longest protein per gene. Works for both the
# single-isoform GN annotation and multi-isoform NCBI-style GFFs.
species=()
while IFS=$'\t' read -r sp label genome gff _; do
  [[ "${sp}" == "species" || -z "${sp}" ]] && continue
  species+=("${sp}")
  faa="${OUT}/${sp}.genes.faa"
  [[ -s "${faa}" ]] && continue
  echo "[$(date)] ${sp} (${label}) proteins"
  for f in "${genome}" "${gff}"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
  done
  [[ -s "${genome}.fai" ]] || samtools faidx "${genome}"

  # -S: '*' for stop codons (DIAMOND reads '.' as a residue)
  gffread "${gff}" -g "${genome}" -y "${OUT}/${sp}.transcripts.faa" -S

  awk -F'\t' '$3 == "mRNA" {
      id = ""; par = ""
      n = split($9, kv, ";")
      for (i = 1; i <= n; i++) {
        if (kv[i] ~ /^ID=/)     id  = substr(kv[i], 4)
        if (kv[i] ~ /^Parent=/) par = substr(kv[i], 8)
      }
      if (id != "" && par != "") print id "\t" par
    }' "${gff}" > "${OUT}/${sp}.tx2gene.tsv"

  awk -v map="${OUT}/${sp}.tx2gene.tsv" '
    BEGIN { while ((getline l < map) > 0) { split(l, a, "\t"); gene[a[1]] = a[2] } }
    function flush() {
      if (tx != "" && (tx in gene)) {
        g = gene[tx]; s = seq; sub(/\*$/, "", s)
        if (length(s) > len[g]) { len[g] = length(s); best[g] = s }
      }
    }
    /^>/ { flush(); tx = substr($1, 2); seq = ""; next }
    { seq = seq $0 }
    END {
      flush()
      for (g in best) printf(">%s\n%s\n", g, best[g])
    }' "${OUT}/${sp}.transcripts.faa" > "${faa}"

  echo "  $(grep -c '>' "${OUT}/${sp}.transcripts.faa") transcripts -> $(grep -c '>' "${faa}") genes"
done < "${SPECIES_TSV}"

(( ${#species[@]} == 2 )) || { echo "ERROR: expected 2 species in ${SPECIES_TSV}" >&2; exit 1; }
a="${species[0]}"; b="${species[1]}"

# ---- 2. reciprocal best hits --------------------------------------------------------
# outfmt 6 column order: qseqid sseqid pident length evalue bitscore qcovhsp scovhsp
fmt=(6 qseqid sseqid pident length evalue bitscore qcovhsp scovhsp)
for pair in "${a} ${b}" "${b} ${a}"; do
  read -r q s <<< "${pair}"
  hits="${OUT}/${q}_vs_${s}.diamond.tsv"
  [[ -s "${hits}" ]] && continue
  [[ -s "${OUT}/${s}.dmnd" ]] || diamond makedb --in "${OUT}/${s}.genes.faa" -d "${OUT}/${s}" -p "${THREADS}"
  echo "[$(date)] diamond ${q} -> ${s}"
  diamond blastp \
    -q "${OUT}/${q}.genes.faa" -d "${OUT}/${s}" \
    --more-sensitive --evalue "${EVALUE}" \
    --max-target-seqs 1 --max-hsps 1 \
    -p "${THREADS}" \
    --outfmt "${fmt[@]}" \
    -o "${hits}"
done

# A pair is kept when each gene is the other's top hit (first line per query).
rbh="${OUT}/${a}_${b}_rbh.tsv"
awk -F'\t' -v OFS='\t' -v a="${a}" -v b="${b}" '
  BEGIN { print a "_gene", b "_gene", "pident", "aln_length", "evalue", "bitscore", a "_cov", b "_cov" }
  NR == FNR { if (!($1 in rev)) rev[$1] = $2; next }
  !($1 in seen) {
    seen[$1] = 1
    if (rev[$2] == $1) print $1, $2, $3, $4, $5, $6, $7, $8
  }
' "${OUT}/${b}_vs_${a}.diamond.tsv" "${OUT}/${a}_vs_${b}.diamond.tsv" > "${rbh}"

echo "[$(date)] $(( $(wc -l < "${rbh}") - 1 )) one-to-one orthologs -> ${rbh}"
