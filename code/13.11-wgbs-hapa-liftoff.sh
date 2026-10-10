#!/usr/bin/env bash
#SBATCH --job-name=13.11-hapa-liftoff
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=output/13-wgbs-hapa/logs/%x_%j.out

# Gene models on HapA coordinates, so gene-level steps (13.6) can run on the
# WGBS_REF=hapa alignment. Och_HapA_assembly.fa has no annotation; GN.gene.gff3
# is on merged_out.fasta, which shares HapA's 2A/4A/6A/10A but uses HapB's
# 1B/5B/6B/7B/8B/9B. Liftoff maps each GN gene onto HapA (genes on the shared
# chromosomes land in place; the rest move to their HapA homologs), keeping
# the GN gene IDs so the 13.5 Swiss-Prot/GO table still applies.
#
# Needs liftoff + minimap2, which are not in `myflow`. One-time setup:
#   conda create -n wgbs-liftoff -c conda-forge -c bioconda liftoff minimap2
#
# Run from the repo root after 13.1 has staged both genomes
# (WGBS_REF=hapa for HapA). Idempotent: skips the lift if the GFF exists.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

REF_FA="data/13-wgbs-genome/merged_out.fasta"
REF_GFF="data/13-wgbs-genome/GN.gene.gff3"
HAPA_DIR="data/13-wgbs-genome-hapa"
HAPA_FA="${HAPA_DIR}/Och_HapA_assembly.fa"
OUT_GFF="${HAPA_DIR}/GN.gene.hapa.gff3"
WORK="output/13-wgbs-hapa/13.11-liftoff"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

for f in "${REF_FA}" "${REF_GFF}" "${HAPA_FA}"; do
  [[ -s "${f}" ]] || { echo "ERROR: ${f} not found - run 13.1 (both references) first" >&2; exit 1; }
done
mkdir -p "${WORK}" output/13-wgbs-hapa/logs

source /mmfs1/gscratch/srlab/sr320/miniforge3/etc/profile.d/conda.sh
conda activate wgbs-liftoff

if [[ ! -s "${OUT_GFF}" ]]; then
  # -exclude_partial: genes below the coverage/identity thresholds (-a/-s, 0.5)
  # go to the unmapped list instead of into the GFF with partial flags.
  liftoff -g "${REF_GFF}" -o "${WORK}/GN.gene.hapa.gff3" -u "${WORK}/unmapped_features.txt" \
    -exclude_partial -p "${THREADS}" -dir "${WORK}/intermediate" \
    "${HAPA_FA}" "${REF_FA}"
  mv "${WORK}/GN.gene.hapa.gff3" "${OUT_GFF}"
  rm -rf "${WORK}/intermediate" "${REF_GFF}_db" "${REF_GFF}.db"   # re-creatable
fi

# Summary: genes lifted overall and per HapA chromosome, and how many kept
# their coordinates on the four chromosomes shared with merged_out.fasta.
n_ref=$(awk -F'\t' '$3 == "gene"' "${REF_GFF}" | wc -l)
n_lift=$(awk -F'\t' '$3 == "gene"' "${OUT_GFF}" | wc -l)
n_unmapped=$(grep -c . "${WORK}/unmapped_features.txt" || true)
same=$(awk -F'\t' -v OFS='\t' '
  $3 != "gene" { next }
  { match($9, /ID=[^;]+/); id = substr($9, RSTART + 3, RLENGTH - 3) }
  NR == FNR { ref[id] = $1 ":" $4 "-" $5; next }
  ($1 ~ /^Chromosome_(2|4|6|10)A$/) { n++; if (ref[id] == $1 ":" $4 "-" $5) s++ }
  END { printf "%d of %d", s, n }' "${REF_GFF}" "${OUT_GFF}")
{
  printf 'reference_genes\t%s\nlifted_genes\t%s\nunmapped_features\t%s\n' "${n_ref}" "${n_lift}" "${n_unmapped}"
  printf 'shared_chrom_genes_same_coords\t%s\n' "${same}"
  awk -F'\t' '$3 == "gene" { n[$1]++ } END { for (c in n) printf "genes_%s\t%d\n", c, n[c] }' "${OUT_GFF}" | sort -V
} | tee "${WORK}/liftoff_summary.tsv"
echo "[$(date)] done: ${OUT_GFF}"
