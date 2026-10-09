#!/usr/bin/env bash
#SBATCH --job-name=13.1-wgbs-genome
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# Download the reference + annotation and build the Bismark (bowtie2) index.
# Idempotent: skips any step whose output already exists.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh
activate_env

mkdir -p "${GENOME_DIR}" "${LOGS}"

for f in merged_out.fasta GN.gene.gff3; do
  if [[ ! -s "${GENOME_DIR}/${f}" ]]; then
    echo "Downloading ${f}"
    curl -fL --retry 3 -o "${GENOME_DIR}/${f}.part" "${GANNET}/${f}"
    # curl 8.14 can exit 0 after a failed --retry, so check the file itself.
    [[ -s "${GENOME_DIR}/${f}.part" ]] || { echo "ERROR: download of ${f} failed" >&2; exit 1; }
    mv "${GENOME_DIR}/${f}.part" "${GENOME_DIR}/${f}"
  fi
done

[[ -s "${GENOME_FA}.fai" ]] || samtools faidx "${GENOME_FA}"
cut -f1,2 "${GENOME_FA}.fai" > "${GENOME_DIR}/genome.chrom.sizes"

if [[ ! -d "${GENOME_DIR}/Bisulfite_Genome" ]]; then
  # --parallel N launches two bowtie2-build jobs (CT and GA) with N threads each
  bismark_genome_preparation --bowtie2 --parallel 8 --verbose "${GENOME_DIR}"
fi

echo "Genome prep done: $(date)"
