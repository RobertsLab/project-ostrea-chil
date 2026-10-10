#!/usr/bin/env bash
# Bisulfite-aware SNP check for one gene or region: is a methylation difference
# epigenetic, or a C->T SNP / divergent allele that only looks like one?
#
# Usage (from the repo root):
#   bash code/13.7-wgbs-bs-snp-check.sh GN020540              # gene ID from the GFF
#   bash code/13.7-wgbs-bs-snp-check.sh Chromosome_9B:7221497-7227996
#   bash code/13.7-wgbs-bs-snp-check.sh GN020540 5000         # flank in bp (default 2000)
#
# Splits each oyster's deduplicated BAM by Bismark strand (XG:CT / XG:GA),
# piles up the target +/- flank, and runs 13.7-wgbs-bs-snp-check.py, which
# counts each base only from the strand where bisulfite can't change it.
# Output: output/13-wgbs/13.7-bs-snp-check/<label>/{variant_sites,cpg_summary}.tsv
# plus summary.txt. A single gene takes seconds, so this is fine on a login node;
# loop over many regions inside an salloc session.

set -euo pipefail
cd "$(dirname "$0")/.."
source code/13.0-wgbs-config.sh
activate_env

target="${1:?give a gene ID (GN...) or region chr:start-end}"
flank="${2:-2000}"

if [[ "${target}" == *:*-* ]]; then
  chrom="${target%%:*}"; range="${target#*:}"
  start="${range%-*}";   end="${range#*-}"
  label="${chrom}_${start}-${end}"
else
  read -r chrom start end < <(awk -F'\t' -v id="${target}" \
    '$3 == "gene" && $9 ~ ("ID=" id ";") {print $1, $4, $5; exit}' "${GFF}")
  [[ -n "${chrom:-}" ]] || { echo "ERROR: ${target} not found as a gene in ${GFF}" >&2; exit 1; }
  label="${target}"
fi

chrom_len=$(awk -v c="${chrom}" '$1 == c {print $2}' "${GENOME_FA}.fai")
[[ -n "${chrom_len}" ]] || { echo "ERROR: ${chrom} not in ${GENOME_FA}.fai" >&2; exit 1; }
pad_start=$(( start > flank ? start - flank : 1 ))
pad_end=$(( end + flank < chrom_len ? end + flank : chrom_len ))
region="${chrom}:${pad_start}-${pad_end}"

WD="${OUT}/13.7-bs-snp-check/${label}"
mkdir -p "${WD}"
echo "Target ${chrom}:${start}-${end}; piling up ${region}"

samtools faidx "${GENOME_FA}" "${region}" | tail -n +2 | tr -d '\n' > "${WD}/ref.seq"

while IFS=$'\t' read -r sample _; do
  bam="${OUT}/dedup/${sample}.sorted.bam"
  [[ -s "${bam}" ]] || { echo "ERROR: missing ${bam}" >&2; exit 1; }
  for xg in CT GA; do
    samtools view -u -d "XG:${xg}" "${bam}" "${region}" |
      samtools mpileup -B -Q 20 -q 10 -d 10000 -f "${GENOME_FA}" - 2>/dev/null \
      > "${WD}/${sample}.${xg}.pileup"
  done
done < <(tail -n +2 "${SAMPLES}")

python3 -I code/13.7-wgbs-bs-snp-check.py \
  --workdir "${WD}" --samples "${SAMPLES}" \
  --region-start "${pad_start}" --target-start "${start}" --target-end "${end}" \
  | tee "${WD}/summary.txt"
