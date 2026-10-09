#!/usr/bin/env bash
#SBATCH --job-name=13.3-wgbs-summary
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# Per-sample Bismark HTML reports, a cross-sample Bismark summary, MultiQC,
# and a flat table of the numbers worth checking before 13.4.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh
activate_env

REPORTS="${OUT}/reports"
mkdir -p "${REPORTS}"

# bismark2report pairs files by basename, so give it each sample's reports.
# Samples without a dedup report yet (e.g. after a --test run) are skipped.
while IFS=$'\t' read -r sample _; do
  [[ -s "${OUT}/dedup/${sample}_pe.deduplication_report.txt" ]] || continue
  bismark2report \
    --alignment_report "${OUT}/bismark/${sample}_PE_report.txt" \
    --dedup_report     "${OUT}/dedup/${sample}_pe.deduplication_report.txt" \
    --splitting_report "${OUT}/methylation/${sample}_pe.deduplicated_splitting_report.txt" \
    --mbias_report     "${OUT}/methylation/${sample}_pe.deduplicated.M-bias.txt" \
    --nucleotide_report "${OUT}/dedup/${sample}_pe.deduplicated.nucleotide_stats.txt" \
    --dir "${REPORTS}" \
    -o "${sample}_bismark_report.html"
done < <(tail -n +2 "${SAMPLES}")

# Name the BAMs explicitly: left alone, bismark2summary picks up any *_pe.bam
# in the folder, including stray files from an interrupted run.
mapfile -t bams < <(tail -n +2 "${SAMPLES}" | cut -f1 | sed 's/$/_pe.bam/' |
                    while read -r b; do [[ -s "${OUT}/bismark/${b}" ]] && echo "${b}"; done)
( cd "${OUT}/bismark" && bismark2summary --basename ../reports/bismark_summary_report "${bams[@]}" )

multiqc --force --filename multiqc_13-wgbs --outdir "${REPORTS}" \
  "${OUT}/fastp" "${OUT}/bismark" "${OUT}/dedup" "${OUT}/methylation"

# One row per sample: mapping, duplication, context methylation (CHH ~ 1 - conversion rate).
{
  printf 'sample\tpairs_analysed\tmapping_pct\tdup_pct\tCpG_meth_pct\tCHG_meth_pct\tCHH_meth_pct\n'
  while IFS=$'\t' read -r sample _; do
    aln="${OUT}/bismark/${sample}_PE_report.txt"
    dd="${OUT}/dedup/${sample}_pe.deduplication_report.txt"
    sp="${OUT}/methylation/${sample}_pe.deduplicated_splitting_report.txt"
    [[ -s "${sp}" ]] || continue
    pairs=$(awk -F'\t' '/^Sequence pairs analysed in total/{print $2}' "${aln}")
    map=$(awk -F'\t' '/^Mapping efficiency/{sub("%","",$2); print $2}' "${aln}")
    dup=$(grep -oP 'Total number duplicated alignments removed:\s+\d+ \(\K[0-9.]+' "${dd}")
    # CpG from the post-dedup extraction; CHG/CHH from the alignment report,
    # because --merge_non_CpG leaves only a combined non-CpG line in ${sp}.
    cpg=$(awk -F'\t' '/^C methylated in CpG context/{sub("%","",$2); print $2}' "${sp}")
    chg=$(awk -F'\t' '/^C methylated in CHG context/{sub("%","",$2); print $2}' "${aln}")
    chh=$(awk -F'\t' '/^C methylated in CHH context/{sub("%","",$2); print $2}' "${aln}")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${sample}" "${pairs}" "${map}" "${dup}" "${cpg}" "${chg}" "${chh}"
  done < <(tail -n +2 "${SAMPLES}")
} > "${OUT}/alignment_summary.tsv"

column -t "${OUT}/alignment_summary.tsv"
