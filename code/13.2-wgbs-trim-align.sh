#!/usr/bin/env bash
#SBATCH --job-name=13.2-wgbs-align
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=1-15%5
#SBATCH --output=output/13-wgbs/logs/%x_%A_%a.out

# One array task per sample (row N of the sample sheet, header excluded):
#   fastp trim -> Bismark align -> deduplicate -> methylation extraction
#   -> strand-merged CpG coverage file (input to 13.4) + sorted BAM for IGV.
# Re-running a task skips steps whose outputs already exist.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh
activate_env

THREADS="${SLURM_CPUS_PER_TASK:-8}"
TASK="${SLURM_ARRAY_TASK_ID:?run as a SLURM array, or set SLURM_ARRAY_TASK_ID=N}"

# Raw filenames contain spaces and parentheses ("PUM1 (paired)_R1.fastq.gz"),
# so read the tab-separated row with IFS=tab and keep every expansion quoted.
IFS=$'\t' read -r sample population site r1 r2 \
  < <(awk -v n="$((TASK + 1))" 'NR == n' "${SAMPLES}")
[[ -n "${sample:-}" ]] || { echo "No sample on row ${TASK}" >&2; exit 1; }
echo "[$(date)] ${sample} (${population}/${site})"

TRIM="${OUT}/trimmed";  ALN="${OUT}/bismark";  DEDUP="${OUT}/dedup"
METH="${OUT}/methylation";  QC="${OUT}/fastp"
mkdir -p "${TRIM}" "${ALN}" "${DEDUP}" "${METH}" "${QC}"
TMP="${TMPDIR:-/tmp}/13-wgbs-${sample}-$$"
mkdir -p "${TMP}"; trap 'rm -rf "${TMP}"' EXIT

# ---- 1. trim ------------------------------------------------------------------
t1="${TRIM}/${sample}_R1.fq.gz";  t2="${TRIM}/${sample}_R2.fq.gz"
bam="${ALN}/${sample}_pe.bam"
if [[ ! -s "${bam}" && ! -s "${t2}" ]]; then
  fastp \
    -i "${RAW_DIR}/${r1}" -I "${RAW_DIR}/${r2}" \
    -o "${t1}" -O "${t2}" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --trim_front1 "${TRIM_FRONT1}" --trim_front2 "${TRIM_FRONT2}" \
    --length_required 50 \
    --thread 16 \
    --json "${QC}/${sample}.fastp.json" --html "${QC}/${sample}.fastp.html"
fi

# ---- 2. align -----------------------------------------------------------------
# Directional PE: each --parallel instance runs 2 bowtie2 processes x -p threads,
# so --parallel 8 -p 2 uses ~32 cores.
# Bismark refuses --basename with --parallel, so outputs get its default R1-derived
# names; rename them to ${sample}_pe.bam / ${sample}_PE_report.txt, which dedup,
# 13.3 (bismark2report/bismark2summary) and the report parsing all expect.
if [[ ! -s "${bam}" ]]; then
  bismark \
    --genome "${GENOME_DIR}" \
    -1 "${t1}" -2 "${t2}" \
    --score_min "${SCORE_MIN}" \
    --parallel 8 -p 2 \
    --temp_dir "${TMP}" \
    --output_dir "${ALN}"
  mv "${ALN}/${sample}_R1_bismark_bt2_PE_report.txt" "${ALN}/${sample}_PE_report.txt"
  mv "${ALN}/${sample}_R1_bismark_bt2_pe.bam" "${bam}"
  rm -f "${t1}" "${t2}"   # trimmed reads are re-creatable from raw; ~5 GB each
fi

# ---- 3. deduplicate --------------------------------------------------------------
dbam="${DEDUP}/${sample}_pe.deduplicated.bam"
if [[ ! -s "${dbam}" ]]; then
  deduplicate_bismark --paired --bam --output_dir "${DEDUP}" "${bam}"
fi

# ---- 4. methylation calls -------------------------------------------------------
# --no_overlap: count each CpG once where mates overlap.
# CHG/CHH calls are kept (merged) because the CHH % methylation is the
# bisulfite-conversion-efficiency check — there is no lambda spike-in.
cov="${METH}/${sample}_pe.deduplicated.bismark.cov.gz"
if [[ ! -s "${cov}" ]]; then
  bismark_methylation_extractor \
    --paired-end --no_overlap \
    --comprehensive --merge_non_CpG \
    --bedGraph --gzip \
    --parallel 8 --buffer_size 40G \
    --output "${METH}" \
    "${dbam}"
fi

# ---- 5. merge CpG strands -----------------------------------------------------------
# Combines top/bottom-strand counts per CpG dinucleotide — roughly doubles
# per-site depth, which matters at ~6-9x raw coverage.
merged="${METH}/${sample}.CpG_report.merged_CpG_evidence.cov.gz"
if [[ ! -s "${merged}" ]]; then
  coverage2cytosine \
    --genome_folder "${GENOME_DIR}" \
    --merge_CpG --gzip \
    --dir "${METH}" \
    -o "${sample}" \
    "${cov}"
  rm -f "${METH}/${sample}.CpG_report.txt.gz"   # genome-wide per-C report; large and unused
fi

# ---- 6. sorted/indexed BAM for IGV + nucleotide coverage for the Bismark report -------
sbam="${DEDUP}/${sample}.sorted.bam"
if [[ ! -s "${sbam}.bai" ]]; then
  samtools sort -@ 8 -m 2G -T "${TMP}/sort" -o "${sbam}" "${dbam}"
  samtools index -@ 8 "${sbam}"
fi
if [[ ! -s "${DEDUP}/${sample}_pe.deduplicated.nucleotide_stats.txt" ]]; then
  bam2nuc --genome_folder "${GENOME_DIR}" --dir "${DEDUP}" "${dbam}"
fi

echo "[$(date)] ${sample} done"
