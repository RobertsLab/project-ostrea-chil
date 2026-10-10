#!/usr/bin/env bash
#SBATCH --job-name=13-wgbs-read-overlap
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=6:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# Exact read sharing among the delivered WGBS FASTQs -> output/13-wgbs/raw_read_overlap.tsv
# Several delivered files contain reads of other samples (README Known issues), so
# every R1 read name in all files in RAW_DIR is compared at once: names are pooled with
# their file label, sorted, and each name shared by more than one file adds to every
# pair of files that hold it. Also reports each file's index (barcode) composition.
# All files come from one lane, so the run/flowcell/lane prefix is dropped from names.
# Run from the repository root: sbatch code/13-wgbs-read-overlap.sh

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh

THREADS="${SLURM_CPUS_PER_TASK:-8}"
TMP="${OUT}/tmp-read-overlap"
mkdir -p "${TMP}" "${LOGS}"
trap 'rm -rf "${TMP}"' EXIT

# ---- 1. names and index per file ------------------------------------------------------
shopt -s nullglob
r1s=("${RAW_DIR}"/*_R1.fastq.gz)
echo "[$(date)] ${#r1s[@]} R1 files"
for f in "${r1s[@]}"; do
  label=$(basename "${f}" | sed 's/ (paired)_R1\.fastq\.gz$//')
  (
    zcat "${f}" | awk -v OFS='\t' -v l="${label}" -v ix="${TMP}/${label}.index" '
      NR % 4 == 1 {
        n = split(substr($1, 2), p, ":"); print p[n - 2] ":" p[n - 1] ":" p[n], l
        split($2, a, ":"); b = substr(a[4], 1, 8); gsub(/N/, ".", b); c[b]++; t++
      }
      END { for (k in c) if (c[k] / t >= 0.01) printf "%s\t%s\t%d\t%d\n", l, k, c[k], t > ix }' \
      > "${TMP}/${label}.names"
  ) &
done
wait
echo "[$(date)] names extracted"

# ---- 2. pool, sort, count shared names per file pair -------------------------------------
export LC_ALL=C
cat "${TMP}"/*.names | sort -S 150G --parallel="${THREADS}" -T "${TMP}" -k1,1 \
  | awk -F'\t' -v OFS='\t' '
      function flush(   i, j) {
        for (i = 1; i <= k; i++) { tot[g[i]]++; for (j = 1; j <= k; j++) if (i != j) sh[g[i] SUBSEP g[j]]++ }
      }
      $1 != prev { if (NR > 1) flush(); prev = $1; k = 0 }
      { g[++k] = $2 }
      END {
        flush()
        print "query_file", "found_in_file", "query_reads", "shared_reads", "pct_of_query"
        for (key in sh) { split(key, x, SUBSEP); has[x[1]] = 1
          printf "%s\t%s\t%d\t%d\t%.2f\n", x[1], x[2], tot[x[1]], sh[key], 100 * sh[key] / tot[x[1]] }
        for (f in tot) if (!(f in has)) printf "%s\tnone\t%d\t0\t0.00\n", f, tot[f]
      }' \
  | { read -r header; echo "${header}"; sort -t$'\t' -k5,5gr -k1,1; } > "${OUT}/raw_read_overlap.tsv"

{ printf 'file\tindex\treads\tfile_reads\n'; cat "${TMP}"/*.index | sort -k1,1 -k3,3nr; } \
  > "${OUT}/raw_read_index.tsv"

column -t "${OUT}/raw_read_overlap.tsv"
column -t "${OUT}/raw_read_index.tsv"
echo "[$(date)] done"
