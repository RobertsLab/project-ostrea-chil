#!/usr/bin/env bash
#SBATCH --job-name=13.7-bs-snp-batch
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# Run 13.7 on every region in a table and summarise per region.
# Usage (from the repo root):
#   sbatch code/13.7-wgbs-bs-snp-check-batch.sh output/13-wgbs/13.4-methylation/DMR_Qui_vs_Rio.tsv [flank] [pops]
# The table needs seqnames/start/end columns (13.4 DML/DMR output). Flank
# defaults to 0: for 1 kb windows, the question is the window itself. The
# populations compared default to the two in a <...>_<A>_vs_<B>.tsv filename.
# Writes output/13-wgbs/13.7-bs-snp-check/batch_<table>.tsv with one row per
# region joined to the input table's columns, plus a class per region:
#   repeat-like                 >= 50% of reads MAPQ < 10 (paralog/repeat mapping)
#   genetic: CpG-destroying     a differentiated C>T at C / G>A at G, which reads
#                               as lost methylation in the population carrying it
#   genetic: other              another differentiated variant (divergent allele)
#   no genetic signal           none of the above

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh
activate_env

table="${1:?give a table with seqnames/start/end columns}"
flank="${2:-0}"
compare="${3:-$(basename "${table}" .tsv | sed -nE 's/.*_([^_]+)_vs_([^_]+)$/\1,\2/p')}"
jobs="${SLURM_CPUS_PER_TASK:-4}"
name=$(basename "${table}" .tsv)
summary="${OUT}/13.7-bs-snp-check/batch_${name}.tsv"
mkdir -p "${OUT}/13.7-bs-snp-check"

# One region per line as chr:start-end, from whichever columns are named so.
regions=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
                      {print $c["seqnames"] ":" $c["start"] "-" $c["end"]}' "${table}")
echo "$(wc -l <<< "${regions}") regions from ${table}, flank ${flank}, comparing ${compare:-all populations}, ${jobs} at a time"

# Each region's output folder is named chr_start-end, so runs never collide.
xargs -P "${jobs}" -I{} bash -c \
  'bash code/13.7-wgbs-bs-snp-check.sh "$1" "$2" "$3" > /dev/null 2>&1 || echo "FAILED $1" >&2' _ {} "${flank}" "${compare}" \
  <<< "${regions}"

python3 -I - "${table}" "${OUT}/13.7-bs-snp-check" "${SAMPLES}" > "${summary}" <<'EOF'
import csv, sys
table, root, sheet = sys.argv[1:]
pop_of = {r["sample"]: r["population"] for r in csv.DictReader(open(sheet), delimiter="\t")}
pops = list(dict.fromkeys(pop_of.values()))
rows = list(csv.DictReader(open(table), delimiter="\t"))
cols = list(rows[0]) + ["n_variant_sites", "n_pop_differentiated", "n_cpg_destroying_differentiated",
                        "mapq_lt10_frac", "class"] + \
       [f"cpg_meth_{p}" for p in pops] + [f"CtoT_snp_cpgs_{p}" for p in pops] + [f"depth_{p}" for p in pops]
print("\t".join(cols))
for r in rows:
    wd = f"{root}/{r['seqnames']}_{r['start']}-{r['end']}"
    try:
        var = list(csv.DictReader(open(f"{wd}/variant_sites.tsv"), delimiter="\t"))
        cpg = list(csv.DictReader(open(f"{wd}/cpg_summary.tsv"), delimiter="\t"))
    except FileNotFoundError:
        print("\t".join([r[c] for c in rows[0]] + ["NA"] * (len(cols) - len(rows[0]))))
        continue
    diff = [v for v in var if v["pop_differentiated"] == "True"]
    # Only inside the target: with a flank, a SNP at a neighbouring CpG doesn't
    # explain this region's call (it still counts as "genetic: other").
    ctot = [v for v in diff if v["cpg_destroying"] == "True" and v["in_target"] == "True"]
    mq = list(csv.DictReader(open(f"{wd}/mapq.tsv"), delimiter="\t"))
    reads = sum(int(m["reads"]) for m in mq)
    low = f"{sum(int(m['mapq_lt10']) for m in mq) / reads:.2f}" if reads else ""
    cls = ("repeat-like" if low and float(low) >= 0.5 else
           "genetic: CpG-destroying" if ctot else
           "genetic: other" if diff else "no genetic signal")
    out = [r[c] for c in rows[0]] + [str(len(var)), str(len(diff)), str(len(ctot)), low, cls]
    def agg(pop, key, how):
        vals = [float(x[key]) for x in cpg if x["population"] == pop and x[key] not in ("", None)]
        if not vals: return ""
        return f"{sum(vals) / len(vals):.1f}" if how == "mean" else str(int(sum(vals)))
    out += [agg(p, "cpg_meth_pct_top_strand", "mean") for p in pops]
    out += [agg(p, "cpgs_with_CtoT_snp_evidence", "sum") for p in pops]
    out += [agg(p, "mean_depth", "mean") for p in pops]
    print("\t".join(out))
EOF

echo "Summary: ${summary}"
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
            {n[$c["class"]]++}
            END {for (k in n) printf "%5d  %s\n", n[k], k}' "${summary}" | sort -rn
