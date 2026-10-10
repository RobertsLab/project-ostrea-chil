#!/usr/bin/env bash
#SBATCH --job-name=13.10-hapa-compare
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=output/13-wgbs/logs/%x_%j.out

# Does anything change if the WGBS reads are aligned to Och_HapA alone instead
# of the merged reference (HapA 2A/4A/6A/10A + HapB 1B/5B/6B/7B/8B/9B)?
#
# Needs both runs finished:  bash code/13-wgbs-submit.sh  and
#                            WGBS_REF=hapa bash code/13-wgbs-submit.sh
# plus the 13.9 batch tables for the merged run. Usage, from the repo root:
#   sbatch code/13.10-wgbs-hapa-compare.sh
#
# Lifts every 13.4 DMR window, every DML CpG (+/- 100 bp) and GN020540 from the
# merged reference to HapA with minimap2, then 13.10-wgbs-hapa-compare.py
# compares QC, the four shared chromosomes, and the Qui-vs-Rio differences.
# Output: output/13-wgbs-hapa/13.10-compare/

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source code/13.0-wgbs-config.sh          # merged settings (default WGBS_REF)
activate_env
MERGED_OUT="${OUT}"; MERGED_FA="${GENOME_FA}"; MERGED_GFF="${GFF}"
WGBS_REF=hapa source code/13.0-wgbs-config.sh
HAPA_OUT="${OUT}"; HAPA_FA="${GENOME_FA}"

WD="${HAPA_OUT}/13.10-compare"
mkdir -p "${WD}"
CHECK="${MERGED_OUT}/13.9-bs-snp-check"

# regions.tsv: name kind chrom start end cpg_offset meth_diff_13.4 feature class
# query.fa:    the merged-reference sequence of each region, named by `name`
python3 -I - "${CHECK}" "${MERGED_GFF}" "${WD}" <<'EOF'
import csv, sys
check, gff, wd = sys.argv[1:]
rows = []
for kind, f in (("DMR", "batch_DMR_Qui_vs_Rio.tsv"), ("DML", "batch_DML_Qui_vs_Rio.tsv")):
    for r in csv.DictReader(open(f"{check}/{f}"), delimiter="\t"):
        c, s, e = r["seqnames"], int(r["start"]), int(r["end"])
        if kind == "DML":                       # +/- 100 bp around the CpG; C at offset 100
            qs = max(1, s - 100); off = s - qs; qe = s + 101
        else:
            qs, qe, off = s, e, ""
        rows.append({"name": f"{kind}:{c}:{s}-{e}", "kind": kind, "chrom": c, "start": s, "end": e,
                     "qs": qs, "qe": qe, "cpg_offset": off, "meth_diff_13.4": r["meth.diff"],
                     "feature": r["feature"], "class": r["class"]})
for line in open(gff):
    f = line.split("\t")
    if len(f) > 8 and f[2] == "gene" and "ID=GN020540;" in f[8]:
        rows.append({"name": "gene:GN020540", "kind": "gene", "chrom": f[0], "start": int(f[3]),
                     "end": int(f[4]), "qs": int(f[3]), "qe": int(f[4]), "cpg_offset": "",
                     "meth_diff_13.4": "", "feature": "gene", "class": ""})
with open(f"{wd}/regions.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(f"{wd}/faidx_regions.txt", "w") as fh:
    for r in rows: fh.write(f"{r['chrom']}:{r['qs']}-{r['qe']}\n")
with open(f"{wd}/query_names.txt", "w") as fh:
    for r in rows: fh.write(r["name"] + "\n")
print(f"{len(rows)} regions")
EOF

# samtools writes the regions in the order given; rename headers to region names.
samtools faidx "${MERGED_FA}" -r "${WD}/faidx_regions.txt" |
  awk -v names="${WD}/query_names.txt" 'BEGIN {while ((getline n < names) > 0) N[++i] = n}
                                        /^>/ {print ">" N[++j]; next} {print}' > "${WD}/query.fa"

# Windows and the gene: assembly-to-assembly preset (HapB->HapA can be divergent).
# Single CpGs (202 bp): short-read preset. Both with base-level CIGARs (-c).
awk '/^>/ {keep = ($0 !~ /^>DML:/)} keep' "${WD}/query.fa" > "${WD}/query_long.fa"
awk '/^>/ {keep = ($0 ~ /^>DML:/)} keep' "${WD}/query.fa" > "${WD}/query_dml.fa"
minimap2 -c -x asm20 -t "${SLURM_CPUS_PER_TASK:-8}" "${HAPA_FA}" "${WD}/query_long.fa" > "${WD}/long.paf"
minimap2 -c -x sr    -t "${SLURM_CPUS_PER_TASK:-8}" "${HAPA_FA}" "${WD}/query_dml.fa"  > "${WD}/dml.paf"
cat "${WD}/long.paf" "${WD}/dml.paf" > "${WD}/regions_to_hapa.paf"

python3 -I code/13.10-wgbs-hapa-compare.py \
  --workdir "${WD}" --samples "${SAMPLES}" \
  --merged-out "${MERGED_OUT}" --hapa-out "${HAPA_OUT}" --hapa-fa "${HAPA_FA}" \
  | tee "${WD}/summary.txt"
