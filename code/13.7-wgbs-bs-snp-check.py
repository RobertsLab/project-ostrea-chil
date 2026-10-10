#!/usr/bin/env python3
"""Bisulfite-aware SNP check for one region, from strand-split mpileups.

Called by 13.7-wgbs-bs-snp-check.sh, which writes <sample>.CT.pileup and
<sample>.GA.pileup (Bismark XG:CT = original top strand, XG:GA = original
bottom strand) plus ref.seq for the padded region into --workdir.

A T at a reference C can be a SNP or an unmethylated C, so each base is only
counted from the strand where bisulfite cannot change it:
  XG:CT reads are ignored at reference C (C->T may be conversion)
  XG:GA reads are ignored at reference G (G->A may be conversion)
At a CpG C, bottom-strand reads therefore report the genotype (C, or T if
there is a real C->T SNP), while top-strand reads give methylation.

Writes to --workdir:
  variant_sites.tsv  positions with a non-reference base at >= --min-alt in
                     any oyster (informative depth >= --min-depth), with
                     per-population mean alt frequency and a flag for sites
                     that differ between populations by >= --pop-diff
  cpg_summary.tsv    per-oyster methylation over target CpGs (top-strand
                     reads), CpGs with C->T SNP evidence, and read depth
and prints a short summary.
"""
import argparse
import collections
import csv
import itertools
import re

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--workdir", required=True)
p.add_argument("--samples", required=True, help="data/13-wgbs-samples.tsv")
p.add_argument("--region-start", type=int, required=True, help="1-based start of ref.seq")
p.add_argument("--target-start", type=int, required=True, help="gene/region start (no flank)")
p.add_argument("--target-end", type=int, required=True)
p.add_argument("--min-depth", type=int, default=4)
p.add_argument("--min-alt", type=float, default=0.25)
p.add_argument("--pop-diff", type=float, default=0.5)
a = p.parse_args()

with open(a.samples) as fh:
    sheet = list(csv.DictReader(fh, delimiter="\t"))
samples = [r["sample"] for r in sheet]
pop_of = {r["sample"]: r["population"] for r in sheet}
pops = list(dict.fromkeys(pop_of.values()))
ref = open(f"{a.workdir}/ref.seq").read().strip().upper()


def parse_pileup(path):
    """pos -> (ref base, Counter of read bases)."""
    out = {}
    for line in open(path):
        _, pos, rb, _, bases = line.rstrip("\n").split("\t")[:5]
        rb = rb.upper()
        bases = re.sub(r"\^.", "", bases).replace("$", "")
        counts, i = collections.Counter(), 0
        while i < len(bases):
            ch = bases[i]
            if ch in "+-":  # indel: +3ACG / -2AT
                m = re.match(r"[+-](\d+)", bases[i:])
                i += len(m.group(0)) + int(m.group(1))
                continue
            if ch in ".,":
                counts[rb] += 1
            elif ch.upper() in "ACGT":
                counts[ch.upper()] += 1
            i += 1
        out[int(pos)] = (rb, counts)
    return out


informative = {}   # sample -> pos -> Counter
meth_calls = {}    # sample -> pos -> (C, T) from top-strand reads at reference C
depth = {}         # sample -> summed read depth over the target
for s in samples:
    ct = parse_pileup(f"{a.workdir}/{s}.CT.pileup")
    ga = parse_pileup(f"{a.workdir}/{s}.GA.pileup")
    inf = collections.defaultdict(collections.Counter)
    meth_calls[s] = {}
    for pos, (rb, c) in ct.items():
        if rb == "C":
            meth_calls[s][pos] = (c["C"], c["T"])
        else:
            inf[pos].update(c)
    for pos, (rb, c) in ga.items():
        if rb != "G":
            inf[pos].update(c)
    informative[s] = inf
    depth[s] = sum(sum(c.values()) for d in (ct, ga) for pos, (_, c) in d.items()
                   if a.target_start <= pos <= a.target_end)


def pop_mean(freqs, pop):
    v = [f for s, (n, f) in freqs.items() if pop_of[s] == pop and n >= a.min_depth]
    return sum(v) / len(v) if v else None


# ---- variant sites ----------------------------------------------------------------
variant_rows = []
for off, rb in enumerate(ref):
    pos = a.region_start + off
    counts = {s: informative[s].get(pos, collections.Counter()) for s in samples}
    total = sum(counts.values(), collections.Counter())
    alts = [b for b in "ACGT" if b != rb and total[b] > 0]
    if not alts:
        continue
    alt = max(alts, key=lambda b: total[b])
    freqs = {}
    for s, c in counts.items():
        n = sum(c.values())
        freqs[s] = (n, c[alt] / n if n else None)
    if not any(n >= a.min_depth and f >= a.min_alt for n, f in freqs.values() if f is not None):
        continue
    means = {pop: pop_mean(freqs, pop) for pop in pops}
    diffs = [abs(means[x] - means[y]) for x, y in itertools.combinations(pops, 2)
             if means[x] is not None and means[y] is not None]
    variant_rows.append({
        "pos": pos, "ref": rb, "alt": alt,
        "in_target": a.target_start <= pos <= a.target_end,
        "n_oysters_with_alt": sum(1 for n, f in freqs.values() if f and n >= a.min_depth and f >= a.min_alt),
        **{f"mean_alt_{pop}": "" if means[pop] is None else round(means[pop], 3) for pop in pops},
        "pop_differentiated": bool(diffs) and max(diffs) >= a.pop_diff,
        **{s: "" if f is None else f"{f:.2f}/{n}" for s, (n, f) in freqs.items()},
    })

with open(f"{a.workdir}/variant_sites.tsv", "w", newline="") as fh:
    cols = ["pos", "ref", "alt", "in_target", "n_oysters_with_alt"] + \
           [f"mean_alt_{pop}" for pop in pops] + ["pop_differentiated"] + samples
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(variant_rows)

# ---- target CpGs: methylation and C->T SNP evidence -----------------------------------
cpgs = [a.region_start + i for i in range(len(ref) - 1)
        if ref[i] == "C" and ref[i + 1] == "G"
        and a.target_start <= a.region_start + i <= a.target_end]
target_len = a.target_end - a.target_start + 1
cpg_rows = []
for s in samples:
    meth = calls = snp = 0
    for pos in cpgs:
        c, t = meth_calls[s].get(pos, (0, 0))
        meth += c
        calls += c + t
        g = informative[s].get(pos, collections.Counter())
        n = sum(g.values())
        if n >= 3 and g["T"] / n >= a.min_alt:
            snp += 1
    cpg_rows.append({
        "sample": s, "population": pop_of[s],
        "cpg_meth_pct_top_strand": "" if not calls else round(100 * meth / calls, 1),
        "meth_calls": calls,
        "cpgs_with_CtoT_snp_evidence": snp,
        "mean_depth": round(depth[s] / target_len, 1),
    })

with open(f"{a.workdir}/cpg_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(cpg_rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(cpg_rows)

# ---- printed summary ---------------------------------------------------------------------
diffd = [r for r in variant_rows if r["pop_differentiated"]]
shared = [r for r in variant_rows if r["n_oysters_with_alt"] >= 2]
print(f"Variant sites (alt >= {a.min_alt:.0%} at informative depth >= {a.min_depth} in any oyster): "
      f"{len(variant_rows)}")
print(f"  seen in >= 2 oysters: {len(shared)}")
print(f"  differ between populations (mean alt difference >= {a.pop_diff}): {len(diffd)}")
for r in diffd:
    print(f"    {r['pos']} {r['ref']}>{r['alt']}  " +
          "  ".join(f"{pop}={r[f'mean_alt_{pop}']}" for pop in pops))
print(f"\n{len(cpgs)} CpGs in target")
print(f"{'sample':8} {'pop':4} {'meth%':>6} {'calls':>6} {'CtoT_snp':>8} {'depth':>6}")
for r in cpg_rows:
    print(f"{r['sample']:8} {r['population']:4} {r['cpg_meth_pct_top_strand']!s:>6} "
          f"{r['meth_calls']:>6} {r['cpgs_with_CtoT_snp_evidence']:>8} {r['mean_depth']:>6}")
