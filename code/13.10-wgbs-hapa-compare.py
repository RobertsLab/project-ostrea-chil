#!/usr/bin/env python3
"""Compare the 13-wgbs results on the merged reference with a re-run on Och_HapA.

Called by 13.10-wgbs-hapa-compare.sh after it has lifted the merged-reference
regions to HapA with minimap2 (PAF with cg tags, in --workdir).

The merged reference is HapA's 2A/4A/6A/10A plus HapB's 1B/5B/6B/7B/8B/9B, so:
  1. qc_compare.tsv            per oyster: mapping, multi-mapping, duplication,
                               CpG methylation, CpGs >= 5x, on each reference
  2. shared_chrom_concordance.tsv
                               per oyster, on the four chromosomes identical in
                               both references: CpGs >= 5x in both, depth
                               change, and correlation of per-CpG methylation.
                               Differences here come only from reads that the
                               B chromosomes drew away (or didn't)
  3. lift_dmr.tsv / lift_dml.tsv
                               every 13.4 DMR window / DML CpG with its 13.9
                               class, where it lands on HapA, and the
                               Quilhua-vs-Rio difference measured the same way
                               (pooled reads per population, Rio - Qui) on both
                               references
  4. gn020540.tsv              per-oyster gene-body methylation on both references
and prints a summary.
"""
import argparse
import csv
import gzip
import math
import re
import statistics
from bisect import bisect_right
from collections import defaultdict

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--workdir", required=True)
p.add_argument("--samples", required=True)
p.add_argument("--merged-out", required=True, help="output/13-wgbs")
p.add_argument("--hapa-out", required=True, help="output/13-wgbs-hapa")
p.add_argument("--hapa-fa", required=True)
p.add_argument("--compare", default="Qui,Rio")
p.add_argument("--min-pop-depth", type=int, default=10,
               help="pooled CpG calls per population needed to measure a region")
a = p.parse_args()

SHARED = ["Chromosome_2A", "Chromosome_4A", "Chromosome_6A", "Chromosome_10A"]
pop_of = {r["sample"]: r["population"] for r in csv.DictReader(open(a.samples), delimiter="\t")}
samples = list(pop_of)
qpop, rpop = a.compare.split(",")
cov_path = lambda root, s: f"{root}/methylation/{s}.CpG_report.merged_CpG_evidence.cov.gz"


def read_report(path, key):
    for line in open(path):
        if line.startswith(key):
            return line.rstrip("\n").split("\t")[1].strip().rstrip("%")
    return ""


# ---- 1. QC --------------------------------------------------------------------------
qc_rows = []
for s in samples:
    row = {"sample": s, "population": pop_of[s]}
    for tag, root in (("merged", a.merged_out), ("hapa", a.hapa_out)):
        rep = f"{root}/bismark/{s}_PE_report.txt"
        pairs = int(read_report(rep, "Sequence pairs analysed in total"))
        uniq = int(read_report(rep, "Number of paired-end alignments with a unique best hit"))
        ambig = int(read_report(rep, "Sequence pairs did not map uniquely"))
        dd = open(f"{root}/dedup/{s}_pe.deduplication_report.txt").read()
        dup = re.search(r"removed:\s+\d+ \(([0-9.]+)%\)", dd).group(1)
        sp = f"{root}/methylation/{s}_pe.deduplicated_splitting_report.txt"
        n5 = 0
        with gzip.open(cov_path(root, s), "rt") as fh:
            for line in fh:
                f = line.split("\t")
                if int(f[4]) + int(f[5]) >= 5:
                    n5 += 1
        row.update({f"{tag}_mapping_pct": round(100 * uniq / pairs, 1),
                    f"{tag}_multimap_pct": round(100 * ambig / pairs, 1),
                    f"{tag}_dup_pct": dup,
                    f"{tag}_cpg_meth_pct": read_report(sp, "C methylated in CpG context"),
                    f"{tag}_cpgs_5x": n5})
    qc_rows.append(row)


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


write(f"{a.workdir}/qc_compare.tsv", qc_rows)

# ---- 2. shared chromosomes ------------------------------------------------------------
def load_shared(path):
    d = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[0] in SHARED:
                d[(f[0], int(f[1]))] = (int(f[4]), int(f[4]) + int(f[5]))
    return d


def pearson(x, y):
    if len(x) < 3:
        return float("nan")
    mx, my = statistics.fmean(x), statistics.fmean(y)
    sxy = sum((u - mx) * (v - my) for u, v in zip(x, y))
    sx = math.sqrt(sum((u - mx) ** 2 for u in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    return sxy / (sx * sy) if sx and sy else float("nan")


shared_rows = []
for s in samples:
    m, h = load_shared(cov_path(a.merged_out, s)), load_shared(cov_path(a.hapa_out, s))
    both = [k for k in m.keys() & h.keys() if m[k][1] >= 5 and h[k][1] >= 5]
    xm = [100 * m[k][0] / m[k][1] for k in both]
    xh = [100 * h[k][0] / h[k][1] for k in both]
    dm = sum(v[1] for v in m.values())
    dh = sum(v[1] for v in h.values())
    shared_rows.append({"sample": s, "population": pop_of[s],
                        "cpgs_5x_merged": sum(1 for v in m.values() if v[1] >= 5),
                        "cpgs_5x_hapa": sum(1 for v in h.values() if v[1] >= 5),
                        "cpgs_5x_both": len(both),
                        "depth_ratio_hapa_over_merged": round(dh / dm, 3) if dm else "",
                        "meth_pearson_r": round(pearson(xm, xh), 4),
                        "mean_abs_meth_diff": round(statistics.fmean(abs(u - v) for u, v in zip(xm, xh)), 2)
                        if both else ""})
write(f"{a.workdir}/shared_chrom_concordance.tsv", shared_rows)

# ---- 3./4. lifted regions ------------------------------------------------------------------
def load_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def parse_paf(path):
    """Best primary hit per query (highest matching bases)."""
    best = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        tags = dict(t.split(":", 2)[::2] for t in f[12:])
        if tags.get("tp") != "P":
            continue
        hit = {"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
               "t": f[5], "ts": int(f[7]), "te": int(f[8]), "match": int(f[9]), "alen": int(f[10]),
               "mapq": int(f[11]), "cg": tags.get("cg", "")}
        if f[0] not in best or hit["match"] > best[f[0]]["match"]:
            best[f[0]] = hit
    return best


def qpos_to_t(hit, q):
    """0-based query position -> 0-based target position, None if in an indel/outside."""
    if not (hit["qs"] <= q < hit["qe"]):
        return None
    # walk the CIGAR in query orientation; for '-' the alignment runs backwards on the target
    qi = hit["qs"] if hit["strand"] == "+" else hit["qlen"] - hit["qe"]
    qq = q if hit["strand"] == "+" else hit["qlen"] - 1 - q
    ti = hit["ts"]
    for n, op in re.findall(r"(\d+)([MIDNX=])", hit["cg"]):
        n = int(n)
        if op in "M=X":
            if qi <= qq < qi + n:
                return ti + (qq - qi)
            qi += n; ti += n
        elif op == "I":
            if qi <= qq < qi + n:
                return None
            qi += n
        elif op in "DN":
            ti += n
    return None


class CovIndex:
    """CpG calls per sample on the chromosomes/regions of interest, queryable by interval."""

    def __init__(self, root, wanted):
        self.d = {}  # sample -> chrom -> (sorted starts, [(meth, total)])
        for s in samples:
            per = defaultdict(list)
            with gzip.open(cov_path(root, s), "rt") as fh:
                for line in fh:
                    f = line.split("\t", 6)
                    if f[0] in wanted:
                        per[f[0]].append((int(f[1]), int(f[4]), int(f[4]) + int(f[5])))
            self.d[s] = {c: ([x[0] for x in v], [(x[1], x[2]) for x in v]) for c, v in per.items()}

    def pooled(self, chrom, start, end, pop):
        """Pooled (meth, total) for a population over CpG starts in [start, end] (1-based)."""
        m = t = 0
        for s in samples:
            if pop_of[s] != pop or chrom not in self.d[s]:
                continue
            starts, vals = self.d[s][chrom]
            i = bisect_right(starts, start - 1)
            while i < len(starts) and starts[i] <= end:
                m += vals[i][0]; t += vals[i][1]; i += 1
        return m, t

    def per_sample(self, chrom, start, end):
        out = {}
        for s in samples:
            m = t = 0
            if chrom in self.d[s]:
                starts, vals = self.d[s][chrom]
                i = bisect_right(starts, start - 1)
                while i < len(starts) and starts[i] <= end:
                    m += vals[i][0]; t += vals[i][1]; i += 1
            out[s] = (m, t)
        return out


def diff(idx, chrom, start, end):
    mq, tq = idx.pooled(chrom, start, end, qpop)
    mr, tr = idx.pooled(chrom, start, end, rpop)
    if tq < a.min_pop_depth or tr < a.min_pop_depth:
        return None, tq, tr
    return 100 * mr / tr - 100 * mq / tq, tq, tr


regions = {}  # name -> dict (from the .sh: regions.tsv)
for r in csv.DictReader(open(f"{a.workdir}/regions.tsv"), delimiter="\t"):
    regions[r["name"]] = r
hits = parse_paf(f"{a.workdir}/regions_to_hapa.paf")
hapa = load_fasta(a.hapa_fa)

lifted = {}
for name, r in regions.items():
    h = hits.get(name)
    if not h:
        continue
    qcov = (h["qe"] - h["qs"]) / h["qlen"]
    ident = h["match"] / h["alen"] if h["alen"] else 0
    if r["kind"] == "DML":
        q = int(r["cpg_offset"])           # 0-based offset of the CpG's C in the query
        t1, t2 = qpos_to_t(h, q), qpos_to_t(h, q + 1)
        if t1 is None or t2 is None:
            lifted[name] = dict(h, qcov=qcov, ident=ident, ok=False, why="CpG in an indel")
            continue
        tc = min(t1, t2)                   # 0-based C of the CpG on HapA's + strand
        cg = hapa[h["t"]][tc:tc + 2]
        lifted[name] = dict(h, qcov=qcov, ident=ident, ok=True, tstart=tc + 1, tend=tc + 2, cg=cg)
    else:
        lifted[name] = dict(h, qcov=qcov, ident=ident, ok=qcov >= 0.5,
                            why="" if qcov >= 0.5 else "< 50% of query aligned",
                            tstart=h["ts"] + 1, tend=h["te"])

wanted_m = {r["chrom"] for r in regions.values()}
wanted_h = {h["t"] for h in lifted.values()}
idx_m = CovIndex(a.merged_out, wanted_m)
idx_h = CovIndex(a.hapa_out, wanted_h)


def lift_table(kind):
    rows = []
    for name, r in regions.items():
        if r["kind"] != kind:
            continue
        s, e = int(r["start"]), int(r["end"])
        dm, tqm, trm = diff(idx_m, r["chrom"], s, e)
        row = {k: r[k] for k in ("name", "chrom", "start", "end", "meth_diff_13.4", "feature", "class")}
        row.update({"merged_diff_rio_minus_qui": "" if dm is None else round(dm, 1)})
        L = lifted.get(name)
        if not L:
            row.update({"hapa_chrom": "", "hapa_start": "", "hapa_end": "", "identity": "", "query_cov": "",
                        "hapa_cpg": "", "lift": "unmapped", "hapa_diff_rio_minus_qui": ""})
        elif not L["ok"]:
            row.update({"hapa_chrom": L["t"], "hapa_start": "", "hapa_end": "",
                        "identity": round(L["ident"], 3), "query_cov": round(L["qcov"], 2),
                        "hapa_cpg": "", "lift": L["why"], "hapa_diff_rio_minus_qui": ""})
        else:
            dh, _, _ = diff(idx_h, L["t"], L["tstart"], L["tend"])
            row.update({"hapa_chrom": L["t"], "hapa_start": L["tstart"], "hapa_end": L["tend"],
                        "identity": round(L["ident"], 3), "query_cov": round(L["qcov"], 2),
                        "hapa_cpg": L.get("cg", ""),
                        "lift": "ok" if kind != "DML" or L.get("cg") == "CG" else "no CpG on HapA",
                        "hapa_diff_rio_minus_qui": "" if dh is None else round(dh, 1)})
        rows.append(row)
    write(f"{a.workdir}/lift_{kind.lower()}.tsv", rows)
    return rows


def summarise(kind, rows, threshold):
    print(f"\n{kind}: {len(rows)} regions")
    lift_counts = defaultdict(int)
    for r in rows:
        lift_counts[r["lift"]] += 1
    print("  lift to HapA: " + ", ".join(f"{k} {v}" for k, v in sorted(lift_counts.items(), key=lambda x: -x[1])))
    print(f"  {'class':26s} {'n':>5} {'both':>5} {'same sign':>9} {'replicated':>10} {'r':>6}")
    by = defaultdict(list)
    for r in rows:
        by[r["class"]].append(r)
    for cls in ["no genetic signal", "untested", "genetic: CpG-destroying", "repeat-like", "genetic: other"]:
        rs = by.get(cls, [])
        pairs = [(float(r["merged_diff_rio_minus_qui"]), float(r["hapa_diff_rio_minus_qui"]))
                 for r in rs if r["merged_diff_rio_minus_qui"] != "" and r["hapa_diff_rio_minus_qui"] != ""]
        same = sum(1 for m, h in pairs if m * h > 0)
        rep = sum(1 for m, h in pairs if m * h > 0 and abs(h) >= threshold)
        r_ = pearson([m for m, _ in pairs], [h for _, h in pairs])
        pct = lambda k: f"{k} ({100 * k / len(pairs):.0f}%)" if pairs else "-"
        print(f"  {cls:26s} {len(rs):>5} {len(pairs):>5} {pct(same):>9} {pct(rep):>10} {r_:6.2f}")
    print(f"  (both = measurable on both references; replicated = same sign and |HapA diff| >= {threshold}%)")


dmr_rows = lift_table("DMR")
dml_rows = lift_table("DML")

# GN020540
g = next(r for r in regions.values() if r["kind"] == "gene")
L = lifted.get(g["name"])
gm = idx_m.per_sample(g["chrom"], int(g["start"]), int(g["end"]))
gh = idx_h.per_sample(L["t"], L["tstart"], L["tend"]) if L and L["ok"] else {}
gene_rows = [{"sample": s, "population": pop_of[s],
              "merged_meth_pct": round(100 * gm[s][0] / gm[s][1], 1) if gm[s][1] else "",
              "merged_calls": gm[s][1],
              "hapa_meth_pct": round(100 * gh[s][0] / gh[s][1], 1) if gh.get(s, (0, 0))[1] else "",
              "hapa_calls": gh.get(s, (0, 0))[1]} for s in samples]
write(f"{a.workdir}/gn020540.tsv", gene_rows)

# ---- printed summary ---------------------------------------------------------------------
print("QC (merged -> HapA), mean over oysters:")
for k in ("mapping_pct", "multimap_pct", "cpgs_5x"):
    mv = statistics.fmean(float(r[f"merged_{k}"]) for r in qc_rows)
    hv = statistics.fmean(float(r[f"hapa_{k}"]) for r in qc_rows)
    print(f"  {k:14s} {mv:14,.1f} -> {hv:14,.1f}")
print("\nShared chromosomes (2A/4A/6A/10A, identical in both):")
for r in shared_rows:
    print(f"  {r['sample']:5s} CpGs>=5x both {r['cpgs_5x_both']:>9,}  depth ratio {r['depth_ratio_hapa_over_merged']}"
          f"  r {r['meth_pearson_r']}  mean |diff| {r['mean_abs_meth_diff']}%")
summarise("DMR windows", dmr_rows, 10)
summarise("DML CpGs", dml_rows, 20)
if L:
    print(f"\nGN020540 -> {L['t']}:{L.get('tstart')}-{L.get('tend')} ({L['strand']}), identity {L['ident']:.3f}, "
          f"query covered {L['qcov']:.2f}")
for r in gene_rows:
    print(f"  {r['sample']:5s} merged {r['merged_meth_pct']!s:>5}% ({r['merged_calls']:>4} calls)   "
          f"HapA {r['hapa_meth_pct']!s:>5}% ({r['hapa_calls']:>4} calls)")
