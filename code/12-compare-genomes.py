#!/usr/bin/env python3
"""Download and compare the three Ostrea chilensis genome assemblies.

The script is designed to be run from anywhere. Inputs are stored in the
repository's data directory and results are written to a directory whose name
matches this script's stem: output/12-compare-genomes/.

Only Python's standard library is required for downloads and summaries.
minimap2 is required for the whole-genome comparisons.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class Assembly:
    key: str
    label: str
    filename: str
    fasta_url: str

    @property
    def fai_url(self) -> str:
        return self.fasta_url + ".fai"


ASSEMBLIES = (
    Assembly(
        "assembly1",
        "Assembly 1.0",
        "merged_out.fasta.TBtools.fa",
        "https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/"
        "OCEAN/docs/jbrowse/data/v1/merged_out.fasta.TBtools.fa",
    ),
    Assembly(
        "hapa",
        "Och_HapA",
        "Och_HapA_assembly.fa",
        "https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/"
        "OCEAN/docs/jbrowse/data/HapA/Och_HapA_assembly.fa",
    ),
    Assembly(
        "hapb",
        "Och_HapB",
        "Och_HapB_assembly.fa",
        "https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/"
        "OCEAN/docs/jbrowse/data/HapB/Och_HapB_assembly.fa",
    ),
)

# minimap2 target, query, and an output-safe comparison name.
COMPARISONS = (
    ("assembly1", "hapa", "assembly1_vs_hapa"),
    ("assembly1", "hapb", "assembly1_vs_hapb"),
    ("hapa", "hapb", "hapa_vs_hapb"),
)


def parse_args() -> argparse.Namespace:
    script_stem = Path(__file__).stem
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Download, summarize, align, and plot three genome assemblies."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=repo_root / "data",
        help="FASTA download directory (default: repository data/)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "output" / script_stem,
        help=f"result directory (default: output/{script_stem}/)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, min(16, os.cpu_count() or 1)),
        help="threads supplied to minimap2 (default: up to 16)",
    )
    parser.add_argument(
        "--minimap2",
        default="minimap2",
        help="minimap2 executable or path (default: minimap2)",
    )
    parser.add_argument(
        "--stats-only",
        action="store_true",
        help="download only the small .fai indexes and write assembly statistics",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="use existing local FASTA and .fai files without network access",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="rerun alignments even when completed .paf.gz files exist",
    )
    parser.add_argument(
        "--plot-min-block",
        type=int,
        default=100_000,
        help="minimum PAF block length drawn in dot plots (default: 100000)",
    )
    parser.add_argument(
        "--max-plot-alignments",
        type=int,
        default=20_000,
        help="maximum longest alignments drawn per dot plot (default: 20000)",
    )
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.plot_min_block < 0 or args.max_plot_alignments < 1:
        parser.error("plot limits must be positive")
    return args


def run_curl(url: str, destination: Path) -> None:
    """Download to a partial file, resuming when possible, then rename atomically."""
    curl = shutil.which("curl")
    if curl is None:
        raise RuntimeError("curl is required to download the genome files")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    print(f"Downloading {url}\n        -> {destination}", flush=True)
    command = [
        curl,
        "--fail",
        "--location",
        "--retry",
        "3",
        "--retry-delay",
        "5",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        url,
    ]
    subprocess.run(command, check=True)
    partial.replace(destination)


def ensure_file(url: str, destination: Path, skip_download: bool) -> None:
    if destination.is_file() and destination.stat().st_size > 0:
        print(f"Using existing {destination}")
        return
    if skip_download:
        raise FileNotFoundError(f"required local file is missing: {destination}")
    run_curl(url, destination)


def read_fai(path: Path) -> List[Tuple[str, int, int, int, int]]:
    records: List[Tuple[str, int, int, int, int]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"invalid FAI record at {path}:{line_number}")
            records.append(
                (fields[0], int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4]))
            )
    if not records:
        raise ValueError(f"FAI has no records: {path}")
    return records


def validate_fasta_against_fai(fasta: Path, fai_records: Sequence[Tuple[str, int, int, int, int]]) -> None:
    """Catch interrupted downloads using the final byte described by the FAI."""
    _name, length, offset, line_bases, line_width = fai_records[-1]
    if length == 0:
        minimum_bytes = offset
    else:
        minimum_bytes = (
            offset
            + ((length - 1) // line_bases) * line_width
            + ((length - 1) % line_bases)
            + 1
        )
    observed_bytes = fasta.stat().st_size
    if observed_bytes < minimum_bytes:
        raise ValueError(
            f"{fasta} is incomplete: {observed_bytes:,} bytes; "
            f"the FAI requires at least {minimum_bytes:,}"
        )


def nx(lengths: Sequence[int], fraction: float) -> Tuple[int, int]:
    threshold = sum(lengths) * fraction
    cumulative = 0
    for rank, length in enumerate(sorted(lengths, reverse=True), 1):
        cumulative += length
        if cumulative >= threshold:
            return length, rank
    raise ValueError("cannot calculate Nx for an empty assembly")


def write_assembly_tables(
    output_dir: Path,
    fai_by_key: Dict[str, List[Tuple[str, int, int, int, int]]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    stats_path = output_dir / "assembly_stats.tsv"
    lengths_path = output_dir / "sequence_lengths.tsv"
    with stats_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "assembly",
                "sequence_count",
                "total_bp",
                "shortest_bp",
                "longest_bp",
                "N50_bp",
                "L50",
                "N90_bp",
                "L90",
            ]
        )
        for assembly in ASSEMBLIES:
            lengths = [record[1] for record in fai_by_key[assembly.key]]
            n50, l50 = nx(lengths, 0.50)
            n90, l90 = nx(lengths, 0.90)
            writer.writerow(
                [
                    assembly.label,
                    len(lengths),
                    sum(lengths),
                    min(lengths),
                    max(lengths),
                    n50,
                    l50,
                    n90,
                    l90,
                ]
            )
    with lengths_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["assembly", "sequence", "length_bp"])
        for assembly in ASSEMBLIES:
            for name, length, _offset, _line_bases, _line_width in fai_by_key[assembly.key]:
                writer.writerow([assembly.label, name, length])
    print(f"Wrote {stats_path}")
    print(f"Wrote {lengths_path}")


def run_alignment(
    minimap2: str,
    threads: int,
    target: Path,
    query: Path,
    paf_gz: Path,
    log_path: Path,
    force: bool,
) -> None:
    if paf_gz.is_file() and paf_gz.stat().st_size > 0 and not force:
        print(f"Using existing alignment {paf_gz}")
        return
    executable = shutil.which(minimap2) if os.sep not in minimap2 else minimap2
    if not executable or not Path(executable).is_file():
        raise RuntimeError(
            f"minimap2 executable not found: {minimap2!r}. "
            "Install minimap2 or pass --minimap2 /full/path/to/minimap2."
        )
    paf_gz.parent.mkdir(parents=True, exist_ok=True)
    temporary = paf_gz.with_name(paf_gz.name + ".part")
    command = [
        str(executable),
        "-x",
        "asm5",
        "--secondary=no",
        "-t",
        str(threads),
        str(target),
        str(query),
    ]
    print("Running " + " ".join(command), flush=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log_handle)
        assert process.stdout is not None
        with gzip.open(temporary, "wb", compresslevel=6) as paf_handle:
            shutil.copyfileobj(process.stdout, paf_handle, length=1024 * 1024)
        return_code = process.wait()
    if return_code != 0:
        temporary.unlink(missing_ok=True)
        raise subprocess.CalledProcessError(return_code, command)
    temporary.replace(paf_gz)


def merge_interval_length(intervals: Iterable[Tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def parse_paf(
    path: Path,
) -> Tuple[
    List[Tuple[str, int, int, int, str, str, int, int, int, int, int, int]],
    Dict[str, object],
    Dict[Tuple[str, str, str], List[int]],
]:
    alignments = []
    query_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    target_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    correspondence: Dict[Tuple[str, str, str], List[int]] = defaultdict(
        lambda: [0, 0, 0]
    )
    matching_bp = 0
    block_bp = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"invalid PAF record at {path}:{line_number}")
            record = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                int(fields[3]),
                fields[4],
                fields[5],
                int(fields[6]),
                int(fields[7]),
                int(fields[8]),
                int(fields[9]),
                int(fields[10]),
                int(fields[11]),
            )
            alignments.append(record)
            qname, _qlen, qstart, qend, strand, tname, _tlen, tstart, tend, matches, block, _mapq = record
            query_intervals[qname].append((qstart, qend))
            target_intervals[tname].append((tstart, tend))
            matching_bp += matches
            block_bp += block
            values = correspondence[(qname, tname, strand)]
            values[0] += 1
            values[1] += matches
            values[2] += block
    summary: Dict[str, object] = {
        "alignment_count": len(alignments),
        "matching_bp": matching_bp,
        "alignment_block_bp": block_bp,
        "identity_pct": (100.0 * matching_bp / block_bp) if block_bp else 0.0,
        "query_aligned_bp": sum(merge_interval_length(v) for v in query_intervals.values()),
        "target_aligned_bp": sum(merge_interval_length(v) for v in target_intervals.values()),
    }
    return alignments, summary, correspondence


def cumulative_offsets(
    records: Sequence[Tuple[str, int, int, int, int]]
) -> Tuple[Dict[str, int], int]:
    offsets: Dict[str, int] = {}
    position = 0
    for name, length, _offset, _line_bases, _line_width in records:
        offsets[name] = position
        position += length
    return offsets, position


def write_dotplot(
    output_path: Path,
    title: str,
    alignments: Sequence[Tuple[str, int, int, int, str, str, int, int, int, int, int, int]],
    target_records: Sequence[Tuple[str, int, int, int, int]],
    query_records: Sequence[Tuple[str, int, int, int, int]],
    minimum_block: int,
    maximum_alignments: int,
) -> None:
    width, height = 1200, 1000
    left, right, top, bottom = 135, 40, 65, 135
    plot_width = width - left - right
    plot_height = height - top - bottom
    target_offsets, target_total = cumulative_offsets(target_records)
    query_offsets, query_total = cumulative_offsets(query_records)
    selected = [record for record in alignments if record[10] >= minimum_block]
    selected.sort(key=lambda record: record[10], reverse=True)
    selected = selected[:maximum_alignments]

    def xcoord(value: int) -> float:
        return left + plot_width * value / target_total

    def ycoord(value: int) -> float:
        return top + plot_height * value / query_total

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="20">{html.escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="#fafafa" stroke="#333"/>',
    ]
    for record in selected:
        qname, _qlen, qstart, qend, strand, tname, _tlen, tstart, tend, _matches, _block, _mapq = record
        if qname not in query_offsets or tname not in target_offsets:
            continue
        x1 = xcoord(target_offsets[tname] + tstart)
        x2 = xcoord(target_offsets[tname] + tend)
        if strand == "+":
            y1 = ycoord(query_offsets[qname] + qstart)
            y2 = ycoord(query_offsets[qname] + qend)
            color = "#1769aa"
        else:
            y1 = ycoord(query_offsets[qname] + qend)
            y2 = ycoord(query_offsets[qname] + qstart)
            color = "#d1495b"
        elements.append(
            f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="1" stroke-opacity="0.55"/>'
        )
    for name, length, _offset, _line_bases, _line_width in target_records:
        start = target_offsets[name]
        x = xcoord(start)
        midpoint = xcoord(start + length // 2)
        elements.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" stroke="#d8d8d8"/>')
        elements.append(
            f'<text x="{midpoint:.2f}" y="{top + plot_height + 10}" transform="rotate(55 {midpoint:.2f} {top + plot_height + 10})" text-anchor="start" font-family="sans-serif" font-size="11">{html.escape(name)}</text>'
        )
    for name, length, _offset, _line_bases, _line_width in query_records:
        start = query_offsets[name]
        y = ycoord(start)
        midpoint = ycoord(start + length // 2)
        elements.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#d8d8d8"/>')
        elements.append(
            f'<text x="{left - 8}" y="{midpoint + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="11">{html.escape(name)}</text>'
        )
    elements.extend(
        [
            f'<text x="{left + plot_width / 2}" y="{height - 18}" text-anchor="middle" font-family="sans-serif" font-size="15">Target genome</text>',
            f'<text x="22" y="{top + plot_height / 2}" transform="rotate(-90 22 {top + plot_height / 2})" text-anchor="middle" font-family="sans-serif" font-size="15">Query genome</text>',
            f'<line x1="{width - 235}" y1="43" x2="{width - 205}" y2="43" stroke="#1769aa" stroke-width="3"/><text x="{width - 198}" y="47" font-family="sans-serif" font-size="12">forward</text>',
            f'<line x1="{width - 120}" y1="43" x2="{width - 90}" y2="43" stroke="#d1495b" stroke-width="3"/><text x="{width - 83}" y="47" font-family="sans-serif" font-size="12">reverse</text>',
            f'<text x="{left}" y="{height - 2}" font-family="sans-serif" font-size="10" fill="#555">Showing {len(selected):,} alignments with block length ≥ {minimum_block:,} bp</text>',
            "</svg>",
        ]
    )
    output_path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_comparison_tables(
    output_dir: Path,
    summaries: Sequence[Dict[str, object]],
    correspondences: Sequence[Tuple[str, Dict[Tuple[str, str, str], List[int]]]],
) -> None:
    summary_path = output_dir / "pairwise_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "comparison",
            "target",
            "query",
            "alignment_count",
            "matching_bp",
            "alignment_block_bp",
            "identity_pct",
            "target_aligned_bp",
            "target_total_bp",
            "target_coverage_pct",
            "query_aligned_bp",
            "query_total_bp",
            "query_coverage_pct",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for summary in summaries:
            formatted = dict(summary)
            for key in ("identity_pct", "target_coverage_pct", "query_coverage_pct"):
                formatted[key] = f"{float(formatted[key]):.4f}"
            writer.writerow(formatted)

    correspondence_path = output_dir / "sequence_correspondence.tsv"
    with correspondence_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "comparison",
                "query_sequence",
                "target_sequence",
                "strand",
                "alignment_count",
                "matching_bp",
                "alignment_block_bp",
                "identity_pct",
            ]
        )
        for comparison_name, correspondence in correspondences:
            rows = sorted(
                correspondence.items(), key=lambda item: item[1][2], reverse=True
            )
            for (query_name, target_name, strand), (count, matches, block) in rows:
                identity = 100.0 * matches / block if block else 0.0
                writer.writerow(
                    [
                        comparison_name,
                        query_name,
                        target_name,
                        strand,
                        count,
                        matches,
                        block,
                        f"{identity:.4f}",
                    ]
                )
    print(f"Wrote {summary_path}")
    print(f"Wrote {correspondence_path}")


def write_manifest(output_dir: Path, data_dir: Path) -> None:
    manifest_path = output_dir / "input_manifest.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["assembly", "local_fasta", "fasta_bytes", "source_url"])
        for assembly in ASSEMBLIES:
            fasta = data_dir / assembly.filename
            writer.writerow(
                [
                    assembly.label,
                    fasta,
                    fasta.stat().st_size if fasta.is_file() else "not_downloaded",
                    assembly.fasta_url,
                ]
            )


def main() -> int:
    args = parse_args()
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fai_by_key: Dict[str, List[Tuple[str, int, int, int, int]]] = {}
    assembly_by_key = {assembly.key: assembly for assembly in ASSEMBLIES}
    for assembly in ASSEMBLIES:
        fai_path = args.data_dir / (assembly.filename + ".fai")
        ensure_file(assembly.fai_url, fai_path, args.skip_download)
        fai_by_key[assembly.key] = read_fai(fai_path)

    write_assembly_tables(args.output_dir, fai_by_key)
    if args.stats_only:
        write_manifest(args.output_dir, args.data_dir)
        print("Stats-only run complete; FASTA files and alignments were not requested.")
        return 0

    for assembly in ASSEMBLIES:
        fasta_path = args.data_dir / assembly.filename
        ensure_file(assembly.fasta_url, fasta_path, args.skip_download)
        validate_fasta_against_fai(fasta_path, fai_by_key[assembly.key])

    summaries: List[Dict[str, object]] = []
    correspondences = []
    alignments_dir = args.output_dir / "alignments"
    logs_dir = args.output_dir / "logs"
    plots_dir = args.output_dir / "plots"
    for directory in (alignments_dir, logs_dir, plots_dir):
        directory.mkdir(parents=True, exist_ok=True)

    for target_key, query_key, comparison_name in COMPARISONS:
        target = assembly_by_key[target_key]
        query = assembly_by_key[query_key]
        paf_gz = alignments_dir / f"{comparison_name}.paf.gz"
        run_alignment(
            args.minimap2,
            args.threads,
            args.data_dir / target.filename,
            args.data_dir / query.filename,
            paf_gz,
            logs_dir / f"{comparison_name}.minimap2.log",
            args.force,
        )
        alignments, summary, correspondence = parse_paf(paf_gz)
        target_total = sum(record[1] for record in fai_by_key[target_key])
        query_total = sum(record[1] for record in fai_by_key[query_key])
        summary.update(
            {
                "comparison": comparison_name,
                "target": target.label,
                "query": query.label,
                "target_total_bp": target_total,
                "query_total_bp": query_total,
                "target_coverage_pct": 100.0
                * int(summary["target_aligned_bp"])
                / target_total,
                "query_coverage_pct": 100.0
                * int(summary["query_aligned_bp"])
                / query_total,
            }
        )
        summaries.append(summary)
        correspondences.append((comparison_name, correspondence))
        write_dotplot(
            plots_dir / f"{comparison_name}.svg",
            f"{target.label} (target) vs {query.label} (query)",
            alignments,
            fai_by_key[target_key],
            fai_by_key[query_key],
            args.plot_min_block,
            args.max_plot_alignments,
        )

    write_comparison_tables(args.output_dir, summaries, correspondences)
    write_manifest(args.output_dir, args.data_dir)
    print(f"Genome comparison complete: {args.output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
