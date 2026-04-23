#!/usr/bin/env python3
"""Generate a concise aggregate report for a pipeline run."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def load_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    with path.open() as fh:
        return json.load(fh)


def load_csv_rows(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh))


def count_files(path: Path, pattern: str) -> int:
    return sum(1 for _ in path.glob(pattern)) if path.is_dir() else 0


def maybe_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def summarize_bin_analysis(rows: list[dict]) -> dict:
    by_binner = Counter()
    by_assembler = Counter()
    samples = set()
    for row in rows:
        by_binner[row.get("binner", "")] += 1
        by_assembler[row.get("assembler", "")] += 1
        samples.add(row.get("sample_base") or row.get("sample") or "")
    return {
        "rows": len(rows),
        "samples": len(samples - {""}),
        "by_binner": dict(sorted((k, v) for k, v in by_binner.items() if k)),
        "by_assembler": dict(sorted((k, v) for k, v in by_assembler.items() if k)),
    }


def summarize_das_qc(rows: list[dict]) -> dict:
    completeness = []
    qc_hits = 0
    for row in rows:
        if row.get("qc_ref_name"):
            qc_hits += 1
        val = maybe_float(row.get("busco_complete_pct", ""))
        if val is not None:
            completeness.append(val)

    avg_complete = round(sum(completeness) / len(completeness), 2) if completeness else None
    return {
        "rows": len(rows),
        "quickclade_hits": qc_hits,
        "busco_rows": len(completeness),
        "avg_busco_complete_pct": avg_complete,
    }


def build_summary(run_root: Path, metadata: dict) -> dict:
    filtered_dir = run_root / "02_filtered_reads"
    assembly_dir = run_root / "03_assemblies"
    binning_dir = run_root / "04_binning"
    gtdb_dir = run_root / "05_gtdbtk"

    bin_analysis_rows = load_csv_rows(binning_dir / "bin_analysis_results.csv")
    das_qc_rows = load_csv_rows(binning_dir / "das_bins_qc_summary.csv")

    gtdb_reports = sorted(gtdb_dir.glob("**/*summary*.tsv")) if gtdb_dir.is_dir() else []

    return {
        "run_root": str(run_root),
        "metadata": metadata,
        "artifacts": {
            "filtered_reads": count_files(filtered_dir, "*_interleaved_filtered.fastq.gz"),
            "spades_fastas": count_files(assembly_dir, "*_spades.fasta"),
            "megahit_sensitive_fastas": count_files(assembly_dir, "*_megahit_sensitive.fasta"),
            "megahit_large_fastas": count_files(assembly_dir, "*_megahit_large.fasta"),
            "binning_sample_dirs": sum(1 for p in binning_dir.iterdir() if p.is_dir()) if binning_dir.is_dir() else 0,
            "gtdb_summary_files": len(gtdb_reports),
        },
        "bin_analysis": summarize_bin_analysis(bin_analysis_rows),
        "das_bins_qc": summarize_das_qc(das_qc_rows),
    }


def render_markdown(summary: dict) -> str:
    md = []
    meta = summary.get("metadata", {})
    artifacts = summary["artifacts"]
    bin_analysis = summary["bin_analysis"]
    das_qc = summary["das_bins_qc"]

    md.append("# Pipeline Run Report")
    md.append("")
    md.append("## Run")
    md.append(f"- Pipeline: {meta.get('pipeline_name', '')}")
    md.append(f"- Run root: `{summary['run_root']}`")
    md.append(f"- Started (UTC): {meta.get('started_at_utc', 'unknown')}")
    md.append(f"- User: {meta.get('user', '')}")
    md.append(f"- Host: {meta.get('hostname', '')}")
    if meta.get("slurm_job_id"):
        md.append(f"- SLURM job ID: {meta['slurm_job_id']}")
    if meta.get("sample_filter"):
        md.append(f"- Sample filter: `{meta['sample_filter']}`")

    md.append("")
    md.append("## Artifacts")
    md.append(f"- Filtered read files: {artifacts['filtered_reads']}")
    md.append(f"- SPAdes assemblies: {artifacts['spades_fastas']}")
    md.append(f"- MEGAHIT sensitive assemblies: {artifacts['megahit_sensitive_fastas']}")
    md.append(f"- MEGAHIT large assemblies: {artifacts['megahit_large_fastas']}")
    md.append(f"- Binning sample directories: {artifacts['binning_sample_dirs']}")
    md.append(f"- GTDB summary files: {artifacts['gtdb_summary_files']}")

    md.append("")
    md.append("## Bin Analysis")
    md.append(f"- Total rows: {bin_analysis['rows']}")
    md.append(f"- Samples represented: {bin_analysis['samples']}")
    if bin_analysis["by_binner"]:
        md.append(f"- Rows by binner: {json.dumps(bin_analysis['by_binner'], sort_keys=True)}")
    if bin_analysis["by_assembler"]:
        md.append(f"- Rows by assembler: {json.dumps(bin_analysis['by_assembler'], sort_keys=True)}")

    md.append("")
    md.append("## DAS Tool QC")
    md.append(f"- Total rows: {das_qc['rows']}")
    md.append(f"- QuickClade hits: {das_qc['quickclade_hits']}")
    md.append(f"- BUSCO rows: {das_qc['busco_rows']}")
    if das_qc["avg_busco_complete_pct"] is not None:
        md.append(f"- Mean BUSCO completeness: {das_qc['avg_busco_complete_pct']:.2f}%")

    return "\n".join(md) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--metadata", required=False, default="")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    run_root = Path(args.run_root).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata = load_json(Path(args.metadata)) if args.metadata else {}
    summary = build_summary(run_root, metadata)

    summary_path = out_dir / "run_summary.json"
    report_path = out_dir / "run_report.md"

    with summary_path.open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    report_path.write_text(render_markdown(summary))

    print(f"Wrote summary JSON -> {summary_path}")
    print(f"Wrote report Markdown -> {report_path}")


if __name__ == "__main__":
    main()
