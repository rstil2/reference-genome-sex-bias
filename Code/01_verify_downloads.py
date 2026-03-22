#!/usr/bin/env python3
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

EXPECTED_CHR = [str(i) for i in range(1, 23)] + ["X"]


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def expected_files_t2t():
    vcfs = [f"1KGP.CHM13v2.0.chr{chrom}.recalibrated.snp_indel.pass.vcf.gz" for chrom in EXPECTED_CHR]
    tbis = [f"{name}.tbi" for name in vcfs]
    return vcfs, tbis


def expected_files_grch38():
    vcfs = [f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz" for chrom in EXPECTED_CHR]
    tbis = [f"{name}.tbi" for name in vcfs]
    return vcfs, tbis


def evaluate_directory(directory: Path, expected_vcfs, expected_tbis):
    present = {p.name: p for p in directory.glob("*") if p.is_file()}

    missing_vcfs = [f for f in expected_vcfs if f not in present]
    missing_tbis = [f for f in expected_tbis if f not in present]

    partial_files = sorted([p.name for p in directory.glob("*.part")])
    empty_files = sorted([name for name, p in present.items() if p.stat().st_size == 0])

    completed_vcfs = [f for f in expected_vcfs if f in present]
    completed_tbis = [f for f in expected_tbis if f in present]

    return {
        "directory": str(directory),
        "completed_vcfs": len(completed_vcfs),
        "completed_tbis": len(completed_tbis),
        "expected_vcfs": len(expected_vcfs),
        "expected_tbis": len(expected_tbis),
        "missing_vcfs": missing_vcfs,
        "missing_tbis": missing_tbis,
        "partial_files": partial_files,
        "empty_files": empty_files,
        "complete": len(missing_vcfs) == 0 and len(missing_tbis) == 0 and len(partial_files) == 0 and len(empty_files) == 0,
    }


def main():
    parser = argparse.ArgumentParser(description="Verify completion state of Project 33 VCF downloads.")
    parser.add_argument("--t2t-dir", required=True)
    parser.add_argument("--grch38-dir", required=True)
    parser.add_argument("--out-json", required=True)
    args = parser.parse_args()

    t2t_dir = Path(args.t2t_dir)
    grch38_dir = Path(args.grch38_dir)

    t2t_vcfs, t2t_tbis = expected_files_t2t()
    g_vcfs, g_tbis = expected_files_grch38()

    t2t_result = evaluate_directory(t2t_dir, t2t_vcfs, t2t_tbis)
    grch38_result = evaluate_directory(grch38_dir, g_vcfs, g_tbis)

    report = {
        "generated_utc": now_utc_iso(),
        "t2t": t2t_result,
        "grch38": grch38_result,
        "all_complete": t2t_result["complete"] and grch38_result["complete"],
    }

    out_path = Path(args.out_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n")

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
