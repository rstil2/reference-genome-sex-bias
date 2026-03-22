#!/usr/bin/env python3
"""
02f_merge_chromosome_features.py
---------------------------------
Merge per-chromosome feature CSVs (output of 02b_extract_variant_features.py)
into a single genome-wide feature matrix for one reference.

Each per-chromosome CSV has columns:
  patient_id, sex, superpopulation, <ref>__<variant_id>, ...

This script:
  1. Loads all per-chromosome CSVs
  2. Verifies consistent metadata (patient_id, sex, superpopulation)
  3. Merges variant feature columns on patient_id
  4. Writes the merged matrix and a JSON summary
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


BASE_COLS = ["patient_id", "sex", "superpopulation"]


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-chromosome feature CSVs into a genome-wide matrix."
    )
    parser.add_argument(
        "--input-csvs",
        nargs="+",
        required=True,
        help="Ordered list of per-chromosome feature CSVs (e.g. chr1.csv chr2.csv ...)",
    )
    parser.add_argument("--out-features-csv", required=True, help="Output merged feature matrix")
    parser.add_argument("--out-summary-json", required=True, help="Output summary JSON")
    args = parser.parse_args()

    input_paths = [Path(p) for p in args.input_csvs]
    missing = [str(p) for p in input_paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Input CSVs not found: {missing}")

    # Load first file to establish the base patient roster + metadata
    print(f"Loading {len(input_paths)} chromosome feature files...")
    base_df = pd.read_csv(input_paths[0])
    for col in BASE_COLS:
        if col not in base_df.columns:
            raise ValueError(f"Missing base column '{col}' in {input_paths[0]}")

    merged = base_df[BASE_COLS].copy()
    feature_counts = []
    total_features = 0

    for path in input_paths:
        df = pd.read_csv(path)
        for col in BASE_COLS:
            if col not in df.columns:
                raise ValueError(f"Missing base column '{col}' in {path}")

        # Verify patient roster is identical
        if set(df["patient_id"]) != set(merged["patient_id"]):
            shared = set(df["patient_id"]) & set(merged["patient_id"])
            print(
                f"  WARNING: {path.name} has {len(df['patient_id'])} patients; "
                f"{len(shared)} shared with base roster. Using inner join."
            )

        variant_cols = [c for c in df.columns if c not in BASE_COLS]
        if not variant_cols:
            print(f"  WARNING: {path.name} has no variant columns — skipping.")
            continue

        feature_counts.append({"file": path.name, "n_features": len(variant_cols)})
        total_features += len(variant_cols)
        merged = merged.merge(df[["patient_id"] + variant_cols], on="patient_id", how="inner")
        print(f"  {path.name}: +{len(variant_cols)} features (running total: {total_features})")

    out_csv = Path(args.out_features_csv)
    out_json = Path(args.out_summary_json)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(out_csv, index=False)

    summary = {
        "generated_utc": now_utc_iso(),
        "n_input_files": len(input_paths),
        "n_samples": int(len(merged)),
        "n_total_variant_features": int(merged.shape[1] - len(BASE_COLS)),
        "per_chromosome_feature_counts": feature_counts,
        "out_features_csv": str(out_csv),
    }
    out_json.write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps({"status": "ok", **summary}, indent=2))


if __name__ == "__main__":
    main()
