#!/usr/bin/env python3
"""
02c_align_reference_features.py
---------------------------------
Align the GRCh38 and T2T-CHM13 feature matrices to an identical patient roster.

The two feature matrices come from the same samples but mapped to different
reference genomes, so their variant coordinates are NOT comparable across
references. Each reference uses its own best informative variants.
This script only ensures both matrices contain exactly the same patients in
the same order — a prerequisite for the paired sex-stratified comparison.

Variant columns are prefixed with the reference label and remain independent.
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
        description="Align GRCh38 and T2T feature matrices to a common patient roster."
    )
    parser.add_argument("--grch38-features", required=True)
    parser.add_argument("--t2t-features", required=True)
    parser.add_argument("--min-common-patients", type=int, default=100)
    parser.add_argument("--min-features-per-ref", type=int, default=20)
    parser.add_argument("--out-grch38-aligned", required=True)
    parser.add_argument("--out-t2t-aligned", required=True)
    parser.add_argument("--out-summary-json", required=True)
    # Legacy arg kept for backwards compatibility — no longer used
    parser.add_argument("--min-common-features", type=int, default=20)
    args = parser.parse_args()

    g = pd.read_csv(args.grch38_features)
    t = pd.read_csv(args.t2t_features)

    for col in BASE_COLS:
        if col not in g.columns:
            raise ValueError(f"Missing required column '{col}' in GRCh38 features")
        if col not in t.columns:
            raise ValueError(f"Missing required column '{col}' in T2T features")

    g_variants = [c for c in g.columns if c not in BASE_COLS]
    t_variants = [c for c in t.columns if c not in BASE_COLS]

    if len(g_variants) < args.min_features_per_ref:
        raise ValueError(f"GRCh38 feature count too low: {len(g_variants)} < {args.min_features_per_ref}")
    if len(t_variants) < args.min_features_per_ref:
        raise ValueError(f"T2T feature count too low: {len(t_variants)} < {args.min_features_per_ref}")

    # Align to common patient set
    common_patients = sorted(set(g["patient_id"]) & set(t["patient_id"]))
    if len(common_patients) < args.min_common_patients:
        raise ValueError(
            f"Too few shared patients: {len(common_patients)} < {args.min_common_patients}"
        )

    g2 = g[g["patient_id"].isin(common_patients)].sort_values("patient_id").reset_index(drop=True)
    t2 = t[t["patient_id"].isin(common_patients)].sort_values("patient_id").reset_index(drop=True)

    out_g = Path(args.out_grch38_aligned)
    out_t = Path(args.out_t2t_aligned)
    out_s = Path(args.out_summary_json)
    out_g.parent.mkdir(parents=True, exist_ok=True)
    out_t.parent.mkdir(parents=True, exist_ok=True)
    out_s.parent.mkdir(parents=True, exist_ok=True)

    g2.to_csv(out_g, index=False)
    t2.to_csv(out_t, index=False)

    summary = {
        "generated_utc": now_utc_iso(),
        "n_common_patients": int(len(common_patients)),
        "n_grch38_variant_features": int(len(g_variants)),
        "n_t2t_variant_features": int(len(t_variants)),
        "note": "Features are reference-specific and not cross-aligned by position. Patient roster is aligned.",
        "grch38_input": args.grch38_features,
        "t2t_input": args.t2t_features,
        "grch38_aligned_output": str(out_g),
        "t2t_aligned_output": str(out_t),
    }
    out_s.write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps({"status": "ok", **summary}, indent=2))


if __name__ == "__main__":
    main()

