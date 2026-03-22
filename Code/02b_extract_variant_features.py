#!/usr/bin/env python3
import argparse
import csv
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_gt_to_dosage(gt: str):
    if gt in {".", "./.", ".|."}:
        return np.nan
    sep = "/" if "/" in gt else "|" if "|" in gt else None
    if sep is None:
        return np.nan
    alleles = gt.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return np.nan
    try:
        return float(int(alleles[0]) + int(alleles[1]))
    except ValueError:
        return np.nan


def load_metadata(metadata_csv: Path):
    meta = pd.read_csv(metadata_csv)
    required = {"patient_id", "sex", "superpopulation"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing columns: {sorted(missing)}")
    return meta


def get_samples_from_vcf(vcf_path: Path):
    cmd = ["bcftools", "query", "-l", str(vcf_path)]
    out = subprocess.check_output(cmd, text=True)
    samples = [x.strip() for x in out.splitlines() if x.strip()]
    if not samples:
        raise ValueError(f"No samples found in VCF: {vcf_path}")
    return samples


def stream_variants(vcf_path: Path, max_variants: int, maf_threshold: float, call_rate_threshold: float):
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%ID[\t%GT]\n", str(vcf_path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None

    accepted = []
    seen_ids = set()
    stopped_early = False

    try:
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, pos, vid = parts[0], parts[1], parts[2]
            gt_values = parts[3:]

            dosages = np.array([parse_gt_to_dosage(gt) for gt in gt_values], dtype=float)
            non_missing = ~np.isnan(dosages)
            call_rate = float(np.mean(non_missing))
            if call_rate < call_rate_threshold:
                continue

            called = dosages[non_missing]
            if called.size == 0:
                continue

            af = float(np.mean(called) / 2.0)
            maf = min(af, 1.0 - af)
            if maf < maf_threshold:
                continue

            if vid and vid != ".":
                variant_id = vid
            else:
                variant_id = f"{chrom}:{pos}"

            if variant_id in seen_ids:
                continue
            seen_ids.add(variant_id)

            accepted.append((variant_id, dosages))
            if len(accepted) >= max_variants:
                stopped_early = True
                break
    finally:
        # Always close stdout so bcftools gets SIGPIPE and exits cleanly.
        # Without this, proc.wait() deadlocks when we break early because
        # the pipe buffer fills and bcftools blocks on write.
        proc.stdout.close()

    if stopped_early:
        proc.kill()

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()

    return_code = proc.returncode
    # A non-zero return code is only an error when we didn't stop early
    # (kill/SIGPIPE produce -9 / -13 which are expected for early exit).
    if not stopped_early and return_code != 0:
        raise RuntimeError(f"bcftools query failed (code {return_code}): {stderr}")

    return accepted


def build_matrix(samples, metadata_df, variants, reference_label):
    # Build via concat instead of repeated column assignment to avoid fragmentation
    base = pd.DataFrame({"patient_id": samples})
    feature_cols = {
        f"{reference_label}__{variant_id}": dosages
        for variant_id, dosages in variants
    }
    feature_df = pd.concat([base, pd.DataFrame(feature_cols)], axis=1)

    merged = feature_df.merge(metadata_df[["patient_id", "sex", "superpopulation"]], on="patient_id", how="inner")
    base_cols = ["patient_id", "sex", "superpopulation"]
    variant_cols = [c for c in merged.columns if c not in base_cols]
    merged = merged[base_cols + variant_cols]
    return merged


def main():
    parser = argparse.ArgumentParser(description="Extract first-pass variant feature matrix from a VCF using bcftools query.")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--reference-label", required=True, choices=["GRCh38", "T2T_CHM13"])
    parser.add_argument("--max-variants", type=int, default=2000)
    parser.add_argument("--maf-threshold", type=float, default=0.01)
    parser.add_argument("--call-rate-threshold", type=float, default=0.95)
    parser.add_argument("--out-features-csv", required=True)
    parser.add_argument("--out-summary-json", required=True)
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    metadata_path = Path(args.metadata)

    metadata_df = load_metadata(metadata_path)
    samples = get_samples_from_vcf(vcf_path)
    variants = stream_variants(vcf_path, args.max_variants, args.maf_threshold, args.call_rate_threshold)

    matrix = build_matrix(samples, metadata_df, variants, args.reference_label)

    out_csv = Path(args.out_features_csv)
    out_json = Path(args.out_summary_json)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    matrix.to_csv(out_csv, index=False)

    summary = {
        "generated_utc": now_utc_iso(),
        "reference_label": args.reference_label,
        "vcf": str(vcf_path),
        "metadata": str(metadata_path),
        "n_samples_in_vcf": len(samples),
        "n_samples_after_metadata_merge": int(matrix.shape[0]),
        "n_variant_features": int(matrix.shape[1] - 3),
        "max_variants_requested": args.max_variants,
        "maf_threshold": args.maf_threshold,
        "call_rate_threshold": args.call_rate_threshold,
        "out_features_csv": str(out_csv),
    }
    out_json.write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps({"status": "ok", **summary}, indent=2))


if __name__ == "__main__":
    main()
