#!/usr/bin/env python3
"""Analyze Project 33 sex-gap differential effect between references.

Supports both single-model and multi-model (Holm-Bonferroni corrected) analysis
per Gold Standard Protocol §2, §5.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def summarize_by_repeat(metrics_df: pd.DataFrame) -> pd.DataFrame:
    required = {"repeat", "group", "balanced_accuracy"}
    missing = required - set(metrics_df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    pivot = metrics_df.pivot_table(
        index="repeat",
        columns="group",
        values="balanced_accuracy",
        aggfunc="mean",
    )

    if "sex_M" not in pivot.columns or "sex_F" not in pivot.columns:
        raise ValueError("Metrics must include groups 'sex_M' and 'sex_F'.")

    out = pivot.reset_index().copy()
    out["sex_gap_male_minus_female"] = out["sex_M"] - out["sex_F"]
    return out[["repeat", "sex_M", "sex_F", "sex_gap_male_minus_female"]]


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int):
    rng = np.random.default_rng(seed)
    n = len(values)
    boots = []
    for _ in range(n_boot):
        sample = rng.choice(values, size=n, replace=True)
        boots.append(float(np.mean(sample)))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def paired_sign_permutation_test(values: np.ndarray, n_perm: int, seed: int):
    rng = np.random.default_rng(seed)
    observed = float(np.mean(values))
    abs_obs = abs(observed)

    null_stats = []
    for _ in range(n_perm):
        signs = rng.choice([-1.0, 1.0], size=len(values), replace=True)
        null_stats.append(float(np.mean(values * signs)))

    null_stats = np.array(null_stats)
    p_value = float((np.sum(np.abs(null_stats) >= abs_obs) + 1) / (len(null_stats) + 1))
    return observed, p_value


def holm_bonferroni(p_values: list[float]) -> list[float]:
    """Apply Holm-Bonferroni step-down correction to a list of p-values."""
    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n
    cummax = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        adj = p * (n - rank)
        adj = min(adj, 1.0)
        cummax = max(cummax, adj)
        adjusted[orig_idx] = cummax
    return adjusted


def analyze_single_model(grch38_path, t2t_path, bootstrap_iters, perm_iters, seed):
    """Analyze a single model pair and return summary dict + merged DataFrame."""
    grch38 = pd.read_csv(grch38_path)
    t2t = pd.read_csv(t2t_path)

    g = summarize_by_repeat(grch38).rename(
        columns={
            "sex_M": "grch38_sex_M",
            "sex_F": "grch38_sex_F",
            "sex_gap_male_minus_female": "grch38_gap",
        }
    )
    t = summarize_by_repeat(t2t).rename(
        columns={
            "sex_M": "t2t_sex_M",
            "sex_F": "t2t_sex_F",
            "sex_gap_male_minus_female": "t2t_gap",
        }
    )

    merged = pd.merge(g, t, on="repeat", how="inner")
    if merged.empty:
        raise ValueError("No overlapping repeats found between GRCh38 and T2T metrics.")

    merged["delta_delta"] = merged["grch38_gap"] - merged["t2t_gap"]

    values = merged["delta_delta"].to_numpy(dtype=float)
    mean_dd = float(np.mean(values))
    ci_lo, ci_hi = bootstrap_ci(values, bootstrap_iters, seed)
    observed, p_value = paired_sign_permutation_test(values, perm_iters, seed + 1)

    summary = {
        "n_repeats": int(len(merged)),
        "delta_delta_mean": mean_dd,
        "delta_delta_bootstrap_ci_95": [ci_lo, ci_hi],
        "paired_sign_permutation_p_value": p_value,
        "direction_supports_hypothesis": bool(mean_dd > 0),
        "observed_statistic": observed,
        "effect_size_threshold_met": bool(abs(mean_dd) >= 0.01),
    }
    return summary, merged


def main():
    parser = argparse.ArgumentParser(
        description="Analyze Project 33 sex-gap differential effect between references."
    )
    parser.add_argument("--grch38-metrics", required=True,
                        help="Single CSV or comma-separated list of model CSVs")
    parser.add_argument("--t2t-metrics", required=True,
                        help="Single CSV or comma-separated list of model CSVs (same order)")
    parser.add_argument("--model-names",
                        help="Comma-separated model names (required if multiple metrics)")
    parser.add_argument("--bootstrap-iters", type=int, default=5000)
    parser.add_argument("--permutation-iters", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260222)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    grch38_paths = [p.strip() for p in args.grch38_metrics.split(",")]
    t2t_paths = [p.strip() for p in args.t2t_metrics.split(",")]

    if len(grch38_paths) != len(t2t_paths):
        raise ValueError("Must provide same number of GRCh38 and T2T metric files.")

    if len(grch38_paths) > 1 and not args.model_names:
        raise ValueError("--model-names required when analyzing multiple models.")

    model_names = (
        [n.strip() for n in args.model_names.split(",")]
        if args.model_names
        else ["model"]
    )
    if len(model_names) != len(grch38_paths):
        raise ValueError("Number of model names must match number of metric file pairs.")

    # Analyze each model
    all_summaries = {}
    all_merged = {}
    raw_p_values = []

    for i, (gpath, tpath, mname) in enumerate(zip(grch38_paths, t2t_paths, model_names)):
        summary, merged = analyze_single_model(
            gpath, tpath, args.bootstrap_iters, args.permutation_iters,
            args.seed + i * 100,
        )
        summary["inputs"] = {"grch38_metrics": gpath, "t2t_metrics": tpath}
        all_summaries[mname] = summary
        all_merged[mname] = merged
        raw_p_values.append(summary["paired_sign_permutation_p_value"])

    # Apply Holm-Bonferroni correction across model families (Protocol §5.5)
    adjusted_p = holm_bonferroni(raw_p_values)
    for i, mname in enumerate(model_names):
        all_summaries[mname]["holm_adjusted_p_value"] = adjusted_p[i]
        all_summaries[mname]["holm_significant_at_005"] = bool(adjusted_p[i] < 0.05)

    # Determine confirmatory decision (Protocol §2)
    dd_means = [all_summaries[m]["delta_delta_mean"] for m in model_names]
    median_dd = float(np.median(dd_means))
    ci_excludes_zero = []
    for m in model_names:
        lo, hi = all_summaries[m]["delta_delta_bootstrap_ci_95"]
        ci_excludes_zero.append(lo > 0 or hi < 0)

    any_holm_sig = any(all_summaries[m]["holm_significant_at_005"] for m in model_names)
    any_effect_threshold = any(all_summaries[m]["effect_size_threshold_met"] for m in model_names)

    confirmatory_criteria = {
        "criterion_1_directional_consistency": bool(median_dd > 0),
        "criterion_2_any_ci_excludes_zero": any(ci_excludes_zero),
        "criterion_3_holm_significant": any_holm_sig,
        "criterion_4_effect_size_threshold": any_effect_threshold,
        "median_delta_delta_across_models": median_dd,
    }

    all_pass = all(confirmatory_criteria[f"criterion_{i}_{k}"]
                   for i, k in [(1, "directional_consistency"),
                                (2, "any_ci_excludes_zero"),
                                (3, "holm_significant"),
                                (4, "effect_size_threshold")])
    suggestive = (confirmatory_criteria["criterion_1_directional_consistency"]
                  and confirmatory_criteria["criterion_2_any_ci_excludes_zero"])

    if all_pass:
        evidence_tier = "A_confirmatory"
    elif suggestive:
        evidence_tier = "B_suggestive"
    else:
        evidence_tier = "C_exploratory"

    confirmatory_criteria["evidence_tier"] = evidence_tier

    # Save outputs
    out_csv = Path(args.out_csv)
    out_json = Path(args.out_json)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    # Combine all merged DataFrames with model column
    combined_dfs = []
    for mname, mdf in all_merged.items():
        mdf = mdf.copy()
        mdf["model"] = mname
        combined_dfs.append(mdf)
    combined = pd.concat(combined_dfs, ignore_index=True)
    combined.to_csv(out_csv, index=False)

    output = {
        "generated_utc": now_utc_iso(),
        "estimand": "delta_delta = (male-female gap in GRCh38) - (male-female gap in T2T)",
        "n_models": len(model_names),
        "model_results": all_summaries,
        "confirmatory_decision": confirmatory_criteria,
    }
    out_json.write_text(json.dumps(output, indent=2) + "\n")

    print(json.dumps({"status": "ok", "out_csv": str(out_csv), "out_json": str(out_json),
                       "evidence_tier": evidence_tier}, indent=2))


if __name__ == "__main__":
    main()
