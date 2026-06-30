#!/usr/bin/env python3
"""
Redesigned downstream panel battery (v2).

Runs reference-sensitive and reference-insensitive downstream tasks,
computes ΔΔ per panel, and writes a combined summary.

See DOWNSTREAM_REDESIGN_v2.md for panel definitions.
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent


PANELS = {
    "panel_A_chrx_summary_superpop": {
        "description": "Superpopulation classification from chrX per-sample summary statistics",
        "reference_sensitive": True,
        "label_column": "superpopulation",
        "stratify_column": "sex",
        "feature_builder": "chrx_summary",
    },
    "panel_B_chrx_summary_sex": {
        "description": "Sex classification from chrX per-sample summary statistics",
        "reference_sensitive": True,
        "label_column": "sex",
        "stratify_column": "sex",
        "feature_builder": "chrx_summary",
    },
    "panel_C_genome_wide_snp": {
        "description": "Superpopulation from genome-wide common SNPs (v1 design, negative control)",
        "reference_sensitive": False,
        "label_column": "superpopulation",
        "stratify_column": "sex",
        "feature_builder": "genome_wide",
    },
    "panel_D_autosome_only": {
        "description": "Superpopulation from autosomal common SNPs only (negative control)",
        "reference_sensitive": False,
        "label_column": "superpopulation",
        "stratify_column": "sex",
        "feature_builder": "autosome_only",
    },
}

MODELS = ["logistic_regression", "random_forest", "xgboost", "svm", "mlp"]
BASE_COLS = ["patient_id", "sex", "superpopulation"]
CHRX_SUMMARY_COLS = ["het_rate", "call_rate", "hom_rate", "n_het", "n_called", "n_variants"]


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_chrx_summary_features(per_sample_chrX: Path, reference: str) -> pd.DataFrame:
    df = pd.read_csv(per_sample_chrX)
    sub = df[df["reference"] == reference].copy()
    cols = BASE_COLS + [c for c in CHRX_SUMMARY_COLS if c in sub.columns]
    out = sub[cols].copy()
    out = out.rename(columns={c: f"chrx_{c}" for c in CHRX_SUMMARY_COLS if c in out.columns})
    return out.sort_values("patient_id").reset_index(drop=True)


def build_genome_wide_features(aligned_csv: Path) -> pd.DataFrame:
    return pd.read_csv(aligned_csv)


def build_autosome_features(aligned_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(aligned_csv)
    feature_cols = [
        c for c in df.columns
        if c not in BASE_COLS and "chrX" not in c and "__chrX:" not in c
    ]
    return df[BASE_COLS + feature_cols]


def ensure_features(panel_id: str, builder: str, paths: dict, features_dir: Path) -> tuple[Path, Path]:
    features_dir.mkdir(parents=True, exist_ok=True)
    g_out = features_dir / f"{panel_id}_grch38_features.csv"
    t_out = features_dir / f"{panel_id}_t2t_features.csv"

    if builder == "chrx_summary":
        build_chrx_summary_features(paths["per_sample_chrX"], "GRCh38").to_csv(g_out, index=False)
        build_chrx_summary_features(paths["per_sample_chrX"], "T2T").to_csv(t_out, index=False)
    elif builder == "genome_wide":
        build_genome_wide_features(paths["grch38_genome"]).to_csv(g_out, index=False)
        build_genome_wide_features(paths["t2t_genome"]).to_csv(t_out, index=False)
    elif builder == "autosome_only":
        build_autosome_features(paths["grch38_genome"]).to_csv(g_out, index=False)
        build_autosome_features(paths["t2t_genome"]).to_csv(t_out, index=False)
    else:
        raise ValueError(f"Unknown feature builder: {builder}")

    return g_out, t_out


def run_training(python: str, panel_id: str, panel_cfg: dict, g_feat: Path, t_feat: Path,
                 split_manifest: Path, out_dir: Path, seed: int, max_repeats: int, models: list[str]):
    panel_dir = out_dir / panel_id
    panel_dir.mkdir(parents=True, exist_ok=True)

    grch38_metrics = []
    t2t_metrics = []

    train_script = CODE_DIR / "03_train_models.py"
    for model in models:
        for ref, feat, metric_list in [
            ("grch38", g_feat, grch38_metrics),
            ("t2t", t_feat, t2t_metrics),
        ]:
            metrics_out = panel_dir / f"{ref}_{model}_metrics.csv"
            preds_out = panel_dir / f"{ref}_{model}_predictions.csv"
            summary_out = panel_dir / f"{ref}_{model}_summary.json"
            cmd = [
                python, str(train_script),
                "--features-csv", str(feat),
                "--split-manifest", str(split_manifest),
                "--model", model,
                "--label-column", panel_cfg["label_column"],
                "--stratify-column", panel_cfg["stratify_column"],
                "--seed", str(seed),
                "--max-repeats", str(max_repeats),
                "--metrics-out", str(metrics_out),
                "--predictions-out", str(preds_out),
                "--summary-out", str(summary_out),
            ]
            subprocess.run(cmd, check=True)
            metric_list.append(str(metrics_out))

    return grch38_metrics, t2t_metrics, panel_dir


def analyze_panel(panel_id: str, panel_cfg: dict, grch38_metrics: list[str],
                  t2t_metrics: list[str], panel_dir: Path, seed: int) -> dict:
    grch38_arg = ",".join(grch38_metrics)
    t2t_arg = ",".join(t2t_metrics)
    models_arg = ",".join(MODELS)

    out_json = panel_dir / "delta_delta_summary.json"
    out_csv = panel_dir / "delta_delta_by_repeat.csv"

    cmd = [
        sys.executable, str(CODE_DIR / "04_analyze_results.py"),
        "--grch38-metrics", grch38_arg,
        "--t2t-metrics", t2t_arg,
        "--model-names", models_arg,
        "--seed", str(seed),
        "--out-json", str(out_json),
        "--out-csv", str(out_csv),
    ]
    subprocess.run(cmd, check=True)
    result = json.loads(out_json.read_text())
    result["panel_id"] = panel_id
    result["description"] = panel_cfg["description"]
    result["reference_sensitive"] = panel_cfg["reference_sensitive"]
    return result


def main():
    parser = argparse.ArgumentParser(description="Run redesigned downstream panel battery.")
    parser.add_argument("--study-root", default=str(CODE_DIR.parent))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--seed", type=int, default=20260222)
    parser.add_argument("--max-repeats", type=int, default=20)
    parser.add_argument("--models", nargs="+", default=MODELS)
    parser.add_argument("--panels", nargs="+", default=list(PANELS.keys()))
    parser.add_argument("--skip-training", action="store_true")
    args = parser.parse_args()

    root = Path(args.study_root)
    results = root / "Results"
    paths = {
        "per_sample_chrX": results / "variant_analysis" / "per_sample_chrX.csv",
        "grch38_genome": results / "features" / "grch38_genome_features_aligned.csv",
        "t2t_genome": results / "features" / "t2t_genome_features_aligned.csv",
        "split_manifest": results / "splits" / "split_manifest_v1.csv",
        "grch38_chrX": results / "features" / "grch38_chrX_features.csv",
        "t2t_chrX": results / "features" / "t2t_chrX_features.csv",
    }

    out_dir = results / "metrics" / "downstream_v2"
    features_dir = results / "features" / "downstream_v2"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Panel E: site overlap
    overlap_json = out_dir / "site_overlap_chrX.json"
    subprocess.run([
        args.python, str(CODE_DIR / "02g_compare_reference_variant_sets.py"),
        "--grch38-features", str(paths["grch38_chrX"]),
        "--t2t-features", str(paths["t2t_chrX"]),
        "--out-json", str(overlap_json),
    ], check=True)

    panel_results = {}
    for panel_id in args.panels:
        if panel_id not in PANELS:
            raise ValueError(f"Unknown panel: {panel_id}")
        cfg = PANELS[panel_id]
        g_feat, t_feat = ensure_features(panel_id, cfg["feature_builder"], paths, features_dir)

        panel_dir = out_dir / panel_id
        if args.skip_training:
            grch38_metrics = [str(panel_dir / f"grch38_{m}_metrics.csv") for m in args.models]
            t2t_metrics = [str(panel_dir / f"t2t_{m}_metrics.csv") for m in args.models]
        else:
            grch38_metrics, t2t_metrics, panel_dir = run_training(
                args.python, panel_id, cfg, g_feat, t_feat,
                paths["split_manifest"], out_dir, args.seed, args.max_repeats, args.models,
            )

        panel_results[panel_id] = analyze_panel(
            panel_id, cfg, grch38_metrics, t2t_metrics, panel_dir, args.seed,
        )

    # Cross-panel summary
    summary_rows = []
    for pid, res in panel_results.items():
        median_dd = res["confirmatory_decision"]["median_delta_delta_across_models"]
        any_sig = res["confirmatory_decision"]["criterion_3_holm_significant"]
        for model, mres in res["model_results"].items():
            summary_rows.append({
                "panel": pid,
                "reference_sensitive": PANELS[pid]["reference_sensitive"],
                "model": model,
                "delta_delta_mean": mres["delta_delta_mean"],
                "ci_lo": mres["delta_delta_bootstrap_ci_95"][0],
                "ci_hi": mres["delta_delta_bootstrap_ci_95"][1],
                "holm_p": mres["holm_adjusted_p_value"],
                "holm_significant": mres["holm_significant_at_005"],
            })
        summary_rows.append({
            "panel": pid,
            "reference_sensitive": PANELS[pid]["reference_sensitive"],
            "model": "_panel_median",
            "delta_delta_mean": median_dd,
            "ci_lo": None,
            "ci_hi": None,
            "holm_p": None,
            "holm_significant": any_sig,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = out_dir / "panel_battery_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    battery_summary = {
        "generated_utc": now_utc_iso(),
        "panels_run": args.panels,
        "site_overlap": json.loads(overlap_json.read_text()),
        "panel_results": panel_results,
        "summary_csv": str(summary_csv),
    }
    battery_json = out_dir / "panel_battery_summary.json"
    battery_json.write_text(json.dumps(battery_summary, indent=2) + "\n")

    print(json.dumps({
        "status": "ok",
        "out_dir": str(out_dir),
        "battery_summary": str(battery_json),
    }, indent=2))


if __name__ == "__main__":
    main()
