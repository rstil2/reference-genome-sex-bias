#!/usr/bin/env python3
"""Generate revision figures for downstream v2 panel battery."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

MODEL_LABELS = {
    "logistic_regression": "Logistic Reg.",
    "random_forest": "Random Forest",
    "xgboost": "XGBoost",
    "svm": "SVM",
    "mlp": "MLP",
}
PANEL_LABELS = {
    "panel_A_chrx_summary_superpop": "Panel A\n(chrX summaries)",
    "panel_B_chrx_summary_sex": "Panel B\n(sex clf)",
    "panel_C_genome_wide_snp": "Panel C\n(genome SNPs)",
    "panel_D_autosome_only": "Panel D\n(autosomes)",
}
PANEL_COLORS = {
    "panel_A_chrx_summary_superpop": "#E64B35",
    "panel_B_chrx_summary_sex": "#F39B7F",
    "panel_C_genome_wide_snp": "#4DBBD5",
    "panel_D_autosome_only": "#91D1DB",
}
TIER_COLORS = {"A_confirmatory": "#00A087", "B_suggestive": "#FFC107", "C_exploratory": "#8491B4"}


def load_panel_dd(study_root: Path) -> pd.DataFrame:
    rows = []
    base = study_root / "Results" / "metrics" / "downstream_v2"
    for panel_dir in sorted(base.glob("panel_*")):
        summary_path = panel_dir / "delta_delta_summary.json"
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text())
        for model, res in summary["model_results"].items():
            lo, hi = res["delta_delta_bootstrap_ci_95"]
            rows.append({
                "panel": panel_dir.name,
                "model": model,
                "delta_delta": res["delta_delta_mean"],
                "ci_lo": lo,
                "ci_hi": hi,
                "holm_p": res["holm_adjusted_p_value"],
                "holm_sig": res["holm_significant_at_005"],
                "tier": summary["confirmatory_decision"]["evidence_tier"],
            })
    return pd.DataFrame(rows)


def load_panel_metrics(study_root: Path, panel: str, ref: str) -> pd.DataFrame:
    base = study_root / "Results" / "metrics" / "downstream_v2" / panel
    records = []
    for csv in sorted(base.glob(f"{ref}_*_metrics.csv")):
        model = csv.stem.replace(f"{ref}_", "").replace("_metrics", "")
        df = pd.read_csv(csv)
        records.append({"model": model, "metrics": df, "reference": ref})
    return records


def figure5_panel_battery_dashboard(study_root: Path, out_dir: Path):
    """Fixed Figure 5: caption matches panel layout."""
    dd = load_panel_dd(study_root)
    models = [m for m in MODEL_LABELS if m in dd["model"].unique()]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    panel_a_dir = study_root / "Results/metrics/downstream_v2/panel_A_chrx_summary_superpop"

    # A: ΔΔ for Panel A vs C with CIs
    ax = axes[0, 0]
    sub = dd[dd["panel"].isin(["panel_A_chrx_summary_superpop", "panel_C_genome_wide_snp"])]
    y = np.arange(len(models))
    width = 0.35
    for i, panel in enumerate(["panel_A_chrx_summary_superpop", "panel_C_genome_wide_snp"]):
        pdf = sub[sub["panel"] == panel].set_index("model").reindex(models)
        offset = (i - 0.5) * width
        ax.barh(y + offset, pdf["delta_delta"], height=width * 0.9,
                xerr=[pdf["delta_delta"] - pdf["ci_lo"], pdf["ci_hi"] - pdf["delta_delta"]],
                color=PANEL_COLORS[panel], capsize=3, label=PANEL_LABELS[panel], alpha=0.85)
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([MODEL_LABELS[m] for m in models])
    ax.set_xlabel(r"$\Delta\Delta$ (95% bootstrap CI)")
    ax.set_title("A) Reference-sensitive vs insensitive $\Delta\Delta$")
    ax.legend(loc="lower right", fontsize=8)
    ax.invert_yaxis()

    # B: Sex gap by reference for Panel A
    ax = axes[0, 1]
    dd_a = pd.read_csv(panel_a_dir / "delta_delta_by_repeat.csv")
    gap_rows = []
    for model in models:
        mdf = dd_a[dd_a["model"] == model]
        for ref_col, ref_label in [("grch38_gap", "GRCh38"), ("t2t_gap", "T2T-CHM13")]:
            gap_rows.append({
                "model": MODEL_LABELS[model],
                "reference": ref_label,
                "gap": mdf[ref_col].mean(),
            })
    gdf = pd.DataFrame(gap_rows)
    pivot = gdf.pivot(index="model", columns="reference", values="gap")
    pivot = pivot.reindex([MODEL_LABELS[m] for m in models])
    pivot.plot(kind="barh", ax=ax, color=["#E64B35", "#4DBBD5"], width=0.7)
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel(r"Sex gap $\Delta$ (BA$_M$ − BA$_F$)")
    ax.set_title("B) Panel A sex disparity by reference")
    ax.legend(title="Reference", fontsize=8)
    ax.invert_yaxis()

    # C: Overall BA Panel C by reference
    ax = axes[1, 0]
    ba_records = []
    for ref_tag, ref_label in [("grch38", "GRCh38"), ("t2t", "T2T-CHM13")]:
        for csv in (study_root / "Results/metrics/downstream_v2/panel_C_genome_wide_snp").glob(f"{ref_tag}_*_metrics.csv"):
            model = csv.stem.replace(f"{ref_tag}_", "").replace("_metrics", "")
            df = pd.read_csv(csv)
            overall = df[df["group"] == "overall"]["balanced_accuracy"]
            ba_records.append({"model": MODEL_LABELS.get(model, model), "reference": ref_label,
                               "ba": overall.mean()})
    badf = pd.DataFrame(ba_records)
    if not badf.empty:
        pivot = badf.pivot(index="model", columns="reference", values="ba")
        pivot.plot(kind="barh", ax=ax, color=["#E64B35", "#4DBBD5"], width=0.7)
        ax.set_xlim(0.9, 1.0)
        ax.set_xlabel("Overall balanced accuracy")
        ax.set_title("C) Panel C performance by reference")
        ax.legend(title="Reference", fontsize=8)
    ax.invert_yaxis()

    # D: Evidence tier heatmap by panel
    ax = axes[1, 1]
    tiers = []
    for panel in ["panel_A_chrx_summary_superpop", "panel_B_chrx_summary_sex",
                  "panel_C_genome_wide_snp", "panel_D_autosome_only"]:
        sp = study_root / "Results/metrics/downstream_v2" / panel / "delta_delta_summary.json"
        if sp.exists():
            tier = json.loads(sp.read_text())["confirmatory_decision"]["evidence_tier"]
            tiers.append({"panel": PANEL_LABELS[panel].replace("\n", " "), "tier": tier.split("_")[0]})
    tdf = pd.DataFrame(tiers)
    tier_map = {"A": 3, "B": 2, "C": 1}
    colors_list = [TIER_COLORS.get(f"{t}_confirmatory" if t == "A" else f"{t}_suggestive" if t == "B" else f"{t}_exploratory", "#ccc")
                   for t in tdf["tier"]]
    ax.barh(tdf["panel"], [tier_map[t] for t in tdf["tier"]], color=colors_list, edgecolor="black", linewidth=0.5)
    ax.set_xlim(0, 4)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["C (exploratory)", "B (suggestive)", "A (confirmatory)"])
    ax.set_xlabel("Evidence tier")
    ax.set_title("D) Confirmatory tier by panel")
    ax.invert_yaxis()

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"Figure5_Panel_Battery_Dashboard.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Figure 5 saved.")


def figure6_panel_a_vs_c(study_root: Path, out_dir: Path):
    """Primary revision figure: Panel A vs Panel C comparison."""
    dd = load_panel_dd(study_root)
    models = [m for m in MODEL_LABELS if m in dd["model"].unique()]
    sub = dd[dd["panel"].isin(["panel_A_chrx_summary_superpop", "panel_C_genome_wide_snp"])]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    ax = axes[0]
    y = np.arange(len(models))
    for i, panel in enumerate(["panel_A_chrx_summary_superpop", "panel_C_genome_wide_snp"]):
        pdf = sub[sub["panel"] == panel].set_index("model").reindex(models)
        ax.errorbar(pdf["delta_delta"], y + (0.15 if i else -0.15),
                    xerr=[pdf["delta_delta"] - pdf["ci_lo"], pdf["ci_hi"] - pdf["delta_delta"]],
                    fmt="o" if i else "s", capsize=4, markersize=8,
                    color=PANEL_COLORS[panel], label=PANEL_LABELS[panel].replace("\n", " "))
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(0.01, color="green", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.axvline(-0.01, color="green", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([MODEL_LABELS[m] for m in models])
    ax.set_xlabel(r"$\Delta\Delta$ (95% bootstrap CI)")
    ax.set_title("Reference effect on sex-disparity shift")
    ax.legend(loc="lower right")
    ax.invert_yaxis()

    ax2 = axes[1]
    a_vals = sub[sub["panel"] == "panel_A_chrx_summary_superpop"].set_index("model").reindex(models)["delta_delta"].abs()
    c_vals = sub[sub["panel"] == "panel_C_genome_wide_snp"].set_index("model").reindex(models)["delta_delta"].abs()
    ax2.scatter(c_vals, a_vals, s=80, c="#3C5488", edgecolors="black", linewidth=0.5, zorder=3)
    for m, x, yv in zip(models, c_vals, a_vals):
        ax2.annotate(MODEL_LABELS[m], (x, yv), textcoords="offset points", xytext=(5, 3), fontsize=8)
    lim = max(a_vals.max(), c_vals.max()) * 1.2
    ax2.plot([0, lim], [0, lim], "--", color="grey", linewidth=0.8, label="Equal effect")
    ax2.set_xlabel(r"|$\Delta\Delta$| Panel C (genome-wide SNPs)")
    ax2.set_ylabel(r"|$\Delta\Delta$| Panel A (chrX summaries)")
    ax2.set_title("Task–signal alignment")
    ax2.set_xlim(-0.002, lim)
    ax2.set_ylim(-0.002, lim)
    ax2.legend(fontsize=8)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out_dir / f"Figure6_PanelA_vs_PanelC.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Figure 6 saved.")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", default="Study_v2_Real_Data")
    parser.add_argument("--out-dir", default="Revisions/Figures")
    args = parser.parse_args()

    study_root = Path(args.study_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    figure5_panel_battery_dashboard(study_root, out_dir)
    figure6_panel_a_vs_c(study_root, out_dir)
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
