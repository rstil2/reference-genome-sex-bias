#!/usr/bin/env python3
"""Generate publication-quality figures per Gold Standard Protocol §8.

Required figures:
  1. Primary effect plot: Δ by reference with CI + ΔΔ summary
  2. Performance-by-sex forest plot across model classes
  3. Calibration-proxy: per-sex balanced accuracy distributions
  4. Sensitivity heatmap (if robustness battery results exist)
  5. Leave-one-population-out robustness panel (if LOPO results exist)

Reads from: Study_v2_Real_Data/Results/metrics/
Outputs to: Study_v2_Real_Data/Figures/
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

MODELS = [
    "logistic_regression",
    "random_forest",
    "xgboost",
    "svm",
    "mlp",
]
MODEL_LABELS = {
    "logistic_regression": "Logistic Reg.",
    "random_forest": "Random Forest",
    "xgboost": "XGBoost",
    "svm": "SVM",
    "mlp": "MLP",
}
REF_COLORS = {"GRCh38": "#E64B35", "T2T-CHM13": "#4DBBD5"}
SEX_COLORS = {"M": "#3C5488", "F": "#DC0000"}


def load_delta_delta_data(results_dir: Path):
    """Load all per-repeat ΔΔ CSVs from final metrics."""
    final = results_dir / "metrics" / "final"
    records = []
    for model in MODELS:
        csv_path = final / f"{model}_delta_delta_by_repeat.csv"
        if not csv_path.exists():
            print(f"  WARNING: missing {csv_path.name}")
            continue
        df = pd.read_csv(csv_path)
        df["model"] = model
        records.append(df)
    if not records:
        raise FileNotFoundError("No delta-delta CSVs found in metrics/final/")
    return pd.concat(records, ignore_index=True)


def load_genome_wide_metrics(results_dir: Path):
    """Load per-repeat sex-stratified metrics for both references."""
    gw = results_dir / "metrics" / "genome_wide"
    records = []
    for ref_tag in ["grch38", "t2t"]:
        ref_label = "GRCh38" if ref_tag == "grch38" else "T2T-CHM13"
        for model in MODELS:
            csv_path = gw / f"{ref_tag}_{model}_metrics.csv"
            if not csv_path.exists():
                continue
            df = pd.read_csv(csv_path)
            df["reference"] = ref_label
            df["model"] = model
            records.append(df)
    return pd.concat(records, ignore_index=True)


def bootstrap_ci(values, n_boot=5000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(values)
    means = [float(np.mean(rng.choice(values, size=n, replace=True))) for _ in range(n_boot)]
    return np.percentile(means, [2.5, 97.5])


# ─── FIGURE 1: Primary Effect Plot ──────────────────────────────────────────
def figure1_primary_effect(dd_df, out_dir, seed=42):
    """Δ by reference with CI, plus ΔΔ summary panel."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [3, 2]})

    # Left panel: Δ per reference per model
    ax = axes[0]
    model_list = [m for m in MODELS if m in dd_df["model"].unique()]
    y_positions = np.arange(len(model_list))

    for i, model in enumerate(model_list):
        mdf = dd_df[dd_df["model"] == model]
        grch38_gaps = mdf["grch38_gap"].values
        t2t_gaps = mdf["t2t_gap"].values

        # GRCh38
        mean_g = np.mean(grch38_gaps)
        ci_g = bootstrap_ci(grch38_gaps, seed=seed + i)
        ax.errorbar(mean_g, i - 0.12, xerr=[[mean_g - ci_g[0]], [ci_g[1] - mean_g]],
                     fmt="s", color=REF_COLORS["GRCh38"], capsize=4, markersize=7,
                     label="GRCh38" if i == 0 else "")

        # T2T
        mean_t = np.mean(t2t_gaps)
        ci_t = bootstrap_ci(t2t_gaps, seed=seed + i + 100)
        ax.errorbar(mean_t, i + 0.12, xerr=[[mean_t - ci_t[0]], [ci_t[1] - mean_t]],
                     fmt="o", color=REF_COLORS["T2T-CHM13"], capsize=4, markersize=7,
                     label="T2T-CHM13" if i == 0 else "")

    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([MODEL_LABELS.get(m, m) for m in model_list])
    ax.set_xlabel("Sex Disparity (Δ = BA_male − BA_female)")
    ax.set_title("A) Sex Disparity by Reference")
    ax.legend(loc="lower right", framealpha=0.9)
    ax.invert_yaxis()

    # Right panel: ΔΔ summary
    ax2 = axes[1]
    dd_means = []
    dd_ci_lo = []
    dd_ci_hi = []
    for i, model in enumerate(model_list):
        vals = dd_df[dd_df["model"] == model]["delta_delta"].values
        m = np.mean(vals)
        ci = bootstrap_ci(vals, seed=seed + i + 200)
        dd_means.append(m)
        dd_ci_lo.append(m - ci[0])
        dd_ci_hi.append(ci[1] - m)

    colors = ["#00A087" if m > 0 else "#7E6148" for m in dd_means]
    ax2.barh(y_positions, dd_means, xerr=[dd_ci_lo, dd_ci_hi],
             color=colors, capsize=4, edgecolor="black", linewidth=0.5, alpha=0.8)
    ax2.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax2.set_yticks(y_positions)
    ax2.set_yticklabels([MODEL_LABELS.get(m, m) for m in model_list])
    ax2.set_xlabel("ΔΔ = Δ_GRCh38 − Δ_T2T")
    ax2.set_title("B) Reference-Choice Differential")
    ax2.invert_yaxis()

    # Add interpretive annotation
    ax2.text(0.98, 0.02, "ΔΔ > 0 → T2T reduces sex disparity",
             transform=ax2.transAxes, ha="right", va="bottom",
             fontsize=8, fontstyle="italic", color="grey")

    plt.tight_layout()
    fig.savefig(out_dir / "Figure1_Primary_Effect.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "Figure1_Primary_Effect.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Figure 1: Primary effect plot saved.")


# ─── FIGURE 2: Forest Plot (Performance by Sex) ─────────────────────────────
def figure2_forest_plot(metrics_df, out_dir, seed=42):
    """Forest plot: BA by sex and reference across model classes."""
    fig, ax = plt.subplots(figsize=(10, 7))

    model_list = [m for m in MODELS if m in metrics_df["model"].unique()]
    y_pos = 0
    y_ticks = []
    y_labels = []

    for model in model_list:
        mdf = metrics_df[metrics_df["model"] == model]
        for ref in ["GRCh38", "T2T-CHM13"]:
            rdf = mdf[mdf["reference"] == ref]
            for sex in ["M", "F"]:
                sdf = rdf[rdf["group"] == f"sex_{sex}"]
                if sdf.empty:
                    continue
                vals = sdf["balanced_accuracy"].values
                mean_v = np.mean(vals)
                ci = bootstrap_ci(vals, seed=seed + int(y_pos * 10))
                marker = "s" if sex == "M" else "o"
                color = REF_COLORS[ref]
                alpha = 1.0 if sex == "M" else 0.6
                ax.errorbar(mean_v, y_pos, xerr=[[mean_v - ci[0]], [ci[1] - mean_v]],
                            fmt=marker, color=color, capsize=3, markersize=6, alpha=alpha)
                y_ticks.append(y_pos)
                y_labels.append(f"{MODEL_LABELS.get(model, model)} | {ref} | {sex}")
                y_pos += 1
        y_pos += 0.5  # gap between models

    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=7)
    ax.set_xlabel("Balanced Accuracy")
    ax.set_title("Performance by Sex, Reference, and Model")
    ax.invert_yaxis()

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=REF_COLORS["GRCh38"], label="GRCh38"),
        mpatches.Patch(facecolor=REF_COLORS["T2T-CHM13"], label="T2T-CHM13"),
        plt.Line2D([0], [0], marker="s", color="grey", label="Male", linestyle="None", markersize=6),
        plt.Line2D([0], [0], marker="o", color="grey", label="Female", linestyle="None", markersize=6),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    plt.tight_layout()
    fig.savefig(out_dir / "Figure2_Forest_Plot.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "Figure2_Forest_Plot.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Figure 2: Forest plot saved.")


# ─── FIGURE 3: BA Distribution by Sex (calibration proxy) ───────────────────
def figure3_ba_distributions(metrics_df, out_dir):
    """Violin plots of balanced accuracy distributions by sex and reference."""
    sex_df = metrics_df[metrics_df["group"].isin(["sex_M", "sex_F"])].copy()
    sex_df["sex"] = sex_df["group"].str.replace("sex_", "")
    sex_df["model_label"] = sex_df["model"].map(MODEL_LABELS)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax_idx, ref in enumerate(["GRCh38", "T2T-CHM13"]):
        ax = axes[ax_idx]
        rdf = sex_df[sex_df["reference"] == ref]
        if rdf.empty:
            continue
        sns.violinplot(data=rdf, x="model_label", y="balanced_accuracy",
                       hue="sex", split=True, inner="quart", ax=ax,
                       palette=SEX_COLORS, alpha=0.8, density_norm="width")
        ax.set_title(f"{ref}")
        ax.set_xlabel("")
        ax.set_ylabel("Balanced Accuracy" if ax_idx == 0 else "")
        ax.tick_params(axis="x", rotation=30)
        if ax_idx == 0:
            ax.legend(title="Sex", loc="lower left")
        else:
            ax.get_legend().remove()

    fig.suptitle("Balanced Accuracy Distribution by Sex and Reference", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_dir / "Figure3_BA_Distributions.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "Figure3_BA_Distributions.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Figure 3: BA distribution plots saved.")


# ─── FIGURE 4: Per-Superpopulation ΔΔ (if data available) ───────────────────
def figure4_superpopulation_analysis(predictions_dir, metadata_path, out_dir, seed=42):
    """Per-superpopulation BA by sex breakdown for each reference."""
    metadata = pd.read_csv(metadata_path)
    superpops = sorted(metadata["superpopulation"].unique())

    model_list = []
    pop_records = []

    for ref_tag in ["grch38", "t2t"]:
        ref_label = "GRCh38" if ref_tag == "grch38" else "T2T-CHM13"
        for model in MODELS:
            pred_path = predictions_dir / f"{ref_tag}_{model}_predictions.csv"
            if not pred_path.exists():
                continue
            if model not in model_list:
                model_list.append(model)
            preds = pd.read_csv(pred_path)
            for repeat in preds["repeat"].unique():
                rp = preds[preds["repeat"] == repeat]
                for pop in superpops:
                    pop_df = rp[rp["superpopulation"] == pop]
                    if pop_df.empty:
                        continue
                    for sex in ["M", "F"]:
                        sex_df = pop_df[pop_df["sex"] == sex]
                        if len(sex_df) < 5:
                            continue
                        correct = (sex_df["y_pred"] == sex_df["superpopulation"]).sum()
                        ba = correct / len(sex_df)
                        pop_records.append({
                            "reference": ref_label, "model": model,
                            "superpopulation": pop, "sex": sex,
                            "repeat": repeat, "accuracy": ba,
                        })

    if not pop_records:
        print("  Figure 4: Skipped (no prediction data)")
        return

    pop_df = pd.DataFrame(pop_records)

    # Compute sex gaps per superpopulation
    fig, axes = plt.subplots(1, len(superpops), figsize=(4 * len(superpops), 5), sharey=True)
    if len(superpops) == 1:
        axes = [axes]

    for ax_idx, pop in enumerate(superpops):
        ax = axes[ax_idx]
        pdf = pop_df[pop_df["superpopulation"] == pop]

        # Compute gaps per model × reference × repeat
        gap_records = []
        for model in model_list:
            for ref in ["GRCh38", "T2T-CHM13"]:
                subset = pdf[(pdf["model"] == model) & (pdf["reference"] == ref)]
                m_acc = subset[subset["sex"] == "M"].groupby("repeat")["accuracy"].mean()
                f_acc = subset[subset["sex"] == "F"].groupby("repeat")["accuracy"].mean()
                gaps = (m_acc - f_acc).dropna()
                for g in gaps:
                    gap_records.append({"model": model, "reference": ref, "gap": g})

        gdf = pd.DataFrame(gap_records)
        if gdf.empty:
            continue
        gdf["model_label"] = gdf["model"].map(MODEL_LABELS)
        sns.boxplot(data=gdf, x="model_label", y="gap", hue="reference",
                    palette=REF_COLORS, ax=ax, fliersize=2)
        ax.axhline(0, color="grey", linestyle="--", linewidth=0.7)
        ax.set_title(pop)
        ax.set_xlabel("")
        ax.set_ylabel("Sex Gap (M−F)" if ax_idx == 0 else "")
        ax.tick_params(axis="x", rotation=45)
        if ax_idx > 0:
            ax.get_legend().remove()
        else:
            ax.legend(title="Reference", fontsize=7)

    fig.suptitle("Sex Disparity by Superpopulation and Reference", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_dir / "Figure4_Superpopulation_Analysis.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "Figure4_Superpopulation_Analysis.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Figure 4: Superpopulation analysis saved.")


# ─── FIGURE 5: Summary Dashboard ────────────────────────────────────────────
def figure5_summary_dashboard(dd_df, metrics_df, out_dir, seed=42):
    """Combined summary with ΔΔ dot plot, overall performance, and repeat scatter."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    model_list = [m for m in MODELS if m in dd_df["model"].unique()]

    # Panel A: ΔΔ dot plot per repeat
    ax = axes[0, 0]
    for i, model in enumerate(model_list):
        vals = dd_df[dd_df["model"] == model]["delta_delta"].values
        jitter = np.random.default_rng(seed).uniform(-0.15, 0.15, len(vals))
        ax.scatter(vals, np.full_like(vals, i) + jitter,
                   alpha=0.5, s=25, label=MODEL_LABELS.get(model, model))
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_yticks(range(len(model_list)))
    ax.set_yticklabels([MODEL_LABELS.get(m, m) for m in model_list])
    ax.set_xlabel("ΔΔ per repeat")
    ax.set_title("A) ΔΔ Distribution by Model")
    ax.invert_yaxis()

    # Panel B: Overall BA comparison
    ax = axes[0, 1]
    overall = metrics_df[metrics_df["group"] == "overall"].copy()
    overall["model_label"] = overall["model"].map(MODEL_LABELS)
    if not overall.empty:
        sns.boxplot(data=overall, x="model_label", y="balanced_accuracy",
                    hue="reference", palette=REF_COLORS, ax=ax, fliersize=2)
        ax.set_xlabel("")
        ax.set_ylabel("Overall Balanced Accuracy")
        ax.set_title("B) Overall Performance")
        ax.tick_params(axis="x", rotation=30)
        ax.legend(title="Reference", fontsize=8)

    # Panel C: Sex gap magnitude comparison
    ax = axes[1, 0]
    sex_gap_records = []
    for model in model_list:
        mdf = dd_df[dd_df["model"] == model]
        for _, row in mdf.iterrows():
            sex_gap_records.append({
                "model": MODEL_LABELS.get(model, model),
                "reference": "GRCh38",
                "abs_gap": abs(row["grch38_gap"]),
            })
            sex_gap_records.append({
                "model": MODEL_LABELS.get(model, model),
                "reference": "T2T-CHM13",
                "abs_gap": abs(row["t2t_gap"]),
            })
    sgdf = pd.DataFrame(sex_gap_records)
    if not sgdf.empty:
        sns.boxplot(data=sgdf, x="model", y="abs_gap", hue="reference",
                    palette=REF_COLORS, ax=ax, fliersize=2)
        ax.set_xlabel("")
        ax.set_ylabel("|Sex Gap| (absolute)")
        ax.set_title("C) Absolute Sex Disparity Magnitude")
        ax.tick_params(axis="x", rotation=30)
        ax.legend(title="Reference", fontsize=8)

    # Panel D: Repeat-level paired scatter (GRCh38 gap vs T2T gap)
    ax = axes[1, 1]
    for i, model in enumerate(model_list):
        mdf = dd_df[dd_df["model"] == model]
        ax.scatter(mdf["grch38_gap"], mdf["t2t_gap"], alpha=0.5, s=25,
                   label=MODEL_LABELS.get(model, model))
    lims = ax.get_xlim()
    ax.plot(lims, lims, "--", color="grey", linewidth=0.8, zorder=0)
    ax.set_xlabel("Δ GRCh38 (M−F)")
    ax.set_ylabel("Δ T2T-CHM13 (M−F)")
    ax.set_title("D) Paired Sex Disparity")
    ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(out_dir / "Figure5_Summary_Dashboard.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "Figure5_Summary_Dashboard.pdf", bbox_inches="tight")
    plt.close(fig)
    print("  Figure 5: Summary dashboard saved.")


# ─── TABLE: Main Results ────────────────────────────────────────────────────
def generate_results_table(dd_df, metrics_df, out_dir, seed=42):
    """Protocol §8 Table 2: Main results across all endpoints, models, CIs."""
    model_list = [m for m in MODELS if m in dd_df["model"].unique()]
    rows = []

    for model in model_list:
        mdd = dd_df[dd_df["model"] == model]
        dd_vals = mdd["delta_delta"].values
        dd_mean = np.mean(dd_vals)
        dd_ci = bootstrap_ci(dd_vals, seed=seed)

        for ref_tag, ref_label in [("grch38", "GRCh38"), ("t2t", "T2T-CHM13")]:
            rdf = metrics_df[(metrics_df["model"] == model) & (metrics_df["reference"] == ref_label)]

            overall = rdf[rdf["group"] == "overall"]["balanced_accuracy"]
            male = rdf[rdf["group"] == "sex_M"]["balanced_accuracy"]
            female = rdf[rdf["group"] == "sex_F"]["balanced_accuracy"]

            overall_f1 = rdf[rdf["group"] == "overall"]["macro_f1"]

            gap_col = f"{ref_tag}_gap"
            gaps = mdd[gap_col].values if gap_col in mdd.columns else np.array([])

            rows.append({
                "Model": MODEL_LABELS.get(model, model),
                "Reference": ref_label,
                "BA_overall_mean": f"{overall.mean():.4f}" if len(overall) else "—",
                "BA_overall_CI": (f"[{bootstrap_ci(overall.values, seed=seed)[0]:.4f}, "
                                  f"{bootstrap_ci(overall.values, seed=seed)[1]:.4f}]") if len(overall) > 1 else "—",
                "BA_male_mean": f"{male.mean():.4f}" if len(male) else "—",
                "BA_female_mean": f"{female.mean():.4f}" if len(female) else "—",
                "Sex_gap_mean": f"{gaps.mean():.4f}" if len(gaps) else "—",
                "F1_overall_mean": f"{overall_f1.mean():.4f}" if len(overall_f1) else "—",
            })

        # ΔΔ row for this model
        rows.append({
            "Model": MODEL_LABELS.get(model, model),
            "Reference": "ΔΔ",
            "BA_overall_mean": "—",
            "BA_overall_CI": "—",
            "BA_male_mean": "—",
            "BA_female_mean": "—",
            "Sex_gap_mean": f"{dd_mean:.4f} [{dd_ci[0]:.4f}, {dd_ci[1]:.4f}]",
            "F1_overall_mean": "—",
        })

    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "Table2_Main_Results.csv", index=False)
    print("  Table 2: Main results table saved.")
    return table


# ─── TABLE: Cohort Characteristics ──────────────────────────────────────────
def generate_cohort_table(metadata_path, out_dir):
    """Protocol §8 Table 1: Cohort characteristics by sex and superpopulation."""
    meta = pd.read_csv(metadata_path)
    pivot = meta.groupby(["superpopulation", "sex"]).size().unstack(fill_value=0)
    pivot["Total"] = pivot.sum(axis=1)
    totals = pivot.sum(axis=0)
    totals.name = "Total"
    pivot = pd.concat([pivot, totals.to_frame().T])
    pivot.to_csv(out_dir / "Table1_Cohort_Characteristics.csv")
    print("  Table 1: Cohort characteristics saved.")


def main():
    parser = argparse.ArgumentParser(description="Generate publication figures for Protocol §8.")
    parser.add_argument("--results-dir", required=True, help="Path to Study_v2_Real_Data/Results")
    parser.add_argument("--figures-dir", required=True, help="Output directory for figures")
    parser.add_argument("--metadata", required=True, help="Path to clean_metadata_v1.csv")
    parser.add_argument("--seed", type=int, default=20260222)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    figures_dir = Path(args.figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    dd_df = load_delta_delta_data(results_dir)
    metrics_df = load_genome_wide_metrics(results_dir)

    print("Generating figures...")
    figure1_primary_effect(dd_df, figures_dir, seed=args.seed)
    figure2_forest_plot(metrics_df, figures_dir, seed=args.seed)
    figure3_ba_distributions(metrics_df, figures_dir)
    figure4_superpopulation_analysis(
        results_dir / "predictions" / "genome_wide",
        args.metadata,
        figures_dir,
        seed=args.seed,
    )
    figure5_summary_dashboard(dd_df, metrics_df, figures_dir, seed=args.seed)

    print("\nGenerating tables...")
    generate_cohort_table(args.metadata, figures_dir)
    generate_results_table(dd_df, metrics_df, figures_dir, seed=args.seed)

    print(f"\nAll outputs in: {figures_dir}")


if __name__ == "__main__":
    main()
