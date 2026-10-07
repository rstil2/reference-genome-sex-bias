#!/usr/bin/env python3
"""Regenerate Figures 1-6 for the revision with metric names that match the
stored quantities. Reads existing per-sample PSC tables and panel-battery
outputs, plus revision_count_decomposition.csv and revision_panelA_relatedness.csv.
Writes Figures/Revised_Fig*.png.
"""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
STUDY = ROOT
OUT = HERE / "Figures"
OUT.mkdir(exist_ok=True)

C38, CT2T = "#3B6FB6", "#C0392B"
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
MODELS = ["logistic_regression", "random_forest", "xgboost", "svm", "mlp"]
LABELS = {
    "logistic_regression": "Logistic Reg.",
    "random_forest": "Random Forest",
    "xgboost": "XGBoost",
    "svm": "SVM",
    "mlp": "MLP",
}

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

chrx = pd.read_csv(STUDY / "Results/variant_analysis/per_sample_chrX.csv")
decomp = pd.read_csv(HERE / "revision_count_decomposition.csv")
SHARE = "Heterozygous share of\nnon-reference genotypes"


def save(fig, name):
    fig.savefig(OUT / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def fig1():
    f = chrx[chrx["sex"] == "F"]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1, 1.6]})
    data = [f[f["reference"] == r]["het_rate"].to_numpy() for r in ("GRCh38", "T2T")]
    parts = ax[0].violinplot(data, showmedians=True, widths=0.8)
    for body, col in zip(parts["bodies"], (C38, CT2T)):
        body.set_facecolor(col)
        body.set_alpha(0.55)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[key].set_color("#444444")
    ax[0].set_xticks([1, 2])
    ax[0].set_xticklabels(["GRCh38 callset", "T2T-CHM13 callset"])
    ax[0].set_ylabel(SHARE)
    ax[0].set_title("a   Females, chrX (n = 1,603)", loc="left")

    width = 0.38
    x = np.arange(len(SUPERPOPS))
    for i, (ref, col) in enumerate((("GRCh38", C38), ("T2T", CT2T))):
        means, errs = [], []
        for sp in SUPERPOPS:
            v = f[(f["reference"] == ref) & (f["superpopulation"] == sp)]["het_rate"].to_numpy()
            means.append(v.mean())
            errs.append(1.96 * v.std(ddof=1) / np.sqrt(len(v)))
        ax[1].bar(x + (i - 0.5) * width, means, width, yerr=errs, capsize=3, color=col,
                  label="GRCh38 callset" if ref == "GRCh38" else "T2T-CHM13 callset")
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(SUPERPOPS)
    ax[1].set_ylim(0.4, 0.72)
    ax[1].set_ylabel(SHARE)
    ax[1].set_title("b   By superpopulation (mean, 95% CI)", loc="left")
    ax[1].legend(frameon=False, loc="upper right")
    fig.tight_layout()
    save(fig, "Revised_Fig1_female_chrX_het_share.png")


def count_fig(chrom, groups, group_labels, filename, title, sex="F"):
    sub = decomp[(decomp["chrom"] == chrom)]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    width = 0.38
    x = np.arange(len(groups))
    for i, (q, col, lab) in enumerate((("nHet", "#8E44AD", "Heterozygous genotypes"),
                                       ("nNonRefHom", "#E39B2D", "Non-reference homozygous genotypes"))):
        vals, lo, hi = [], [], []
        for g, s in groups:
            r = sub[(sub["quantity"] == q) & (sub["sex"] == s) & (sub["group"] == g)].iloc[0]
            vals.append(r["delta_mean"])
            lo.append(r["delta_mean"] - r["ci_lo"])
            hi.append(r["ci_hi"] - r["delta_mean"])
        ax.bar(x + (i - 0.5) * width, vals, width, yerr=[lo, hi], capsize=3, color=col, label=lab)
    ax.axhline(0, color="#444444", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(group_labels)
    ax.set_ylabel("Paired change per person,\nT2T minus GRCh38 (count)")
    ax.set_title(title, loc="left")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9)
    fig.tight_layout()
    save(fig, filename)


def fig2():
    groups = [(sp, "F") for sp in SUPERPOPS] + [("all", "F")]
    count_fig("chrX", groups, SUPERPOPS + ["All"], "Revised_Fig2_female_chrX_count_change.png",
              "Females, chrX: change in genotype counts (mean, 95% bootstrap CI)")


def fig3():
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cats = [("chrX", "F", "chrX\nfemales"), ("chrX", "M", "chrX\nmales"),
            ("chr22", "F", "chr22\nfemales"), ("chr22", "M", "chr22\nmales")]
    width = 0.38
    x = np.arange(len(cats))
    g, t = [], []
    for chrom, sex, _ in cats:
        r = decomp[(decomp["chrom"] == chrom) & (decomp["sex"] == sex) &
                   (decomp["group"] == "all") & (decomp["quantity"] == "nMissing")].iloc[0]
        g.append(r["grch38_mean"])
        t.append(r["t2t_mean"])
    ax.bar(x - width / 2, g, width, color=C38, label="GRCh38 callset")
    ax.bar(x + width / 2, t, width, color=CT2T, label="T2T-CHM13 callset")
    ax.set_xticks(x)
    ax.set_xticklabels([c[2] for c in cats])
    ax.set_ylabel("Missing genotypes per person\n(mean count among VCF records)")
    ax.set_title("Missing genotypes in each callset (GRCh38: none in any person)", loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "Revised_Fig3_missing_genotypes.png")


def fig4():
    groups = [("all", "F"), ("all", "M")]
    count_fig("chr22", groups, ["Females", "Males"], "Revised_Fig4_chr22_count_change.png",
              "chr22: change in genotype counts (mean, 95% bootstrap CI)")


def panel_dd():
    rows = []
    base = STUDY / "Results/metrics/downstream_v2"
    for pdir in sorted(base.glob("panel_*")):
        js = pdir / "delta_delta_summary.json"
        if not js.exists():
            continue
        summ = json.loads(js.read_text())
        for model, res in summ["model_results"].items():
            lo, hi = res["delta_delta_bootstrap_ci_95"]
            rows.append({"panel": pdir.name, "model": model, "dd": res["delta_delta_mean"],
                         "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


def fig5():
    dd = panel_dd()
    a = dd[dd["panel"] == "panel_A_chrx_summary_superpop"].set_index("model").reindex(MODELS)
    c = dd[dd["panel"] == "panel_C_genome_wide_snp"].set_index("model").reindex(MODELS)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    y = np.arange(len(MODELS))

    for off, frame, col, lab in ((-0.18, a, "#E64B35", "Panel A (chrX summaries)"),
                                 (0.18, c, "#4DBBD5", "Panel C (VCF-order SNPs)")):
        ax[0, 0].barh(y + off, frame["dd"], 0.34, color=col, label=lab,
                      xerr=[frame["dd"] - frame["lo"], frame["hi"] - frame["dd"]], capsize=3)
    ax[0, 0].axvline(0, color="#444444", linewidth=0.8)
    ax[0, 0].set_yticks(y)
    ax[0, 0].set_yticklabels([LABELS[m] for m in MODELS])
    ax[0, 0].invert_yaxis()
    ax[0, 0].set_xlabel("ΔΔ (95% bootstrap CI)")
    ax[0, 0].set_title("a   Change in male–female gap", loc="left")
    ax[0, 0].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

    by_rep = pd.read_csv(STUDY / "Results/metrics/downstream_v2/panel_A_chrx_summary_superpop/delta_delta_by_repeat.csv")
    for off, col, lab, key in ((-0.18, C38, "GRCh38 callset", "grch38_gap"),
                               (0.18, CT2T, "T2T-CHM13 callset", "t2t_gap")):
        vals = [by_rep[by_rep["model"] == m][key].mean() for m in MODELS]
        ax[0, 1].barh(y + off, vals, 0.34, color=col, label=lab)
    ax[0, 1].axvline(0, color="#444444", linewidth=0.8)
    ax[0, 1].set_yticks(y)
    ax[0, 1].set_yticklabels([LABELS[m] for m in MODELS])
    ax[0, 1].invert_yaxis()
    ax[0, 1].set_xlabel("Male minus female balanced accuracy")
    ax[0, 1].set_title("b   Panel A sex gap by callset", loc="left")
    ax[0, 1].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

    base = STUDY / "Results/metrics/downstream_v2/panel_C_genome_wide_snp"
    for off, col, lab, tag in ((-0.18, C38, "GRCh38 callset", "grch38"),
                               (0.18, CT2T, "T2T-CHM13 callset", "t2t")):
        vals = []
        for m in MODELS:
            df = pd.read_csv(base / f"{tag}_{m}_metrics.csv")
            vals.append(df[df["group"] == "overall"]["balanced_accuracy"].mean())
        ax[1, 0].barh(y + off, vals, 0.34, color=col, label=lab)
    ax[1, 0].set_xlim(0.9, 1.0)
    ax[1, 0].set_yticks(y)
    ax[1, 0].set_yticklabels([LABELS[m] for m in MODELS])
    ax[1, 0].invert_yaxis()
    ax[1, 0].set_xlabel("Overall balanced accuracy")
    ax[1, 0].set_title("c   Panel C accuracy by callset", loc="left")
    ax[1, 0].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

    rel = pd.read_csv(HERE / "revision_panelA_relatedness.csv")
    series = [("Submitted splits", a["dd"].to_numpy(), "#444444"),
              ("Phase 3 panel (2,504)", None, "#7FB3D5"),
              ("Family-blocked (3,202)", None, "#F5B041"),
              ("One per family (1,879)", None, "#58D68D")]
    keys = {"Phase 3 panel (2,504)": "phase3_2504", "Family-blocked (3,202)": "family_blocked_3202",
            "One per family (1,879)": "one_per_family"}
    h = 0.19
    for i, (lab, vals, col) in enumerate(series):
        if vals is None:
            sub = rel[rel["analysis"] == keys[lab]].set_index("model").reindex(MODELS)
            vals = sub["delta_delta_mean"].to_numpy()
        ax[1, 1].barh(y + (i - 1.5) * h, vals, h, color=col, label=lab)
    ax[1, 1].axvline(0, color="#444444", linewidth=0.8)
    ax[1, 1].set_yticks(y)
    ax[1, 1].set_yticklabels([LABELS[m] for m in MODELS])
    ax[1, 1].invert_yaxis()
    ax[1, 1].set_xlabel("Mean ΔΔ over 20 splits")
    ax[1, 1].set_title("d   Panel A with relatedness constrained", loc="left")
    ax[1, 1].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    fig.tight_layout()
    save(fig, "Revised_Fig5_panel_battery.png")


def fig6():
    dd = panel_dd()
    a = dd[dd["panel"] == "panel_A_chrx_summary_superpop"].set_index("model").reindex(MODELS)
    c = dd[dd["panel"] == "panel_C_genome_wide_snp"].set_index("model").reindex(MODELS)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
    y = np.arange(len(MODELS))
    for off, frame, col, mk, lab in ((-0.15, a, "#E64B35", "s", "Panel A (chrX summaries)"),
                                     (0.15, c, "#4DBBD5", "o", "Panel C (VCF-order SNPs)")):
        ax[0].errorbar(frame["dd"], y + off,
                       xerr=[frame["dd"] - frame["lo"], frame["hi"] - frame["dd"]],
                       fmt=mk, color=col, capsize=4, markersize=7, label=lab)
    ax[0].axvline(0, color="#444444", linewidth=0.8)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([LABELS[m] for m in MODELS])
    ax[0].invert_yaxis()
    ax[0].set_xlabel("ΔΔ (95% bootstrap CI)")
    ax[0].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(0.0, 0.62))
    ax[0].set_title("a   Change in male–female gap", loc="left")

    ax[1].scatter(c["dd"].abs(), a["dd"].abs(), s=60, c="#3C5488", edgecolor="black", zorder=3)
    for m in MODELS:
        ax[1].annotate(LABELS[m], (abs(c.loc[m, "dd"]), abs(a.loc[m, "dd"])),
                       xytext=(5, 3), textcoords="offset points", fontsize=8)
    lim = max(a["dd"].abs().max(), c["dd"].abs().max()) * 1.2
    ax[1].plot([0, lim], [0, lim], "--", color="grey", linewidth=0.8)
    ax[1].set_xlim(-0.002, lim)
    ax[1].set_ylim(-0.002, lim)
    ax[1].set_xlabel("|ΔΔ|, Panel C")
    ax[1].set_ylabel("|ΔΔ|, Panel A")
    ax[1].set_title("b   Size of the change, by panel", loc="left")
    fig.tight_layout()
    save(fig, "Revised_Fig6_panelA_vs_panelC.png")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
