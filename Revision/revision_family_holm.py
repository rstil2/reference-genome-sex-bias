#!/usr/bin/env python3
"""Holm-Bonferroni across the whole panel x model family.

Uses the raw paired sign-permutation p-values stored in the panel-battery
summaries (04_analyze_results.py). Reports Holm-adjusted p within each panel
(as in the original analysis) and across all 20 tests (4 panels x 5 models).
The permutation unit is a random split, and splits share people, so these p
values describe stability across splits, not sampling error in a population.
"""

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "Results/metrics/downstream_v2"
PANELS = [
    "panel_A_chrx_summary_superpop",
    "panel_B_chrx_summary_sex",
    "panel_C_genome_wide_snp",
    "panel_D_autosome_only",
]


def holm(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    adj = [0.0] * len(pvals)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, pvals[i] * (len(pvals) - rank)))
        adj[i] = running
    return adj


def main():
    rows = []
    for panel in PANELS:
        summ = json.loads((BASE / panel / "delta_delta_summary.json").read_text())
        for model, r in summ["model_results"].items():
            rows.append(
                {
                    "panel": panel,
                    "model": model,
                    "delta_delta": r["delta_delta_mean"],
                    "raw_p": r["paired_sign_permutation_p_value"],
                    "holm_within_panel": r["holm_adjusted_p_value"],
                }
            )
    df = pd.DataFrame(rows)
    df["holm_all_20"] = holm(df["raw_p"].tolist())
    subset = df["panel"] != "panel_B_chrx_summary_sex"
    df["holm_all_15_excl_B"] = float("nan")
    df.loc[subset, "holm_all_15_excl_B"] = holm(df.loc[subset, "raw_p"].tolist())
    df.to_csv(HERE / "revision_family_holm.csv", index=False)
    pd.set_option("display.width", 200)
    print(df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
