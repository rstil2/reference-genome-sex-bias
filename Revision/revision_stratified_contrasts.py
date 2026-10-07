#!/usr/bin/env python3
"""Revision analyses on existing per-sample bcftools PSC tables.

Does not re-call reads. Writes revision_stratified_contrasts.csv with:
  - female chrX heterozygosity, paired T2T minus GRCh38, by superpopulation
  - Kruskal-Wallis test of those paired differences across superpopulations
  - chromosome 22 paired heterozygosity differences by sex
  - mean PSC record counts by sex and reference on chrX
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
CHR_X = ROOT / "Results/variant_analysis/per_sample_chrX.csv"
CHR_22 = ROOT / "Results/variant_analysis/per_sample_chr22.csv"
OUT = Path(__file__).resolve().parent / "revision_stratified_contrasts.csv"

RNG = np.random.default_rng(20260222)
N_BOOT = 2000


def boot_mean(x, n=N_BOOT):
    x = np.asarray(x, dtype=float)
    idx = RNG.integers(0, len(x), size=(n, len(x)))
    means = x[idx].mean(axis=1)
    return float(np.mean(x)), float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def paired(df, value="het_rate"):
    wide = df.pivot(index="patient_id", columns="reference", values=value).dropna()
    return wide["GRCh38"], wide["T2T"], wide["T2T"] - wide["GRCh38"]


def main():
    chr_x = pd.read_csv(CHR_X)
    chr_22 = pd.read_csv(CHR_22)
    rows = []

    females = chr_x[chr_x["sex"] == "F"].copy()
    pieces = []
    for sp, g in females.groupby("superpopulation"):
        g38, t2t, delta = paired(g)
        mean, lo, hi = boot_mean(delta.to_numpy())
        _, p = stats.wilcoxon(t2t, g38)
        rows.append(
            {
                "analysis": "chrX_female_paired_het",
                "group": sp,
                "n": int(len(delta)),
                "estimate": mean,
                "ci_lo": lo,
                "ci_hi": hi,
                "p_value": p,
                "note": "estimate is mean (T2T - GRCh38); bootstrap 95% CI; Wilcoxon signed-rank p",
            }
        )
        pieces.append(delta.to_numpy())

    h, kp = stats.kruskal(*pieces)
    rows.append(
        {
            "analysis": "chrX_female_het_kruskal",
            "group": "five_superpopulations",
            "n": int(sum(len(p) for p in pieces)),
            "estimate": h,
            "ci_lo": np.nan,
            "ci_hi": np.nan,
            "p_value": kp,
            "note": "estimate is the Kruskal-Wallis H statistic on the paired differences",
        }
    )

    for sex, label in [("F", "female"), ("M", "male")]:
        g38, t2t, delta = paired(chr_22[chr_22["sex"] == sex])
        mean, lo, hi = boot_mean(delta.to_numpy())
        rows.append(
            {
                "analysis": "chr22_paired_het",
                "group": label,
                "n": int(len(delta)),
                "estimate": mean,
                "ci_lo": lo,
                "ci_hi": hi,
                "p_value": np.nan,
                "note": "estimate is mean (T2T - GRCh38); bootstrap 95% CI",
            }
        )

    counts = (
        chr_x.groupby(["sex", "reference"])[["n_total", "n_called", "n_missing", "call_rate"]]
        .mean()
        .reset_index()
    )
    for rec in counts.to_dict(orient="records"):
        rows.append(
            {
                "analysis": "chrX_mean_n_total",
                "group": f"{rec['sex']}_{rec['reference']}",
                "n": np.nan,
                "estimate": rec["n_total"],
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "p_value": np.nan,
                "note": f"mean n_called={rec['n_called']:.4f}; mean n_missing={rec['n_missing']:.4f}; mean call_rate={rec['call_rate']:.6f}",
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
