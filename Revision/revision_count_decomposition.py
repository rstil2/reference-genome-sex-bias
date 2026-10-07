#!/usr/bin/env python3
"""Decompose the stored heterozygosity ratio into its genotype counts.

The per-sample PSC tables were written by 06_variant_level_bias_analysis.py.
That parser mapped bcftools columns nRefHom and nNonRefHom onto the same name
(n_ref_hom), so the column called n_ref_hom holds nNonRefHom (non-reference
homozygous genotypes) and the true reference-homozygous count is not stored.
n_hom_alt is zero in every row. What the tables do hold correctly: nHet,
nNonRefHom (under the name n_ref_hom), nHapRef, nHapAlt, nMissing.

This script therefore reports:
  het_share = nHet / (nHet + nNonRefHom) for diploid samples (the stored
              het_rate, which excludes reference-homozygous genotypes)
  and the paired T2T-minus-GRCh38 change in each count, inside each group.
It does not recompute call rate, which would need the reference-homozygous
counts.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CHR_X = ROOT / "Results/variant_analysis/per_sample_chrX.csv"
CHR_22 = ROOT / "Results/variant_analysis/per_sample_chr22.csv"
OUT = Path(__file__).resolve().parent / "revision_count_decomposition.csv"

RNG = np.random.default_rng(20260222)
N_BOOT = 2000


def boot(x):
    x = np.asarray(x, dtype=float)
    idx = RNG.integers(0, len(x), size=(N_BOOT, len(x)))
    m = x[idx].mean(axis=1)
    return float(x.mean()), float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


def paired_rows(df, label, chrom, sex, group):
    rows = []
    for col, name in [
        ("n_het", "nHet"),
        ("n_ref_hom", "nNonRefHom"),
        ("n_hap_ref", "nHapRef"),
        ("n_hap_alt", "nHapAlt"),
        ("n_missing", "nMissing"),
        ("het_rate", "stored_het_rate"),
    ]:
        wide = df.pivot(index="patient_id", columns="reference", values=col).dropna()
        if wide["GRCh38"].abs().sum() + wide["T2T"].abs().sum() == 0:
            continue
        d = (wide["T2T"] - wide["GRCh38"]).to_numpy()
        mean, lo, hi = boot(d)
        rows.append(
            {
                "chrom": chrom,
                "sex": sex,
                "group": group,
                "quantity": name,
                "n": int(len(d)),
                "grch38_mean": float(wide["GRCh38"].mean()),
                "t2t_mean": float(wide["T2T"].mean()),
                "delta_mean": mean,
                "ci_lo": lo,
                "ci_hi": hi,
            }
        )
    return rows


def main():
    rows = []
    for chrom, path in (("chrX", CHR_X), ("chr22", CHR_22)):
        df = pd.read_csv(path)
        for sex in ("F", "M"):
            sub = df[df["sex"] == sex]
            rows += paired_rows(sub, "all", chrom, sex, "all")
            if chrom == "chrX" and sex == "F":
                for sp, g in sub.groupby("superpopulation"):
                    rows += paired_rows(g, sp, chrom, sex, sp)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    pd.set_option("display.width", 200)
    print(out.round(4).to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
