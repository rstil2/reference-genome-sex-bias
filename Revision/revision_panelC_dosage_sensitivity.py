#!/usr/bin/env python3
"""Panel C sensitivity to multiallelic dosage coding in the T2T feature matrix.

02b_extract_variant_features.py coded a genotype as the SUM of the two allele
indices. That equals the alternate-allele count only at biallelic sites. The
T2T callset keeps multiallelic records, so genotypes such as 1|2 and 2|2 were
stored as 3 and 4. The GRCh38 matrix contains only 0, 1, 2.

This script (1) counts the affected cells, (2) caps T2T dosages at 2 (1|2 and
2|2 become 2, which is correct; 0|2 stays 2 and is miscoded as 2 rather than 1,
which cannot be repaired from the stored sum), and (3) refits the five Panel C
models on the capped T2T matrix with the original splits and seeds, reusing the
stored GRCh38 results.

Usage: python3 revision_panelC_dosage_sensitivity.py [--skip-training]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent  # repository root
STUDY = ROOT
CODE = STUDY / "Code"
FEAT = STUDY / "Results/features"
METRICS = STUDY / "Results/metrics/downstream_v2/panel_C_genome_wide_snp"
SPLITS = STUDY / "Results/splits/split_manifest_v1.csv"
WORK = HERE / "revision_data" / "panelC_capped"
MODELS = ["logistic_regression", "random_forest", "xgboost", "svm", "mlp"]
BASE = ["patient_id", "sex", "superpopulation"]
SEED = 20260222


def census(df: pd.DataFrame, label: str) -> dict:
    cols = [c for c in df.columns if c not in BASE]
    x = df[cols].to_numpy(dtype=float)
    v = x[~np.isnan(x)]
    out = {
        "callset": label,
        "n_snps": len(cols),
        "cells": int(x.size),
        "missing_cells": int(np.isnan(x).sum()),
        "dosage_0": int((v == 0).sum()),
        "dosage_1": int((v == 1).sum()),
        "dosage_2": int((v == 2).sum()),
        "dosage_3": int((v == 3).sum()),
        "dosage_4": int((v == 4).sum()),
        "snps_with_dosage_gt2": int((np.nan_to_num(x) > 2).any(axis=0).sum()),
    }
    out["share_cells_gt2"] = (out["dosage_3"] + out["dosage_4"]) / max(len(v), 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-training", action="store_true")
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)

    g = pd.read_csv(FEAT / "grch38_genome_features_aligned.csv")
    t = pd.read_csv(FEAT / "t2t_genome_features_aligned.csv")
    rows = [census(g, "GRCh38"), census(t, "T2T_CHM13")]
    pd.DataFrame(rows).to_csv(HERE / "revision_panelC_dosage_census.csv", index=False)
    print(pd.DataFrame(rows).T)

    # chrX columns in Panel C: positions, to document that they are pseudoautosomal
    pos = {}
    for lab, d in (("GRCh38", g), ("T2T_CHM13", t)):
        xs = [int(c.split(":")[1]) for c in d.columns if "chrX:" in c]
        pos[lab] = {"n": len(xs), "min": min(xs), "max": max(xs)}
    (HERE / "revision_panelC_chrX_positions.json").write_text(json.dumps(pos, indent=2))
    print(pos)

    cols = [c for c in t.columns if c not in BASE]
    t_cap = t.copy()
    t_cap[cols] = t_cap[cols].clip(upper=2)
    t_path = WORK / "panel_C_t2t_capped_features.csv"
    t_cap.to_csv(t_path, index=False)

    if args.skip_training:
        return

    t2t_metrics = []
    for m in MODELS:
        out = WORK / f"t2t_capped_{m}_metrics.csv"
        cmd = [
            sys.executable, str(CODE / "03_train_models.py"),
            "--features-csv", str(t_path),
            "--split-manifest", str(SPLITS),
            "--model", m,
            "--seed", str(SEED),
            "--max-repeats", "20",
            "--metrics-out", str(out),
            "--predictions-out", str(WORK / f"t2t_capped_{m}_predictions.csv"),
            "--summary-out", str(WORK / f"t2t_capped_{m}_summary.json"),
        ]
        print("fitting", m, flush=True)
        subprocess.run(cmd, check=True)
        t2t_metrics.append(str(out))

    g_metrics = [str(METRICS / f"grch38_{m}_metrics.csv") for m in MODELS]
    cmd = [
        sys.executable, str(CODE / "04_analyze_results.py"),
        "--grch38-metrics", ",".join(g_metrics),
        "--t2t-metrics", ",".join(t2t_metrics),
        "--model-names", ",".join(MODELS),
        "--seed", str(SEED),
        "--out-json", str(HERE / "revision_panelC_capped_delta_delta.json"),
        "--out-csv", str(HERE / "revision_panelC_capped_delta_delta_by_repeat.csv"),
    ]
    subprocess.run(cmd, check=True)
    print("done")


if __name__ == "__main__":
    main()
