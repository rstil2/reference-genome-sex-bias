#!/usr/bin/env python3
"""Panel A sensitivity to relatedness, using existing chrX PSC features.

Relatedness sensitivity for Panel A, in three cuts:
  1. phase3_2504: the 2,504-sample phase 3 panel (still contains relatives)
  2. family_blocked: all 3,202 samples, families kept inside one split
  3. one_per_family: one person from each pedigree family

Does not retrain the genome-wide SNP panels. Models and the delta-delta
estimand match Code/03_train_models.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Code"))
from importlib.machinery import SourceFileLoader

train = SourceFileLoader(
    "train_models", str(ROOT / "Code" / "03_train_models.py")
).load_module()

CHR_X = ROOT / "Results/variant_analysis/per_sample_chrX.csv"
PANEL = Path(__file__).resolve().parent / "data" / "integrated_call_samples.panel"
PED = Path(__file__).resolve().parent / "data" / "20130606_g1k.ped"
OUT = Path(__file__).resolve().parent / "revision_panelA_relatedness.csv"

FEATURES = ["het_rate", "call_rate", "hom_rate", "n_het", "n_called", "n_variants"]
MODELS = ["logistic_regression", "random_forest", "xgboost", "svm", "mlp"]
SEED = 20260222
N_REPEATS = 20
ENCODED = {"xgboost", "mlp"}


def load_features(reference: str) -> pd.DataFrame:
    df = pd.read_csv(CHR_X)
    sub = df[df["reference"] == reference][
        ["patient_id", "sex", "superpopulation"] + FEATURES
    ].copy()
    return sub.sort_values("patient_id").reset_index(drop=True)


def pedigree_tables():
    panel = pd.read_csv(PANEL, sep="\t")
    ped = pd.read_csv(PED, sep="\t")
    ped = ped.rename(
        columns={"Family ID": "family_id", "Individual ID": "patient_id", "Relationship": "relationship"}
    )
    return set(panel["sample"]), ped[["family_id", "patient_id", "relationship"]]


def one_per_family_ids(ped: pd.DataFrame, allowed: set[str]) -> set[str]:
    sub = ped[ped["patient_id"].isin(allowed)].copy()
    chosen = []
    for _, g in sub.groupby("family_id"):
        unrel = g[g["relationship"].astype(str).str.startswith("unrel")]
        pick = unrel if len(unrel) else g
        chosen.append(sorted(pick["patient_id"])[0])
    return set(chosen)


def split_ids(meta: pd.DataFrame, family: pd.Series | None, seed: int):
    """Return train and test patient_id sets. Validation rows are unused, as in the original gap."""
    rng_seed = seed
    if family is None:
        strat = meta["superpopulation"].astype(str) + "_" + meta["sex"].astype(str)
        train_val, test = train_test_split(
            meta["patient_id"], test_size=0.15, random_state=rng_seed, stratify=strat
        )
        meta_tv = meta[meta["patient_id"].isin(set(train_val))]
        strat_tv = meta_tv["superpopulation"].astype(str) + "_" + meta_tv["sex"].astype(str)
        train, _val = train_test_split(
            meta_tv["patient_id"], test_size=0.15 / 0.85, random_state=rng_seed, stratify=strat_tv
        )
        return set(train), set(test)

    fam = (
        meta.assign(family_id=meta["patient_id"].map(family))
        .groupby("family_id")
        .agg(superpopulation=("superpopulation", "first"), n=("patient_id", "size"))
        .reset_index()
    )
    train_fams = []
    test_fams = []
    for sp, g in fam.groupby("superpopulation"):
        ids = g["family_id"].tolist()
        if len(ids) < 3:
            raise RuntimeError(f"Too few families in {sp}")
        tr, te = train_test_split(ids, test_size=0.15, random_state=rng_seed)
        # peel validation off train so test stays ~15% and train ~70%
        val_size = 0.15 / 0.85
        if len(tr) >= 2:
            tr, _va = train_test_split(tr, test_size=val_size, random_state=rng_seed)
        train_fams.extend(tr)
        test_fams.extend(te)
    train_fams, test_fams = set(train_fams), set(test_fams)
    mapped = meta["patient_id"].map(family)
    train_ids = set(meta.loc[mapped.isin(train_fams), "patient_id"])
    test_ids = set(meta.loc[mapped.isin(test_fams), "patient_id"])
    return train_ids, test_ids


def fit_predict(model_name: str, train_df: pd.DataFrame, test_df: pd.DataFrame, seed: int):
    model = train.build_model(model_name, seed)
    y_train = train_df["superpopulation"]
    X_train = train_df[FEATURES]
    X_test = test_df[FEATURES]
    if model_name in ENCODED:
        le = LabelEncoder()
        model.fit(X_train, le.fit_transform(y_train))
        preds = le.inverse_transform(model.predict(X_test))
    else:
        model.fit(X_train, y_train)
        preds = model.predict(X_test)
    return preds


def sex_gap(test_df: pd.DataFrame, preds: np.ndarray) -> float:
    y = test_df["superpopulation"].to_numpy()
    sex = test_df["sex"].to_numpy()
    bas = {}
    for label in ("M", "F"):
        mask = sex == label
        bas[label] = balanced_accuracy_score(y[mask], preds[mask])
    return float(bas["M"] - bas["F"])


def run_analysis(name: str, g38: pd.DataFrame, t2t: pd.DataFrame, family: pd.Series | None):
    rows = []
    meta = g38[["patient_id", "sex", "superpopulation"]].copy()
    for model_name in MODELS:
        deltas = []
        for repeat in range(N_REPEATS):
            train_ids, test_ids = split_ids(meta, family, SEED + repeat)
            gaps = {}
            for ref, frame in (("GRCh38", g38), ("T2T", t2t)):
                tr = frame[frame["patient_id"].isin(train_ids)]
                te = frame[frame["patient_id"].isin(test_ids)]
                preds = fit_predict(model_name, tr, te, SEED + repeat)
                gaps[ref] = sex_gap(te, preds)
            deltas.append(gaps["GRCh38"] - gaps["T2T"])
        arr = np.asarray(deltas, dtype=float)
        rows.append(
            {
                "analysis": name,
                "model": model_name,
                "n_samples": int(len(meta)),
                "n_repeats": N_REPEATS,
                "delta_delta_mean": float(arr.mean()),
                "delta_delta_min": float(arr.min()),
                "delta_delta_max": float(arr.max()),
                "n_positive": int((arr > 0).sum()),
            }
        )
        print(name, model_name, round(float(arr.mean()), 4), flush=True)
    return rows


def main():
    phase3, ped = pedigree_tables()
    g38 = load_features("GRCh38")
    t2t = load_features("T2T")
    family_map = ped.set_index("patient_id")["family_id"]
    # every sample should map
    missing = set(g38["patient_id"]) - set(family_map.index)
    if missing:
        raise SystemExit(f"{len(missing)} samples lack a family id")

    analyses = {
        "phase3_2504": (g38[g38.patient_id.isin(phase3)].copy(), None),
        "family_blocked_3202": (g38.copy(), family_map),
        "one_per_family": (
            g38[g38.patient_id.isin(one_per_family_ids(ped, set(g38.patient_id)))].copy(),
            None,
        ),
    }
    rows = []
    for name, (base, fam) in analyses.items():
        t2 = t2t[t2t.patient_id.isin(set(base.patient_id))].copy()
        rows.extend(run_analysis(name, base.reset_index(drop=True), t2.reset_index(drop=True), fam))
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
