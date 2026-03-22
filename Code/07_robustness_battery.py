#!/usr/bin/env python3
"""Robustness and sensitivity battery per Gold Standard Protocol §6.

Runs the following checks against existing results and, where needed,
invokes the training pipeline with different parameters.

Required checks:
  1. Label permutation sanity check
  2. Leave-one-superpopulation-out analysis
  3. Sex-balance stress test (downsample to parity)

NOTE: Sensitivity grid over MAF/call-rate/k (§6.2-3) requires re-running
feature extraction + training, which is computationally expensive.
This script handles the fast checks (1, 2, 3) using existing features.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC

try:
    from xgboost import XGBClassifier
    _XGBOOST_AVAILABLE = True
except Exception:
    _XGBOOST_AVAILABLE = False


REQUIRED_COLUMNS = {"patient_id", "sex", "superpopulation"}


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_model(name: str, seed: int):
    if name == "logistic_regression":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler(with_mean=False)),
            ("clf", LogisticRegression(max_iter=2000, solver="lbfgs", random_state=seed)),
        ])
    if name == "random_forest":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("clf", RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=-1,
                                           class_weight="balanced_subsample")),
        ])
    if name == "xgboost":
        if not _XGBOOST_AVAILABLE:
            raise ImportError("xgboost not installed")
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("clf", XGBClassifier(n_estimators=300, max_depth=6, learning_rate=0.1,
                                   eval_metric="mlogloss", random_state=seed, n_jobs=-1, verbosity=0)),
        ])
    if name == "svm":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler(with_mean=False)),
            ("clf", SVC(kernel="rbf", C=1.0, class_weight="balanced", random_state=seed)),
        ])
    if name == "mlp":
        return Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler(with_mean=False)),
            ("clf", MLPClassifier(hidden_layer_sizes=(256, 128), activation="relu",
                                   max_iter=500, early_stopping=True,
                                   validation_fraction=0.1, random_state=seed)),
        ])
    raise ValueError(f"Unknown model: {name}")


def train_evaluate(X_train, y_train, X_test, y_test, model_name, seed):
    """Train and evaluate, return balanced accuracy."""
    model = build_model(model_name, seed)
    if model_name in ("xgboost", "mlp"):
        le = LabelEncoder()
        le.fit(y_train)
        model.fit(X_train, le.transform(y_train))
        preds = le.inverse_transform(model.predict(X_test))
    else:
        model.fit(X_train, y_train)
        preds = model.predict(X_test)
    return balanced_accuracy_score(y_test, preds)


# ─── CHECK 1: Label Permutation Sanity (§6.7) ───────────────────────────────
def label_permutation_check(data_df, split_df, feature_cols, models, seed, n_perms=5):
    """Permute superpopulation labels; performance should collapse to chance (~0.20 for 5 classes)."""
    print("\n=== CHECK 1: Label Permutation Sanity ===")
    rng = np.random.default_rng(seed)
    results = []

    # Use only repeat 0
    rsplit = split_df[split_df["repeat"] == 0]
    train_ids = set(rsplit.loc[rsplit["split"] == "train", "patient_id"])
    test_ids = set(rsplit.loc[rsplit["split"] == "test", "patient_id"])
    train_df = data_df[data_df["patient_id"].isin(train_ids)]
    test_df = data_df[data_df["patient_id"].isin(test_ids)]

    X_train = train_df[feature_cols].values
    X_test = test_df[feature_cols].values
    y_test = test_df["superpopulation"].values

    for model_name in models:
        # True labels
        y_train_true = train_df["superpopulation"].values
        true_ba = train_evaluate(X_train, y_train_true, X_test, y_test, model_name, seed)

        # Permuted labels
        perm_bas = []
        for p in range(n_perms):
            y_perm = rng.permutation(y_train_true)
            pba = train_evaluate(X_train, y_perm, X_test, y_test, model_name, seed + p + 1)
            perm_bas.append(pba)

        mean_perm = np.mean(perm_bas)
        results.append({
            "model": model_name,
            "true_ba": true_ba,
            "mean_permuted_ba": mean_perm,
            "chance_level": 1.0 / len(np.unique(y_train_true)),
            "collapsed_to_chance": bool(mean_perm < 0.30),
        })
        print(f"  {model_name}: true={true_ba:.4f}, permuted={mean_perm:.4f}, "
              f"chance={results[-1]['chance_level']:.2f}, "
              f"collapsed={'✓' if results[-1]['collapsed_to_chance'] else '✗'}")

    return results


# ─── CHECK 2: Leave-One-Superpopulation-Out (§6.4) ──────────────────────────
def lopo_analysis(data_df, split_df, feature_cols, models, seed):
    """Train with one superpopulation removed; evaluate on all test samples."""
    print("\n=== CHECK 2: Leave-One-Superpopulation-Out ===")
    superpops = sorted(data_df["superpopulation"].unique())
    results = []

    rsplit = split_df[split_df["repeat"] == 0]
    train_ids = set(rsplit.loc[rsplit["split"] == "train", "patient_id"])
    test_ids = set(rsplit.loc[rsplit["split"] == "test", "patient_id"])

    train_df = data_df[data_df["patient_id"].isin(train_ids)]
    test_df = data_df[data_df["patient_id"].isin(test_ids)]

    X_test = test_df[feature_cols].values
    y_test = test_df["superpopulation"].values

    for model_name in models:
        for held_out_pop in superpops:
            train_subset = train_df[train_df["superpopulation"] != held_out_pop]
            X_train = train_subset[feature_cols].values
            y_train = train_subset["superpopulation"].values

            ba = train_evaluate(X_train, y_train, X_test, y_test, model_name, seed)

            # Sex-stratified
            for sex in ["M", "F"]:
                mask = test_df["sex"] == sex
                if mask.sum() == 0:
                    continue
                sex_ba = train_evaluate(
                    X_train, y_train,
                    test_df[mask][feature_cols].values,
                    test_df[mask]["superpopulation"].values,
                    model_name, seed,
                )
                results.append({
                    "model": model_name,
                    "held_out_pop": held_out_pop,
                    "sex": sex,
                    "balanced_accuracy": sex_ba,
                })

            results.append({
                "model": model_name,
                "held_out_pop": held_out_pop,
                "sex": "overall",
                "balanced_accuracy": ba,
            })
            print(f"  {model_name} (drop {held_out_pop}): BA={ba:.4f}")

    return results


# ─── CHECK 3: Sex-Balance Stress Test (§6.5) ────────────────────────────────
def sex_balance_stress_test(data_df, split_df, feature_cols, models, seed):
    """Downsample to exact sex parity and re-run."""
    print("\n=== CHECK 3: Sex-Balance Stress Test ===")
    rng = np.random.default_rng(seed + 999)
    results = []

    rsplit = split_df[split_df["repeat"] == 0]
    train_ids = set(rsplit.loc[rsplit["split"] == "train", "patient_id"])
    test_ids = set(rsplit.loc[rsplit["split"] == "test", "patient_id"])

    train_df = data_df[data_df["patient_id"].isin(train_ids)]
    test_df = data_df[data_df["patient_id"].isin(test_ids)]

    # Downsample training to sex parity per superpopulation
    balanced_dfs = []
    for pop in train_df["superpopulation"].unique():
        pop_df = train_df[train_df["superpopulation"] == pop]
        n_male = (pop_df["sex"] == "M").sum()
        n_female = (pop_df["sex"] == "F").sum()
        n_min = min(n_male, n_female)
        males = pop_df[pop_df["sex"] == "M"].sample(n=n_min, random_state=rng.integers(1e9))
        females = pop_df[pop_df["sex"] == "F"].sample(n=n_min, random_state=rng.integers(1e9))
        balanced_dfs.append(pd.concat([males, females]))

    train_balanced = pd.concat(balanced_dfs)
    X_train = train_balanced[feature_cols].values
    y_train = train_balanced["superpopulation"].values

    print(f"  Downsampled training: {len(train_balanced)} samples "
          f"(M={sum(train_balanced['sex']=='M')}, F={sum(train_balanced['sex']=='F')})")

    for model_name in models:
        for sex in ["M", "F", "overall"]:
            if sex == "overall":
                X_test = test_df[feature_cols].values
                y_test = test_df["superpopulation"].values
            else:
                mask = test_df["sex"] == sex
                X_test = test_df[mask][feature_cols].values
                y_test = test_df[mask]["superpopulation"].values

            ba = train_evaluate(X_train, y_train, X_test, y_test, model_name, seed)
            results.append({
                "model": model_name,
                "sex": sex,
                "balanced_accuracy": ba,
                "train_n": len(train_balanced),
            })
            if sex == "overall":
                print(f"  {model_name} (balanced): BA={ba:.4f}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Robustness battery per Protocol §6")
    parser.add_argument("--features-csv", required=True, help="Aligned feature matrix (either ref)")
    parser.add_argument("--split-manifest", required=True)
    parser.add_argument("--out-dir", required=True, help="Output directory for robustness results")
    parser.add_argument("--models", default="logistic_regression,random_forest",
                        help="Comma-separated model list (default: LR,RF for speed)")
    parser.add_argument("--seed", type=int, default=20260222)
    parser.add_argument("--reference-label", default="grch38", help="Label for this reference")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    models = [m.strip() for m in args.models.split(",")]

    data_df = pd.read_csv(args.features_csv)
    split_df = pd.read_csv(args.split_manifest)
    feature_cols = [c for c in data_df.columns if c not in REQUIRED_COLUMNS]

    ref = args.reference_label

    # Check 1: Label permutation
    perm_results = label_permutation_check(data_df, split_df, feature_cols, models, args.seed)
    pd.DataFrame(perm_results).to_csv(out_dir / f"{ref}_label_permutation.csv", index=False)

    # Check 2: Leave-one-superpopulation-out
    lopo_results = lopo_analysis(data_df, split_df, feature_cols, models, args.seed)
    pd.DataFrame(lopo_results).to_csv(out_dir / f"{ref}_lopo_analysis.csv", index=False)

    # Check 3: Sex-balance stress test
    balance_results = sex_balance_stress_test(data_df, split_df, feature_cols, models, args.seed)
    pd.DataFrame(balance_results).to_csv(out_dir / f"{ref}_sex_balance_stress.csv", index=False)

    # Summary
    summary = {
        "generated_utc": now_utc_iso(),
        "reference": ref,
        "models": models,
        "label_permutation_all_collapsed": all(r["collapsed_to_chance"] for r in perm_results),
        "n_lopo_tests": len(lopo_results),
        "n_balance_tests": len(balance_results),
    }
    (out_dir / f"{ref}_robustness_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nRobustness results saved to: {out_dir}")


if __name__ == "__main__":
    main()
