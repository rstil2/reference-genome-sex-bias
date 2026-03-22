#!/usr/bin/env python3
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
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler(with_mean=False)),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=2000,
                        solver="lbfgs",
                        random_state=seed,
                    ),
                ),
            ]
        )
    if name == "random_forest":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "clf",
                    RandomForestClassifier(
                        n_estimators=500,
                        max_depth=None,
                        min_samples_leaf=1,
                        random_state=seed,
                        n_jobs=-1,
                        class_weight="balanced_subsample",
                    ),
                ),
            ]
        )
    if name == "xgboost":
        if not _XGBOOST_AVAILABLE:
            raise ImportError("xgboost is not installed. Run: pip install xgboost")
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "clf",
                    XGBClassifier(
                        n_estimators=300,
                        max_depth=6,
                        learning_rate=0.1,
                        eval_metric="mlogloss",
                        random_state=seed,
                        n_jobs=-1,
                        verbosity=0,
                    ),
                ),
            ]
        )
    if name == "svm":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler(with_mean=False)),
                (
                    "clf",
                    SVC(
                        kernel="rbf",
                        C=1.0,
                        class_weight="balanced",
                        random_state=seed,
                        decision_function_shape="ovr",
                    ),
                ),
            ]
        )
    if name == "mlp":
        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler(with_mean=False)),
                (
                    "clf",
                    MLPClassifier(
                        hidden_layer_sizes=(256, 128),
                        activation="relu",
                        max_iter=500,
                        early_stopping=True,
                        validation_fraction=0.1,
                        random_state=seed,
                    ),
                ),
            ]
        )
    raise ValueError(f"Unsupported model: {name}")


def compute_metrics(y_true, y_pred):
    return {
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro")),
    }


def run_repeat(data_df, split_df, repeat, feature_cols, model_name, seed):
    rsplit = split_df[split_df["repeat"] == repeat]
    train_ids = set(rsplit.loc[rsplit["split"] == "train", "patient_id"])
    test_ids = set(rsplit.loc[rsplit["split"] == "test", "patient_id"])

    train_df = data_df[data_df["patient_id"].isin(train_ids)].copy()
    test_df = data_df[data_df["patient_id"].isin(test_ids)].copy()

    if train_df.empty or test_df.empty:
        raise ValueError(f"Repeat {repeat}: train/test split empty")

    X_train = train_df[feature_cols]
    y_train = train_df["superpopulation"]
    X_test = test_df[feature_cols]
    y_test = test_df["superpopulation"]

    model = build_model(model_name, seed + repeat)
    # XGBoost requires integer labels [0, n_classes-1]; encode transparently.
    # MLP also needs encoding because sklearn 1.5+ calls np.isnan(y_pred) in
    # early-stopping validation, which fails for string label arrays.
    if model_name in ("xgboost", "mlp"):
        le = LabelEncoder()
        le.fit(y_train)
        model.fit(X_train, le.transform(y_train))
        preds = le.inverse_transform(model.predict(X_test))
    else:
        model.fit(X_train, y_train)
        preds = model.predict(X_test)

    overall = compute_metrics(y_test, preds)

    sex_metrics = {}
    for sex in sorted(test_df["sex"].unique()):
        mask = test_df["sex"] == sex
        if mask.sum() == 0:
            continue
        sex_metrics[sex] = compute_metrics(y_test[mask], preds[mask])

    out_pred = test_df[["patient_id", "sex", "superpopulation"]].copy()
    out_pred["repeat"] = repeat
    out_pred["y_pred"] = preds
    out_pred["model"] = model_name

    return overall, sex_metrics, out_pred


def main():
    parser = argparse.ArgumentParser(description="Leakage-safe repeated-split training/evaluation for Project 33.")
    parser.add_argument("--features-csv", required=True)
    parser.add_argument("--split-manifest", required=True)
    parser.add_argument("--model", default="logistic_regression",
                        choices=["logistic_regression", "random_forest", "xgboost", "svm", "mlp"])
    parser.add_argument("--seed", type=int, default=20260222)
    parser.add_argument("--max-repeats", type=int, default=20)
    parser.add_argument("--metrics-out", required=True)
    parser.add_argument("--predictions-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    data_df = pd.read_csv(args.features_csv)
    split_df = pd.read_csv(args.split_manifest)

    missing = REQUIRED_COLUMNS - set(data_df.columns)
    if missing:
        raise ValueError(f"Features CSV missing required columns: {sorted(missing)}")

    feature_cols = [c for c in data_df.columns if c not in REQUIRED_COLUMNS]
    if not feature_cols:
        raise ValueError("No numeric feature columns found.")

    split_df["repeat"] = split_df["repeat"].astype(int)
    repeats = sorted(split_df["repeat"].unique())[: args.max_repeats]

    metric_rows = []
    pred_dfs = []

    for repeat in repeats:
        overall, sex_metrics, pred_df = run_repeat(
            data_df, split_df, repeat, feature_cols, args.model, args.seed
        )

        metric_rows.append(
            {
                "repeat": repeat,
                "model": args.model,
                "group": "overall",
                **overall,
            }
        )
        for sex, values in sex_metrics.items():
            metric_rows.append(
                {
                    "repeat": repeat,
                    "model": args.model,
                    "group": f"sex_{sex}",
                    **values,
                }
            )

        pred_dfs.append(pred_df)

    metrics_df = pd.DataFrame(metric_rows)
    preds_df = pd.concat(pred_dfs, ignore_index=True)

    metrics_out = Path(args.metrics_out)
    preds_out = Path(args.predictions_out)
    summary_out = Path(args.summary_out)
    metrics_out.parent.mkdir(parents=True, exist_ok=True)
    preds_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    metrics_df.to_csv(metrics_out, index=False)
    preds_df.to_csv(preds_out, index=False)

    overall_df = metrics_df[metrics_df["group"] == "overall"]
    male_df = metrics_df[metrics_df["group"] == "sex_M"]
    female_df = metrics_df[metrics_df["group"] == "sex_F"]

    summary = {
        "generated_utc": now_utc_iso(),
        "model": args.model,
        "n_repeats": int(len(repeats)),
        "n_samples": int(len(data_df)),
        "n_features": int(len(feature_cols)),
        "overall_balanced_accuracy_mean": float(overall_df["balanced_accuracy"].mean()),
        "overall_balanced_accuracy_std": float(overall_df["balanced_accuracy"].std(ddof=1) if len(overall_df) > 1 else 0.0),
        "male_balanced_accuracy_mean": float(male_df["balanced_accuracy"].mean()) if not male_df.empty else None,
        "female_balanced_accuracy_mean": float(female_df["balanced_accuracy"].mean()) if not female_df.empty else None,
    }

    if summary["male_balanced_accuracy_mean"] is not None and summary["female_balanced_accuracy_mean"] is not None:
        summary["sex_gap_male_minus_female"] = (
            summary["male_balanced_accuracy_mean"] - summary["female_balanced_accuracy_mean"]
        )

    summary_out.write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps({
        "status": "ok",
        "metrics_out": str(metrics_out),
        "predictions_out": str(preds_out),
        "summary_out": str(summary_out),
    }, indent=2))


if __name__ == "__main__":
    main()
