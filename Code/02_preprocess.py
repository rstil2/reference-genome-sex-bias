#!/usr/bin/env python3
import argparse
import csv
import json
import random
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_metadata(metadata_path: Path):
    rows = []
    with metadata_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=" ")
        for raw in reader:
            row = {k: v for k, v in raw.items() if k is not None and k != ""}
            if not row:
                continue
            sex_value = row.get("Sex", "")
            if sex_value == "1":
                sex_label = "M"
            elif sex_value == "2":
                sex_label = "F"
            else:
                sex_label = "UNK"
            rows.append(
                {
                    "patient_id": row.get("SampleID", ""),
                    "sex": sex_label,
                    "population": row.get("Population", "UNK"),
                    "superpopulation": row.get("Superpopulation", "UNK"),
                }
            )
    return [r for r in rows if r["patient_id"]]


def write_clean_metadata(rows, out_csv: Path):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["patient_id", "sex", "population", "superpopulation"],
        )
        writer.writeheader()
        writer.writerows(rows)


def stratified_split(rows, train_frac, val_frac, test_frac, seed, repeats):
    if round(train_frac + val_frac + test_frac, 6) != 1.0:
        raise ValueError("train/val/test fractions must sum to 1.0")

    by_stratum = defaultdict(list)
    for r in rows:
        key = f"{r['superpopulation']}__{r['sex']}"
        by_stratum[key].append(r["patient_id"])

    all_splits = []
    for repeat in range(repeats):
        rng = random.Random(seed + repeat)
        split_rows = []

        for stratum, sample_ids in by_stratum.items():
            ids = sample_ids[:]
            rng.shuffle(ids)

            n = len(ids)
            n_train = int(n * train_frac)
            n_val = int(n * val_frac)
            n_test = n - n_train - n_val

            if n_train == 0 and n >= 3:
                n_train = 1
                n_test = max(1, n_test - 1)
            if n_val == 0 and n >= 6:
                n_val = 1
                n_test = max(1, n_test - 1)

            train_ids = ids[:n_train]
            val_ids = ids[n_train : n_train + n_val]
            test_ids = ids[n_train + n_val :]

            assert len(train_ids) + len(val_ids) + len(test_ids) == n

            for sid in train_ids:
                split_rows.append({"repeat": repeat, "patient_id": sid, "split": "train", "stratum": stratum})
            for sid in val_ids:
                split_rows.append({"repeat": repeat, "patient_id": sid, "split": "val", "stratum": stratum})
            for sid in test_ids:
                split_rows.append({"repeat": repeat, "patient_id": sid, "split": "test", "stratum": stratum})

        all_splits.extend(split_rows)

    return all_splits


def write_splits(splits, split_manifest_csv: Path):
    split_manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    with split_manifest_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["repeat", "patient_id", "split", "stratum"],
        )
        writer.writeheader()
        writer.writerows(splits)


def build_summary(rows, splits, repeats):
    sex_counts = defaultdict(int)
    superpop_counts = defaultdict(int)
    for r in rows:
        sex_counts[r["sex"]] += 1
        superpop_counts[r["superpopulation"]] += 1

    by_repeat_split = defaultdict(int)
    for s in splits:
        by_repeat_split[f"repeat_{s['repeat']}__{s['split']}"] += 1

    return {
        "generated_utc": now_utc_iso(),
        "n_samples": len(rows),
        "n_repeats": repeats,
        "sex_counts": dict(sorted(sex_counts.items())),
        "superpopulation_counts": dict(sorted(superpop_counts.items())),
        "split_counts": dict(sorted(by_repeat_split.items())),
    }


def main():
    parser = argparse.ArgumentParser(description="Project 33 preprocessing phase 1: metadata normalization and deterministic split manifest generation.")
    parser.add_argument("--metadata", required=True, help="Path to 1000G metadata file (ped_population.txt)")
    parser.add_argument("--clean-metadata-out", required=True, help="Output cleaned metadata CSV")
    parser.add_argument("--split-manifest-out", required=True, help="Output split manifest CSV")
    parser.add_argument("--summary-out", required=True, help="Output preprocessing summary JSON")
    parser.add_argument("--train-frac", type=float, default=0.70)
    parser.add_argument("--val-frac", type=float, default=0.15)
    parser.add_argument("--test-frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=20260222)
    parser.add_argument("--repeats", type=int, default=20)
    args = parser.parse_args()

    metadata_path = Path(args.metadata)
    rows = read_metadata(metadata_path)
    if not rows:
        raise ValueError(f"No samples parsed from metadata: {metadata_path}")

    clean_metadata_out = Path(args.clean_metadata_out)
    split_manifest_out = Path(args.split_manifest_out)
    summary_out = Path(args.summary_out)

    write_clean_metadata(rows, clean_metadata_out)
    splits = stratified_split(rows, args.train_frac, args.val_frac, args.test_frac, args.seed, args.repeats)
    write_splits(splits, split_manifest_out)

    summary = build_summary(rows, splits, args.repeats)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps({
        "status": "ok",
        "generated_utc": now_utc_iso(),
        "n_samples": len(rows),
        "split_rows": len(splits),
        "clean_metadata_out": str(clean_metadata_out),
        "split_manifest_out": str(split_manifest_out),
        "summary_out": str(summary_out),
    }, indent=2))


if __name__ == "__main__":
    main()
