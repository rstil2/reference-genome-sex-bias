#!/usr/bin/env python3
"""Compare variant position overlap between GRCh38 and T2T feature matrices."""

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_positions(columns: list[str]) -> set[int]:
    positions = set()
    for col in columns:
        m = re.search(r"chrX:(\d+)", col, re.IGNORECASE)
        if m:
            positions.add(int(m.group(1)))
    return positions


def genotype_concordance(grch38_df: pd.DataFrame, t2t_df: pd.DataFrame, shared_positions: set[int]) -> dict:
    """For shared positions (by coordinate), compare dosages if column naming allows."""
    # Position-aligned comparison requires matching patient order
    g = grch38_df.sort_values("patient_id").reset_index(drop=True)
    t = t2t_df.sort_values("patient_id").reset_index(drop=True)
    if not (g["patient_id"].values == t["patient_id"].values).all():
        return {"n_shared_positions": len(shared_positions), "concordance_computed": False}

    concordances = []
    for pos in shared_positions:
        gcol = next((c for c in g.columns if f"chrX:{pos}" in c), None)
        tcol = next((c for c in t.columns if f"chrX:{pos}" in c), None)
        if gcol is None or tcol is None:
            continue
        gv, tv = g[gcol].to_numpy(), t[tcol].to_numpy()
        mask = ~(pd.isna(gv) | pd.isna(tv))
        if mask.sum() == 0:
            continue
        concordances.append(float((gv[mask] == tv[mask]).mean()))

    if not concordances:
        return {
            "n_shared_positions": len(shared_positions),
            "concordance_computed": False,
            "mean_genotype_concordance": None,
        }

    return {
        "n_shared_positions": len(shared_positions),
        "concordance_computed": True,
        "mean_genotype_concordance": float(sum(concordances) / len(concordances)),
        "n_positions_with_concordance": len(concordances),
    }


def main():
    parser = argparse.ArgumentParser(description="Audit chrX variant position overlap across references.")
    parser.add_argument("--grch38-features", required=True)
    parser.add_argument("--t2t-features", required=True)
    parser.add_argument("--out-json", required=True)
    args = parser.parse_args()

    base_cols = {"patient_id", "sex", "superpopulation"}
    g = pd.read_csv(args.grch38_features, nrows=0)
    t = pd.read_csv(args.t2t_features, nrows=0)
    g_cols = [c for c in g.columns if c not in base_cols]
    t_cols = [c for c in t.columns if c not in base_cols]

    g_pos = parse_positions(g_cols)
    t_pos = parse_positions(t_cols)
    shared = g_pos & t_pos

    g_full = pd.read_csv(args.grch38_features)
    t_full = pd.read_csv(args.t2t_features)
    conc = genotype_concordance(g_full, t_full, shared)

    summary = {
        "generated_utc": now_utc_iso(),
        "grch38_features": args.grch38_features,
        "t2t_features": args.t2t_features,
        "n_grch38_sites": len(g_pos),
        "n_t2t_sites": len(t_pos),
        "n_shared_positions": len(shared),
        "n_grch38_only": len(g_pos - t_pos),
        "n_t2t_only": len(t_pos - g_pos),
        "fraction_shared": float(len(shared) / max(len(g_pos | t_pos), 1)),
        **conc,
    }

    out = Path(args.out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": "ok", **summary}, indent=2))


if __name__ == "__main__":
    main()
