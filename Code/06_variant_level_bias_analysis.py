#!/usr/bin/env python3
"""06_variant_level_bias_analysis.py

Variant-level bias analysis: GRCh38 vs T2T-CHM13 sex chromosome calling
statistics stratified by biological sex and superpopulation.

Primary question: does reference genome choice differentially affect chrX
variant calling in males vs females?

Key metrics per sample:
  - call_rate     : fraction of sites with non-missing GT
  - het_rate      : fraction of called sites that are heterozygous
  - var_rate      : fraction of called sites that are non-reference
  - norm_het_rate : chrX het_rate / chr22 het_rate (auto-normalised)

Analyses performed:
  1. chrX whole-chromosome  (GRCh38 and T2T)
  2. chr22 autosomal control (GRCh38 and T2T)
  3. chrX PAR1 and PAR2 sub-regions (GRCh38 and T2T, correct coordinates each)
  4. Group comparisons: male vs female × reference; Δhet (T2T − GRCh38) by sex
  5. Superpopulation-stratified summary

Outputs (all in Study_v2_Real_Data/Results/variant_analysis/):
  - per_sample_chrX.csv        : per-sample chrX metrics, both refs
  - per_sample_chr22.csv       : per-sample chr22 control, both refs
  - per_sample_par.csv         : per-sample PAR1/PAR2 metrics
  - group_comparisons.csv      : group means, SDs, p-values, effect sizes
  - superpop_summary.csv       : superpopulation-stratified metrics
  - summary.json               : machine-readable summary for manuscript

Usage:
    python 06_variant_level_bias_analysis.py
    python 06_variant_level_bias_analysis.py --test-mode      # 2 Mb pilot
    python 06_variant_level_bias_analysis.py --skip-chr22     # skip autosomal ctrl
    python 06_variant_level_bias_analysis.py --skip-par       # skip PAR sub-regions
"""

import argparse
import json
import logging
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
STUDY_DIR    = PROJECT_ROOT / "Study_v2_Real_Data"
DATA_DIR     = STUDY_DIR / "Data"
RESULTS_DIR  = STUDY_DIR / "Results"
METADATA_CSV = RESULTS_DIR / "manifests" / "clean_metadata_v1.csv"
OUTPUT_DIR   = RESULTS_DIR / "variant_analysis"
BCFTOOLS     = "/Users/stillwell/miniconda3/bin/bcftools"

# GRCh38 VCFs
_GRCH38 = DATA_DIR / "GRCh38"
_T2T    = DATA_DIR / "T2T_CHM13"

VCF = {
    "GRCh38": {
        **{
            str(c): _GRCH38 / f"CCDG_14151_B01_GRM_WGS_2020-08-05_chr{c}.filtered.shapeit2-duohmm-phased.vcf.gz"
            for c in range(1, 23)
        },
        "X": _GRCH38 / "CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz",
    },
    "T2T": {
        **{
            str(c): _T2T / f"1KGP.CHM13v2.0.chr{c}.recalibrated.snp_indel.pass.vcf.gz"
            for c in range(1, 23)
        },
        "X": _T2T / "1KGP.CHM13v2.0.chrX.recalibrated.snp_indel.pass.vcf.gz",
    },
}

# PAR coordinates per reference (chrX)
PAR = {
    "GRCh38": {
        "PAR1": "chrX:60001-2699520",
        "PAR2": "chrX:154931044-155260560",
    },
    "T2T": {
        "PAR1": "chrX:10001-2781479",
        "PAR2": "chrX:155701383-156030895",
    },
}

# Test-mode regions (small window; validates parsing quickly)
TEST_REGIONS = {
    "GRCh38": {"X": "chrX:2700000-4700000", "22": "chr22:20000000-22000000"},
    "T2T":    {"X": "chrX:2800000-4800000", "22": "chr22:20000000-22000000"},
}


# ── bcftools stats / PSC parsing ──────────────────────────────────────────────

def _read_psc_header(text: str) -> dict[int, str]:
    """Extract 0-indexed column position → name map from bcftools stats header.

    The stats output contains a comment like:
        # PSC\\t[2]id\\t[3]sample\\t[4]nRefHom\\t[5]nNonRefHom\\t[6]nHet\\t...
    where [N] is 1-based.  We return {0-based_index: canonical_name}.
    """
    for line in text.splitlines():
        if line.startswith("# PSC\t"):
            parts = line.split("\t")
            col_map = {}
            for i, part in enumerate(parts):
                m = re.match(r"\[(\d+)\](\S+)", part.strip())
                if m:
                    col_map[int(m.group(1)) - 1] = m.group(2)
            return col_map
    return {}


def _canonical_col(raw: str) -> str:
    # bcftools 1.23 PSC columns (from "# PSC [2]id [3]sample [4]nRefHom …"):
    #   nRefHom, nNonRefHom, nHets, nTransitions, nTransversions,
    #   nIndels, average depth, nSingletons, nHapRef, nHapAlt, nMissing
    lc = raw.lower()
    if lc == "sample":
        return "sample_id"
    if "nrefhom" in lc:
        return "n_ref_hom"
    if "nnonrefhom" in lc:
        return "n_hom_alt"
    if lc in ("nhet", "nhets"):        # bcftools ≤1.19: nHet, ≥1.20: nHets
        return "n_het"
    if "nhapref" in lc:
        return "n_hap_ref"
    if "nhapalt" in lc:
        return "n_hap_alt"
    if "nmissing" in lc:
        return "n_missing"
    if "nsnp" in lc:
        return "n_snps"
    if "nindel" in lc:
        return "n_indels"
    return lc


def parse_psc(text: str) -> pd.DataFrame:
    """Parse PSC lines from bcftools stats stdout into a DataFrame."""
    col_map = _read_psc_header(text)

    psc_lines = [l for l in text.splitlines() if l.startswith("PSC\t")]
    if not psc_lines:
        log.warning("No PSC lines found in bcftools stats output.")
        return pd.DataFrame()

    records = []
    for line in psc_lines:
        parts = line.split("\t")
        r = {"_raw": line}
        for idx, name in col_map.items():
            if idx < len(parts):
                r[_canonical_col(name)] = parts[idx]
        records.append(r)

    df = pd.DataFrame(records).drop(columns=["_raw"], errors="ignore")

    int_cols = ["n_ref_hom", "n_hom_alt", "n_het", "n_hap_ref", "n_hap_alt",
                "n_missing", "n_snps", "n_indels"]
    for c in int_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    return df


def run_bcftools_stats(
    vcf_path: Path, region: str | None = None, timeout: int = 7200
) -> pd.DataFrame:
    """Run `bcftools stats -s - [-r region] vcf` and return per-sample PSC df."""
    cmd = [BCFTOOLS, "stats", "-s", "-"]
    if region:
        cmd += ["-r", region]
    cmd.append(str(vcf_path))

    log.info("  bcftools stats %s %s", f"-r {region}" if region else "", vcf_path.name)
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        log.error("bcftools stats timed out (>%ds) for %s", timeout, vcf_path.name)
        return pd.DataFrame()
    except subprocess.CalledProcessError as exc:
        log.error(
            "bcftools stats returned %d for %s:\n%s",
            exc.returncode, vcf_path.name, exc.stderr[:800],
        )
        return pd.DataFrame()

    df = parse_psc(result.stdout)
    log.info("    → %d samples parsed", len(df))
    return df


# ── Derived metrics ───────────────────────────────────────────────────────────

def add_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add call_rate, het_rate, var_rate, hom_alt_rate to a PSC dataframe."""
    if df.empty:
        return df
    df = df.copy()

    # Diploid genotype counts
    hom_ref = df.get("n_ref_hom", pd.Series(0, index=df.index))
    hom_alt = df.get("n_hom_alt", pd.Series(0, index=df.index))
    het     = df.get("n_het",     pd.Series(0, index=df.index))
    missing = df.get("n_missing", pd.Series(0, index=df.index))

    # Haploid counts (males on chrX in some pipelines)
    hap_ref = df.get("n_hap_ref", pd.Series(0, index=df.index))
    hap_alt = df.get("n_hap_alt", pd.Series(0, index=df.index))

    called  = hom_ref + hom_alt + het + hap_ref + hap_alt
    total   = called + missing

    df["n_total"]    = total
    df["n_called"]   = called
    df["n_variants"] = hom_alt + het + hap_alt  # non-ref called genotypes
    # Preserve raw counts as explicit columns (used in downstream selection)
    df["n_het"]     = het
    df["n_hom_alt"] = hom_alt
    df["n_hap_alt"] = hap_alt

    df["call_rate"] = np.where(total  > 0, called / total,  np.nan)
    df["het_rate"]  = np.where(called > 0, het    / called, np.nan)
    df["var_rate"]  = np.where(called > 0, (hom_alt + het + hap_alt) / called, np.nan)

    # Hom-alt fraction among variant sites (high in males on chrX → hemizygous)
    n_var = hom_alt + het + hap_alt
    df["hom_rate"] = np.where(n_var > 0, (hom_alt + hap_alt) / n_var, np.nan)

    return df


# ── Statistics helpers ────────────────────────────────────────────────────────

def mannwhitney(a: np.ndarray, b: np.ndarray) -> dict:
    """Two-sided Mann-Whitney U test + effect size (rank-biserial r)."""
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return {"stat": np.nan, "p": np.nan, "effect_r": np.nan, "n_a": len(a), "n_b": len(b)}
    u, p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
    r = 1 - (2 * u) / (len(a) * len(b))  # rank-biserial correlation
    return {"stat": float(u), "p": float(p), "effect_r": float(r),
            "n_a": len(a), "n_b": len(b)}


def wilcoxon_paired(a: np.ndarray, b: np.ndarray) -> dict:
    """Wilcoxon signed-rank test for paired samples (a vs b, same ordering)."""
    diff = a - b
    valid = ~np.isnan(diff)
    diff = diff[valid]
    if len(diff) < 10:
        return {"stat": np.nan, "p": np.nan, "n": len(diff)}
    try:
        stat, p = scipy_stats.wilcoxon(diff, alternative="two-sided")
        return {"stat": float(stat), "p": float(p), "n": int(len(diff))}
    except ValueError as e:
        return {"stat": np.nan, "p": np.nan, "n": int(len(diff)), "error": str(e)}


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt((np.std(a, ddof=1) ** 2 + np.std(b, ddof=1) ** 2) / 2)
    if pooled == 0:
        return np.nan
    return float((np.mean(a) - np.mean(b)) / pooled)


def ci95_bootstrap(x: np.ndarray, n_boot: int = 2000) -> tuple[float, float]:
    """Bootstrap 95% CI for the mean."""
    x = x[~np.isnan(x)]
    if len(x) < 4:
        return (np.nan, np.nan)
    rng = np.random.default_rng(42)
    boots = rng.choice(x, size=(n_boot, len(x)), replace=True).mean(axis=1)
    return (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))


# ── Comparison table builder ──────────────────────────────────────────────────

def make_group_row(
    label: str,
    ref: str,
    chrom: str,
    region: str,
    metric: str,
    grp_a_name: str,
    grp_b_name: str,
    a: np.ndarray,
    b: np.ndarray,
) -> dict:
    mw = mannwhitney(a, b)
    lo_a, hi_a = ci95_bootstrap(a)
    lo_b, hi_b = ci95_bootstrap(b)
    return {
        "label":       label,
        "reference":   ref,
        "chromosome":  chrom,
        "region":      region,
        "metric":      metric,
        "group_a":     grp_a_name,
        "group_b":     grp_b_name,
        "mean_a":      float(np.nanmean(a)),
        "sd_a":        float(np.nanstd(a, ddof=1)),
        "ci95_lo_a":   lo_a,
        "ci95_hi_a":   hi_a,
        "n_a":         mw["n_a"],
        "mean_b":      float(np.nanmean(b)),
        "sd_b":        float(np.nanstd(b, ddof=1)),
        "ci95_lo_b":   lo_b,
        "ci95_hi_b":   hi_b,
        "n_b":         mw["n_b"],
        "mean_diff":   float(np.nanmean(a) - np.nanmean(b)),
        "mw_u":        mw["stat"],
        "mw_p":        mw["p"],
        "effect_r":    mw["effect_r"],
        "cohens_d":    cohens_d(a, b),
    }


# ── Figures ───────────────────────────────────────────────────────────────────

def violin_by_sex_ref(data: pd.DataFrame, metric: str, title: str, out_path: Path):
    """Violin plot of `metric` grouped by sex × reference."""
    fig, ax = plt.subplots(figsize=(8, 5))
    positions = []
    xs = []
    labels = []
    colors = []
    colour_map = {"GRCh38": "#2196F3", "T2T": "#F44336"}
    pos = 1
    for ref in ["GRCh38", "T2T"]:
        for sex in ["F", "M"]:
            sub = data[(data["reference"] == ref) & (data["sex"] == sex)][metric].dropna().values
            if len(sub) == 0:
                continue
            vp = ax.violinplot(sub, positions=[pos], showmedians=True, widths=0.7)
            for part in vp["bodies"]:
                part.set_facecolor(colour_map[ref])
                part.set_alpha(0.6)
            for part in ("cmedians", "cmins", "cmaxes", "cbars"):
                if part in vp:
                    vp[part].set_color(colour_map[ref])
            positions.append(pos)
            xs.append(sub)
            labels.append(f"{ref}\n{'Female' if sex == 'F' else 'Male'}")
            colors.append(colour_map[ref])
            pos += 1
        pos += 0.4  # gap between refs

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(metric.replace("_", " ").title(), fontsize=11)
    ax.set_title(title, fontsize=12)
    patches = [mpatches.Patch(color=c, label=r) for r, c in colour_map.items()]
    ax.legend(handles=patches, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved figure: %s", out_path)


def scatter_norm_het(data: pd.DataFrame, out_path: Path):
    """Scatter: normalised chrX het rate (chrX/chr22) coloured by sex × reference."""
    fig, ax = plt.subplots(figsize=(6, 5))
    marker_map = {"F": "o", "M": "s"}
    colour_map = {"GRCh38": "#2196F3", "T2T": "#F44336"}
    for ref in ["GRCh38", "T2T"]:
        for sex in ["F", "M"]:
            sub = data[(data["reference"] == ref) & (data["sex"] == sex)]
            ax.scatter(
                sub["chr22_het_rate"],
                sub["chrX_het_rate"],
                marker=marker_map[sex],
                color=colour_map[ref],
                alpha=0.20,
                s=10,
                label=f"{ref} {'♀' if sex == 'F' else '♂'}",
                rasterized=True,
            )
    ax.set_xlabel("Chr22 het rate (autosomal baseline)", fontsize=11)
    ax.set_ylabel("ChrX het rate", fontsize=11)
    ax.set_title("ChrX vs Chr22 heterozygosity by sex × reference", fontsize=11)
    # Identity line
    lim = (0, max(data["chr22_het_rate"].max(), data["chrX_het_rate"].max()) * 1.05)
    ax.plot(lim, lim, "k--", lw=0.8, alpha=0.5, label="x = y")
    ax.set_xlim(lim)
    ax.set_ylim(0, lim[1])
    ax.legend(markerscale=2, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved figure: %s", out_path)


def bar_superpop(data: pd.DataFrame, metric: str, title: str, out_path: Path):
    """Grouped bar chart of metric by superpopulation × sex × reference."""
    superpops = sorted(data["superpopulation"].dropna().unique())
    refs = ["GRCh38", "T2T"]
    sexes = ["F", "M"]
    colour_map = {
        ("GRCh38", "F"): "#1565C0",
        ("GRCh38", "M"): "#90CAF9",
        ("T2T", "F"):    "#B71C1C",
        ("T2T", "M"):    "#EF9A9A",
    }

    n_sp = len(superpops)
    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.18
    x = np.arange(n_sp)

    for i, (ref, sex) in enumerate([(r, s) for r in refs for s in sexes]):
        means, errs = [], []
        for sp in superpops:
            sub = data[
                (data["reference"] == ref)
                & (data["sex"] == sex)
                & (data["superpopulation"] == sp)
            ][metric].dropna()
            means.append(sub.mean())
            errs.append(sub.std(ddof=1) / np.sqrt(len(sub)) if len(sub) > 1 else 0)
        offset = (i - 1.5) * width
        ax.bar(
            x + offset, means, width,
            yerr=errs, capsize=3,
            color=colour_map[(ref, sex)],
            label=f"{ref} {'Female' if sex == 'F' else 'Male'}",
            error_kw={"linewidth": 0.8},
        )

    ax.set_xticks(x)
    ax.set_xticklabels(superpops, fontsize=10)
    ax.set_ylabel(metric.replace("_", " ").title(), fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log.info("Saved figure: %s", out_path)


# ── Main pipeline ─────────────────────────────────────────────────────────────

def collect_stats(ref: str, chrom_key: str, region: str | None, test_mode: bool) -> pd.DataFrame:
    """Run bcftools stats for one ref × chromosome (× optional region)."""
    vcf = VCF[ref][chrom_key]
    if not vcf.exists():
        log.error("VCF not found: %s", vcf)
        return pd.DataFrame()

    # In test-mode, override region to a 2 Mb window
    if test_mode and region is None:
        chrom_label = "X" if chrom_key == "X" else chrom_key
        tr = TEST_REGIONS.get(ref, {}).get(chrom_label)
        if tr:
            region = tr
            log.info("TEST MODE: restricting to %s", region)

    df = run_bcftools_stats(vcf, region=region)
    if df.empty:
        return pd.DataFrame()
    df = add_metrics(df)
    df["reference"] = ref
    df["chrom_key"] = chrom_key
    df["region"]    = region or f"chr{chrom_key}"
    return df


def main():
    parser = argparse.ArgumentParser(description="Variant-level bias analysis")
    parser.add_argument("--test-mode",   action="store_true",
                        help="Run on 2 Mb pilot region for quick validation")
    parser.add_argument("--skip-chr22",  action="store_true",
                        help="Skip chr22 autosomal control")
    parser.add_argument("--skip-par",    action="store_true",
                        help="Skip PAR1/PAR2 sub-region analysis")
    parser.add_argument("--output-dir",  type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    log.info("=== Variant-Level Bias Analysis ===")
    log.info("Output directory: %s", out)
    if args.test_mode:
        log.info("*** TEST MODE: results cover 2 Mb window only ***")

    # ── 1. Load metadata ──────────────────────────────────────────────────────
    log.info("Loading metadata from %s", METADATA_CSV)
    meta = pd.read_csv(METADATA_CSV)
    # normalise sex column: M/F
    meta["sex"] = meta["sex"].str.strip().str.upper().str[0]

    # ── 2. chrX stats (primary analysis) ─────────────────────────────────────
    log.info("\n--- chrX whole-chromosome stats ---")
    chrx_parts = []
    for ref in ["GRCh38", "T2T"]:
        df = collect_stats(ref, "X", region=None, test_mode=args.test_mode)
        if not df.empty:
            chrx_parts.append(df)

    if not chrx_parts:
        log.error("No chrX data collected. Aborting.")
        sys.exit(1)

    chrx_df = pd.concat(chrx_parts, ignore_index=True)
    # Merge metadata
    chrx_df = chrx_df.rename(columns={"sample_id": "patient_id"})
    chrx_df = chrx_df.merge(
        meta[["patient_id", "sex", "superpopulation", "population"]],
        on="patient_id", how="inner"
    )
    log.info("chrX: %d sample-reference records after metadata merge", len(chrx_df))

    # Save chrX per-sample
    chrx_out = out / "per_sample_chrX.csv"
    chrx_df.to_csv(chrx_out, index=False)
    log.info("Saved: %s", chrx_out)

    # ── 3. chr22 autosomal control ────────────────────────────────────────────
    chr22_df = pd.DataFrame()
    if not args.skip_chr22:
        log.info("\n--- chr22 autosomal control ---")
        chr22_parts = []
        for ref in ["GRCh38", "T2T"]:
            df = collect_stats(ref, "22", region=None, test_mode=args.test_mode)
            if not df.empty:
                chr22_parts.append(df)

        if chr22_parts:
            chr22_df = pd.concat(chr22_parts, ignore_index=True)
            chr22_df = chr22_df.rename(columns={"sample_id": "patient_id"})
            chr22_df = chr22_df.merge(
                meta[["patient_id", "sex", "superpopulation"]],
                on="patient_id", how="inner"
            )
            chr22_out = out / "per_sample_chr22.csv"
            chr22_df.to_csv(chr22_out, index=False)
            log.info("Saved: %s", chr22_out)

    # ── 4. PAR sub-region stats ───────────────────────────────────────────────
    par_df = pd.DataFrame()
    if not args.skip_par:
        log.info("\n--- PAR1 / PAR2 sub-region stats ---")
        par_parts = []
        for ref in ["GRCh38", "T2T"]:
            for par_name, region in PAR[ref].items():
                log.info("  %s %s: %s", ref, par_name, region)
                df = run_bcftools_stats(VCF[ref]["X"], region=region)
                if df.empty:
                    continue
                df = add_metrics(df)
                df["reference"] = ref
                df["chrom_key"] = "X"
                df["region"]    = region
                df["par_region"] = par_name
                par_parts.append(df)

        if par_parts:
            par_df = pd.concat(par_parts, ignore_index=True)
            par_df = par_df.rename(columns={"sample_id": "patient_id"})
            par_df = par_df.merge(
                meta[["patient_id", "sex", "superpopulation"]],
                on="patient_id", how="inner"
            )
            par_out = out / "per_sample_par.csv"
            par_df.to_csv(par_out, index=False)
            log.info("Saved: %s", par_out)

    # ── 5. Build normalised het rate (chrX / chr22) ───────────────────────────
    combined_df = None
    if not chr22_df.empty:
        log.info("\n--- Computing normalised het rate (chrX / chr22) ---")
        chrx_sub = chrx_df[["patient_id", "reference", "sex", "superpopulation",
                              "het_rate", "var_rate", "call_rate", "hom_rate",
                              "n_variants", "n_het", "n_called"]].copy()
        chrx_sub.columns = ["patient_id", "reference", "sex", "superpopulation",
                             "chrX_het_rate", "chrX_var_rate", "chrX_call_rate",
                             "chrX_hom_rate", "chrX_n_variants", "chrX_n_het",
                             "chrX_n_called"]

        chr22_sub = chr22_df[["patient_id", "reference", "het_rate", "var_rate",
                               "call_rate"]].copy()
        chr22_sub.columns = ["patient_id", "reference", "chr22_het_rate",
                              "chr22_var_rate", "chr22_call_rate"]

        combined_df = chrx_sub.merge(chr22_sub, on=["patient_id", "reference"], how="inner")
        combined_df["norm_het_rate"] = np.where(
            combined_df["chr22_het_rate"] > 0,
            combined_df["chrX_het_rate"] / combined_df["chr22_het_rate"],
            np.nan,
        )
        norm_out = out / "per_sample_normalised.csv"
        combined_df.to_csv(norm_out, index=False)
        log.info("Saved: %s", norm_out)

    # ── 6. Group comparisons ──────────────────────────────────────────────────
    log.info("\n--- Group comparisons ---")
    rows = []
    metrics = ["het_rate", "call_rate", "var_rate", "hom_rate"]

    for metric in metrics:
        for ref in ["GRCh38", "T2T"]:
            sub = chrx_df[chrx_df["reference"] == ref]
            f_vals = sub[sub["sex"] == "F"][metric].dropna().values
            m_vals = sub[sub["sex"] == "M"][metric].dropna().values
            rows.append(make_group_row(
                "F vs M on chrX", ref, "X", "chrX_whole", metric,
                "Female", "Male", f_vals, m_vals
            ))

    # GRCh38 vs T2T within sex (paired Wilcoxon)
    chrx_wide = chrx_df.pivot_table(
        index=["patient_id", "sex", "superpopulation"],
        columns="reference",
        values=metrics,
    ).reset_index()
    chrx_wide.columns = ["_".join(c).strip("_") for c in chrx_wide.columns.values]

    for metric in metrics:
        g38_col = f"{metric}_GRCh38"
        t2t_col = f"{metric}_T2T"
        if g38_col not in chrx_wide.columns or t2t_col not in chrx_wide.columns:
            continue
        for sex in ["F", "M"]:
            sub = chrx_wide[chrx_wide["sex"] == sex]
            a = sub[t2t_col].values
            b = sub[g38_col].values
            wc = wilcoxon_paired(a - np.nan, b - np.nan)  # diff = T2T − GRCh38
            diff = a - b
            valid = ~np.isnan(diff)
            wc = wilcoxon_paired(a[valid], b[valid]) if valid.any() else {"stat": np.nan, "p": np.nan, "n": 0}
            lo, hi = ci95_bootstrap(diff[valid])
            rows.append({
                "label":      f"T2T vs GRCh38 on chrX ({sex})",
                "reference":  "T2T_vs_GRCh38",
                "chromosome": "X",
                "region":     "chrX_whole",
                "metric":     metric,
                "group_a":    "T2T",
                "group_b":    "GRCh38",
                "mean_a":     float(np.nanmean(a)),
                "sd_a":       float(np.nanstd(a, ddof=1)),
                "ci95_lo_a":  np.nan,
                "ci95_hi_a":  np.nan,
                "n_a":        int(valid.sum()),
                "mean_b":     float(np.nanmean(b)),
                "sd_b":       float(np.nanstd(b, ddof=1)),
                "ci95_lo_b":  np.nan,
                "ci95_hi_b":  np.nan,
                "n_b":        int(valid.sum()),
                "mean_diff":  float(np.nanmean(diff)),
                "ci95_lo_diff": lo,
                "ci95_hi_diff": hi,
                "wilcoxon_stat": wc["stat"],
                "wilcoxon_p":    wc["p"],
                "mw_u":       np.nan,
                "mw_p":       np.nan,
                "effect_r":   np.nan,
                "cohens_d":   cohens_d(a[valid], b[valid]),
                "sex":        sex,
            })

    if not par_df.empty:
        for par_name in par_df["par_region"].unique():
            par_sub = par_df[par_df["par_region"] == par_name]
            for ref in ["GRCh38", "T2T"]:
                sub = par_sub[par_sub["reference"] == ref]
                f_vals = sub[sub["sex"] == "F"]["het_rate"].dropna().values
                m_vals = sub[sub["sex"] == "M"]["het_rate"].dropna().values
                rows.append(make_group_row(
                    f"F vs M on {par_name}", ref, "X",
                    PAR[ref].get(par_name, par_name), "het_rate",
                    "Female", "Male", f_vals, m_vals
                ))

    if combined_df is not None:
        for ref in ["GRCh38", "T2T"]:
            sub = combined_df[combined_df["reference"] == ref]
            f_vals = sub[sub["sex"] == "F"]["norm_het_rate"].dropna().values
            m_vals = sub[sub["sex"] == "M"]["norm_het_rate"].dropna().values
            rows.append(make_group_row(
                "F vs M: norm_het_rate (chrX/chr22)", ref, "X", "normalised", "norm_het_rate",
                "Female", "Male", f_vals, m_vals
            ))

    cmp_df = pd.DataFrame(rows)
    # Bonferroni correction on mw_p (primary tests only)
    mask = cmp_df["mw_p"].notna()
    n_tests = mask.sum()
    if n_tests > 0:
        cmp_df.loc[mask, "mw_p_bonferroni"] = (cmp_df.loc[mask, "mw_p"] * n_tests).clip(upper=1.0)
    cmp_out = out / "group_comparisons.csv"
    cmp_df.to_csv(cmp_out, index=False)
    log.info("Saved: %s", cmp_out)

    # ── 7. Superpopulation summary ────────────────────────────────────────────
    log.info("\n--- Superpopulation-stratified summary ---")
    sp_rows = []
    for ref in ["GRCh38", "T2T"]:
        for sex in ["F", "M"]:
            for sp in sorted(chrx_df["superpopulation"].dropna().unique()):
                sub = chrx_df[
                    (chrx_df["reference"] == ref)
                    & (chrx_df["sex"] == sex)
                    & (chrx_df["superpopulation"] == sp)
                ]
                if sub.empty:
                    continue
                for metric in metrics:
                    vals = sub[metric].dropna()
                    sp_rows.append({
                        "reference":      ref,
                        "sex":            sex,
                        "superpopulation": sp,
                        "metric":         metric,
                        "n":              len(vals),
                        "mean":           float(vals.mean()),
                        "sd":             float(vals.std(ddof=1)) if len(vals) > 1 else np.nan,
                        "median":         float(vals.median()),
                        "q25":            float(vals.quantile(0.25)),
                        "q75":            float(vals.quantile(0.75)),
                    })

    sp_df = pd.DataFrame(sp_rows)
    sp_out = out / "superpop_summary.csv"
    sp_df.to_csv(sp_out, index=False)
    log.info("Saved: %s", sp_out)

    # ── 8. Figures ────────────────────────────────────────────────────────────
    log.info("\n--- Generating figures ---")
    try:
        violin_by_sex_ref(
            chrx_df, "het_rate",
            "ChrX heterozygosity rate by sex and reference genome",
            out / "fig_chrX_het_rate.png",
        )
        violin_by_sex_ref(
            chrx_df, "call_rate",
            "ChrX genotype call rate by sex and reference genome",
            out / "fig_chrX_call_rate.png",
        )
        violin_by_sex_ref(
            chrx_df, "hom_rate",
            "ChrX homozygous fraction among variant sites (sex × reference)",
            out / "fig_chrX_hom_rate.png",
        )
        bar_superpop(
            chrx_df, "het_rate",
            "ChrX het rate by superpopulation, sex, reference",
            out / "fig_chrX_het_by_superpop.png",
        )
        if combined_df is not None:
            scatter_norm_het(combined_df, out / "fig_chrX_vs_chr22_het.png")
            violin_by_sex_ref(
                combined_df.rename(columns={"chrX_het_rate": "het_rate",
                                            "chrX_call_rate": "call_rate",
                                            "reference": "reference",
                                            "sex": "sex"}),
                "norm_het_rate" if "norm_het_rate" in combined_df.columns else "chrX_het_rate",
                "Normalised chrX het rate (chrX / chr22) by sex × reference",
                out / "fig_norm_het_rate.png",
            )
    except Exception as exc:
        log.warning("Figure generation failed (non-fatal): %s", exc)

    # ── 9. JSON summary ───────────────────────────────────────────────────────
    log.info("\n--- Writing JSON summary ---")

    def group_stats_json(df, ref, sex, metric):
        sub = df[(df["reference"] == ref) & (df["sex"] == sex)][metric].dropna()
        lo, hi = ci95_bootstrap(sub.values)
        return {
            "n": int(len(sub)),
            "mean": float(sub.mean()),
            "sd":   float(sub.std(ddof=1)) if len(sub) > 1 else None,
            "median": float(sub.median()),
            "ci95": [lo, hi],
        }

    # Key comparisons for manuscript text
    key = {}
    for ref in ["GRCh38", "T2T"]:
        key[ref] = {}
        for metric in ["het_rate", "call_rate", "hom_rate"]:
            key[ref][metric] = {
                "female_chrX": group_stats_json(chrx_df, ref, "F", metric),
                "male_chrX":   group_stats_json(chrx_df, ref, "M", metric),
            }

    # Delta het (T2T − GRCh38) for females
    wide_f = chrx_wide[chrx_wide["sex"] == "F"]
    for metric in ["het_rate", "call_rate"]:
        g38c = f"{metric}_GRCh38"
        t2tc = f"{metric}_T2T"
        if g38c in wide_f.columns and t2tc in wide_f.columns:
            diff = (wide_f[t2tc] - wide_f[g38c]).dropna()
            lo, hi = ci95_bootstrap(diff.values)
            wc = wilcoxon_paired(wide_f[t2tc].dropna().values,
                                 wide_f[g38c].dropna().values)
            key[f"delta_female_{metric}_T2T_minus_GRCh38"] = {
                "mean_delta": float(diff.mean()),
                "sd_delta":   float(diff.std(ddof=1)),
                "ci95":       [lo, hi],
                "wilcoxon_p": wc["p"],
                "n":          wc["n"],
            }

    # Delta het (T2T − GRCh38) for males
    wide_m = chrx_wide[chrx_wide["sex"] == "M"]
    for metric in ["het_rate"]:
        g38c = f"{metric}_GRCh38"
        t2tc = f"{metric}_T2T"
        if g38c in wide_m.columns and t2tc in wide_m.columns:
            diff = (wide_m[t2tc] - wide_m[g38c]).dropna()
            lo, hi = ci95_bootstrap(diff.values)
            wc = wilcoxon_paired(wide_m[t2tc].dropna().values,
                                 wide_m[g38c].dropna().values)
            key[f"delta_male_{metric}_T2T_minus_GRCh38"] = {
                "mean_delta": float(diff.mean()),
                "sd_delta":   float(diff.std(ddof=1)),
                "ci95":       [lo, hi],
                "wilcoxon_p": wc["p"],
                "n":          wc["n"],
            }

    summary = {
        "run_timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "test_mode":     args.test_mode,
        "n_samples":     int(meta.shape[0]),
        "n_female":      int((meta["sex"] == "F").sum()),
        "n_male":        int((meta["sex"] == "M").sum()),
        "chrX_n_records": int(len(chrx_df)),
        "key_statistics": key,
        "outputs": {
            "per_sample_chrX":      str(chrx_out),
            "group_comparisons":    str(cmp_out),
            "superpop_summary":     str(sp_out),
        },
    }

    summary_path = out / "summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log.info("Saved: %s", summary_path)

    # ── 10. Console report ────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("VARIANT-LEVEL BIAS ANALYSIS — RESULTS")
    print("=" * 70)

    for ref in ["GRCh38", "T2T"]:
        print(f"\n[{ref}] chrX het rate:")
        for sex, label in [("F", "Females"), ("M", "Males")]:
            sub = chrx_df[(chrx_df["reference"] == ref) & (chrx_df["sex"] == sex)]["het_rate"].dropna()
            if len(sub) > 0:
                print(f"  {label:8s}: mean={sub.mean():.4f}  sd={sub.std(ddof=1):.4f}  n={len(sub)}")

    print("\n[Δ het rate (T2T − GRCh38)]")
    for sex, label in [("F", "Females"), ("M", "Males")]:
        g38c = "het_rate_GRCh38"
        t2tc = "het_rate_T2T"
        if g38c in chrx_wide.columns and t2tc in chrx_wide.columns:
            sub = chrx_wide[chrx_wide["sex"] == sex]
            diff = (sub[t2tc] - sub[g38c]).dropna()
            wc = wilcoxon_paired(sub[t2tc].dropna().values, sub[g38c].dropna().values)
            lo, hi = ci95_bootstrap(diff.values)
            print(f"  {label:8s}: Δ={diff.mean():+.5f}  95%CI=[{lo:+.5f}, {hi:+.5f}]  "
                  f"Wilcoxon p={wc['p']:.4g}  n={wc['n']}")

    print(f"\nAll outputs in: {out}")
    print("=" * 70)


if __name__ == "__main__":
    main()
