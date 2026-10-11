# Reference callsets and chrX sex disparity

Reanalysis of 3,202 genomes from the 1000 Genomes Project, comparing the published GRCh38 and T2T-CHM13 callsets. The question is whether genotype summaries on chromosome X, and the male–female gap in classifiers trained on those summaries, differ between the two callsets.

**Paper.** *The Reference Shapes the Signal: chrX Genotype Summaries and Classifier Sex Disparity Between GRCh38 and T2T-CHM13 Callsets.* Accepted at *Genomics & Informatics*, 10 October 2026.

The two callsets were produced by different alignment, calling, phasing, and filtering workflows. Differences belong to the reference and its pipeline together. Women are classified better than men under both callsets. Whether that gap widens or narrows on T2T depends on the model. The gap is a difference in classifier performance, not a biological sex difference.

This repository is Python. The original analysis used Python 3.12.13, scikit-learn 1.8.0, XGBoost 3.0.3, and SciPy 1.16.0. The revision analyses used Python 3.13.1, scikit-learn 1.8.0, XGBoost 3.3.0, and SciPy 1.14.1. Random seeds are 20260222.

## What a reader should open first

| Path | What it is |
|---|---|
| `Revision/` | Analyses and figures that match the accepted paper. Start here. |
| `Results/metrics/downstream_v2/` | Classifier metrics and ΔΔ for the four panels |
| `Results/variant_analysis/per_sample_chrX.csv` | Stored per-sample chrX counts. Read the caveat below before using `het_rate` or `call_rate` |
| `Code/` | Pipeline that built those results from the VCFs |
| `Manuscript/` | An earlier draft. It does not match the accepted paper |

## Repository map

```
Code/           download, genotype counts, SNP matrices, model training, ΔΔ
Data/Metadata/  sample population labels (public). VCFs are not stored
Results/        metrics, splits, per-sample tables, environment record
Revision/       post-review analyses, corrected figures, pedigree file
Figures/        figures from the earlier draft, not the accepted paper
Manuscript/     earlier draft
```

## Pipeline

The source VCFs are public and large (the T2T chrX file is about 54 GB). They are listed in `.gitignore` and are not in this repository. Feature matrices under `Results/features/` are also omitted. The metrics, per-sample tables, and revision outputs are included, so the reported numbers can be checked without rebuilding the VCFs.

1. **Download.** `Code/download_grch38_vcfs.sh` and `Code/download_t2t_chm13_vcfs.sh` fetch the public callsets. Several other `download_*.sh` scripts are earlier attempts at the same downloads.
2. **Per-sample counts.** `Code/06_variant_level_bias_analysis.py` runs `bcftools stats` and writes `Results/variant_analysis/`.
3. **SNP matrices.** `Code/02b_extract_variant_features.py` takes the first 500 SNPs per chromosome in VCF order (minor-allele frequency ≥ 0.01, call rate ≥ 0.95).
4. **Classifiers.** `Code/10_downstream_panel_battery.py` trains five models on four panels over 20 stratified 70/15/15 splits. `Code/03_train_models.py` fits one model. `Code/04_analyze_results.py` computes ΔΔ.
5. **Revision.** From `Revision/`, the scripts `revision_count_decomposition.py`, `revision_stratified_contrasts.py`, `revision_panelA_relatedness.py`, `revision_family_holm.py`, `revision_panelC_dosage_sensitivity.py`, and `revision_figures.py` rebuild the tables and figures used in the accepted paper. See `Revision/README.md`.

To rerun the panel battery once the feature files and per-sample tables are in place:

```bash
python3 Code/10_downstream_panel_battery.py --study-root . --seed 20260222 --max-repeats 20
```

`Code/05_run_redesigned_downstream.sh` is the same run, with a log file under `Results/logs/`. It uses `python3` on `PATH`, or `.venv/bin/python` if that environment exists in the repository root.

Dependencies for a fresh environment are in `requirements.txt` (versions from the original run).

## Panels

| Panel | Task | Features |
|---|---|---|
| A | Superpopulation | Six stored chrX summaries |
| B | Sex | The same six summaries |
| C | Superpopulation | 500 SNPs × 23 chromosomes, VCF order |
| D | Superpopulation | The same rule, autosomes only |

ΔΔ is the male–female gap in balanced accuracy on GRCh38 minus the same gap on T2T. In Panel A the gap is negative for every model, because women are classified better than men. A positive ΔΔ means the female advantage is larger on T2T.

## Known limitations of the stored files

These are the same limitations reported in the paper.

**Per-sample tables.** An earlier version of `06_variant_level_bias_analysis.py` matched the column name `nRefHom` before `nNonRefHom`. The substring `nrefhom` occurs inside `nnonrefhom`, so both bcftools columns were written to `n_ref_hom`. In the stored CSVs:

- `n_ref_hom` holds the non-reference homozygous count, and `n_hom_alt` is 0 in every row
- in women, `het_rate` is nHet / (nHet + nNonRefHom), not heterozygosity among all called genotypes
- the male and female values of that ratio are not comparable, because haploid calls enter the male denominator
- `call_rate` omits diploid reference-homozygous genotypes

The script is fixed. The stored tables were not regenerated. Panel A and Panel B were fit on these stored columns.

**Panels C and D.** Genotypes were coded as the sum of the two allele indices. That equals the alternate-allele count only at biallelic sites. The T2T matrix contains dosages 3 and 4 (about 157,000 of 36.8 million cells). Capping those values at 2 and refitting Panel C did not change its conclusion (`Revision/revision_panelC_capped_delta_delta.json`). The chrX SNPs in Panel C lie in the pseudoautosomal region.

**Relatedness.** The 3,202 people include families. `Revision/revision_panelA_relatedness.py` refits Panel A on the phase 3 panel, with family-blocked splits, and with one person per family. The pedigree file is `Revision/data/20130606_g1k.ped`.

## Earlier draft

`Manuscript/Manuscript.md` and `Figures/Figure1`–`Figure5` predate peer review. They describe the heterozygosity ratio incorrectly and they do not match the accepted results. Use `Revision/` for the numbers and figures in the paper.
