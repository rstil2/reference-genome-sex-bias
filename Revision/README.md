# Revision analyses (October 2026)

Scripts, outputs and figures added when the study was revised. They sit beside the original pipeline in `../Code` and read the stored tables in `../Results`. Run each script from this folder.

## Correction to the stored per-sample tables

`Results/variant_analysis/per_sample_chrX.csv` and `per_sample_chr22.csv` were written by a version of `Code/06_variant_level_bias_analysis.py` that tested for `nRefHom` before `nNonRefHom`. Because "nrefhom" is a substring of "nnonrefhom", both bcftools PSC columns were mapped to `n_ref_hom`. In the stored files:

- `n_ref_hom` holds nNonRefHom, and `n_hom_alt` is 0 in every row. The reference-homozygous count was lost.
- `het_rate` in women is nHet / (nHet + nNonRefHom), the heterozygous fraction of non-reference genotypes. In men the denominator includes haploid calls, so the male and female values are not comparable.
- `call_rate` leaves out diploid reference-homozygous genotypes.

`Code/06_variant_level_bias_analysis.py` is now fixed. The stored tables were **not** regenerated, because the source VCFs are not in the repository (the T2T chrX file alone is about 54 GB). The Panel A and B results were computed from the stored values and are unchanged.

## Genotype coding in Panels C and D

`Code/02b_extract_variant_features.py` codes a genotype as the sum of the two allele indices. That is the alternate-allele count only at biallelic sites. The GRCh38 feature matrix contains 0, 1 and 2. The T2T matrix also contains 3 and 4 at multiallelic records (about 157,000 of 36.8 million cells). `revision_panelC_dosage_census.csv` counts them. The first 500 chrX SNPs per callset fall in the pseudoautosomal region (`revision_panelC_chrX_positions.json`).

## Files

| File | What it does |
|:---|:---|
| `revision_stratified_contrasts.py` / `.csv` | Paired T2T − GRCh38 change in the female chrX heterozygous fraction by superpopulation, bootstrap CIs, Kruskal–Wallis; chr22 contrasts |
| `revision_count_decomposition.py` / `.csv` | Paired changes in nHet, nNonRefHom, nHapRef, nHapAlt and nMissing, chrX and chr22, by sex and superpopulation |
| `revision_panelA_relatedness.py` / `.csv` | Panel A refit on the 2,504 phase 3 panel, with family-blocked splits, and with one person per family |
| `revision_family_holm.py` / `.csv` | Holm–Bonferroni across all 20 panel × model tests |
| `revision_panelC_dosage_sensitivity.py`, `revision_panelC_dosage_census.csv`, `revision_panelC_capped_delta_delta*.{json,csv}`, `revision_panelC_chrX_positions.json` | Dosage census and Panel C refit with T2T dosages capped at 2 |
| `revision_figures.py`, `Figures/Revised_Fig1…6` | Figures for the revised manuscript |
| `data/` | 1000 Genomes `20130606_g1k.ped` and `integrated_call_samples.panel` (public) |

`revision_panelC_dosage_sensitivity.py` needs `Results/features/*_genome_features_aligned.csv`, which are not stored here (see `.gitignore`); regenerate them with the pipeline in `../Code`. The original analyses ran under Python 3.12.13, scikit-learn 1.8.0, XGBoost 3.0.3 and SciPy 1.16.0. The revision analyses ran under Python 3.13.1, scikit-learn 1.8.0, XGBoost 3.3.0 and SciPy 1.14.1. Bootstrap intervals use seed 20260222; last-digit differences in p-values across SciPy versions are expected.
