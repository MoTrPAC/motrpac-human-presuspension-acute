# figures/adipose

Code for the adipose tissue manuscript. Each R script covers a specific analysis area; the `Files/` folder contains external reference data from the literature, with each file's source and usage documented in the script that loads it.

## Scripts

| Script | Description |
|--------|-------------|
| `precawg_adi_da.R` | DA results and volcano plots; cross-references published adipose exercise studies (Fabre et al. 2018, Ahn et al. 2025) for external validation |
| `precawg_adi_pathway.R` | Pathway enrichment analysis; loads PTM-SEA results |
| `precawg_adi_secretome.R` | Exercise-induced secretome; uses HPA secretome annotations |
| `precawg_adi_clinic_decon_wgcna.R` | Clinical × omics, cell type deconvolution, and WGCNA integration |
| `precawg_adi_bcka.R` | Branched-chain keto acid analysis; cross-references Wlejko et al. phosphosite tables |
| `precawg_adi_celltype.R` | Cell type analysis |
| `precawg_adi_gdcat.R` | GD-CAT analysis |
| `precawg_adi_ssec.R` | Sex-stratified exercise comparisons |
| `baseline_sex_analyses_figures.R` | Baseline sex difference figures; reads from `Files/sex_da/` and the all-omes enrichment file |
| `figure_2_cmeans.Rmd` | C-means clustering figure |

The `sex_difference_analysis/` subfolder contains sensitivity analyses for baseline sex differences:
- `baseline_sex_differences.Rmd`
- `sex_diff_helpers.R`
- `sex_differences_sens_analysis.Rmd`

## Files/

External reference data from the literature. Each file's source and usage is documented in the R script that loads it.

| File | Source | Used in |
|------|--------|---------|
| `Ahn_2025_acute_ex.csv` | Ahn et al. 2025 | `precawg_adi_da.R` |
| `Fabre_2018_untrained.xlsx` | Fabre et al. 2018 | `precawg_adi_da.R` |
| `wlejko_sup2.xlsx` | Wlejko et al. | `precawg_adi_bcka.R` |
| `sa_location_Secreted.tsv` | Human Protein Atlas (HPA) | `precawg_adi_secretome.R` |
| `n3_ptm-sea-results-combined.gct` | PTM-SEA output - coming to `MotrpacHumanPreSuspensionAnalysis` soon | `precawg_adi_pathway.R` |
| `human-precovid-sed-adu_t11-adipose_all-omes_enrichment_sex_diff_baseline_simple_v1.31.txt` | Pre-computed sex difference enrichments | `baseline_sex_analyses_figures.R` |
| `sex_da/` | Per-platform sex DA results (proteomics, phosphoproteomics, metabolomics, RNA-seq) | `baseline_sex_analyses_figures.R` |

## Data Access

Like all analyses in this repository, most data is loaded via `MotrpacHumanPreSuspensionAnalysis` and `MotrpacHumanPreSuspensionData`. Individual-level data requires a formal data access request to the MoTrPAC consortium.
