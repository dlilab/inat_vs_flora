# Data and code: Flora-based vs. iNaturalist-derived flowering phenology

Data and R code to reproduce the statistical analyses and figures in
*Flora-based phenological information is not a reliable proxy of actual observations* (Li & Guralnick).

Only the data used directly in the statistical analyses and the analysis code are included.
Code for downloading, cleaning, and exploring the raw data (including estimating phenology from
the Phenobase/iNaturalist observations) is not included.

## Usage

Run from the root of this folder:

```r
source("run_all.R")   # or: Rscript run_all.R
```

Outputs are written to `data_output/` (statistics and intermediate objects) and `figures/`.
Required R packages: `tidyverse`, `lubridate`, `sf`, `phylolm`, `broom`, `phytools`, `rtrees`,
`units`, `cowplot`, `patchwork`.

## Code (`code/`)

| Script | Purpose | Main outputs |
|---|---|---|
| `01_compare_wildflower.R` | National comparison: iNaturalist-derived vs. Wildflower Center bloom periods (paired tests, RMSE, predictive R², phylogenetic signal, range-size and trait PGLS) | `data_output/stats_wildflower.RData`, `pheno_all_wildflower.rds`, `est_pheno_df_wildflower.rds`, `range_area_wildflower.rds`, figures in `figures/wildflower/` |
| `02_compare_calflora.R` | Regional comparison with the Jepson eFlora (via CalFlora) | `data_output/stats_calflora.RData`, `pheno_all_calflora.rds`, etc., figures in `figures/calflora/` |
| `03_combined_figures.R` | Combined national + regional main-text figures | `figures/combined/` |
| `04_flora_jitter_sensitivity.R` | Monte Carlo sensitivity of the results to the monthly resolution of flora data | `data_output/flora_jitter_sensitivity.rds` |
| `05_range_size_effort_robustness.R` | Range-size PGLS controlling for number of observations / grid cells | `data_output/range_effort_robustness.rds` |

## Data (`data/`)

| File | Description |
|---|---|
| `est_pheno_df_phenobase.rds` | Grid-cell × year iNaturalist-derived flowering onset, peak and offset (`phenesse` 5th/50th/95th percentiles; 25 km × 25 km Albers equal-area cells), national species set. Columns include `scientific_name`, `observed_year`, `id_cells`, `onset_doy`, `peak_doy`, `offset_doy`. |
| `dat_used_phenobase.rds` | Flowering observations (sf points; `scientific_name`, `observed_on`, `observed_year`, `id_cells`, `geometry`) used for the national estimates; used to compute range size (convex hull) and number of observations. |
| `wildflower_time_phenobase.csv` | Lady Bird Johnson Wildflower Center bloom periods (`bloom_time`, `start_month`, `end_month`) and traits (`life_form`, `growth_form`) for the species. |
| `est_pheno_df_calflora.rds`, `dat_used_calflora.rds` | Same as above for the regional (California, Oregon, Nevada, Arizona) comparison. |
| `calflora_phenobase.csv` | Jepson eFlora (via CalFlora) bloom periods (`start_month`, `end_month`) and `habit_info` text. |
| `wildflower_info_df.rds` | Wildflower Center trait table (`life_form`, `growth_form`) used to assign growth forms to CalFlora species. |

Observation data originate from Phenobase (https://phenobase.org; iNaturalist observations annotated
with PhenoVision). Bloom periods were obtained from the Lady Bird Johnson Wildflower Center
(https://www.wildflower.org/plants/) and the Jepson eFlora via CalFlora (https://ucjeps.berkeley.edu/eflora/).

## Notes

- Reproducibility: all estimates (means, RMSE, R², PGLS coefficients, Blomberg's K and Pagel's lambda) are
  reproduced exactly. The p-values for Blomberg's K come from randomization tests (`phytools::phylosig`)
  with no seed set, so they vary slightly between runs (e.g., 0.040-0.046 for flora-based onset K).
- Other random steps (figure subsampling, jitter sensitivity) use fixed seeds.
- Running `run_all.R` takes about 3 minutes.
