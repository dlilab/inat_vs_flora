library(tidyverse)
library(sf)
library(phylolm)
library(broom)

# Is the negative range-size ~ discrepancy relationship (Reviewer 1, comment 2) an
# artifact of observation effort? Wide-ranging species have more iNaturalist records,
# so (i) their observation-based estimates are more reliable, and (ii) convex-hull
# area grows with n. We (a) add log10(n observations) and log10(n grid cells) as
# covariates in the range-size PGLS, and (b) refit on species with high observation
# effort only (n_obs >= median).

run_scale <- function(scale, pheno_path, range_path, obs_path, est_path) {
  pheno_all <- readRDS(pheno_path)
  range_df <- readRDS(range_path)
  obs <- st_drop_geometry(readRDS(obs_path))
  est <- readRDS(est_path)

  effort <- obs |>
    count(scientific_name, name = "n_obs") |>
    left_join(count(distinct(est, scientific_name, id_cells), scientific_name, name = "n_cells"),
      by = "scientific_name"
    ) |>
    transmute(sp = scientific_name, log_nobs = log10(n_obs), log_ncells = log10(n_cells), n_obs)

  sp_list <- gsub(" ", "_", pheno_all$scientific_name)
  phy <- rtrees::get_tree(rtrees::sp_list_df(sp_list, taxon = "plant"), taxon = "plant")
  phy$root.edge <- 0

  d <- range_df |>
    select(sp, area_km2_log10, Onset, Peak, Offset, Duration) |>
    left_join(effort, by = "sp") |>
    drop_na(log_nobs) |>
    mutate(sp = gsub(" ", "_", sp)) |>
    filter(sp %in% phy$tip.label) |>
    as.data.frame() |>
    column_to_rownames("sp")

  cat(scale, ": n species =", nrow(d),
    "; cor(log range, log n_obs) =", round(cor(d$area_km2_log10, d$log_nobs), 2), "\n")

  fit <- function(dat, rhs, model_name) {
    map_dfr(c("Onset", "Peak", "Offset", "Duration"), \(ph) {
      dd <- dat[!is.na(dat[[ph]]), ]
      m <- phylolm(as.formula(paste(ph, "~", rhs)),
        data = dd, phy = drop.tip(phy, setdiff(phy$tip.label, rownames(dd))), model = "lambda"
      )
      class(m) <- c(class(m), "lm")
      tidy(m) |> mutate(phenostage = ph, model = model_name, n = nrow(dd))
    })
  }

  bind_rows(
    fit(d, "area_km2_log10", "range only"),
    fit(d, "area_km2_log10 + log_nobs", "range + log10(n obs)"),
    fit(d, "area_km2_log10 + log_ncells", "range + log10(n cells)"),
    fit(filter(d, n_obs >= median(n_obs)), "area_km2_log10", "range only, n obs >= median")
  ) |> mutate(scale = scale)
}

res <- bind_rows(
  run_scale("National", "data_output/pheno_all_wildflower.rds",
    "data_output/range_area_wildflower.rds", "data/dat_used_phenobase.rds",
    "data/est_pheno_df_phenobase.rds"),
  run_scale("Regional", "data_output/pheno_all_calflora.rds",
    "data_output/range_area_calflora.rds", "data/dat_used_calflora.rds",
    "data/est_pheno_df_calflora.rds")
)

range_effort_robustness <- res |>
  filter(term == "area_km2_log10") |>
  select(scale, model, phenostage, n, estimate, std.error, p.value) |>
  mutate(phenostage = factor(phenostage, c("Onset", "Peak", "Offset", "Duration"))) |>
  arrange(scale, model, phenostage)

saveRDS(range_effort_robustness, "data_output/range_effort_robustness.rds")
saveRDS(filter(res, term %in% c("log_nobs", "log_ncells")), "data_output/range_effort_covariate_terms.rds")
print(as.data.frame(mutate(range_effort_robustness, across(where(is.numeric), \(x) round(x, 3)))))
cat("\nEffort covariate terms:\n")
print(as.data.frame(mutate(filter(res, term %in% c("log_nobs", "log_ncells")) |>
  select(scale, model, phenostage, term, estimate, p.value), across(where(is.numeric), \(x) round(x, 3)))))
