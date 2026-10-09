library(tidyverse)

# Sensitivity of flora-vs-iNat discrepancies to the coarse monthly resolution of
# flora records (Reviewer 1, comment 1). Flora onset/offset/peak are redrawn
# uniformly within +/- w days of the monthly-boundary / midpoint values used in the
# main analysis; the comparison statistics are recomputed for each draw.

set.seed(2026)
n_rep <- 1000
windows <- c(0, 5, 10, 15)

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
r2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

jitter_stats <- function(d, w) {
  n <- nrow(d)
  on_j <- runif(n, -w, w)
  off_j <- runif(n, -w, w)
  clip <- function(x) pmin(pmax(x, 1), 366)
  onset <- clip(d$onset_doy_flora + on_j)
  offset <- clip(d$offset_doy_flora + off_j)
  peak <- clip(d$peak_doy_flora + runif(n, -w, w))
  dur <- offset - onset

  pairs <- list(
    Onset = list(d$onset_doy_inat, onset),
    Peak = list(d$peak_doy_inat, peak),
    Offset = list(d$offset_doy_inat, offset),
    Duration = list(d$dur_inat, dur)
  )
  imap_dfr(pairs, \(p, nm) {
    # use the same complete pairs for every statistic
    ok <- complete.cases(p[[1]], p[[2]])
    x <- p[[1]][ok]
    y <- p[[2]][ok]
    tibble(
      phenophase = nm,
      n = sum(ok),
      mean_diff = mean(x - y),
      RMSE = rmse(x, y),
      R2 = r2(x, y),
      p_wilcox = suppressWarnings(wilcox.test(x, y, paired = TRUE)$p.value)
    )
  })
}

run_one <- function(path, scale) {
  d <- readRDS(path) |> drop_na(onset_doy_flora, offset_doy_flora)
  map_dfr(windows, \(w) {
    reps <- if (w == 0) 1 else n_rep
    map_dfr(seq_len(reps), \(i) jitter_stats(d, w)) |>
      mutate(w = w)
  }) |>
    mutate(scale = scale)
}

draws <- bind_rows(
  run_one("data_output/pheno_all_wildflower.rds", "National"),
  run_one("data_output/pheno_all_calflora.rds", "Regional")
)

jitter_summary <- draws |>
  group_by(scale, w, phenophase) |>
  summarise(
    across(c(mean_diff, RMSE, R2), list(
      med = median,
      lo = \(x) quantile(x, 0.025),
      hi = \(x) quantile(x, 0.975)
    )),
    prop_sig = mean(p_wilcox < 0.05),
    .groups = "drop"
  ) |>
  mutate(phenophase = factor(phenophase, c("Onset", "Peak", "Offset", "Duration"))) |>
  arrange(scale, w, phenophase)

saveRDS(jitter_summary, "data_output/flora_jitter_sensitivity.rds")
print(jitter_summary, n = Inf, width = Inf)
