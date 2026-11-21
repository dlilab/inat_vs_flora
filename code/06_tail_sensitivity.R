library(tidyverse)
library(patchwork)
library(cowplot)

# Sensitivity of flora-vs-iNat discrepancies to the tail proportion used to aggregate
# grid-cell onset/offset estimates to the species level (Figure S1).
# Main analysis: onset = mean of the earliest 10% of grid cells, offset = mean of the
# latest 10%; here 5%, 10%, and 15% are compared. Peak (median of all grid cells) does
# not depend on the tail proportion. Run from the root of this folder after scripts 01 and 02.

props <- c(0.05, 0.10, 0.15)

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
r2 <- function(obs_inat, pred_flora) {
  ss_res <- sum((obs_inat - pred_flora)^2, na.rm = TRUE)
  ss_null <- sum((obs_inat - mean(obs_inat, na.rm = TRUE))^2, na.rm = TRUE)
  1 - ss_res / ss_null
}

tail_stats <- function(est, flora, p) {
  onset <- est |>
    group_by(scientific_name) |>
    slice_min(prop = p, order_by = onset_doy) |>
    summarise(onset_doy_inat = round(mean(onset_doy, na.rm = TRUE)), .groups = "drop")
  offset <- est |>
    group_by(scientific_name) |>
    slice_max(prop = p, order_by = offset_doy) |>
    summarise(offset_doy_inat = round(mean(offset_doy, na.rm = TRUE)), .groups = "drop")
  peak <- est |>
    group_by(scientific_name) |>
    summarise(peak_doy_inat = round(median(peak_doy, na.rm = TRUE)), .groups = "drop")

  d <- onset |>
    left_join(offset, by = "scientific_name") |>
    left_join(peak, by = "scientific_name") |>
    mutate(dur_inat = offset_doy_inat - onset_doy_inat) |>
    left_join(flora, by = "scientific_name")

  pairs <- list(
    Onset = c("onset_doy_inat", "onset_doy_flora"),
    Peak = c("peak_doy_inat", "peak_doy_flora"),
    Offset = c("offset_doy_inat", "offset_doy_flora"),
    Duration = c("dur_inat", "dur_flora")
  )
  imap_dfr(pairs, \(v, nm) {
    x <- d[[v[1]]]
    y <- d[[v[2]]]
    tibble(
      phenophase = nm, n = sum(complete.cases(x, y)),
      mean_diff = mean(x - y, na.rm = TRUE), RMSE = rmse(x, y), R2 = r2(x, y),
      wilcox_p = wilcox.test(x, y, paired = TRUE)$p.value
    )
  }) |> mutate(prop = p)
}

run_one <- function(scale, dataset) {
  est <- readRDS(paste0("data_output/est_pheno_df_", scale, ".rds"))
  flora <- readRDS(paste0("data_output/pheno_all_", scale, ".rds")) |>
    select(scientific_name, onset_doy_flora, offset_doy_flora, peak_doy_flora, dur_flora)
  map_dfr(props, \(p) tail_stats(est, flora, p)) |> mutate(dataset = dataset)
}

sensitivity_results <- bind_rows(
  run_one("wildflower", "National (Wildflower Center)"),
  run_one("calflora", "Regional (Jepson eFlora)")
) |>
  mutate(
    phenophase = factor(phenophase, levels = c("Onset", "Peak", "Offset", "Duration")),
    prop_label = factor(paste0(prop * 100, "% tail"), levels = paste0(props * 100, "% tail")),
    sig = case_when(wilcox_p < 0.001 ~ "***", wilcox_p < 0.01 ~ "**", wilcox_p < 0.05 ~ "*", TRUE ~ "ns")
  ) |>
  relocate(dataset, prop, phenophase)

print(as.data.frame(sensitivity_results), digits = 3)
saveRDS(sensitivity_results, "data_output/sensitivity_results.rds")

# ── Figure S1 ─────────────────────────────────────────────────────────────────
cols <- c("National (Wildflower Center)" = "#1b7837", "Regional (Jepson eFlora)" = "#762a83")
shps <- c("National (Wildflower Center)" = 16, "Regional (Jepson eFlora)" = 17)

panel <- function(yvar, ylab) {
  ggplot(sensitivity_results, aes(x = prop_label, y = .data[[yvar]],
                                  colour = dataset, shape = dataset, group = dataset)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    facet_wrap(~phenophase, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = cols) +
    scale_shape_manual(values = shps) +
    labs(x = "Tail proportion for onset/offset aggregation", y = ylab, colour = NULL, shape = NULL) +
    theme_cowplot(font_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey85", colour = NA),
          legend.position = "bottom")
}

fig_s1 <- panel("mean_diff", "iNat-based − flora-based (days)") /
  panel("R2", expression(Predictive ~ R^2 ~ "(NSE)")) +
  plot_annotation(tag_levels = "A")

dir.create("figures/sensitivity", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/sensitivity/fig_sensitivity.pdf", fig_s1, width = 11, height = 8, device = cairo_pdf)
ggsave("figures/sensitivity/fig_sensitivity.png", fig_s1, width = 11, height = 8, dpi = 300, bg = "white")
