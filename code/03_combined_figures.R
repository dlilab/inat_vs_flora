library(tidyverse)
library(patchwork)
library(cowplot)

dir.create("figures/combined", recursive = TRUE, showWarnings = FALSE)

# ── Load saved data ───────────────────────────────────────────────────────────
pheno_all_wf    <- readRDS("data_output/pheno_all_wildflower.rds")
est_pheno_df_wf <- readRDS("data_output/est_pheno_df_wildflower.rds")
range_area_wf   <- readRDS("data_output/range_area_wildflower.rds")
pgls_annot_wf   <- readRDS("data_output/pgls_range_annot_wildflower.rds")

pheno_all_cf    <- readRDS("data_output/pheno_all_calflora.rds")
est_pheno_df_cf <- readRDS("data_output/est_pheno_df_calflora.rds")
range_area_cf   <- readRDS("data_output/range_area_calflora.rds")
pgls_annot_cf   <- readRDS("data_output/pgls_range_annot_calflora.rds")

# ── Helper functions ──────────────────────────────────────────────────────────

make_overview_panel <- function(pheno_all, est_pheno_df, n_sp = 15, seed = 7, ncols = 5) {
  set.seed(seed)
  fig_dat <- slice_sample(pheno_all, n = n_sp)

  d_fig <- filter(est_pheno_df, scientific_name %in% fig_dat$scientific_name) |>
    arrange(scientific_name, peak_doy) |>
    group_by(scientific_name) |>
    mutate(idx = row_number()) |>
    ungroup() |>
    mutate(scientific_name = fct_reorder(scientific_name, onset_doy))

  fig_dat <- mutate(fig_dat,
    scientific_name = factor(scientific_name, levels = levels(d_fig$scientific_name))
  )

  ggplot(d_fig) +
    geom_rect(
      data = fig_dat,
      aes(xmin = onset_doy_flora, xmax = offset_doy_flora, ymin = -Inf, ymax = Inf),
      fill = "lightblue", alpha = 1
    ) +
    geom_segment(
      aes(x = onset_doy, y = idx, xend = offset_doy, yend = idx),
      linewidth = 0.7, alpha = 0.7
    ) +
    geom_segment(
      data = fig_dat,
      aes(x = onset_doy_inat, y = Inf, xend = offset_doy_inat, yend = Inf),
      linewidth = 2, color = "red"
    ) +
    facet_wrap(~scientific_name, scales = "free", ncol = ncols) +
    labs(x = "Day of year", y = "") +
    theme_bw() +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text   = element_text(face = "italic", size = 8)
    )
}

make_lollipop_panel <- function(pheno_all) {
  sp_order_peak <- pheno_all |>
    mutate(Peak = peak_doy_inat - peak_doy_flora) |>
    arrange(Peak) |>
    pull(scientific_name)

  pheno_all |>
    mutate(
      Onset    = onset_doy_inat  - onset_doy_flora,
      Peak     = peak_doy_inat   - peak_doy_flora,
      Offset   = offset_doy_inat - offset_doy_flora,
      Duration = dur_inat        - dur_flora,
      scientific_name = factor(scientific_name, levels = sp_order_peak)
    ) |>
    select(scientific_name, Onset, Peak, Offset, Duration) |>
    pivot_longer(-scientific_name, names_to = "variable") |>
    mutate(variable = factor(variable, levels = c("Onset", "Peak", "Offset", "Duration"))) |>
    ggplot(aes(x = value, y = scientific_name)) +
    geom_point(size = 0.8) +
    geom_vline(xintercept = 0, linetype = 1) +
    geom_segment(aes(x = 0, xend = value, yend = scientific_name), linewidth = 0.3) +
    facet_wrap(~variable, scales = "free_x", nrow = 1) +
    labs(x = "iNat-based − flora-based (days)", y = "Species") +
    theme_cowplot(font_size = 14) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

make_rangesize_panel <- function(range_area, pgls_annot) {
  pivot_longer(range_area,
    cols = c(Onset, Peak, Offset, Duration),
    names_to = "type", values_to = "delta"
  ) |>
    mutate(type = factor(type, levels = c("Onset", "Peak", "Offset", "Duration"))) |>
    left_join(select(pgls_annot, type, linetype_pgls), by = "type") |>
    ggplot(aes(x = area_km2_log10, y = delta)) +
    geom_point(alpha = 0.2, size = 0.8) +
    geom_smooth(aes(linetype = linetype_pgls), method = "lm", se = TRUE) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_linetype_identity() +
    geom_text(
      data = mutate(pgls_annot, type = factor(type, levels = c("Onset", "Peak", "Offset", "Duration"))),
      aes(x = -Inf, y = -Inf, label = label),
      hjust = -0.05, vjust = -0.5, size = 2.8, inherit.aes = FALSE
    ) +
    facet_wrap(~type, scales = "free_y", nrow = 1) +
    labs(
      x = expression(Log[10] ~ range ~ area ~ (km^2)),
      y = "iNat-based − flora-based (days)"
    ) +
    theme_cowplot(font_size = 14)
}

# ── Fig 1: combined species overview ─────────────────────────────────────────
pA <- make_overview_panel(pheno_all_wf, est_pheno_df_wf, n_sp = 10, seed = 10)
pB <- make_overview_panel(pheno_all_cf, est_pheno_df_cf, n_sp = 10, seed = 8)

fig1 <- pA / pB + plot_annotation(tag_levels = "A")
ggsave("figures/combined/fig_1_overview.pdf", plot = fig1, width = 11, height = 10)
ggsave("figures/combined/fig_1_overview.png", plot = fig1, width = 11, height = 10)
cat("Saved figures/combined/fig_1_overview\n")

# ── Fig 2: combined lollipop differences ─────────────────────────────────────
p2A <- make_lollipop_panel(pheno_all_wf)
p2B <- make_lollipop_panel(pheno_all_cf)

fig2 <- p2A / p2B + plot_annotation(tag_levels = "A")
ggsave("figures/combined/fig_2_differences.pdf", plot = fig2, width = 10, height = 10)
ggsave("figures/combined/fig_2_differences.png", plot = fig2, width = 10, height = 10)
cat("Saved figures/combined/fig_2_differences\n")

# ── Fig 4: combined range size ────────────────────────────────────────────────
p4A <- make_rangesize_panel(range_area_wf, pgls_annot_wf)
p4B <- make_rangesize_panel(range_area_cf, pgls_annot_cf)

fig4 <- p4A / p4B + plot_annotation(tag_levels = "A")
ggsave("figures/combined/fig_4_rangesize.pdf", plot = fig4, width = 10, height = 8)
ggsave("figures/combined/fig_4_rangesize.png", plot = fig4, width = 10, height = 8)
cat("Saved figures/combined/fig_4_rangesize\n")

cat("\nDone. Combined figures in figures/combined/\n")
