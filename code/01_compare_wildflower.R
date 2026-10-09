library(tidyverse)
library(lubridate)
library(sf)
library(phylolm)
library(broom)

dir.create("figures/wildflower", recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
est_pheno_df <- readRDS("data/est_pheno_df_phenobase.rds")

wildflower_raw <- read_csv("data/wildflower_time_phenobase.csv", show_col_types = FALSE)
# 7 exact duplicates — keep one row per species
wildflower <- distinct(wildflower_raw, sp, .keep_all = TRUE)

# ── 2. Filter: ≥10 grid cells per species ────────────────────────────────────
sp_keep <- est_pheno_df |>
  group_by(scientific_name) |>
  tally() |>
  filter(n >= 10) |>
  pull(scientific_name)

est_pheno_df <- filter(est_pheno_df, scientific_name %in% sp_keep)
cat("Species with ≥10 grid cells:", length(sp_keep), "\n")

# ── 2b. Grid-cell outlier filtering ──────────────────────────────────────────
# Within each species, flag cells with onset too early, offset too late, or
# duration too long using Tukey's k=2 IQR fences (catches extreme outliers only).
n_before <- nrow(est_pheno_df)

est_pheno_df <- est_pheno_df |>
  mutate(dur_doy = offset_doy - onset_doy) |>
  group_by(scientific_name) |>
  mutate(
    onset_lo  = quantile(onset_doy, 0.25, na.rm = TRUE) - 2 * IQR(onset_doy, na.rm = TRUE),
    offset_hi = quantile(offset_doy, 0.75, na.rm = TRUE) + 2 * IQR(offset_doy, na.rm = TRUE),
    dur_hi    = quantile(dur_doy, 0.75, na.rm = TRUE) + 2 * IQR(dur_doy, na.rm = TRUE)
  ) |>
  filter(
    onset_doy >= onset_lo,
    offset_doy <= offset_hi,
    dur_doy <= dur_hi,
    dur_doy >= 0 # drop impossible cells (offset before onset)
  ) |>
  select(-onset_lo, -offset_hi, -dur_hi, -dur_doy) |>
  ungroup()

cat("Grid-cell rows removed as outliers:", n_before - nrow(est_pheno_df), "\n")

# ── 3. Prepare flora data ─────────────────────────────────────────────────────
# start_month / end_month already numeric in wildflower CSV
flora <- wildflower |>
  filter(!is.na(start_month), !is.na(end_month)) |>
  mutate(
    start_doy = yday(make_date(2020, start_month, 1)),
    end_doy = yday(make_date(2020, end_month + 1, 1) - 1),
    peak_doy_flora = start_doy + (end_doy - start_doy) / 2,
    dur_flora = end_doy - start_doy,

    # bloom_group: spring = onset before June, summer = June onward
    bloom_group = if_else(start_month < 6, "Spring", "Summer"),

    # growth_form_3: consolidate into Woody_perennial / Herb_perennial / Herb_annual/biennial
    growth_form_3 = case_when(
      str_detect(growth_form, "Tree|Shrub|Subshrub") ~ "Woody_perennial",
      str_detect(life_form, "Annual|Biennial") ~ "Herb_annual/biennial",
      TRUE ~ "Herb_perennial"
    )
  )

# ── 4. Common species ─────────────────────────────────────────────────────────
common_sp <- intersect(unique(est_pheno_df$scientific_name), flora$sp)
cat("Common species:", length(common_sp), "\n")

est_pheno_df <- filter(est_pheno_df, scientific_name %in% common_sp)

# ── 5. Species-level aggregation  ────────────────
# onset  = mean of earliest 10% of grid-cell onset values
# peak   = median across all grid cells
# offset = mean of latest  10% of grid-cell offset values

onset_inat <- est_pheno_df |>
  group_by(scientific_name) |>
  slice_min(prop = 0.1, order_by = onset_doy) |>
  summarise(onset_doy_inat = round(mean(onset_doy, na.rm = TRUE)), .groups = "drop")

offset_inat <- est_pheno_df |>
  group_by(scientific_name) |>
  slice_max(prop = 0.1, order_by = offset_doy) |>
  summarise(offset_doy_inat = round(mean(offset_doy, na.rm = TRUE)), .groups = "drop")

peak_inat <- est_pheno_df |>
  group_by(scientific_name) |>
  summarise(peak_doy_inat = round(median(peak_doy, na.rm = TRUE)), .groups = "drop")

flora_join <- flora |>
  filter(sp %in% common_sp) |>
  select(
    scientific_name = sp, onset_doy_flora = start_doy, offset_doy_flora = end_doy,
    peak_doy_flora, dur_flora, bloom_group, growth_form_3
  )

pheno_all <- onset_inat |>
  left_join(offset_inat, by = "scientific_name") |>
  left_join(peak_inat, by = "scientific_name") |>
  mutate(dur_inat = offset_doy_inat - onset_doy_inat) |>
  left_join(flora_join, by = "scientific_name")

cat("\nSpecies in final comparison:", nrow(pheno_all), "\n")

# ── 6. Paired tests ───────────────────────────────────────────────────────────
ttest_results_wildflower <- bind_rows(
  tidy(t.test(pheno_all$onset_doy_inat, pheno_all$onset_doy_flora, paired = TRUE)),
  tidy(t.test(pheno_all$peak_doy_inat, pheno_all$peak_doy_flora, paired = TRUE)),
  tidy(t.test(pheno_all$offset_doy_inat, pheno_all$offset_doy_flora, paired = TRUE)),
  tidy(t.test(pheno_all$dur_inat, pheno_all$dur_flora, paired = TRUE))
) |> mutate(phenophase = c("Onset", "Peak", "Offset", "Duration"))

wilcox_results_wildflower <- bind_rows(
  tidy(wilcox.test(pheno_all$onset_doy_inat, pheno_all$onset_doy_flora, paired = TRUE)),
  tidy(wilcox.test(pheno_all$peak_doy_inat, pheno_all$peak_doy_flora, paired = TRUE)),
  tidy(wilcox.test(pheno_all$offset_doy_inat, pheno_all$offset_doy_flora, paired = TRUE)),
  tidy(wilcox.test(pheno_all$dur_inat, pheno_all$dur_flora, paired = TRUE))
) |>
  mutate(phenophase = c("Onset", "Peak", "Offset", "Duration")) |>
  select(phenophase, statistic, p.value)

cat("\n── Paired t-test (mean diff iNat − flora) ──\n")
print(select(ttest_results_wildflower, phenophase, estimate, p.value))

cat("\n── Wilcoxon signed-rank ──\n")
print(wilcox_results_wildflower)

# ── 7. RMSE and R² (vs. null model) ──────────────────────────────────────────
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
r2 <- function(obs_inat, pred_flora) {
  ss_res <- sum((obs_inat - pred_flora)^2, na.rm = TRUE)
  ss_null <- sum((obs_inat - mean(obs_inat, na.rm = TRUE))^2, na.rm = TRUE)
  1 - ss_res / ss_null
}

precision_stats_wildflower <- tibble(
  phenophase = c("Onset", "Peak", "Offset", "Duration"),
  mean_diff = c(
    mean(pheno_all$onset_doy_inat - pheno_all$onset_doy_flora, na.rm = TRUE),
    mean(pheno_all$peak_doy_inat - pheno_all$peak_doy_flora, na.rm = TRUE),
    mean(pheno_all$offset_doy_inat - pheno_all$offset_doy_flora, na.rm = TRUE),
    mean(pheno_all$dur_inat - pheno_all$dur_flora, na.rm = TRUE)
  ),
  RMSE = c(
    rmse(pheno_all$onset_doy_inat, pheno_all$onset_doy_flora),
    rmse(pheno_all$peak_doy_inat, pheno_all$peak_doy_flora),
    rmse(pheno_all$offset_doy_inat, pheno_all$offset_doy_flora),
    rmse(pheno_all$dur_inat, pheno_all$dur_flora)
  ),
  R2 = c(
    r2(pheno_all$onset_doy_inat, pheno_all$onset_doy_flora),
    r2(pheno_all$peak_doy_inat, pheno_all$peak_doy_flora),
    r2(pheno_all$offset_doy_inat, pheno_all$offset_doy_flora),
    r2(pheno_all$dur_inat, pheno_all$dur_flora)
  )
) |> mutate(across(where(is.numeric), \(x) round(x, 2)))

cat("\n── Precision statistics ──\n")
print(precision_stats_wildflower)

# ── 8. Phylogeny ──────────────────────────────────────────────────────────────
sp_list <- gsub(" ", "_", pheno_all$scientific_name)
sp_list_df <- rtrees::sp_list_df(sp_list, taxon = "plant")
filter(sp_list_df, is.na(family))
phy <- rtrees::get_tree(sp_list_df, taxon = "plant")
phy$root.edge <- 0
setdiff(sp_list, phy$tip.label)

pheno_diff <- pheno_all |>
  mutate(
    Onset    = onset_doy_inat - onset_doy_flora,
    Peak     = peak_doy_inat - peak_doy_flora,
    Offset   = offset_doy_inat - offset_doy_flora,
    Duration = dur_inat - dur_flora,
    sp       = gsub(" ", "_", scientific_name)
  ) |>
  select(sp, Onset, Peak, Offset, Duration) |>
  as.data.frame()
row.names(pheno_diff) <- pheno_diff$sp
pheno_diff <- pheno_diff[phy$tip.label, ]

# ── 9. Phylogenetic signal ────────────────────────────────────────────────────
get_phy_signal <- function(v, phy, names_v) {
  names(v) <- names_v
  k <- phytools::phylosig(phy, v, method = "K", test = TRUE)
  l <- phytools::phylosig(phy, v, method = "lambda", test = TRUE)
  tibble(
    method = c("K", "lambda"),
    statistic = c(k$K, l$lambda),
    p = c(k$P, l$P)
  )
}

physig_all_wildflower <- bind_rows(
  get_phy_signal(pheno_all$onset_doy_inat, phy, sp_list) |> mutate(source = "iNat", type = "Onset"),
  get_phy_signal(pheno_all$onset_doy_flora, phy, sp_list) |> mutate(source = "flora", type = "Onset"),
  get_phy_signal(pheno_all$peak_doy_inat, phy, sp_list) |> mutate(source = "iNat", type = "Peak"),
  get_phy_signal(pheno_all$peak_doy_flora, phy, sp_list) |> mutate(source = "flora", type = "Peak"),
  get_phy_signal(pheno_all$offset_doy_inat, phy, sp_list) |> mutate(source = "iNat", type = "Offset"),
  get_phy_signal(pheno_all$offset_doy_flora, phy, sp_list) |> mutate(source = "flora", type = "Offset"),
  get_phy_signal(pheno_all$dur_inat, phy, sp_list) |> mutate(source = "iNat", type = "Duration"),
  get_phy_signal(pheno_all$dur_flora, phy, sp_list) |> mutate(source = "flora", type = "Duration")
) |> arrange(type, method)

cat("\n── Phylogenetic signal ──\n")
print(physig_all_wildflower)

# ── 10. Trait PGLS ───────────────────────────────────────────────────────────
# bloom_group
bg_df <- pheno_diff |>
  mutate(sp = row.names(pheno_diff)) |>
  left_join(
    mutate(flora_join, sp = gsub(" ", "_", scientific_name)),
    by = "sp"
  ) |>
  as.data.frame()
row.names(bg_df) <- bg_df$sp
bg_df <- bg_df[phy$tip.label, ]

# growth_form_3: reference = Herb_annual/biennial
bg_df$bloom_group <- factor(bg_df$bloom_group, levels = c("Summer", "Spring"))
bg_df$growth_form_3 <- factor(bg_df$growth_form_3,
  levels = c("Herb_annual/biennial", "Herb_perennial", "Woody_perennial")
)

run_pgls <- function(df, response, predictor, phy) {
  f <- as.formula(paste(response, "~", predictor))
  m <- phylolm(f, data = df, phy = phy, model = "lambda")
  class(m) <- c(class(m), "lm")
  tidy(m)
}

trait_pgls_all_wildflower <- bind_rows(
  run_pgls(bg_df, "Onset", "bloom_group", phy) |> mutate(trait = "Blooming Group", phenostage = "Onset"),
  run_pgls(bg_df, "Peak", "bloom_group", phy) |> mutate(trait = "Blooming Group", phenostage = "Peak"),
  run_pgls(bg_df, "Offset", "bloom_group", phy) |> mutate(trait = "Blooming Group", phenostage = "Offset"),
  run_pgls(bg_df, "Duration", "bloom_group", phy) |> mutate(trait = "Blooming Group", phenostage = "Duration"),
  run_pgls(bg_df, "Onset", "growth_form_3", phy) |> mutate(trait = "Growth Form", phenostage = "Onset"),
  run_pgls(bg_df, "Peak", "growth_form_3", phy) |> mutate(trait = "Growth Form", phenostage = "Peak"),
  run_pgls(bg_df, "Offset", "growth_form_3", phy) |> mutate(trait = "Growth Form", phenostage = "Offset"),
  run_pgls(bg_df, "Duration", "growth_form_3", phy) |> mutate(trait = "Growth Form", phenostage = "Duration")
) |>
  mutate(
    term = str_remove(term, "bloom_group|growth_form_3"),
    phenostage = factor(phenostage, levels = c("Onset", "Peak", "Offset", "Duration"))
  ) |>
  select(trait, phenostage, term, estimate, std.error, statistic, p.value) |>
  arrange(trait, phenostage)

cat("\n── Trait PGLS ──\n")
print(as.data.frame(mutate_if(trait_pgls_all_wildflower, is.numeric, round, 3)))

# ── 11. Range size ────────────────────────────────────────────────────────────
all_d <- readRDS("data/dat_used_phenobase.rds")

# convex hull polygon per species
fc <- function(x) {
  select(x, scientific_name) |>
    mutate(geometry = st_convex_hull(st_combine(x))) |>
    unique()
}
all_d_polygon <- all_d |>
  group_by(scientific_name) |>
  do(fc(.)) |>
  ungroup() |>
  st_as_sf()

range_area <- tibble(
  sp       = all_d_polygon$scientific_name,
  area_km2 = as.vector(units::set_units(st_area(all_d_polygon), km^2))
)

all_sp_range_area <- range_area |>
  left_join(
    pheno_all |>
      mutate(
        Onset    = onset_doy_inat - onset_doy_flora,
        Peak     = peak_doy_inat - peak_doy_flora,
        Offset   = offset_doy_inat - offset_doy_flora,
        Duration = dur_inat - dur_flora
      ) |>
      select(sp = scientific_name, Onset, Peak, Offset, Duration),
    by = "sp"
  ) |>
  mutate(area_km2_log10 = log10(area_km2))

cat("\n── Range size ~ delta PGLS ──\n")
range_pgls_df <- all_sp_range_area |>
  mutate(sp = gsub(" ", "_", sp)) |>
  as.data.frame() |>
  column_to_rownames("sp")
range_pgls_df <- range_pgls_df[phy$tip.label[phy$tip.label %in% rownames(range_pgls_df)], ]

range_delta_pgls_all_coef_wildflower <- lapply(
  c("Onset", "Peak", "Offset", "Duration"),
  function(ph) {
    df <- range_pgls_df[!is.na(range_pgls_df[[ph]]), ]
    m <- phylolm(as.formula(paste(ph, "~ area_km2_log10")),
      data = df, phy = drop.tip(phy, setdiff(phy$tip.label, rownames(df))),
      model = "lambda"
    )
    class(m) <- c(class(m), "lm")
    tidy(m) |> mutate(phenostage = ph)
  }
) |> bind_rows()

print(as.data.frame(mutate_if(
  select(range_delta_pgls_all_coef_wildflower, phenostage, term, estimate, std.error, p.value),
  is.numeric, round, 3
)))

# ── 12. Figures ──────────────────────────────────────────────────────────────

# Fig 1 — species-level overview (onset/offset bars + flora window), 30 random species
set.seed(42)
fig1_dat <- pheno_all |>
  slice_sample(n = 30)

d_fig1 <- filter(est_pheno_df, scientific_name %in% fig1_dat$scientific_name) |>
  arrange(scientific_name, peak_doy) %>%
  group_by(scientific_name) %>%
  mutate(idx = 1:n()) %>%
  ungroup() %>%
  mutate(scientific_name = fct_reorder(scientific_name, onset_doy))

fig1_dat <- fig1_dat %>%
  mutate(scientific_name = factor(scientific_name, levels = levels(d_fig1$scientific_name)))


p1 <- ggplot(d_fig1) +
  geom_rect(
    data = fig1_dat,
    aes(xmin = onset_doy_flora, xmax = offset_doy_flora, ymin = -Inf, ymax = Inf),
    fill = "lightblue", alpha = 1
  ) +
  # geom_point(aes(x = peak_doy, y = idx)) +
  geom_segment(aes(x = onset_doy, y = idx, xend = offset_doy, yend = idx), size = 0.7, alpha = 0.7) +
  geom_segment(data = fig1_dat, aes(x = onset_doy_inat, y = Inf, xend = offset_doy_inat, yend = Inf), size = 1, color = "red") +
  # coord_cartesian(data = limits, xlim = c(min_doy, max_doy)) +
  # geom_vline(data = dd1, aes(xintercept = ave_peak_doy), linetype = 2) +
  # geom_vline(data = pheno_flora2, aes(xintercept = peak_doy), linetype = 2, color = "blue") +
  facet_wrap(~scientific_name, scales = "free", ncol = 6) +
  labs(x = "Day of year", y = "") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "italic")
  )

ggsave("figures/wildflower/fig_1_overview.pdf", plot = p1, width = 11, height = 9)
ggsave("figures/wildflower/fig_1_overview.png", plot = p1, width = 11, height = 9)

# Fig 2 — differences by species (lollipop), ordered by Peak difference
sp_order_peak <- pheno_all |>
  mutate(Peak = peak_doy_inat - peak_doy_flora) |>
  arrange(Peak) |>
  pull(scientific_name)

p2 <- pheno_all |>
  mutate(
    Onset = onset_doy_inat - onset_doy_flora,
    Peak = peak_doy_inat - peak_doy_flora,
    Offset = offset_doy_inat - offset_doy_flora,
    Duration = dur_inat - dur_flora,
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
  labs(x = "iNat-based - flora-based (days)", y = "Species") +
  cowplot::theme_cowplot(font_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("figures/wildflower/fig_2_differences.pdf", plot = p2, width = 10, height = 7)
ggsave("figures/wildflower/fig_2_differences.png", plot = p2, width = 10, height = 7)

# Fig 3 — trait differences
all_sp_traits_long <- pheno_all |>
  mutate(
    Onset    = onset_doy_inat - onset_doy_flora,
    Peak     = peak_doy_inat - peak_doy_flora,
    Offset   = offset_doy_inat - offset_doy_flora,
    Duration = dur_inat - dur_flora
  ) |>
  select(Onset, Peak, Offset, Duration,
    `Blooming Group` = bloom_group,
    `Growth Form`    = growth_form_3,
    sp               = scientific_name
  ) |>
  pivot_longer(Onset:Duration, names_to = "phenostage", values_to = "diff_value") |>
  pivot_longer(`Blooming Group`:`Growth Form`, names_to = "trait", values_to = "trait_value") |>
  mutate(phenostage = factor(phenostage, levels = c("Onset", "Peak", "Offset", "Duration")))

p3_1 <- ggplot(
  filter(all_sp_traits_long, trait == "Blooming Group"),
  aes(x = trait_value, y = diff_value)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.colour = NA, width = 0.4) +
  geom_jitter(color = "blue", width = 0.1, alpha = 0.4, size = 0.8) +
  facet_wrap(~phenostage, scales = "free_y", nrow = 1) +
  labs(y = "iNat - flora (days)", x = "Blooming Group") +
  cowplot::theme_cowplot()

p3_2 <- ggplot(
  filter(all_sp_traits_long, trait == "Growth Form"),
  aes(x = trait_value, y = diff_value)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.colour = NA, width = 0.4) +
  geom_jitter(color = "blue", width = 0.1, alpha = 0.4, size = 0.8) +
  facet_wrap(~phenostage, scales = "free_y", nrow = 1) +
  labs(y = "iNat - flora (days)", x = "Growth Form") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background.x = element_blank(), strip.text.x = element_blank()
  )

p3 <- cowplot::plot_grid(p3_1, p3_2, nrow = 2, labels = "AUTO") +
  theme(plot.margin = margin(b = 30, r = 5, t = 5, l = 15))
ggsave("figures/wildflower/fig_3_traits.pdf", plot = p3, width = 12, height = 9)
ggsave("figures/wildflower/fig_3_traits.png", plot = p3, width = 12, height = 9)

# Fig 4 — range size vs phenology differences
# annotations derived from PGLS results: line solid if slope p < 0.05
pgls_range_annot <- left_join(
  filter(range_delta_pgls_all_coef_wildflower, term == "(Intercept)") |>
    select(phenostage, intercept = estimate),
  filter(range_delta_pgls_all_coef_wildflower, term == "area_km2_log10") |>
    select(phenostage, slope = estimate, p_slope = p.value),
  by = "phenostage"
) |>
  mutate(
    type = factor(phenostage, levels = c("Onset", "Peak", "Offset", "Duration")),
    linetype_pgls = if_else(p_slope < 0.05, "solid", "dashed"),
    label = paste0(
      "y = ", round(intercept, 1),
      ifelse(slope >= 0, " + ", " − "), abs(round(slope, 2)), "x;",
      "     p = ", ifelse(p_slope < 0.001, "<0.001", round(p_slope, 3))
    )
  )

p4 <- pivot_longer(all_sp_range_area,
  cols = c(Onset, Peak, Offset, Duration),
  names_to = "type", values_to = "delta"
) |>
  mutate(type = factor(type, levels = c("Onset", "Peak", "Offset", "Duration"))) |>
  left_join(select(pgls_range_annot, type, linetype_pgls), by = "type") |>
  ggplot(aes(x = area_km2_log10, y = delta)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(aes(linetype = linetype_pgls), method = "lm", se = TRUE) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_linetype_identity() +
  geom_text(
    data = pgls_range_annot,
    aes(x = -Inf, y = -Inf, label = label),
    hjust = -0.05, vjust = -0.5, size = 2.8, inherit.aes = FALSE
  ) +
  facet_wrap(~type, scales = "free_y") +
  labs(
    x = expression(Log[10] ~ range ~ area ~ (km^2)),
    y = "iNat-based − flora-based (days)"
  ) +
  cowplot::theme_cowplot()
ggsave("figures/wildflower/fig_4_rangesize.pdf", width = 9, height = 7)
ggsave("figures/wildflower/fig_4_rangesize.png", width = 9, height = 7)

# ── 13. Save stats and intermediate data ─────────────────────────────────────
save(wilcox_results_wildflower, ttest_results_wildflower, precision_stats_wildflower,
  physig_all_wildflower, trait_pgls_all_wildflower, range_delta_pgls_all_coef_wildflower,
  file = "data_output/stats_wildflower.RData"
)

saveRDS(pheno_all, "data_output/pheno_all_wildflower.rds")
saveRDS(est_pheno_df, "data_output/est_pheno_df_wildflower.rds")
saveRDS(all_sp_range_area, "data_output/range_area_wildflower.rds")
saveRDS(pgls_range_annot, "data_output/pgls_range_annot_wildflower.rds")

cat("\nDone. Figures in figures/wildflower/, stats in data_output/stats_wildflower.RData\n")
