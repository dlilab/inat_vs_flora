# Reproduce all statistical results and figures for the manuscript
# "Flora-based phenological information is not a reliable proxy of actual observations".
# Run from the root of this folder (the folder containing this file), e.g.
#   Rscript run_all.R
# Scripts are run in order; later scripts use outputs written by earlier ones.

for (d in c("data_output", "figures/wildflower", "figures/calflora", "figures/combined")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

scripts <- c(
  "code/01_compare_wildflower.R",             # national comparison (Wildflower Center)
  "code/02_compare_calflora.R",               # regional comparison (Jepson eFlora via CalFlora)
  "code/03_combined_figures.R",               # combined main-text figures
  "code/04_flora_jitter_sensitivity.R",       # sensitivity to monthly resolution of flora records
  "code/05_range_size_effort_robustness.R"    # range-size effect controlling for observation effort
)

for (s in scripts) {
  message("\n==== Running ", s, " ====")
  source(s, echo = FALSE)
}
