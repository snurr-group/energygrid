# WARNING: DEPRECATED
# Code became unmaintainable, so transformed it into an intermediate version called save_h2_hists.R

# Calculates histograms for H2 using reasonable parameters (which can be overwritten)
# Caches values in hist_vals to avoid long re-computation time when running the code block
# Need to first load "R/load_data.R" and "R/get_energy_stats.R"

R_GAS_KJ <- 8.314 / 1000

if (!exists("OVERRIDE_H2_HIST_PARAMS")) {
  DATA_SPLIT <- c(0.4, 0.4, 0.2)  # Split for hyperparams, training, and test data
  # Minimum histogram parameters, in kJ/mol
  min_bin_width <- 0.01
  hist_range <- c(-8, 1)  # about 1.5kT in the positive direction, at 77 K
}


# Through testing, the default params for a function (tidy_energy_hists) evaluate the global variables at runtime, not "compile"/assignment, so we can reset the appropriate parameters here.
# Also convert J/mol to original K units reported in the energy grid scripts.
BIN_WIDTH <- min_bin_width / R_GAS_KJ  # minimum histogram bin width we'll consider
ENERGY_RANGE <- hist_range / R_GAS_KJ

if (!exists("hist_vals")) {
  hist_vals <- run_energy_stat(ANALYSIS_DIRS, tidy_energy_hists, bin_width = BIN_WIDTH, min_max = ENERGY_RANGE / BIN_WIDTH)
  # Also convert our hist_vals to more convenient units of kJ/mol
  hist_vals <- mutate(hist_vals, lower=lower*R_GAS_KJ, upper=upper*R_GAS_KJ)
}

