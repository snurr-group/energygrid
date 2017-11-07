#!/usr/bin/env Rscript
# Consider running Rscript with the "--vanilla" flag

# Calculates histograms for H2 using reasonable parameters (which can be overwritten).
# Arguments are the output file followed by a list of input analysis directories.
# Could also consider a more explicit `optparse` implementation with command-line config of params

library(dplyr)
library(readr)
source("R/get_energy_stats.R")

# See arg documentation at https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/
args = commandArgs(trailingOnly=TRUE)
OUTFILE = "BigData/hist_vals.Rds"

if (length(args) == 0) {
  stop("Need to at least specify the output file", call.=FALSE)
} else {
  OUTFILE = args[1]
  ANALYSIS_DIRS <- c("BigData/10k-hMOFs/part1/CIF_FILES", "BigData/10k-hMOFs/part2/CIF_FILES")
  if (length(args) > 1) {
    ANALYSIS_DIRS = args[2:length(args)]
  }
}

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

hist_vals <- run_energy_stat(ANALYSIS_DIRS, tidy_energy_hists, bin_width = BIN_WIDTH, min_max = ENERGY_RANGE / BIN_WIDTH)
# Also convert our hist_vals to more convenient units of kJ/mol
hist_vals <- mutate(hist_vals, lower=lower*R_GAS_KJ, upper=upper*R_GAS_KJ)

# See Rds info in the [R for Data Science book](http://r4ds.had.co.nz/data-import.html)
write_rds(hist_vals, OUTFILE)
