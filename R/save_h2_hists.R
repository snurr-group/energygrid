#!/usr/bin/env Rscript
# Consider running Rscript with the "--vanilla" flag

# Calculates histograms for H2 using reasonable parameters (which can be overwritten).
# Arguments are the output file followed by a list of input analysis directories.
# Could also consider a more explicit `optparse` implementation with command-line config of params

# added new feature: auto-tuning the lower-bound (081619)

library(dplyr)
library(readr)
source("R/package_verification.R")
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
  # Minimum histogram parameters, in kJ/mol
  min_bin_width <- 0.25
  hist_range <- c(-25, 5)  # Need -25 kJ/mol for tobmof5885
}

if (args[length(args)] == "autotune") {
  # this will give the code permission to determine the lower_bound by itself
  write("Overriding H2 parameters for the energy histogram", "")
  write("Determining lower bound by itself", "")
  min_bin_width <- 0.10
  GRID_DIR <- ANALYSIS_DIRS[1:(length(ANALYSIS_DIRS)-1)]
  minimum_energy <- determine_lower_bound(path_of_files = GRID_DIR)
  write(paste0("Minimum Energy is: ", minimum_energy, " (Kelvin)"), "")
  # warning message has an additional half of bin width, so minus one to make it lower
  lower_boundary <- floor(minimum_energy*R_GAS_KJ) - 1
  write(paste0("Lower Boundary of the pre-calculated histograms is: ",
               lower_boundary, " (kJ/mol)"), "")
  hist_range <- c(lower_boundary, 30)
  ANALYSIS_DIRS <- ANALYSIS_DIRS[1:(length(ANALYSIS_DIRS)-1)]
}

if (args[length(args)] == "use_CoRE") {
  # Last argument is a flag to change the histogram parameters.
  # This could also be implemented with a double dash flag eventually
  write("Overriding H2 parameters for the energy histogram with an extended range for CoRE MOFs", "")
  min_bin_width <- 0.10
  hist_range <- c(-113, 30)
  ANALYSIS_DIRS <- ANALYSIS_DIRS[1:(length(ANALYSIS_DIRS)-1)]
}

# Through testing, the default params for a function (tidy_energy_hists) evaluate the global variables at runtime, not "compile"/assignment, so we can reset the appropriate parameters here.
# Also convert J/mol to original K units reported in the energy grid scripts.
BIN_WIDTH <- min_bin_width / R_GAS_KJ  # minimum histogram bin width we'll consider
ENERGY_RANGE <- hist_range / R_GAS_KJ

if (!(length(ANALYSIS_DIRS) == 1)){
  stop("Need just one directory!", call.=FALSE)
  }
  hist_vals <- energy_stats(ANALYSIS_DIRS, bin_width = BIN_WIDTH, min_max = ENERGY_RANGE / BIN_WIDTH)

# Also convert our hist_vals to more convenient units of kJ/mol
hist_vals <- mutate(hist_vals, lower=lower*R_GAS_KJ, upper=upper*R_GAS_KJ)


# Warn if the lower bound is not sufficient to capture all of the data
err_filled_lowest_bins <- 
  hist_vals %>% 
  filter(lower < (hist_range[1] + 0.5*min_bin_width) & counts > 0)
if (nrow(err_filled_lowest_bins) > 0) {
  warning(paste(
    "Warning: bin range not sufficient.  More attractive regions captured in lowest bin for",
    nrow(err_filled_lowest_bins), "MOFs:\n",
    err_filled_lowest_bins$id
    ))
}

# See Rds info in the [R for Data Science book](http://r4ds.had.co.nz/data-import.html)
write_rds(hist_vals, OUTFILE)
