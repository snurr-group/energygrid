# Gets statistics on calculated energies from the scotty_map.py code
library(tidyverse)
library(stringr)
library(R.utils)
library(manipulate)

# deprecated cutoff parameters
TEMPERATURE = 77  # K
E_CUTOFF <- 15 * TEMPERATURE  # K

BIN_WIDTH = 5
ENERGY_RANGE = c(-1000, 20)
ANALYSIS_DIRS <- c("BigData/10k-hMOFs/part1/CIF_FILES", "BigData/10k-hMOFs/part2/CIF_FILES")
QUICK_TEST <- FALSE  # Set to true to do a "practice run" instead of all of the files


k_to_kj_mol <- function(energy)  {
  kb <- 1.38064853e-23  # J/K from Wikipedia
  na <- 6.022e23
  energy * kb * na / 1000
}

energy_stats <- function(data_dir, stats_fcn, df_prototype, num_rows = 1) {
  # Run a stats_fcn on all files nestled within data_dir, according to the cif_dir spec
  # Returns a data.frame, where each CIF's stats has num_rows by df_prototype entries
  dirs <- list.files(data_dir)
  if (QUICK_TEST) {
    dirs <- dirs[1:100]  # Debugging trick to only run the first 100 folders
  }
  num_cifs <- length(dirs)
  
  # Preallocate a df using an R.utils function: http://r.789695.n4.nabble.com/idiom-for-constructing-data-frame-td4705353.html
  # Names for the dataframe will automatically be copied from the prototype
  all_stats <- dataFrame(df_prototype, nrow = num_cifs*num_rows)
  
  write(paste0("Compiling statistics from the energy files in ", data_dir, "..."), "")
  pb <- txtProgressBar(min = 1, max = num_cifs, style = 3)  # Add a progress bar, courtesy of https://www.r-bloggers.com/r-monitoring-the-function-progress-with-a-progress-bar/
  current_row = 1
  for (cif_dir in dirs) {
    energy_file <- file.path(data_dir, cif_dir, "Energy_Values.txt")
    energy <- read_tsv(energy_file, col_names = "V1", col_types = "d")$V1
    end_row <- current_row + num_rows - 1
    all_stats[current_row:end_row, ] <- stats_fcn(energy)
    current_row = current_row + num_rows
    setTxtProgressBar(pb, (current_row-1)/num_rows)
  }
  close(pb)
  
  ids <- as.integer(str_sub(dirs, 2, -1))  # strip off leading "h" for hMOF designation
  ids <- rep(ids, each=num_rows)
  all_stats$id <- ids
  all_stats
}


energy_summary <- function(data_dir, upper_cutoff = E_CUTOFF) {
  # Computes the six-number summary for the pretabulated energies, using a max cutoff
  energy_summary_fcn <- function(energy) {
    filtered_energy <- energy[energy < upper_cutoff]
    energy_row <- summary(filtered_energy)
  }
  output_prototype <- rep("numeric", 6)
  names(output_prototype) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  
  energy_stats(data_dir, energy_summary_fcn, output_prototype)
}

energy_metric <- function(data_dir, lower_bound = -200, upper_bound = 0) {
  # Computes the "LJ metric" based on upper and lower cutoffs
  # DEMO function.  Ideally, this should be rapidly done via integrating the histograms
  energy_metric_fcn <- function(energy) {
    filtered_energy <- energy[energy > lower_bound & energy < upper_bound]
    energy_row <- length(filtered_energy) / length(energy)
  }
  output_prototype <- c("numeric")
  names(output_prototype) <- c("LJ.metric")
  
  energy_stats(data_dir, energy_metric_fcn, output_prototype)
}

metric_from_hists <- function(hist_df, lower_bound = -200, upper_bound = 0) {
  # Compute the "LJ metric" based on given cutoffs
  # Set a bound to NA for open intervals (e.g., energy > -200, but no upper bound)
  # TODO: NA code
  if (!(lower_bound %in% hist_df$lower & upper_bound %in% hist_df$upper)) {
    warning("Metric bounds do not exactly line up with a histogram bin")
  }
  
  good_counts <- hist_df %>%
    filter(lower >= lower_bound) %>% filter(upper <= upper_bound) %>% 
    group_by(id) %>% summarize(good = sum(counts))
  
  total_counts <- hist_df %>% group_by(id) %>% summarize(total = sum(counts))
  
  lj_metric <- total_counts %>%
    inner_join(good_counts, by="id") %>% 
    mutate(metric = good / total) %>% 
    select(id, metric)
}


row_energy_hists <- function(data_dir, bin_width = TEMPERATURE, min_max = c(-15, 15)) {
  # Retrieve the histograms from the energy directories
  # min_max: multiplier for minimum and maximum bins, e.g. +/- 15kT
  bins <- seq(from = bin_width * min_max[1],
              to   = bin_width * (min_max[2] + 1),
              by   = bin_width
              )
  hists_fcn <- function(energy) {
    # First cap the maximum energy so that the last bin (energy > max) includes the high end
    energy[energy > bin_width * (min_max[2]) + 1] <- bin_width * (min_max[2] + 0.5)
    energy_row <- hist(energy, breaks = bins, plot = FALSE)$counts
  }
  output_prototype <- rep("integer", length(bins) - 1)
  # TODO: think about how to include the bin information in the data.frame (in the names, perhaps?)
  # It's encouraging that the file reading and IO appears like the most compute-heavy part of this work
  
  energy_stats(data_dir, hists_fcn, output_prototype)
  # Also redone in a tidy format, which will facilitate energy scans
}


tidy_energy_hists <- function(data_dir, bin_width = BIN_WIDTH, min_max = ENERGY_RANGE / BIN_WIDTH) {
  # Retrieve the histograms from the energy directories
  # min_max: multiplier for minimum and maximum bins, e.g. +/- 15kT
  # Returns a tidy data.frame, where the histogram bins are explicitly specified as a column
  # (instead of embedding the energies in the column titles)
  
  bins <- seq(from = bin_width * min_max[1],
              to   = bin_width * (min_max[2] + 1),
              by   = bin_width
  )
  hists_fcn <- function(energy) {
    # First cap the maximum energy so that the last bin (energy > max) includes the high end
    energy[energy > bin_width * (min_max[2]) + 1] <- bin_width * (min_max[2] + 0.5)
    raw_hist <- hist(energy, breaks = bins, plot = FALSE)
    cbind(lower = raw_hist$breaks[1:(length(raw_hist$breaks)-1)],
          upper = raw_hist$breaks[2:length(raw_hist$breaks)],
          counts = raw_hist$counts
          )
  }
  output_prototype <- rep("integer", 3)  # Lower bound, upper bound, count
  names(output_prototype) <- c("lower", "upper", "counts")
  # Could also simplify this with only the lower bound to save 1/3 memory requirements
  # Even so, 100k structures will comfortably fit in 600MB of RAM
  
  energy_stats(data_dir, hists_fcn, output_prototype, length(bins)-1)
  # Need to come back to the summarization command as well
}


bin_labels <- function(bin_width = TEMPERATURE, min_max = c(-15, 15)) {
  # Get the lower boundaries of the histogram "breaks"
  # Uses the same parameters as energy_hists
  lower_bins <- seq(from = bin_width * min_max[1],
                    to   = bin_width * (min_max[2]),
                    by   = bin_width
                    )
}

run_energy_stat <- function(dirs, stat_fcn) {
  # Runs the stat function over multiple directories to reduce copy-paste
  all_stats <- stat_fcn(dirs[1])
  if (length(dirs) > 1) {
    for (x in dirs[2:length(dirs)]) {
      all_stats <- rbind(all_stats, stat_fcn(x))
    }
  }
  all_stats
}

#hist_summary <- run_energy_stat(ANALYSIS_DIRS, energy_summary)
hist_vals <- run_energy_stat(ANALYSIS_DIRS, tidy_energy_hists)
#hist_metric <- run_energy_stat(ANALYSIS_DIRS, energy_metric)
hist_metric <- metric_from_hists(hist_vals)

run_interactive_lj_plot <- function() {
  # Requires that gcmc_data is already loaded into the workspace
  # As shown by the [Rstudio documentation](https://support.rstudio.com/hc/en-us/articles/200551906-Interactive-Plotting-with-Manipulate)
  # manipulate is an easy package to use in this IDE.
  # Could probably be more efficient (fast), but honestly it's impressive it works as well as it does
  print("Click on the gear to play with the upper and lower bounds!")
  print(paste0("Currently, +", (ENERGY_RANGE[2]+1)*BIN_WIDTH, " is a special entry corresponding to >", ENERGY_RANGE[2]*BIN_WIDTH, "K"))
  manipulate(
    hist_vals %>% 
      metric_from_hists(lower=lower_bound, upper=upper_bound) %>%
      left_join(gcmc_data, by="id") %>%
      ggplot(aes(metric, g.L)) + geom_point(),
    lower_bound = slider(ENERGY_RANGE[1], ENERGY_RANGE[2], initial=-200, step=BIN_WIDTH),
    upper_bound = slider(ENERGY_RANGE[1], ENERGY_RANGE[2], initial=-10, step=BIN_WIDTH)
    )
  # hist_vals %>% filter(counts > 0) %>% select(lower) %>% min
  # will return -830 as the lowest occupied energy bin
}

