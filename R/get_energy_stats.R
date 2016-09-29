# Gets statistics on calculated energies from the scotty_map.py code
library(tidyverse)
library(stringr)
library(R.utils)

E_CUTOFF <- 15 * 77  # K

k_to_kj_mol <- function(energy)  {
  kb <- 1.38064853e-23  # J/K from Wikipedia
  na <- 6.022e23
  energy * kb * na / 1000
}

energy_stats <- function(data_dir, stats_fcn, df_prototype) {
  dirs <- list.files(data_dir)
  dirs <- dirs[1:100]  # Debugging trick to only run the first 100 folders
  num_cifs <- length(dirs)
  
  # Preallocate a df using an R.utils function: http://r.789695.n4.nabble.com/idiom-for-constructing-data-frame-td4705353.html
  # Names for the dataframe will automatically be copied from the prototype
  all_stats <- dataFrame(df_prototype, nrow = num_cifs)
  
  write(paste0("Compiling statistics from the energy files in ", data_dir, "..."), "")
  pb <- txtProgressBar(min = 1, max = num_cifs, style = 3)  # Add a progress bar, courtesy of https://www.r-bloggers.com/r-monitoring-the-function-progress-with-a-progress-bar/
  current_row = 1
  for (cif_dir in dirs) {
    energy_file <- file.path(data_dir, cif_dir, "Energy_Values.txt")
    energy <- read_tsv(energy_file, col_names = "V1", col_types = "d")$V1
    all_stats[current_row, ] <- stats_fcn(energy)
    current_row = current_row + 1
    setTxtProgressBar(pb, current_row)
  }
  close(pb)
  
  ids <- str_sub(dirs, 2, -1)  # strip off leading "h" for hMOF designation
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

energy_hists <- function(data_dir, bin_width = 77, min_max = c(-15, 15)) {
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
}

bin_labels <- function(bin_width = 77, min_max = c(-15, 15)) {
  # Get the lower boundaries of the histogram "breaks"
  # Uses the same parameters as energy_hists
  lower_bins <- seq(from = bin_width * min_max[1],
                    to   = bin_width * (min_max[2]),
                    by   = bin_width
                    )
}


all_stats_1 <- energy_summary("BigData/10k-hMOFs/part1/CIF_FILES")
all_stats_2 <- energy_summary("BigData/10k-hMOFs/part2/CIF_FILES")
all_stats_combined <- rbind(all_stats_1, all_stats_2)
temp_hists_1 <- energy_hists("BigData/10k-hMOFs/part1/CIF_FILES")
temp_hists_2 <- energy_hists("BigData/10k-hMOFs/part2/CIF_FILES")
hist_stats <- rbind(temp_hists_1, temp_hists_2)


to_delete <- function(x) {
  bins <- seq(from=-15*77, to=16*77, by=77)
  energy[energy>15*77] <- 15.5*77
  
  # Preallocate array: https://www.r-bloggers.com/pitfall-did-you-really-mean-to-use-matrixnrow-ncol/
  new_hist <- data.frame(matrix(NA_integer_, nrow=134, ncol=31))
  new_hist[2,] <- hist(energy, breaks=bins, plot=FALSE)$counts
}


