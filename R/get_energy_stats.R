# Gets statistics on calculated energies from the scotty_map.py code
library(tidyverse)

E_CUTOFF <- 15 * 77  # K

k_to_kj_mol <- function(energy)  {
  kb <- 1.38064853e-23  # J/K from Wikipedia
  na <- 6.022e23
  energy * kb * na / 1000
}

energy_summary <- function(data_dir, upper_cutoff = E_CUTOFF) {
  dirs <- list.files(data_dir)
  num_cifs <- length(dirs)
  all_stats <- data.frame(numeric(num_cifs), numeric(num_cifs), numeric(num_cifs), numeric(num_cifs), numeric(num_cifs), numeric(num_cifs))
  names(all_stats) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  current_row = 1
  
  write(paste0("Compiling statistics from the energy files in ", data_dir, "..."), "")
  pb <- txtProgressBar(min = 1, max = num_cifs, style = 3)  # Add a progress bar, courtesy of https://www.r-bloggers.com/r-monitoring-the-function-progress-with-a-progress-bar/
  for (cif_dir in dirs) {
    energy_file <- file.path(data_dir, cif_dir, "Energy_Values.txt")
    energy <- read_tsv(energy_file, col_names = "V1", col_types = "d")$V1
    filtered_energy <- energy[energy < upper_cutoff]
    energy_stats <- summary(filtered_energy)
    all_stats[current_row, ] <- energy_stats
    current_row = current_row + 1
    setTxtProgressBar(pb, current_row)
  }
  close(pb)
  
  ids <- str_sub(dirs, 2, -1)  # strip off leading "h" for hMOF designation
  all_stats$id <- ids
  all_stats
}

all_stats_1 <- energy_summary("BigData/10k-hMOFs/part1/CIF_FILES")
all_stats_2 <- energy_summary("BigData/10k-hMOFs/part2/CIF_FILES")
all_stats_combined <- rbind(all_stats_1, all_stats_2)

