# Gets statistics on calculated energies from the scotty_map.py code
library(tidyverse)

DATA_DIR <- "BigData/10k-hMOFs/part1/CIF_FILES"
E_CUTOFF <- 5000  # K

dirs <- list.files(DATA_DIR)
num_cifs <- length(dirs)
all_stats <- data.frame(numeric(num_cifs), numeric(num_cifs), numeric(num_cifs), numeric(num_cifs), numeric(num_cifs), numeric(num_cifs))
names(all_stats) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
current_row = 1

write(paste0("Compiling statistics from the energy files in ", DATA_DIR, "..."), "")
pb <- txtProgressBar(min = 1, max = num_cifs, style = 3)  # Add a progress bar, courtesy of https://www.r-bloggers.com/r-monitoring-the-function-progress-with-a-progress-bar/
for (cif_dir in dirs) {
  energy_file <- file.path(DATA_DIR, cif_dir, "Energy_Values.txt")
  energy <- read_tsv(energy_file, col_names = "V1", col_types = "d")$V1
  filtered_energy <- energy[energy < E_CUTOFF]
  energy_stats <- summary(filtered_energy)
  all_stats[current_row, ] <- energy_stats
  current_row = current_row + 1
  setTxtProgressBar(pb, current_row)
}
close(pb)
all_stats

