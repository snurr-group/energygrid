source("R/package_verification.R")
# compute difference in energy between two energy grids
Grid_1_path_of_files <- "../DATAS/Xe_1A_CoRE/Energies"
file_1 <- list.files(Grid_1_path_of_files)
Grid_2_path_of_files <- "../DATAS/Kr_grid_CORE/Energies"
file_2 <- list.files(Grid_2_path_of_files)
# the directory for writting difference vector
diff_dir <- "../DATAS/Diff_grids_XeKr"
if (length(file_1) != length(file_2)){
  warning("Numbers of files in two folders should be equal.")
}
# loop over all files 
for (grid_file in file_1){
  grid_file_1 <- paste0(Grid_1_path_of_files, "/", grid_file)
  grid_file_2 <- paste0(Grid_2_path_of_files, "/", grid_file)
  energy_1 <- read_table(grid_file_1, col_names = "V1", col_types = "d")$V1 
  energy_2 <- read_table(grid_file_2, col_names = "V1", col_types = "d")$V1 
  difference <- energy_1 - energy_2 # Xe - Kr
  positions_1 <- which(energy_1 == 1.0e+23)
  positions_2 <- which(energy_2 == 1.0e+23)
  # merge these two vectors and keep the unique values.
  positions <- c(positions_1, positions_2) %>% unique()
  difference[positions] <- 1.0e+23
  #for (a in 1:length(energy_1)){
  # if (energy_1[a] == 1.0e23 || energy_2[a] == 1.0e23){
  #    difference[a] = 1.0e23
  #  }
    write.table(difference, paste0(diff_dir, "/", grid_file), 
                row.names = FALSE, 
                col.names = FALSE, 
                sep = "\t")
}

# convert Energy values to exp(deltaE/kBT)
diff_dir_kB <- "../DATAS/Diff_grids_XeKr_kB_Favors_Xe"
Temperature <- 273 # XeKr simulation at 273 Kelvin
for (grid_file in file_1){
  grid_file_1 <- paste0(Grid_1_path_of_files, "/", grid_file)
  grid_file_2 <- paste0(Grid_2_path_of_files, "/", grid_file)
  energy_1 <- read_table(grid_file_1, col_names = "V1", col_types = "d")$V1 
  energy_2 <- read_table(grid_file_2, col_names = "V1", col_types = "d")$V1 
  difference <- energy_1 - energy_2 # Xe - Kr
  positions_1 <- which(energy_1 == 1.0e+23)
  positions_2 <- which(energy_2 == 1.0e+23)
  # merge these two vectors and keep the unique values.
  positions <- c(positions_1, positions_2) %>% unique()
  difference[positions] <- 1.0e+23
  # exp gives huge values, no exponentials
  boltz <- difference/Temperature 
  # one will get a lot of large numbers,consider a threshold, 
  # replace large values with place-holder
  #boltz[which(boltz >= 10)] <- 1.0e+23
  boltz[positions] <- 1.0e+23
  write.table(boltz, paste0(diff_dir_kB, "/", grid_file), 
              row.names = FALSE, 
              col.names = FALSE, 
              sep = "\t")
}
