# combine energy grid rds
grid_file = "All_data/CH3_1A_probe_UFF.rds"
grid2_file = "All_data/CH3_1A_CoRE.rds"

grid_first <- read_rds(grid_file)
grid_second <- read_rds(grid2_file)

all_grid <- rbind(grid_first, grid_second)
# finally, save to a new rds
write_rds(all_grid, "All_data/CH3_1A_combined_tobacco_core.rds")
