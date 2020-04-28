rm(list=ls())
# for checking package availability, if not, install it
source("R/package_verification.R")

source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/read_gcmc_data.R")
source("R/read_tobacco_data.R")
source("R/save_train_test_data.R")
source("R/refined_bins_calc.R")
# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4
# define the input files
grid_file = "All_data/CH3_probe_1A_UFF.rds"
gcmc_file = "All_data/propane_1bar_UFF_2000points.txt"
if (!dir.exists("Results")){
  dir.create("Results")
}
#hmof_h2_grid <- read_rds(grid_file)
gcmc_data <- read_data(gcmc_file)

hmof_h2_grid <- read_rds("All_data/Grids_Propane_298K_1A_GenericMOFs.Rds")

# here we can have two options, feed uniform or non-uniform bins
#binbounds <- bounds_from_params(ch4_binspec)
gcmc_data <- mutate(gcmc_data, id=ID)
new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
binbounds <- automatic_bins(new_grid) # for non-uniform bins
# get the names of files from that path
namess <- list.files(path="Energies/") %>% str_remove(., ".grid")
gcmc_data <- gcmc_data %>% filter(ID %in% namess)
#source("R/get_heat_map.R")
# get the names from gcmc_data$ID

# count = 0
# for (MOF_ID in gcmc_data$ID){
#   abvtf <- heat_maps(MOF_ID, binbounds)
#   count = count + 1
#   if (count == 1) {
#     abvtfs <- abvtf
#   } else{
#     abvtfs <- rbind(abvtfs, abvtf) 
#   }
# }
#saveRDS(abvtfs, "heat_map.rds")
abvtfs <- readRDS("heat_map_more_grads.rds")

abvtfs$rsums <- rowSums(abvtfs %>% select(-ID))
a <- abvtfs %>% select(-ID, -rsums)
a <- a/abvtfs$rsums
a$ID <- abvtfs$ID
gcmc_data <- gcmc_data %>% filter(ID %in% a$ID)
a <- a %>% filter(ID %in% gcmc_data$ID)
aa <-  a %>% select(-ID)
# sometimes gcmc_data are available but grid is not

# knock out columns and normalize
dd1 <- aa[,colSums(aa[, 1:ncol(aa)]) > 0.6]
dd1 <- dd1/rowSums(dd1)
heat_maps <- cbind(dd1, gcmc_data %>% select(ID, Uptake)) # may generate a duplicated colname
# depending on the size of the dataset, data split is going to have two options
if (nrow(gcmc_data) >= 2000){
  DATA_SPLIT <- 1000 # Number of data points used in training data, from setup_data.R
} else {
  DATA_SPLIT <- ceiling(0.5 * (nrow(gcmc_data)))
}
training_rows <- sample(length(heat_maps$ID), DATA_SPLIT)
training_data <- heat_maps[training_rows, ]
testing_data <- heat_maps[-training_rows, ]
trained_model <- training_data %>% select(-ID, -Uptake) %>% fit_glmnet(., training_data$Uptake, lambda = NULL, alpha = 1)

tested <- pred_glmnet(trained_model, testing_data %>% select(-ID, -Uptake))
postResample(tested, testing_data$Uptake)
qplot(x = testing_data$Uptake, y = tested) + geom_abline(slope = 1, intercept = 0)

trained_rf <- randomForest(x = training_data %>% select(-ID, -Uptake), y = training_data$Uptake, ntree = 500)
tested_rf <- predict(trained_rf, testing_data %>% select(-ID, Uptake))
postResample(tested_rf, testing_data$Uptake)
qplot(x = testing_data$Uptake, y = tested_rf) + geom_abline(slope = 1, intercept = 0)