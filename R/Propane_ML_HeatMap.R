source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
library(plotly)
library(randomForest)
source("R/read_tobacco_new_propane.R")
source("R/tobacco_data_for_zhao.R")
source("R/save_train_test_data.R")
source("R/refined_bins_calc.R")
# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4

#ch4_binspec <- c(from=-26 , to=20, step=0.5, width=0.5)

hmof_h2_grid <- read_rds("All_data/Grids_Propane_298K_1A_GenericMOFs.Rds")

# here we can have two options, feed uniform or non-uniform bins
#binbounds <- bounds_from_params(ch4_binspec)
gcmc_data <- mutate(gcmc_data, id=ID)
new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
binbounds <- automatic_bins(new_grid,1000000) # for non-uniform bins
# get the names of files from that path
namess <- list.files(path="Energies/") %>% str_remove(., ".grid")
gcmc_data <- gcmc_data %>% filter(ID %in% namess)
source("R/get_heat_map.R")
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
abvtfs <- readRDS("heat_map.rds")
abvtfs$rsums <- rowSums(abvtfs %>% select(-ID))
a <- abvtfs %>% select(-ID, -rsums)
a <- a/abvtfs$rsums
a$ID <- abvtfs$ID

heat_maps <- cbind(a, gcmc_data %>% select(ID, Molec_cm3overcm3)) # may generate a duplicated colname
# remove the duplicated column name
heat_maps <-heat_maps[, !duplicated(colnames(heat_maps))]

DATA_SPLIT = 1000
training_rows <- sample(length(heat_maps$ID), DATA_SPLIT)
training_data <- heat_maps[training_rows, ]
testing_data <- heat_maps[-training_rows, ]
trained_model <- training_data %>% select(-ID, -Molec_cm3overcm3) %>% fit_glmnet(., training_data$Molec_cm3overcm3, lambda = NULL, alpha = 1)

tested <- pred_glmnet(trained_model, testing_data %>% select(-ID, -Molec_cm3overcm3))
postResample(tested, testing_data$Molec_cm3overcm3)
qplot(x = testing_data$Molec_cm3overcm3, y = tested_rf) + geom_abline(slope = 1, intercept = 0)

trained_rf <- randomForest(x = training_data %>% select(-ID, -Molec_cm3overcm3), y = training_data$Molec_cm3overcm3, ntree = 500)
tested_rf <- predict(trained_rf, testing_data %>% select(-ID, Molec_cm3overcm3))
postResample(tested_rf, testing_data$Molec_cm3overcm3)
qplot(x = testing_data$Molec_cm3overcm3, y = tested_rf) + geom_abline(slope = 1, intercept = 0)
