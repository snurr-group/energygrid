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
gcmc_file = "All_data/Hexane_.1bar_new.txt"
if (!dir.exists("Results")){
  dir.create("Results")
}
hmof_h2_grid <- read_rds(grid_file)
gcmc_data <- read_data(gcmc_file)
# extract the molecule name for picture naming
molecule_name <- sub(".*\\/", "", gcmc_file) %>% sub("\\_.*", "", .)
# here we can have two options, feed uniform or non-uniform bins
#binbounds <- bounds_from_params(ch4_binspec)
gcmc_data <- mutate(gcmc_data, id=ID)
new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
binbounds <- automatic_bins(hmof_h2_grid) # for non-uniform bins
# further lump all repulsive bins together
binbounds <- as.data.frame(binbounds)
binbounds <- binbounds %>% filter(upper <= 0.0)
# get the names of files from that path
namess <- list.files(path="Energies/") %>% str_remove(., ".grid")
gcmc_data <- gcmc_data %>% filter(ID %in% namess)
#source("R/get_heat_map.R")
asdasd <- readRDS("All_data/heat_map_UFF_grads.rds")
a <- asdasd
gcmc_data <- gcmc_data %>% filter(ID %in% a$ID)
a <- a %>% filter(ID %in% gcmc_data$ID)
aa <-  a %>% select(-ID)
# sometimes gcmc_data are available but grid is not

# knock out columns and normalize
dd1 <- aa[,colSums(aa[, 1:ncol(aa)]) > 0.3]
dd1 <- dd1/rowSums(dd1)
dd1$ID <- a$ID
dd1 <- dd1[order(match(dd1$ID, gcmc_data$ID)),] # without this line, it won't be right
heat_maps <- cbind(dd1, gcmc_data %>% select(ID, Molec_cm3overcm3)) # may generate a duplicated colname
# remove the duplicated column name
heat_maps <-heat_maps[, !duplicated(colnames(heat_maps))]

DATA_SPLIT = 1000
training_rows <- sample(length(heat_maps$ID), DATA_SPLIT)
training_data <- heat_maps[training_rows, ]
testing_data <- heat_maps[-training_rows, ]
trained_model <- training_data %>% select(-ID, -Molec_cm3overcm3) %>% fit_glmnet(., training_data$Molec_cm3overcm3, lambda = NULL, alpha = 1)

tested <- pred_glmnet(trained_model, testing_data %>% select(-ID, -Molec_cm3overcm3))
postResample(tested, testing_data$Molec_cm3overcm3)
qplot(x = testing_data$Molec_cm3overcm3, y = tested) + geom_abline(slope = 1, intercept = 0) + scale_x_continuous(limits=c(0, 200))+ scale_y_continuous(limits=c(0, 200))

trained_rf <- randomForest(x = training_data %>% select(-ID, -Molec_cm3overcm3), y = training_data$Molec_cm3overcm3, ntree = 500)
tested_rf <- predict(trained_rf, testing_data %>% select(-ID, Molec_cm3overcm3))
postResample(tested_rf, testing_data$Molec_cm3overcm3)
qplot(x = testing_data$Molec_cm3overcm3, y = tested_rf) + geom_abline(slope = 1, intercept = 0)
