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
 gcmc_file = "All_data/propane_10bar_uff.txt"
 if (!dir.exists("Results")){
   dir.create("Results")
 }
 hmof_h2_grid <- read_rds(grid_file)
 gcmc_data <- read_data(gcmc_file)
 # extract the molecule name, pressure and temperature for picture naming
 molecule_name <- sub(".*\\/", "", gcmc_file) %>% sub("\\_.*", "", .)
 temperature <- gcmc_data$Temp[1]
 pressure <- gcmc_data$Pres[1]
 string_to_paste <- paste(molecule_name, temperature, pressure, sep = "_")
 # create a directory for that condition and molecule
 save_path <- paste0(string_to_paste, "/", sep = "")
 save_path <- paste0("Results/", save_path)
 if (!dir.exists(save_path)){
   dir.create(save_path)
 }
 # here we can have two options, feed uniform or non-uniform bins
 #binbounds <- bounds_from_params(ch4_binspec)
 gcmc_data <- mutate(gcmc_data, id=ID)
 new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
 binbounds <- automatic_bins(new_grid) # for non-uniform bins
 binbounds <- as.data.frame(binbounds)
 binbounds <- binbounds %>% filter(upper <= 0.0) # bind all repulsive bins
 ch4_binspec <- list("from" = head(binbounds$lower, n = 1), "to" = tail(binbounds$upper, n = 1), bounds = binbounds)

 
 # depending on the size of the dataset, data split is going to have two options
 if (nrow(gcmc_data) >= 2000){
   DATA_SPLIT <- 1000 # Number of data points used in training data, from setup_data.R
 } else {
   DATA_SPLIT <- ceiling(0.5 * (nrow(gcmc_data)))
 }
 
p_ch4_sets <- partition_data_subsets(hmof_h2_grid, gcmc_data, DATA_SPLIT)
p_ch4_vol <- gcmc_data %>%
  mutate(g.L = Molec_cm3overcm3) %>%
  run_model_on_partitions(p_ch4_sets, ., ch4_binspec, plot_units="cm\u00B3/cm\u00B3", db_name = "tobacco")

gg_train <- p_ch4_vol$plots$parity_training%>% rescale_ch4_parity() + theme(axis.title.y = element_text(hjust=0.5)) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")
gg_test <- p_ch4_vol$plots$parity_testing%>% rescale_ch4_parity() + theme(axis.title.y = element_text(hjust=0.5)) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")
save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_train.png"), sep = ""),gg_train, base_width = 10, base_height = 8, dpi = 600)
save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_test.png"), sep = ""),gg_test, base_width = 10, base_height = 8, dpi = 600)
# export data for later testing other models, one can also export data to csv files
 train_data <- export_data(p_ch4_sets$training, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec)
 test_data <- export_data(p_ch4_sets$testing, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec)

# a histogram correlation analysis
  train_histo_bins <- train_data %>% select(-y_act)
  test_histo_bins <- test_data %>% select(-y_act)
  cor(as.matrix(t(train_histo_bins[1:2,])))
# Run a random forest model
 rf_model <- randomForest(x = train_data %>% select(-y_act), y = train_data$y_act, ntree = 500)
 # gg <- qplot(x = model$predicted, y = model$y) + geom_abline(slope = 1, intercept = 0)

# # plot error automatically if model is a randomforest
# #plot(model)
# # making predictions/ testing your fit, tell what variables are important
 #varplot <- varImpPlot(rf_model)
 #save_plot(paste("Results/", paste0(string_to_paste, "_Variable_Importance_Plot.png"), sep = ""),varplot, base_width = 10, base_height = 8, dpi = 600)
# 
 tested1 <- predict(rf_model, test_data)

# gg_test <- qplot(x = test_data$y_act, y = tested) + geom_abline(slope = 1, intercept = 0)
 rf_tested_rmse <- postResample(pred = tested1, obs = test_data$y_act) # original model
 gg_rf_train <- qplot(x = rf_model$y, y = rf_model$predicted) %>% rescale_ch4_parity() + geom_abline(slope = 1, intercept = 0, linetype = 2) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")+ geom_point(color='orange')
 save_plot(paste(save_path, paste0(string_to_paste, "_random_forest_train.png"), sep = ""),gg_rf_train, base_width = 10, base_height = 8, dpi = 600)
 gg_rf_test <- qplot(x = test_data$y_act, y = tested1) %>% rescale_ch4_parity() + geom_abline(slope = 1, intercept = 0, linetype = 2)+ xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")+ geom_point(color='lightblue')
 save_plot(paste(save_path, paste0(string_to_paste, "_random_forest_test.png"), sep = ""),gg_rf_test, base_width = 10, base_height = 8, dpi = 600)

# Part below are related to topology, see Paper by Yamil and Diego
# Read the names with topology
topologies <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is topology

# just keep all the structural properties that are numbers
structural_data <- orig_tobacco_data %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)

train_data_with_id <- export_data(p_ch4_sets$training, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec, with_id = TRUE)
test_data_with_id <- export_data(p_ch4_sets$testing, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec, with_id = TRUE)
# calculate the correlation matrix for these data
colnames(train_data_with_id)[colnames(train_data_with_id)=="id"] <- "MOF.ID"
colnames(test_data_with_id)[colnames(test_data_with_id)=="id"] <- "MOF.ID"
correlation_train <- correlation_of_energy_histograms(train_data_with_id, write_to_csv = TRUE)
correlation_test <- correlation_of_energy_histograms(test_data_with_id)
# extract a row of interest
# all_Rs <- subset(correlation_train, rownames(correlation_train) %in% "12616")
# df_all_Rs <- data.frame(t(all_Rs))
# df_all_Rs$y_actual <- train_data_with_id$y_act
# colnames(df_all_Rs)[colnames(df_all_Rs)=="X12616"] <- "R_scores"
# ordered_Rvalues <- df_all_Rs[order(-df_all_Rs$y_actual),]

# add R-score as a predictor
all_data_with_R <- rbind(train_data_with_id, test_data_with_id)
correlation_all <- correlation_of_energy_histograms(all_data_with_R)
all_Rs <- subset(correlation_all, rownames(correlation_all) %in% "13320")
df_all_Rs <- data.frame(t(all_Rs))
colnames(df_all_Rs)[colnames(df_all_Rs)=="X13320"] <- "R_scores"
all_data_with_R$Rscore <- df_all_Rs$R_scores
training_rows <- sample(length(all_data_with_R$MOF.ID), DATA_SPLIT)
training_data_with_R <- all_data_with_R[training_rows, ]
testing_data_with_R <- all_data_with_R[-training_rows, ]

rf_model_Rscore <- randomForest(x = training_data_with_R %>% select(-y_act, -MOF.ID), y = training_data_with_R$y_act, ntree = 700)
tested_Rscore <- predict(rf_model_Rscore, testing_data_with_R)
rf_tested_Rscore_rmse <- postResample(pred = tested_Rscore, obs = testing_data_with_R$y_act) 
gg_rf_train_Rscore <- qplot(x = rf_model_Rscore$y, y = rf_model_Rscore$predicted, label = training_data_with_R$MOF.ID) %>% rescale_ch4_parity() + geom_abline(slope = 1, intercept = 0, linetype = 2) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)") + geom_point(color='orange')
save_plot(paste(save_path, paste0(string_to_paste, "_random_forest_train_Rscore.png"), sep = ""),gg_rf_train_Rscore, base_width = 10, base_height = 8, dpi = 600)
gg_rf_test_Rscore <- qplot(x = testing_data_with_R$y_act, y = tested_Rscore) %>% rescale_ch4_parity() + geom_abline(slope = 1, intercept = 0, linetype = 2)+ xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)") + geom_point(color='lightblue')
save_plot(paste(save_path, paste0(string_to_paste, "_random_forest_test_Rscore.png"), sep = ""),gg_rf_test_Rscore, base_width = 10, base_height = 8, dpi = 600)

# learning using Topology data
topo_train_data<- merge(train_data_with_id, structural_data, by ="MOF.ID")
topo_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")

# perform a random forest model with these structural properties
rf_model_topo <- randomForest(x = topo_train_data %>% select(-y_act, -MOF.ID), y = topo_train_data$y_act, ntree = 500)
tested_topo <- predict(rf_model_topo, topo_test_data)
rf_tested_topo_rmse <- postResample(pred = tested_topo, obs = topo_test_data$y_act) 
gg_rf_train_topo <- qplot(x = rf_model_topo$y, y = rf_model_topo$predicted) %>% rescale_ch4_parity() + geom_abline(slope = 1, intercept = 0, linetype = 2) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)") + geom_point(color='orange')
save_plot(paste(save_path, paste0(string_to_paste, "_random_forest_train_topology.png"), sep = ""),gg_rf_train_topo, base_width = 10, base_height = 8, dpi = 600)
gg_rf_test_topo <- qplot(x = topo_test_data$y_act, y = tested_topo) %>% rescale_ch4_parity() + geom_abline(slope = 1, intercept = 0, linetype = 2)+ xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)") + geom_point(color='lightblue')
save_plot(paste(save_path, paste0(string_to_paste, "_random_forest_test_topology.png"), sep = ""),gg_rf_test_topo, base_width = 10, base_height = 8, dpi = 600)