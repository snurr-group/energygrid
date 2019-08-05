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
 gcmc_file = "All_data/propane_10bar.txt"
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
 ch4_binspec <- list("from" = head(binbounds$lower, n = 1), 
                     "to" = tail(binbounds$upper, n = 1), 
                     bounds = binbounds)

 
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

gg_train <- p_ch4_vol$plots$parity_training%>% rescale_ch4_parity() + 
  theme(axis.title.y = element_text(hjust=0.5)) + 
  xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + 
  ylab("Predicted capacity (cm\u00B3/cm\u00B3)")

gg_test <- p_ch4_vol$plots$parity_testing%>% rescale_ch4_parity() + 
  theme(axis.title.y = element_text(hjust=0.5)) + 
  xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + 
  ylab("Predicted capacity (cm\u00B3/cm\u00B3)")

save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_train.png"), sep = ""), 
          gg_train, base_width = 10, base_height = 8, dpi = 600)
save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_test.png"), sep = ""), 
          gg_test, base_width = 10, base_height = 8, dpi = 600)
# export data for later testing other models, one can also export data to csv files
 train_data <- export_data(p_ch4_sets$training, 
                           rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec)
 test_data <- export_data(p_ch4_sets$testing, 
                          rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec)

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
 make_rf_prediction_plots(condition_name = string_to_paste, 
                          plot_name = "rf", 
                          rf_model= rf_model,  test_data = test_data)

# Part below are related to topology, see Paper by Yamil and Diego
# Read the names with topology
topologies <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is topology

# just keep all the structural properties that are numbers
structural_data <- orig_tobacco_data %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)

train_data_with_id <- export_data(p_ch4_sets$training, 
                                  rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), 
                                  ch4_binspec, with_id = TRUE)
test_data_with_id <- export_data(p_ch4_sets$testing, 
                                 rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), 
                                 ch4_binspec, with_id = TRUE)
# calculate the correlation matrix for these data
colnames(train_data_with_id)[colnames(train_data_with_id)=="id"] <- "MOF.ID"
colnames(test_data_with_id)[colnames(test_data_with_id)=="id"] <- "MOF.ID"
correlation_train <- correlation_of_energy_histograms(train_data_with_id, write_to_csv = TRUE)
correlation_test <- correlation_of_energy_histograms(test_data_with_id)

# add R-score as a predictor
all_data_with_R <- rbind(train_data_with_id, test_data_with_id)
correlation_all <- correlation_of_energy_histograms(all_data_with_R)
# randomly select an ID as the reference for R-scores
random_ID <- sample(all_data_with_R$MOF.ID, size = 1)
all_Rs <- subset(correlation_all, rownames(correlation_all) %in% random_ID)
df_all_Rs <- data.frame(t(all_Rs))
colnames(df_all_Rs)[colnames(df_all_Rs)==paste0("X", random_ID)] <- "R_scores"
all_data_with_R$Rscore <- df_all_Rs$R_scores
training_rows <- sample(length(all_data_with_R$MOF.ID), DATA_SPLIT)
training_data_with_R <- all_data_with_R[training_rows, ]
testing_data_with_R <- all_data_with_R[-training_rows, ]
# show that Rscore and uptake are related
#Rscore_compare <- all_data_with_R %>% select(MOF.ID, y_act, Rscore)
#Rscore_2 <- Rscore_compare[order(-Rscore_compare$Rscore),]


rf_model_Rscore <- randomForest(x = training_data_with_R %>% select(-y_act, -MOF.ID),
                                y = training_data_with_R$y_act, ntree = 500)
make_rf_prediction_plots(condition_name = string_to_paste, 
                         plot_name = "rf_histogram_Rscore", 
                         rf_model= rf_model_Rscore,  test_data = testing_data_with_R)

# what about just using R value?
rf_model_just_Rscore <- randomForest(x = training_data_with_R %>% select(Rscore), 
                                     y = training_data_with_R$y_act, ntree = 500)
make_rf_prediction_plots(condition_name = string_to_paste, 
                         plot_name = "rf_just_Rscore", 
                         rf_model= rf_model_just_Rscore,  test_data = testing_data_with_R)

# learning using Topology data
topo_train_data<- merge(train_data_with_id, structural_data, by ="MOF.ID")
topo_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")

# perform a random forest model with these structural properties
rf_model_topo <- randomForest(x = topo_train_data %>% select(-y_act, -MOF.ID), 
                              y = topo_train_data$y_act, ntree = 500)
make_rf_prediction_plots(condition_name = string_to_paste, 
                         plot_name = "rf_histogram_topology", 
                         rf_model= rf_model_topo,  test_data = topo_test_data)

# generate topology histograms
topology_data <- rbind(topo_train_data, topo_test_data)
make_topology_histograms(condition_name = string_to_paste, topo_data = topology_data)
