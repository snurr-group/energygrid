rm(list=ls())

library(plotly)
library(randomForest)
library(cowplot)
source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/read_gcmc_data.R")
source("R/tobacco_data_for_zhao.R")
source("R/save_train_test_data.R")
source("R/refined_bins_calc.R")
# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4
 # define the input files
 grid_file = "All_data/CH3_probe_1A_UFF.rds"
 gcmc_file = "All_data/ethane_10bar_UFF.txt"
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
 binbounds <- automatic_bins(new_grid) # for non-uniform bins
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
save_plot(paste("Results/", paste0(molecule_name, "_train.png"), sep = ""),gg_train, base_width = 10, base_height = 8, dpi = 600)
save_plot(paste("Results/", paste0(molecule_name, "_test.png"), sep = ""),gg_test, base_width = 10, base_height = 8, dpi = 600)
# export data for later testing other models, one can also export data to csv files
 train_data <- export_data(p_ch4_sets$training, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec)
 test_data <- export_data(p_ch4_sets$testing, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), ch4_binspec)


# ORIGINAL MODEL
 rf_model <- randomForest(x = train_data %>% select(-y_act), y = train_data$y_act, ntree = 500)
 # gg <- qplot(x = model$predicted, y = model$y) + geom_abline(slope = 1, intercept = 0)

# # plot error automatically if model is a randomforest
# #plot(model)
# # making predictions/ testing your fit, tell what variables are important
 #varplot <- varImpPlot(rf_model)
 #save_plot(paste("Results/", paste0(molecule_name, "_Variable_Importance_Plot.png"), sep = ""),varplot, base_width = 10, base_height = 8, dpi = 600)
# 
 tested1 <- predict(rf_model, test_data)

# gg_test <- qplot(x = test_data$y_act, y = tested) + geom_abline(slope = 1, intercept = 0)
 rf_tested_rmse <- postResample(pred = tested1, obs = test_data$y_act) # original model
 gg_rf_train <- qplot(x = rf_model$y, y = rf_model$predicted) + geom_abline(slope = 1, intercept = 0, linetype = 2) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")
 save_plot(paste("Results/", paste0(molecule_name, "_random_forest_train.png"), sep = ""),gg_rf_train, base_width = 10, base_height = 8, dpi = 600)
 gg_rf_test <- qplot(x = test_data$y_act, y = tested1) + geom_abline(slope = 1, intercept = 0, linetype = 2)+ xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")
 save_plot(paste("Results/", paste0(molecule_name, "_random_forest_test.png"), sep = ""),gg_rf_test, base_width = 10, base_height = 8, dpi = 600)
 # # #read the names with topology
topologies <- read.table("All_data/fullnames.txt") # first column is ID, second column is topology
x <- vector(mode="character", length = dim(p_ch4_vol$pred_df)[1])
n1is_sym3 <- vector(mode="character", length = dim(p_ch4_vol$pred_df)[1])
n2is_sym3 <- vector(mode="character", length = dim(p_ch4_vol$pred_df)[1])
n1is_sym3mc0 <- vector(mode="logical", length = dim(p_ch4_vol$pred_df)[1])
n2is_sym3mc0 <- vector(mode="logical", length = dim(p_ch4_vol$pred_df)[1])
for (a in p_ch4_vol$pred_df$ID){
  if (a %in% topologies$V1){
    aa <-as.numeric(a) # row number in the tobacco database topology file
    bb <- which(p_ch4_vol$pred_df$id == a) # row number in this dataframe
    x[bb] <- as.character(topologies$V2[aa])
    n1is_sym3[bb] <- orig_tobacco_data$n1.name[aa]
    n2is_sym3[bb] <- orig_tobacco_data$n2.name[aa]
    n1is_sym3mc0[bb] <- (n1is_sym3[bb] == "sym_3_mc_0")
    n2is_sym3mc0[bb] <- (n2is_sym3[bb] == "sym_3_mc_0")
    }
}




p_ch4_vol$pred_df$Topo <- x
p_ch4_vol$pred_df$n1 <- n1is_sym3
p_ch4_vol$pred_df$n2 <- n2is_sym3
p_ch4_vol$pred_df$n1is_sym3mc0 <- n1is_sym3mc0
p_ch4_vol$pred_df$n2is_sym3mc0 <- n2is_sym3mc0


# extract variables out for datapoints in that cluster of outliers
# delete all columns that doesn't have textual properties
orig_tobacco_data <- orig_tobacco_data[ -c(7:33)]
colnames(p_ch4_vol$pred_df)[colnames(p_ch4_vol$pred_df)=="id"] <- "MOF.ID"

merged_data<- merge(p_ch4_vol$pred_df, orig_tobacco_data, by ="MOF.ID")
merged_data$textual <- merged_data$vf * merged_data$vsa * merged_data$gsa * merged_data$pld * merged_data$lcd/(mean(merged_data$vf) * mean(merged_data$vsa) * mean(merged_data$gsa)* mean(merged_data$pld) * mean(merged_data$lcd))
#gg <- merged_data %>% filter(Topo == "srsb") %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred, size = textual)) + geom_point(aes(text=MOF.ID))
  gg <- merged_data %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred, color = n2.ID)) + geom_point(aes(text=MOF.ID)) + scale_colour_gradientn(colours=rainbow(4))

gg <- ggplotly(gg)
gg

# clean up the columns that have same stuff


gg_heat <- merged_data %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred, color = gsa)) + geom_point(aes(text=MOF.ID)) + scale_colour_gradientn(colours=rainbow(4))
gg_heat
#gg <- p_ch4_vol$pred_df %>% filter(Topo == "srsb") %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred)) + geom_point(aes(text=ID))
#gg <- ggplotly(gg)
# #gg

colnames(merged_data)[colnames(merged_data) == "Host-Guest"] <- "hg"
colnames(merged_data)[colnames(merged_data) == "Guest-Guest"] <- "gg"

# x <- c(2,6,7,8,9,10,11,12,21,22,23,24,25)
# for (valx in x) {
#   for (valy in x) {
#     if (valy > valx){
#       for (valz in x) {
#         if ((valz != valx) & (valz != valy)) {
#           gg_heat <- merged_data %>% ggplot(aes_string(x = names(merged_data[valx]), y = names(merged_data[valy]), color = names(merged_data[valz]))) + geom_point() + scale_colour_gradientn(colours=rainbow(4))
#           ggsave(gg_heat, filename=paste("valx_", valx,"_valy_",valy, "_valz_", valz,".png", sep=""))
#         }
#       }
#     }
#   }
# }



# Idea: train the model for predicting gg interaction based on topology data
export_topology_data(p_ch4_sets$training, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), orig_tobacco_data, "tobacco_train_topo.csv")
export_topology_data(p_ch4_sets$testing, rename(gcmc_data %>% mutate(g.L = Molec_cm3overcm3), y_act=g.L), orig_tobacco_data, "tobacco_test_topo.csv")

topo_train_data <- read_csv("tobacco_train_topo.csv")
topo_test_data <- read_csv("tobacco_test_topo.csv")

# let's think about it in another way: let's predict guest-guest interaction energy!
model_gg <- randomForest(x = topo_train_data %>% select(-`Guest-Guest`), y = topo_train_data$`Guest-Guest`, ntree = 500)
tested_gg <- predict(model_gg,topo_test_data)
tested_rmse_gg <- postResample(pred = tested_gg, obs = topo_test_data$`Guest-Guest`)
# not working: prediction with just topology yields R2 ~ 0.85, not good
# let's do with energy grid
model_gg_grid <- randomForest(x = train_data %>% select(-y_act, -`Host-Guest`, -`Guest-Guest`, -Heat_of_Ads, -Qst_low_loading), y = train_data$`Guest-Guest`, ntree = 500)
tested_gg_grid <- predict(model_gg_grid,test_data)
trained_gg_grid <- predict(model_gg_grid,train_data)
tested_rmse_gg_grid <- postResample(pred = tested_gg_grid, obs = test_data$`Guest-Guest`)

# seems ok though: using energygrid to predict gg interaction gives ~ 0.95 R2
# then can we use the predicted gg interaction to train the model that predicts uptake???
test_data$pred_gg <- tested_gg_grid
train_data$pred_gg <- trained_gg_grid
model_new <- randomForest(x = train_data %>% select(-y_act, -`Host-Guest`, -`Guest-Guest`, -Heat_of_Ads), y = train_data$y_act, ntree = 500)
tested_new <- predict(model_new, test_data)
tested_rmse_new <- postResample(pred = tested_new, obs = test_data$y_act)


model_original<- randomForest(x = train_data %>% select(-y_act, -`Host-Guest`, -`Guest-Guest`, -Heat_of_Ads, -Qst_low_loading, -pred_gg), y = train_data$y_act, ntree = 500)
tested_original <- predict(model_original, test_data)
tested_rmse_original <- postResample(pred = tested_original, obs = test_data$y_act)
