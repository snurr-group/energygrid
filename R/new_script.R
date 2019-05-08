source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
library(plotly)
#source("R/load_data.R")
#source("R/read_tobacco_propane.R")
#source("R/read_tobacco_butane.R")
source("R/read_tobacco_new_propane.R")
#source("R/read_tobacco_ethane.R")
#source("R/read_tobacco_methane.R")
source("R/tobacco_data_for_zhao.R")

# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4

 ch4_binspec <- c(from=-26, to=0, step=0.5, width=0.5)
# p_ch4_grids <- raw_hmof_grids_ch4 %>% 
#   mutate(id = as.integer(str_sub(dirname, 6))) %>% 
#   filter(id %in% p_2bar_data$id)
# simplified_grid_ids_ch4 <- c(
#   simplified_grid_training,
#   sample(simplified_grid_testing, 2250-length(simplified_grid_training))
# )
 
 hmof_h2_grid <- read_rds("whatever.Rds")
 gcmc_data <- mutate(gcmc_data, id=ID)
 DATA_SPLIT <- 1000  # Number of data points used in training data, from setup_data.R
p_ch4_sets <- partition_data_subsets(hmof_h2_grid, gcmc_data, DATA_SPLIT)
p_ch4_vol <- gcmc_data %>%
  mutate(g.L = Molec_cm3overcm3) %>%
  #run_model_on_partitions(p_ch4_sets, ., ch4_binspec, plot_units="cm\u00B3/cm\u00B3", db_name = "tobacco")
  run_model_on_partitions(p_ch4_sets, ., ch4_binspec, plot_units="cm\u00B3/cm\u00B3", db_name = "tobacco")

p_ch4_vol$plots$parity_training%>% rescale_ch4_parity() + theme(axis.title.y = element_text(hjust=0.5)) + xlab("GCMC capacity (cm\u00B3/cm\u00B3)") + ylab("Predicted capacity (cm\u00B3/cm\u00B3)")
# 

# # #read the names with topology
# topologies <- read.table("fullnames.txt") # first column is ID, second column is topology
# x <- vector(mode="character", length = dim(p_ch4_vol$pred_df)[1])
# n1is_sym3 <- vector(mode="character", length = dim(p_ch4_vol$pred_df)[1])
# n2is_sym3 <- vector(mode="character", length = dim(p_ch4_vol$pred_df)[1])
# n1is_sym3mc0 <- vector(mode="logical", length = dim(p_ch4_vol$pred_df)[1])
# n2is_sym3mc0 <- vector(mode="logical", length = dim(p_ch4_vol$pred_df)[1])
# for (a in p_ch4_vol$pred_df$ID){
#   if (a %in% topologies$V1){
#     aa <-as.numeric(a) # row number in the tobacco database topology file
#     bb <- which(p_ch4_vol$pred_df$id == a) # row number in this dataframe
#     x[bb] <- as.character(topologies$V2[aa])
#     n1is_sym3[bb] <- orig_tobacco_data$n1.name[aa]
#     n2is_sym3[bb] <- orig_tobacco_data$n2.name[aa]
#     n1is_sym3mc0[bb] <- (n1is_sym3[bb] == "sym_3_mc_0")
#     n2is_sym3mc0[bb] <- (n2is_sym3[bb] == "sym_3_mc_0")
#     }
# }
# 
# 
# 
# 
# p_ch4_vol$pred_df$Topo <- x
# p_ch4_vol$pred_df$n1 <- n1is_sym3
# p_ch4_vol$pred_df$n2 <- n2is_sym3
# p_ch4_vol$pred_df$n1is_sym3mc0 <- n1is_sym3mc0
# p_ch4_vol$pred_df$n2is_sym3mc0 <- n2is_sym3mc0
# 
# 
# # extract variables out for datapoints in that cluster of outliers
# # delete all columns that doesn't have textual properties
# orig_tobacco_data <- orig_tobacco_data[ -c(7:33)]
# colnames(p_ch4_vol$pred_df)[colnames(p_ch4_vol$pred_df)=="id"] <- "MOF.ID"
# 
# merged_data<- merge(p_ch4_vol$pred_df, orig_tobacco_data, by ="MOF.ID")
# merged_data$textual <- merged_data$vf * merged_data$vsa * merged_data$gsa * merged_data$pld * merged_data$lcd/(mean(merged_data$vf) * mean(merged_data$vsa) * mean(merged_data$gsa)* mean(merged_data$pld) * mean(merged_data$lcd))
# #gg <- merged_data %>% filter(Topo == "srsb") %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred, size = textual)) + geom_point(aes(text=MOF.ID))
#   gg <- merged_data %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred, color = n2.ID)) + geom_point(aes(text=MOF.ID)) + scale_colour_gradientn(colours=rainbow(4))
# 
# gg <- ggplotly(gg)
# gg
# 
# # clean up the columns that have same stuff
# 
# 
# gg_heat <- merged_data %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred, color = gsa)) + geom_point(aes(text=MOF.ID)) + scale_colour_gradientn(colours=rainbow(4))
# gg_heat
# #gg <- p_ch4_vol$pred_df %>% filter(Topo == "srsb") %>% ggplot(aes(x = Molec_cm3overcm3, y = y_pred)) + geom_point(aes(text=ID))
# #gg <- ggplotly(gg)
# # #gg
# 
# 
