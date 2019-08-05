rm(list=ls())
source("R/refined_bins_calc.R")
source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/save_train_test_data.R")
source("R/read_tobacco_data.R")
source("R/package_verification.R")
gcmc_file <- "All_data/XeKr_mix_273_1bar.txt"
gcmc_data <- read_table2(gcmc_file)

gcmc_data <- na.omit(gcmc_data)
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
gcmc_data$Selectivity <- (gcmc_data$Xe_uptake/0.2)/(gcmc_data$Kr_uptake/0.8)
gg_histo_selectivity <- qplot(gcmc_data$Selectivity , geom = "histogram") + 
  xlab("GCMC Selectivity") + ylab("Counts")

save_plot(paste(save_path, paste0(string_to_paste, "_Selectivity_histogram.png"), sep = ""), 
          gg_histo_selectivity, base_width = 10, base_height = 8, dpi = 600)

# other topology stuff: histograms for the population
topologies <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is topology

# just keep all the structural properties that are numbers
structural_data <- orig_tobacco_data %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)
gcmc_with_MOFID <- gcmc_data
colnames(gcmc_with_MOFID)[colnames(gcmc_with_MOFID)=="ID"] <- "MOF.ID"
topology_data <- merge(gcmc_with_MOFID, structural_data, by ="MOF.ID")
make_topology_histograms(condition_name = string_to_paste, topo_data = topology_data)

remove_outliers <- TRUE
focus <- FALSE
limits <- 200
if (remove_outliers){
# get rid of the outliers using 1.5 IQR method, IQR stands for interquantile range: range between 25% and 75%
iqr_factor = 6
iqr = IQR(gcmc_data$Selectivity)
upper_threshold = as.numeric((iqr * iqr_factor) + quantile(gcmc_data$Selectivity)[2]) # the second value is the 25% one
lower_threshold = as.numeric(quantile(gcmc_data$Selectivity)[4] - (iqr * iqr_factor)) # the second value is the 75% one
# then filter out 
gcmc_data <- gcmc_data %>% filter((Selectivity >= lower_threshold) & (Selectivity <= upper_threshold))
string_to_paste <- paste0(string_to_paste, "_without_outliers")
limits <- 30
} else if(focus){
  limits <- 30
  string_to_paste <- paste0(string_to_paste, "less_range")
}
#gcmc_data <- gcmc_data %>% filter(Selectivity < 50)
# convert IDs to characters
gcmc_data$ID <- as.character(gcmc_data$ID)

# 
# h2_types <- paste0("y.h2.", c("g.L", "mol.kg", "wtp"))  # Prefix data with its source (Yamil, Scotty, etc.)
# orig_tobacco_data <- read_xlsx(
#   "Data/CrystGrowthDesign_SI.xlsx",
#   sheet = "data",
#   skip = 3, na = "inf",
#   col_names = c(
#     "MOF.ID",
#     "vf", "vsa", "gsa", "pld", "lcd",
#     paste0(h2_types, ".100.77"),
#     paste0(h2_types, ".100.130"),
#     paste0(h2_types, ".100.200"),
#     paste0(h2_types, ".100.243"),
#     paste0("h2.qst.6.", c(77, 130, 200, 243)),
#     paste0("y.ch4.", rep(c("v.v.", "mg.g."), 2), c("100.298", "100.298", "65.298", "65.298")),
#     "ch4.qst.6.298",
#     paste0("y.xe.kr.1.", c("xe", "kr", "select")),
#     paste0("y.xe.kr.5.", c("xe", "kr", "select")),
#     "topology",
#     paste0("n1.", c("sym", "character", "ID")),
#     paste0("n2.", c("sym", "character", "ID")),
#     "cbb.ID"
#   )
# )
# filtered <-  orig_tobacco_data[orig_tobacco_data$MOF.ID %in% gcmc_data$ID,]
# colnames(filtered)[which(names(filtered) == "MOF.ID")] <- "id"
# filtered$id <- as.character(filtered$id)
# more_select <- inner_join(filtered, gcmc_data, by = "id")
# exported <- more_select %>% select(id, Selectivity, `y.xe.kr.1.select`)
# write.xlsx(exported, "XeKr_selectivity_compare.xlsx")
#tobacco_codes <- read_table2("Data/mofs_map.dat", col_names = c("MOF.ID", "python.id"), col_types="ic")


hmof_h2_grid <- read_rds("Kr_greater_range.Rds")
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
  mutate(g.L = Selectivity) %>%
  run_model_on_partitions(p_ch4_sets, ., ch4_binspec, plot_units="", db_name = "tobacco")

gg_train <- p_ch4_vol$plots$parity_training + 
  scale_x_continuous(limits = c(0,limits)) + scale_y_continuous(limits = c(0,limits)) + 
  xlab("GCMC Selectivity") + 
  ylab("Predicted Selectivity")

save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_train.png"), sep = ""), 
          gg_train, base_width = 10, base_height = 8, dpi = 600)
gg_test <- p_ch4_vol$plots$parity_testing + 
  scale_x_continuous(limits = c(0,limits)) + scale_y_continuous(limits = c(0,limits)) + 
  xlab("GCMC selectivity") + 
  ylab("Predicted selectivity")
save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_test.png"), sep = ""), 
          gg_test, base_width = 10, base_height = 8, dpi = 600)

train_data <- export_data(p_ch4_sets$training, 
                          rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), ch4_binspec)
test_data <- export_data(p_ch4_sets$testing, 
                         rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), ch4_binspec)

model <- randomForest(x = train_data %>% select(-y_act), y = train_data$y_act, ntree = 500)

make_rf_prediction_plots(condition_name = string_to_paste, 
                         plot_name = "rf", 
                         rf_model= model,  test_data = test_data, lim = limits, axis_label = "Selectivity")


# get the histograms and topology, glued together
train_data_with_id <- export_data(p_ch4_sets$training, 
                          rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), 
                          ch4_binspec, with_id = TRUE)
test_data_with_id <- export_data(p_ch4_sets$testing, 
                         rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), 
                         ch4_binspec, with_id = TRUE)
colnames(train_data_with_id)[colnames(train_data_with_id)=="id"] <- "MOF.ID"
colnames(test_data_with_id)[colnames(test_data_with_id)=="id"] <- "MOF.ID"
topo_train_data <- merge(train_data_with_id, structural_data, by ="MOF.ID")
topo_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")
