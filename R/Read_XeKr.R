rm(list=ls())
source("R/refined_bins_calc.R")
source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/save_train_test_data.R")
source("R/read_textural_data.R")
source("R/package_verification.R")
source("R/read_gcmc_data.R")
set.seed(12345) # seems ok for both 1bar and 10bar

unit_for_plot <<- "" # Selectivity is unitless
poster <<- TRUE # save this as a global variable: no need to pass it around
XeKr <<- TRUE
# poster plots should have big dots
if (poster){
  dot_size <<- 3
}else{
  dot_size <<- 1 
}

grid_file <- "All_data/Kr_1A_greater_range.rds"

molecule_name <- "XeKr_Mix"
Temperature <- "273K"
Pressure <- "1Bar"
chosen_unit <- "Molec_cm3overcm3"
gcmc_data <- read_data(gcmc_file, read_SI = TRUE, sheetname = paste(molecule_name, Temperature, Pressure, sep = "_"), 
                       unit_of_ads = chosen_unit, 
                       relax = TRUE, just2k = TRUE, no_low_loading = FALSE)
# extract the probe name, molecule name, pressure and temperature for picture naming
# the naming convention of rds file should be: "probe"_"size of grid"_whatever.rds
if (str_count(grid_file, "_") > 1) {
  probe_name <- sub(".*\\/", "", grid_file) %>% sub("\\_.*", "", .)
  grid_size <- sub(".*\\/", "", grid_file) %>% 
    str_match(., paste0(probe_name, "_(.*?)_"))
  grid_size <- grid_size[2]
  grid_info <- paste(probe_name, grid_size, sep = "_")
} else{
  grid_info <- sub(".*\\/", "", grid_file) %>% sub("\\.rds", "", .)
}
# extract the molecule name, pressure and temperature for picture naming

string_to_paste <- paste(molecule_name, Temperature, Pressure, grid_info, sep = "_")
# create a directory for that condition and molecule
save_path <- paste0(string_to_paste, "/", sep = "")
# if save for poster figures, put poster in the front
if(poster){
  save_path <- paste0("Poster_", save_path)
}
save_path <- paste0("Results/", save_path)
if (!dir.exists(save_path)){
  dir.create(save_path)
}
gcmc_data$Selectivity <- (gcmc_data$Xe_uptake/0.2)/(gcmc_data$Kr_uptake/0.8)
# intentionally filter out that outlier > 1M
#gcmc_data <- gcmc_data %>% filter(Selectivity < 1000)
gg_histo_selectivity <- qplot(gcmc_data$Selectivity , geom = "histogram") + 
  xlab("GCMC Selectivity") + ylab("Counts")

save_plot(paste(save_path, paste0(string_to_paste, "_Selectivity_histogram.png"), sep = ""), 
          gg_histo_selectivity, base_width = 10, base_height = 10, dpi = 600)

# other topology stuff: histograms for the population
topologies <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is topology

# just keep all the structural properties that are numbers
structural_data <- read_textual_data(option = "ToBaCCo")
gcmc_with_MOFID <- gcmc_data
colnames(gcmc_with_MOFID)[colnames(gcmc_with_MOFID)=="ID"] <- "MOF.ID"
topology_data <- merge(gcmc_with_MOFID, structural_data, by ="MOF.ID")
#make_topology_histograms(condition_name = string_to_paste, topo_data = topology_data)

remove_outliers <- FALSE
focus <- FALSE

if (remove_outliers)
{
  # get rid of the outliers using 1.5 IQR method, IQR stands for interquantile range: range between 25% and 75%
  iqr_factor = 30
  iqr = IQR(gcmc_data$Selectivity)
  upper_threshold = as.numeric((iqr * iqr_factor) + quantile(gcmc_data$Selectivity)[2]) # the second value is the 25% one
  lower_threshold = as.numeric(quantile(gcmc_data$Selectivity)[4] - (iqr * iqr_factor)) # the second value is the 75% one
  # then filter out 
  gcmc_data <- gcmc_data %>% filter((Selectivity >= lower_threshold) & (Selectivity <= upper_threshold))
  string_to_paste <- paste0(string_to_paste, "_without_outliers")
  plot_limit <- max(gcmc_data$Selectivity)
  # round to the nearest 10
  plot_limit <- round(plot_limit, digits = -1) + 10
  } else if(focus){
    plot_limit <- 30
    string_to_paste <- paste0(string_to_paste, "less_range")
}

# we can do another filtering: based on loading
# filter out those who has a really low loading
filter_absolute_loading <- FALSE
absolute_loading_threshold <- 0.1
if (filter_absolute_loading)
{
  gcmc_data <- gcmc_data %>% filter(Kr_uptake > absolute_loading_threshold & Xe_uptake > absolute_loading_threshold)
}
plot_limit <- max(gcmc_data$Selectivity)
plot_limit <- plot_limit-mod(plot_limit,50)+50
previous_plot_lim <- plot_limit

#gcmc_data <- gcmc_data %>% filter(Selectivity < 50)
# convert IDs to characters
gcmc_data$ID <- as.character(gcmc_data$ID)

if (nrow(gcmc_data) > 2000)
{
  differ <- nrow(gcmc_data) - 2000
  rows_to_delete <- sample(1:nrow(gcmc_data), differ)
  gcmc_data <- gcmc_data[-rows_to_delete, ]
}
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


hmof_h2_grid <- read_rds(grid_file)
gcmc_data <- mutate(gcmc_data, id=ID)
new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
binbounds <- automatic_bins(new_grid, threshold_ratio = 0.001) # for non-uniform bins
binbounds <- as.data.frame(binbounds)
#binbounds <- binbounds %>% filter(upper <= 0.0) # bind all repulsive bins
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
  run_model_on_partitions(p_ch4_sets, ., ch4_binspec, plot_units="", db_name = "points")

gg_train <- p_ch4_vol$plots$parity_training %>% rescale_ch4_parity(., lims=c(0,plot_limit)) + 
  xlab("GCMC Selectivity") + 
  ylab("Predicted Selectivity")
gg_test <- p_ch4_vol$plots$parity_testing %>% rescale_ch4_parity(., lims=c(0,plot_limit)) + 
  xlab("GCMC selectivity") + 
  ylab("Predicted selectivity")
if (poster)
{
  gg_train <- gg_train + theme(axis.text=element_text(size=30),
                               axis.title=element_text(size=30,face="bold"))
  gg_test <- gg_test + theme(axis.text=element_text(size=30),
                             axis.title=element_text(size=30,face="bold")) 
}
save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_train.png"), sep = ""), 
          gg_train, base_width = 10, base_height = 10, dpi = 600)
save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_test.png"), sep = ""), 
          gg_test, base_width = 10, base_height = 10, dpi = 600)

train_data <- export_data(p_ch4_sets$training, 
                          rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), ch4_binspec)
test_data <- export_data(p_ch4_sets$testing, 
                         rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), ch4_binspec)

model <- randomForest(x = train_data %>% select(-y_act), y = train_data$y_act, ntree = 500)

make_rf_prediction_plots_XeKr(condition_name = string_to_paste, 
                         plot_name = "rf", 
                         rf_model= model,  test_data = test_data, lim = plot_limit, axis_label = "Selectivity")


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

rf_model_topo <- randomForest(x = topo_train_data %>% select(-y_act, -MOF.ID), 
                              y = topo_train_data$y_act, ntree = 500)
make_rf_prediction_plots_XeKr(condition_name = string_to_paste, 
                         plot_name = "rf_histogram_topology", 
                         rf_model= rf_model_topo,  test_data = topo_test_data, axis_label = "Selectivity")

abc <- rbind(topo_train_data, topo_test_data)
# save train and test data
train_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "train_Select.csv"), sep = ""))
test_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "test_Select.csv"), sep = ""))
model %>% saveRDS(., paste(save_path, paste0(string_to_paste, "normal_RF.rds"), sep = ""))
rf_model_topo %>% saveRDS(., paste(save_path, paste0(string_to_paste, "topo_RF.rds"), sep = ""))
# # filter out those mofs with small lcd (between 4.1 and 7)
# # suitable_mofs <- structural_data %>% filter(lcd > 4.1 & lcd < 7)
# # unscreened_suitable_mofs <- suitable_mofs %>% filter(!MOF.ID %in% gcmc_data$ID)
# # mof_xyz <- read_table2("tob_perm_unitcells_no_mc_0.txt")
# # mof_xyz$ID <- as.character(mof_xyz$ID)
# # mof_xyz <- mof_xyz %>% filter(ID %in% unscreened_suitable_mofs$MOF.ID)
# # # convert IDs back to numbers: if write as characters, it will have quotes!
# # mof_xyz$ID <- as.integer(mof_xyz$ID)
# # write.table(mof_xyz, "should_screen_small_lcd_mofs.txt", row.names = FALSE, col.names = FALSE, sep = " ")
# # 
# # # look at how topology data, especially lcd are. Pick selectivity = 20 as threshold
# # topo_data <- rbind(topo_train_data, topo_test_data)
# # outliers <- topo_data %>% filter(y_act > 20)
# 
# # # let's try something special: try changing weights of data points of random forest
# # rf_model_topo_weights <- randomForest(x = topo_train_data %>% select(-y_act, -MOF.ID), 
# #                               y = topo_train_data$y_act, ntree = 500, mtry = 27)
# # make_rf_prediction_plots(condition_name = string_to_paste, 
# #                          plot_name = "rf_histogram_topology_mtry27", 
# #                          rf_model= rf_model_topo_weights,  test_data = topo_test_data)
# # # tuning rf model
# # stay_tuned <- tuneRF(x = topo_train_data %>% select(-y_act, -MOF.ID), 
# #                                                    y = topo_train_data$y_act, ntree = 500)
# # m2 <- tuneRF(
# #   x          = topo_train_data %>% select(-y_act, -MOF.ID),
# #   y          = topo_train_data$y_act,
# #   ntreeTry   = 500,
# #   mtryStart  = 18,
# #   stepFactor = 1.5,
# #   improve    = 0.01,
# #   trace      = FALSE      # to not show real-time progress 
# # )
# 
# # include some other features
# # the smallest number in the grid file
# topo_train_data <- merge(train_data_with_id, structural_data, by ="MOF.ID")
# topo_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")
# read_minimum_energy <- function(dataset, file_dir){
#   # read all the files from the MOF.IDs
#   list_of_files <- dataset$MOF.ID
#   #file_dir <- "../DATAS/Diff_grids_XeKr_kB_Favors_Xe/"
#   y <- data.frame(MOF.ID=character(), min_energy=double(), exp_energy1=double(), 
#                   exp_energy2=double(), exp_energy3=double(), exp_energy4=double())
#   
#   for (grid_file in list_of_files){
#     file_name <- paste0(file_dir, grid_file, ".grid")
#     energies <- read_table(file_name, col_names = "V1", col_types = "d")$V1
#     minimum_energy <- min(energies)
#     exp_min_e <- exp(-minimum_energy)
#     a <- energies %>% unique() %>% sort()
#     exp_a <- exp(-a)
#     y <- rbind(y, data.frame(MOF.ID = grid_file, min_energy = minimum_energy, exp_energy1=exp_a[1], 
#                              exp_energy2=exp_a[2], exp_energy3=exp_a[3], exp_energy4=exp_a[4]))
#   }
#   new_dataset <- merge(dataset, y, by ="MOF.ID")
#   new_dataset
# }
# 
# train_with_min_E <- read_minimum_energy(topo_train_data, 
#                                         file_dir = "../DATAS/Diff_grids_XeKr_kB_Favors_Xe/")
# test_with_min_E <- read_minimum_energy(topo_test_data, 
#                                         file_dir = "../DATAS/Diff_grids_XeKr_kB_Favors_Xe/")
# rf_model_minE <- randomForest(x = train_with_min_E %>% select(-y_act, -MOF.ID), 
#                               y = train_with_min_E$y_act, ntree = 800)
# make_rf_prediction_plots(condition_name = string_to_paste, 
#                          plot_name = "rf_histogram_minE", 
#                          rf_model= rf_model_minE,  test_data = test_with_min_E)
# 
# # stack train and test together
# dataset_with_minE <- rbind(train_with_min_E, test_with_min_E)
# # filter out the large selectivities
# dataset_with_minE %>% filter(y_act > 1000)

#####################################
# we can try fit loading, here.######
#####################################
save_path <- paste0(save_path, "/FitLoading/")
if (!dir.exists(save_path)){
  dir.create(save_path)
}
train_data_Kr <- export_data(p_ch4_sets$training, 
                             rename(gcmc_data %>% mutate(g.L = Kr_uptake), y_act=g.L), ch4_binspec, with_id = TRUE)
test_data_Kr <- export_data(p_ch4_sets$testing, 
                            rename(gcmc_data %>% mutate(g.L = Kr_uptake), y_act=g.L), ch4_binspec, with_id = TRUE)

train_data_Xe <- export_data(p_ch4_sets$training, 
                             rename(gcmc_data %>% mutate(g.L = Xe_uptake), y_act=g.L), ch4_binspec, with_id = TRUE)
test_data_Xe <- export_data(p_ch4_sets$testing, 
                            rename(gcmc_data %>% mutate(g.L = Xe_uptake), y_act=g.L), ch4_binspec, with_id = TRUE)
train_data_Select <- export_data(p_ch4_sets$training, 
                                 rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), ch4_binspec, with_id = TRUE)
test_data_Select <- export_data(p_ch4_sets$testing, 
                                rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), ch4_binspec, with_id = TRUE)
colnames(train_data_Kr)[colnames(train_data_Kr)=="id"] <- "MOF.ID"
colnames(test_data_Kr)[colnames(test_data_Kr)=="id"] <- "MOF.ID"
topo_train_data_Kr <- merge(train_data_Kr, structural_data, by ="MOF.ID")
topo_test_data_Kr <- merge(test_data_Kr, structural_data, by ="MOF.ID")

colnames(train_data_Xe)[colnames(train_data_Xe)=="id"] <- "MOF.ID"
colnames(test_data_Xe)[colnames(test_data_Xe)=="id"] <- "MOF.ID"
topo_train_data_Xe <- merge(train_data_Xe, structural_data, by ="MOF.ID")
topo_test_data_Xe <- merge(test_data_Xe, structural_data, by ="MOF.ID")
model_Kr <- randomForest(x = train_data_Kr %>% select(-y_act, -MOF.ID), y = train_data_Kr$y_act, ntree = 500)
model_Xe <- randomForest(x = train_data_Xe %>% select(-y_act, -MOF.ID), y = train_data_Xe$y_act, ntree = 500)

make_rf_prediction_plots_XeKr(condition_name = string_to_paste, 
                              plot_name = "rf_Kr_Loading", 
                              rf_model= model_Kr,  test_data = test_data_Kr, lim = plot_limit, axis_label = "Loading [cm\u00B3/cm\u00B3]")

make_rf_prediction_plots_XeKr(condition_name = string_to_paste, 
                              plot_name = "rf_Xe_Loading", 
                              rf_model= model_Xe,  test_data = test_data_Xe, lim = plot_limit, axis_label = "Loading [cm\u00B3/cm\u00B3]")
Kr_prediction <- predict(model_Kr, test_data_Kr)
Kr_prediction_train <- model_Kr$predicted
Xe_prediction <- predict(model_Xe, test_data_Xe)
Xe_prediction_train <- model_Xe$predicted
# get the selectivity
XeKr_Select <- (Xe_prediction/0.2)/(Kr_prediction/0.8)
XeKr_Select_train <- (Xe_prediction_train/0.2)/(Kr_prediction_train/0.8)
test_rmse <- postResample(pred = XeKr_Select, obs = test_data_Select$y_act)
train_rmse <- postResample(pred = XeKr_Select_train, obs = train_data_Select$y_act)


makequickplot <- function(condition_name, plot_name, xdata, ydata, axislabel, test = FALSE, postresult, zoomin = FALSE){
  plot_limit <- c(0, max(xdata, ydata))
  if(!zoomin){
    plot_limit[2] <- plot_limit[2]-mod(plot_limit[2],50)+50
    if (previous_plot_lim < plot_limit[2])
    {
      previous_plot_lim <<- plot_limit[2] # use global variable
    }else{
      plot_limit[2] <- previous_plot_lim
    }
  }
  new_string <- paste0(condition_name, "_")
  if(poster)
  {
    new_string <- paste0(new_string, "_poster")
  }
  if(test){
    colorboard <- "#0070C0"
  }else{
    colorboard <- "#CA7C1B"
  }
  gg_rf <- qplot(x = xdata, y = ydata) %>% rescale_ch4_parity_XeKr(., lims=plot_limit) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    xlab(paste0("GCMC ", axislabel)) + 
    ylab(paste0("Predicted ", axislabel)) + 
    geom_point(color=colorboard, size = dot_size) + 
    theme_classic()
  if(test){
    gg_rf <- gg_rf %>% annotate_plot(paste0("Testing data\n", length(xdata)," ", "points"), "top.left", colorboard) %>% 
      label_stats(., postresult, plot_units = "", do_label_r2=TRUE)
  }else{
    gg_rf <- gg_rf %>% annotate_plot(paste0("Training data\n", length(xdata)," ", "points"), "top.left", colorboard) %>% 
      label_stats(., postresult, plot_units = "", do_label_r2=TRUE)
  }
  gg_rf <- gg_rf + theme(axis.text=element_text(size=30),
                         axis.title=element_text(size=30,face="bold"))
  if(test){
    save_plot(paste(save_path, paste0(new_string, paste0(plot_name, "_test.png")), sep = ""), 
              gg_rf, base_width = 10, base_height = 10, dpi = 600)
  }else{
    save_plot(paste(save_path, paste0(new_string, paste0(plot_name, "_train.png")), sep = ""), 
              gg_rf, base_width = 10, base_height = 10, dpi = 600)
  }
}
makequickplot(condition_name = string_to_paste, plot_name = 'Selectivity_train', 
              xdata = train_data_Select$y_act, 
              ydata = XeKr_Select_train, axislabel = "Selectivity",
              test = FALSE, postresult = train_rmse)
makequickplot(condition_name = string_to_paste, plot_name = 'Selectivity_test', 
              xdata = test_data_Select$y_act, 
              ydata = XeKr_Select, axislabel = "Selectivity",
              test = TRUE, postresult = test_rmse)
# try get rid of the outlier when plotting!
test_data_no_outlier <- test_data_Select %>% filter(test_data_Select$y_act < 30)
numbers <- which(test_data_Select$y_act >= 30)
if(length(numbers) > 0){
  XeKr_Select_no_outlier <- XeKr_Select[-numbers]
}else{
  XeKr_Select_no_outlier <- XeKr_Select
}
test_rmse_no_out <- postResample(pred = XeKr_Select_no_outlier, obs = test_data_no_outlier$y_act)
# do the same for train data
train_data_no_outlier <- train_data_Select %>% filter(train_data_Select$y_act < 30)
numbers <- which(train_data_Select$y_act >= 30)
if(length(numbers) > 0){
  XeKr_Select_no_outlier_train <- XeKr_Select_train[-numbers]
}else{
  XeKr_Select_no_outlier_train <- XeKr_Select_train
}
train_rmse_no_out <- postResample(pred = XeKr_Select_no_outlier_train, obs = train_data_no_outlier$y_act)
makequickplot(condition_name = string_to_paste, plot_name = 'Selectivity_test_remove2outliers', 
              xdata = test_data_no_outlier$y_act, 
              ydata = XeKr_Select_no_outlier, axislabel = "Selectivity",
              test = TRUE, postresult = test_rmse_no_out, zoomin = TRUE)
makequickplot(condition_name = string_to_paste, plot_name = 'Selectivity_train_remove2outliers', 
              xdata = train_data_no_outlier$y_act, 
              ydata = XeKr_Select_no_outlier_train, axislabel = "Selectivity",
              test = FALSE, postresult = train_rmse_no_out, zoomin = TRUE)
# print the train and test datasets, also save models!
model_Kr %>% saveRDS(., paste(save_path, paste0(string_to_paste, "Kr_RF.rds"), sep = ""))
model_Xe %>% saveRDS(., paste(save_path, paste0(string_to_paste, "Xe_RF.rds"), sep = ""))
topo_train_data_Kr %>% write.csv(., paste(save_path, paste0(string_to_paste, "Topo_train_Kr.csv"), sep = ""))
topo_test_data_Kr %>% write.csv(., paste(save_path, paste0(string_to_paste, "Topo_test_Kr.csv"), sep = ""))
topo_train_data_Xe %>% write.csv(., paste(save_path, paste0(string_to_paste, "Topo_train_Xe.csv"), sep = ""))
topo_test_data_Xe %>% write.csv(., paste(save_path, paste0(string_to_paste, "Topo_test_Xe.csv"), sep = ""))
train_data_Select %>% write.csv(., paste(save_path, paste0(string_to_paste, "Train_Selectivity.csv"), sep = ""))
test_data_Select %>% write.csv(., paste(save_path, paste0(string_to_paste, "Test_Selectivity.csv"), sep = ""))
