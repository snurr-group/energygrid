rm(list=ls())
# for checking package availability, if not, install it
source("R/package_verification.R")

source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/read_gcmc_data.R")
source("R/read_textural_data.R")
source("R/save_train_test_data.R")
source("R/refined_bins_calc.R")
# to check, use "gcmc_data %>% filter(!(ID %in% hmof_h2_grid$id))"
# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4
# tell it whether to save plots for poster
poster <<- TRUE # save this as a global variable: no need to pass it around
set.seed(12345) #just do it for Xe Kr comparison, disable for alkane fitting: original: 12345
# poster plots should have big dots
if (poster){
  dot_size <<- 3
}else{
  dot_size <<- 1 
}
XeKr <<- FALSE # for distinguishing normal fitting from XeKr selectivity fitting
 # define the input files
 molecule_name <- "Kr"
 Temperature <- "273K"
 Pressure <- "1Bar" # use Bar
 previous_plot_lim <- 300

 grid_file = "All_data/Kr_1A_greater_range.rds" # use for Xe and Kr
 #grid_file = "All_data/CH3_1A_combined_tobacco_core.rds" # use for alkanes
 # extract the directory from the file name
 data_dir <- sub("\\/.*", "", grid_file)
 file_name <- sub(".*\\/", "", grid_file)
 data_dir <- getwd() %>% paste0(., "/", data_dir, "/")
 # find files with that prefix
 # because 0.5 A is too large, we splitted it into 5 rds files
 list_of_files <- list.files(data_dir, pattern = file_name)
 # then check length of the list of files
 if (length(list_of_files) > 1){
   for (i in 1:length(list_of_files)){
     gridi = read_rds(paste0(grid_file, as.character(i), ".rds"))
     if (i == 1) {
       grids = gridi
     } else{
       grids <- rbind(grids, gridi)
     }
   }
   hmof_h2_grid <- grids
 } else{
   hmof_h2_grid <- read_rds(grid_file)
 }

 if (!dir.exists("Results")){
   dir.create("Results")
 }     
       chosen_unit <- "cm3overcm3"
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
       
       string_to_paste <- paste(molecule_name, Temperature, Pressure, grid_info, sep = "_")
       # create a directory for that condition and molecule
       save_path <- paste0(string_to_paste, "/", sep = "")
       # if save for poster figures, put poster in the front
       save_path <- paste0("Results/", save_path)
       if (!dir.exists(save_path)){
         dir.create(save_path)
       }
       # now, we want to consider the effect of units, so add units to save_path
       unit_path <- paste0(save_path, chosen_unit, "/")
       if (!dir.exists(unit_path)){
         dir.create(unit_path)
       }
       save_path <- unit_path
       # here we can have two options, feed uniform or non-uniform bins
       #binbounds <- bounds_from_params(ch4_binspec)
       # # delete the line below, just testing Core MOFs
       # a <- read_table2(gcmc_file)
       # b <- a %>% select(-Heat_of_Ads, -Heat_fluc)
       # bb <- na.omit(b)
       # gcmc_data <- bb %>% select(ID, Temp, Pres, cm3overcm3)
       # gcmc_data <- mutate(gcmc_data, Uptake = cm3overcm3)# delete!
       gcmc_data <- mutate(gcmc_data, id=ID)
       new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
       Uniform_bins <- FALSE
       if (!Uniform_bins){
       binbounds <- automatic_bins(new_grid) # for non-uniform bins
       binbounds <- as.data.frame(binbounds)
       binbounds <- binbounds %>% filter(upper <= 0.0) # bind all repulsive bins
       } else{
       # if we want to use uniform bins
       binbounds <- data.frame()
       # 1 KJ/mol spacing
       lower <- seq(-50, -1, by = 1)
       upper <- seq(-50+1, 0, by = 1)
       binbounds <- data.frame(lower, upper)
       }
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
        mutate(g.L = Uptake) %>%
        run_model_on_partitions(p_ch4_sets, ., ch4_binspec, 
                                plot_units=unit_for_plot, db_name = "points")
      # calculate the plot limit: largest value in the dataset
      plot_limit <- max(gcmc_data$Uptake)
      # round to the nearest 10
      plot_limit <- plot_limit-mod(plot_limit,100)+100
      gg_train <- p_ch4_vol$plots$parity_training + 
        theme(axis.title.y = element_text(hjust=0.5)) + 
        xlab(paste0("GCMC capacity",  "(", unit_for_plot, ")")) + 
        ylab(paste0("Predicted capacity", "(", unit_for_plot , ")"))
      
      gg_test <- p_ch4_vol$plots$parity_testing + 
        theme(axis.title.y = element_text(hjust=0.5)) + 
        xlab(paste0("GCMC capacity",  "(", unit_for_plot, ")")) + 
        ylab(paste0("Predicted capacity", "(", unit_for_plot , ")"))
      if (poster){
        gg_train <- gg_train + theme(axis.text=element_text(size=30),
                         axis.title=element_text(size=30,face="bold"))
        gg_test <- gg_test + theme(axis.text=element_text(size=30),
                                axis.title=element_text(size=30,face="bold")) 
      }
      save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_train.png"), sep = ""), 
                gg_train, base_width = 10, base_height = 10, dpi = 600)
      save_plot(paste(save_path, paste0(string_to_paste, "_LASSO_test.png"), sep = ""), 
                gg_test, base_width = 10, base_height = 10, dpi = 600)
      
      # export data for later testing other models, one can also export data to csv files
      train_data <- export_data(p_ch4_sets$training, 
                               rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), ch4_binspec, with_id = TRUE)
      
      test_data <- export_data(p_ch4_sets$testing, 
                               rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), ch4_binspec, with_id = TRUE)
      # sort by 1st digit
      train_data <- train_data[order(train_data$id),]
      test_data <- test_data[order(test_data$id),]
       # add LASSO predictions to the data
       train_data$lassopred <- predict(p_ch4_vol$trained_mod$mod, as.matrix(p_ch4_vol$trained_mod$x))
       test_data$lassopred <- p_ch4_vol$pred_df$y_pred
       train_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "_LASSO_train.csv"), sep = ""))
       test_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "_LASSO_test.csv"), sep = ""))
       # remove the predictions, and IDs for later
       train_data <- train_data %>% select(-lassopred)
       test_data <- test_data %>% select(-lassopred)
       
##############################
# RANDOM FOREST STARTS HERE###
##############################
       
      # Run a random forest model
       rf_model <- randomForest(x = train_data %>% select(-y_act, -id), 
                                y = train_data$y_act, ntree = 500, 
                                importance = TRUE)
       # gg <- qplot(x = model$predicted, y = model$y) + geom_abline(slope = 1, intercept = 0)
      
       make_rf_prediction_plots(condition_name = string_to_paste, 
                                plot_name = "RF", lim = plot_limit, 
                                rf_model= rf_model,  test_data = test_data)
       textcondition <- paste0("Regression of ", molecule_name, 
                               " at ", toString(gcmc_data$Temp[1]), 
                               " K and ",toString(gcmc_data$Pres[1]/1e5), " Bar")
       modelname <- paste0("RF using ", "Energy Histogram")
       getVarImp(rf_model, modelshort = "RF", 
                 howmany = 10, condition = textcondition, modelname = modelname)
       train_data$predicted <- predict(rf_model, train_data)
       test_data$predicted <- predict(rf_model, test_data)
       train_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "_RF_train.csv"), sep = ""))
       test_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "_RF_test.csv"), sep = ""))
       train_data <- train_data %>% select(-predicted)
       test_data <- test_data %>% select(-predicted)
      # Part below are related to textural properties, see Paper by Yamil and Diego
      textural_prop <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is textural properties
       
      # just keep all the structural properties that are numbers
      # filter out those tobacco data with NAs
      ToBaCCo <- TRUE
      CoRE <- FALSE
      if (ToBaCCo){
        orig_tobacco_data <-  na.omit(read_textual_data())
      }else if(CoRE){
        orig_tobacco_data <-  na.omit(read_textual_data(option = "CoRE"))
      }
      #
      structural_data <- orig_tobacco_data # %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)
      
      
      train_data_with_id <- export_data(p_ch4_sets$training, 
                                        rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), 
                                        ch4_binspec, with_id = TRUE)
      test_data_with_id <- export_data(p_ch4_sets$testing, 
                                       rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), 
                                       ch4_binspec, with_id = TRUE)
      
      # calculate the correlation matrix for these data
      colnames(train_data_with_id)[colnames(train_data_with_id)=="id"] <- "MOF.ID"
      colnames(test_data_with_id)[colnames(test_data_with_id)=="id"] <- "MOF.ID"
########################################
# USE TEXTURAL PROPERTIES AS FEATURES###
########################################
      # learning using Textural properties data
      text_train_data<- merge(train_data_with_id, structural_data, by ="MOF.ID")
      text_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")
      
      # perform a random forest model with these structural properties
      # check if training and testing data has MOF.ID that has NAs in tobacco data
      
      rf_model_text <- randomForest(x = text_train_data %>% select(-y_act, -MOF.ID), 
                                    y = text_train_data$y_act, ntree = 500, 
                                    importance = TRUE)
      make_rf_prediction_plots(condition_name = string_to_paste, 
                               plot_name = "RF_Histogram_Textural", lim = plot_limit, 
                               rf_model= rf_model_text,  test_data = text_test_data)
      
      textcondition <- paste0("Regression of ", molecule_name, 
                              " at ", toString(gcmc_data$Temp[1]), 
                              " K and ",toString(gcmc_data$Pres[1]), " Bar")
      modelname <- paste0("RF using ", "Energy Histogram", " and ", 
                          "Textural Properties")
      getVarImp(rf_model_text, modelshort = "RF-Textural", 
                howmany = 10, condition = textcondition, modelname = modelname)
      text_train_data$predicted <- rf_model_text$predicted
      text_test_data$predicted <- predict(rf_model_text, text_test_data)
      text_train_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "_RF_Textural_train.csv"), sep = ""))
      text_test_data %>% write.csv(., paste(save_path, paste0(string_to_paste, "_RF_Textural_test.csv"), sep = ""))
      text_train_data <- text_train_data %>% select(-predicted)
      text_test_data <- text_test_data %>% select(-predicted)

      # generate textural properties histograms
      textural_data <- rbind(text_train_data, text_test_data)
      make_textural_histograms(condition_name = string_to_paste, text_data = textural_data)
      
###################################
# ENERGY STATISTICS AS FEATURES ###
###################################
      
      # adding new features to the feature space
      # add minimum energy, mean energy, median energy, 1st quartile
      # remember to exclude 1.0e23, which is a place-holder
      read_minimum_energy_iqr <- function(dataset, file_dir){
        # read all the files from the MOF.IDs
        list_of_files <- dataset$MOF.ID
        #file_dir <- "../DATAS/Diff_grids_XeKr_kB_Favors_Xe/"
        y <- data.frame(MOF.ID=character(), min_energy=double(), 
                        #mean_energy=double(),
                        median_energy=double(), 
                        fst_quartile=double())
        
        for (grid_file in list_of_files){
          file_name <- paste0(file_dir, grid_file, ".grid")
          energies <- read_table(file_name, col_names = "V1", col_types = "d")$V1
          # exclude 1.0e23
          energies <- energies[!energies == 1.0e23]
          minimum_energy <- min(energies)
          #mean_E <- mean(energies) Mean value excluded
          median_E <- median(energies)
          fst_iqr_E <- as.double(quantile(energies)[2])
          y <- rbind(y, data.frame(MOF.ID = grid_file, 
                                   Min = minimum_energy, 
                                   #Mean=mean_E, 
                                   Median=median_E, 
                                   Q1=fst_iqr_E))
        }
        new_dataset <- merge(dataset, y, by ="MOF.ID")
        new_dataset
      }

      if (ToBaCCo){
      more_feature_train <- read_minimum_energy_iqr(text_train_data, file_dir = "Raw_energies_data/ToBaCCo_CH3_Energies/")
      more_feature_test <- read_minimum_energy_iqr(text_test_data, file_dir = "Raw_energies_data/ToBaCCo_CH3_Energies/")
      } else if (CoRE){      
      more_feature_train <- read_minimum_energy_iqr(text_train_data, file_dir = "../DATAS/CORE_CH3/Energies/")
      more_feature_test <- read_minimum_energy_iqr(text_test_data, file_dir = "../DATAS/CORE_CH3/Energies/")
      }
      rf_model_stats <- randomForest(x = more_feature_train %>% select(-y_act, -MOF.ID), 
                                    y = more_feature_train$y_act, ntree = 500, 
                                    importance = TRUE)
      make_rf_prediction_plots(condition_name = string_to_paste, 
                               plot_name = "RF_Textural_Stats", lim = plot_limit, 
                               rf_model= rf_model_stats,  test_data = more_feature_test)
      
      textcondition <- paste0("Regression of ", molecule_name, 
                              " at ", toString(gcmc_data$Temp[1]), 
                              " K and ",toString(gcmc_data$Pres[1]/1e5), " Bar")
      modelname <- paste0("RF using ", "Energy Histogram", ", \n", 
                          "Textural Properties," ," and ", "Energy Stats")
      getVarImp(rf_model_stats, modelshort = "RF-Textural-Energy-Stats", 
                howmany = 10, condition = textcondition, modelname = modelname)
      
      more_feature_train$predicted <- predict(rf_model_stats, more_feature_train)
      more_feature_test$predicted <- predict(rf_model_stats, more_feature_test)
      more_feature_train %>% write.csv(., paste(save_path, paste0(string_to_paste, "_RF_Stats_train.csv"), sep = ""))
      more_feature_test %>% write.csv(., paste(save_path, paste0(string_to_paste, "_RF_Stats_test.csv"), sep = ""))
      