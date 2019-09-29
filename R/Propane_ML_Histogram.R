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
# in command prompt: Rscript --vanilla R\save_h2_hists.R whatever.rds Energies\ use_ch4
# tell it whether to save plots for poster
poster <<- TRUE # save this as a global variable: no need to pass it around
# poster plots should have big dots
if (poster){
  dot_size <<- 3
}else{
  dot_size <<- 1 
}
 # define the input files
 gcmc_file = "All_data/propane_1bar.txt"
 grid_file = "All_data/CH3_1A_probe_UFF.rds"
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
 #unit_list <- c("Molec_cm3overcm3", "Mol_kg", "Mill_gram", "Cm3_gram")
 unit_list <- c("Molec_cm3overcm3")
 for (chosen_unit in unit_list){ 
       gcmc_data <- read_data(gcmc_file, unit_of_ads = chosen_unit)
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
       molecule_name <- sub(".*\\/", "", gcmc_file) %>% sub("\\_.*", "", .)
       temperature <- gcmc_data$Temp[1]
       pressure <- gcmc_data$Pres[1]
       string_to_paste <- paste(molecule_name, temperature, pressure, grid_info, sep = "_")
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
       # now, we want to consider the effect of units, so add units to save_path
       unit_path <- paste0(save_path, chosen_unit, "/")
       if (!dir.exists(unit_path)){
         dir.create(unit_path)
       }
       save_path <- unit_path
       # here we can have two options, feed uniform or non-uniform bins
       #binbounds <- bounds_from_params(ch4_binspec)
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
                                plot_units=unit_for_plot, db_name = "tobacco")
      # calculate the plot limit: largest value in the dataset
      plot_limit <- max(gcmc_data$Uptake)
      # round to the nearest 10
      plot_limit <- round(plot_limit, digits = -1) + 10
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
                                 rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), ch4_binspec)
       test_data <- export_data(p_ch4_sets$testing, 
                                rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), ch4_binspec)
       
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
                                plot_name = "rf", lim = plot_limit, 
                                rf_model= rf_model,  test_data = test_data)
      
      # Part below are related to topology, see Paper by Yamil and Diego
      # Read the names with topology
      topologies <- read.table("All_data/fullnames_without_tob_cleaner_py.txt") # first column is ID, second column is topology
      
      # just keep all the structural properties that are numbers
      # filter out those tobacco data with NAs
      orig_tobacco_data <-  na.omit(read_textual_data())
      structural_data <- orig_tobacco_data %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)
      
      
      train_data_with_id <- export_data(p_ch4_sets$training, 
                                        rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), 
                                        ch4_binspec, with_id = TRUE)
      test_data_with_id <- export_data(p_ch4_sets$testing, 
                                       rename(gcmc_data %>% mutate(g.L = Uptake), y_act=g.L), 
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
                               plot_name = "rf_histogram_Rscore", lim = plot_limit, 
                               rf_model= rf_model_Rscore,  test_data = testing_data_with_R)
      
      # learning using Topology data
      topo_train_data<- merge(train_data_with_id, structural_data, by ="MOF.ID")
      topo_test_data <- merge(test_data_with_id, structural_data, by ="MOF.ID")
      
      # perform a random forest model with these structural properties
      # check if training and testing data has MOF.ID that has NAs in tobacco data
      
      rf_model_topo <- randomForest(x = topo_train_data %>% select(-y_act, -MOF.ID), 
                                    y = topo_train_data$y_act, ntree = 500)
      make_rf_prediction_plots(condition_name = string_to_paste, 
                               plot_name = "rf_histogram_topology", lim = plot_limit, 
                               rf_model= rf_model_topo,  test_data = topo_test_data)
      # do this qplot stuff again, to extract what are those outliers
      
      # seems ch4 qst at low loading is interesting...add into feature space
      # data_with_met_qst <- orig_tobacco_data %>% select(MOF.ID, vf, vsa, gsa, pld, lcd)
      # topo_qst_train_data<- merge(train_data_with_id, data_with_met_qst, by ="MOF.ID")
      # topo_qst_test_data <- merge(test_data_with_id, data_with_met_qst, by ="MOF.ID")
      # rf_model_topo_qst <- randomForest(x = topo_qst_train_data %>% select(-y_act, -MOF.ID), 
      #                               y = topo_qst_train_data$y_act, ntree = 500)
      # make_rf_prediction_plots(condition_name = string_to_paste, 
      #                          plot_name = "rf_histogram_topology_qst", 
      #                          rf_model= rf_model_topo_qst,  test_data = topo_qst_test_data)
      # generate topology histograms
      topology_data <- rbind(topo_train_data, topo_test_data)
      make_topology_histograms(condition_name = string_to_paste, topo_data = topology_data)
      
      # for outliers in hexane simulations
      if (grepl("Hexane_298_10000", string_to_paste)){
        # also make histograms for those large deviated points
        # filter out those whose difference between predicted and actual is greater than 40
        tested <- predict(rf_model_topo, topo_test_data)
        train_pred <- predict(rf_model_topo)
        large_deviations_train <- topo_train_data
        large_deviations_test <- topo_test_data
        large_deviations_train$pred <- train_pred
        large_deviations_test$pred <- tested
        large_deviations <- rbind(large_deviations_train, large_deviations_test)
        small_deviations <- large_deviations %>% filter(abs(y_act - pred) < 30)
        large_deviations <- large_deviations %>% filter(abs(y_act - pred) >= 30)
        axis_label <- "capacity (cm\u00B3/cm\u00B3)"
        gg_rf_test_without_outliers <- qplot(x = small_deviations$y_act, y = small_deviations$pred) +
          geom_abline(slope = 1, intercept = 0, linetype = 2) +
          xlab(paste0("GCMC ", axis_label)) +
          ylab(paste0("Predicted ", axis_label))
        # get the gcmc_data for those who has large deviations
        large_deviation_gcmc <- gcmc_data %>% filter(ID %in% large_deviations$MOF.ID)
        # filter the original tobacco data
        large_deviation_tobacco <- orig_tobacco_data %>% filter(MOF.ID %in% large_deviations$MOF.ID)
        new_string_to_paste <- paste0(string_to_paste, "_2")
        make_topology_histograms(condition_name = new_string_to_paste, topo_data = large_deviation_tobacco)


        # outlier ranges for hexane at 0.1 bar
        # for running new simulations, filter out those already simulated
        not_touched <- structural_data %>% filter(!MOF.ID %in% gcmc_data$ID)
        # filter vf between 0.9 and 0.95
        new_mofs_needed <- not_touched %>% filter(vf >= 0.9 & vf <= 0.95)
        # then filter pld: between 15 and 25
        new_mofs_needed <- new_mofs_needed %>% filter(pld >= 15 & pld <= 25)
        # then filter lcd between 20 and 35
        new_mofs_needed <- new_mofs_needed %>% filter(lcd >= 20 & lcd <= 35)
        # then filter gsa between 6250 and 8125
        new_mofs_needed <- new_mofs_needed %>% filter(gsa >= 6250 & gsa <= 8125)
        # finally filter vsa between 800 and 1200
        new_mofs_needed <- new_mofs_needed %>% filter(vsa >= 800 & vsa <= 1200)
        # then filter out those with mc_0 nodes
        mc_0_ids <- read_table("All_data/IDs_with_mc_0.txt")
        new_mofs_needed <- new_mofs_needed %>% filter(!MOF.ID %in% mc_0_ids$IDs_with_mc_0)
        # make a list with unitcell xyz
        mof_xyz <- read_table2("tob_perm_unitcells_no_mc_0.txt")
        mof_xyz$ID <- as.character(mof_xyz$ID)
        mof_xyz <- mof_xyz %>% filter(ID %in% new_mofs_needed$MOF.ID)
        # convert IDs back to numbers: if write as characters, it will have quotes!
        mof_xyz$ID <- as.integer(mof_xyz$ID)
        write.table(mof_xyz, "biasing_list_for_hexane_point_1bar.txt", row.names = FALSE, col.names = FALSE, sep = " ")
      }
      # adding new features to the feature space
      # add minimum energy, mean energy, median energy, 1st quartile
      # remember to exclude 1.0e23, which is a place-holder
      read_minimum_energy_iqr <- function(dataset, file_dir){
        # read all the files from the MOF.IDs
        list_of_files <- dataset$MOF.ID
        #file_dir <- "../DATAS/Diff_grids_XeKr_kB_Favors_Xe/"
        y <- data.frame(MOF.ID=character(), min_energy=double(), 
                        mean_energy=double(),
                        median_energy=double(), 
                        fst_quartile=double())
        
        for (grid_file in list_of_files){
          file_name <- paste0(file_dir, grid_file, ".grid")
          energies <- read_table(file_name, col_names = "V1", col_types = "d")$V1
          # exclude 1.0e23
          energies <- energies[!energies == 1.0e23]
          minimum_energy <- min(energies)
          mean_E <- mean(energies)
          median_E <- median(energies)
          fst_iqr_E <- as.double(quantile(energies)[2])
          y <- rbind(y, data.frame(MOF.ID = grid_file, 
                                   min_energy = minimum_energy, 
                                   mean_energy=mean_E, 
                                   median_energy=median_E, 
                                   fst_quartile=fst_iqr_E))
        }
        new_dataset <- merge(dataset, y, by ="MOF.ID")
        new_dataset
      }
      more_feature_train <- read_minimum_energy_iqr(topo_train_data, file_dir = "Raw_energies_data/Energies/")
      more_feature_test <- read_minimum_energy_iqr(topo_test_data, file_dir = "Raw_energies_data/Energies/")
      rf_model_stats <- randomForest(x = more_feature_train %>% select(-y_act, -MOF.ID), 
                                    y = more_feature_train$y_act, ntree = 500)
      make_rf_prediction_plots(condition_name = string_to_paste, 
                               plot_name = "rf_stats", lim = plot_limit, 
                               rf_model= rf_model_stats,  test_data = more_feature_test)
 }
 
 # make histogram of 2 mofs
 # add the low/up bounds for the last bin: Inf
 binbounds_inf <- binbounds
 inf_bounds <- data.frame(0.0, 2)
 names(inf_bounds) <- c("lower", "upper")
 binbounds_inf <- rbind(binbounds_inf, inf_bounds)
 binbounds_inf$widths <- binbounds_inf$upper - binbounds_inf$lower
 t_binbounds_inf <- as.data.frame(t(binbounds_inf))
 names(t_binbounds_inf) <- 1:(length(binbounds$lower) + 1)

 
 bin_heights <- train_data_with_id %>% filter(MOF.ID == "1713") %>% select(-MOF.ID, -y_act)
 # rename the Inf bin to length(binbounds$lower) + 1
 colnames(bin_heights)[colnames(bin_heights) == "Inf"] = as.character(length(binbounds$lower) + 1)
 # bind t_binbounds_inf with bin_heights
 bin_info <- rbind(t_binbounds_inf, bin_heights)
 # finally, transpose
 t_bin_info <- as.data.frame(t(bin_info))
 names(t_bin_info) <- c("lower", "upper", "range", "height")
 focused <- TRUE
 if (focused){
   t_bin_info <- t_bin_info %>% filter(lower >= -10)
 }
 hist_df <- data.frame(width = t_bin_info$range, height = t_bin_info$height)
 hist_df$w <- cumsum(hist_df$width)
 hist_df$wm <- hist_df$w - hist_df$width
 p <- ggplot(hist_df)
 p1 <- p + geom_rect(aes(xmin = t_bin_info$lower, 
                         xmax = t_bin_info$upper, 
                         ymin = 0, ymax = height)) + 
            xlab("Energies(kJ/mol)") + 
            ylab("Relative Bin Counts") + theme_classic() + 
              theme(axis.text=element_text(size=30),
              axis.title=element_text(size=30,face="bold")) 
 p1 
 save_plot(paste(save_path, paste0(string_to_paste, "1713_histo.png"), sep = ""), 
           p1, base_width = 10, base_height = 10, dpi = 600)
 # analysis on hexane outliers with movies
 out_gcmc <- read_table2("../DATAS/Hexane_.1bar_movies.tar/Hexane_.1bar_movies/file.txt")
 outlier_id <- as.character(out_gcmc$ID)
 # filter out those predictions of outliers
 # use full tobacco data
 tobacco_outlier <- read_textual_data(keep = TRUE) %>% filter(MOF.ID %in% outlier_id)
# lets just compare 1713 and 4510
 tob_1713 <- tobacco_outlier %>% filter(MOF.ID == "1713")
 tob_4510 <- tobacco_outlier %>% filter(MOF.ID == "4510") 
 tob <- rbind(tob_1713, tob_4510)
 
 # new feature: add highest value of derivatives of PSD to the feature space
 # calculate the difference between lcd and diameter where maximum derivative occurs
 psd_data <- read_table2("Data/PSD_tobacco.txt")
 psd_data$MOF.ID <- as.character(psd_data$MOF.ID)
 # filter out those don't have good psd data
 psd_train_data <- merge(topo_train_data, psd_data, by = "MOF.ID")
 psd_test_data <- merge(topo_test_data, psd_data, by = "MOF.ID")
 # calculate the difference between lcd and max_dir
 psd_train_data$distance <- psd_train_data$max_dir_diameter - psd_train_data$lcd
 psd_test_data$distance <- psd_test_data$max_dir_diameter - psd_test_data$lcd
 rf_model_psd <- randomForest(x = psd_train_data %>% select(-y_act, -MOF.ID), y = psd_train_data$y_act, ntree = 500)
 make_rf_prediction_plots(condition_name = string_to_paste, 
                          plot_name = "rf_histogram_psd_max_dir", lim = plot_limit, 
                          rf_model= rf_model_psd,  test_data = psd_test_data)
 
 # try neural net
 # get the names of the features
 # type the formula
 topo_data <- rbind(topo_train_data, topo_test_data)
 V_topo_train_data <- topo_train_data
 V_topo_test_data <- topo_test_data
 colnames(V_topo_train_data) <- paste("V", colnames(V_topo_train_data), sep = "")
 colnames(V_topo_test_data) <- paste("V", colnames(V_topo_test_data), sep = "")
 feature_names <- V_topo_train_data %>% select(-VMOF.ID, -Vy_act, -VInf) %>% colnames()
 f <- as.formula(paste("Vy_act ~ ", paste(feature_names, collapse = " + ")))
 #m <- model.matrix(as.formula(paste("Vy_act ~ ", paste(feature_names, collapse = " + "))), topo_train_data)
 nn <- neuralnet(formula = f, data = V_topo_train_data, hidden = c(5,3), stepmax = 1e6)
 pred <- predict(nn, V_topo_test_data) 
 
 