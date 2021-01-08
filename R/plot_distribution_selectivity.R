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
# poster plots should have big dots
if (poster)
{
  dot_size <<- 3
}else{
  dot_size <<- 1 
}
################################
# Read Data From the SI      ###
################################
gcmc_data_1bar  <- read_data(gcmc_file, read_SI = TRUE, sheetname = "XeKr_Mix_273K_1Bar", 
                             unit_of_ads = chosen_unit, 
                             relax = TRUE, just2k = TRUE, no_low_loading = FALSE)
gcmc_data_1bar$Selectivity <- (gcmc_data_1bar$Xe_uptake/0.2)/(gcmc_data_1bar$Kr_uptake/0.8)

gcmc_data_10bar <- read_data(gcmc_file, read_SI = TRUE, sheetname = "XeKr_Mix_273K_10Bar", 
                             unit_of_ads = chosen_unit, 
                             relax = TRUE, just2k = TRUE, no_low_loading = FALSE)
gcmc_data_10bar$Selectivity <- (gcmc_data_10bar$Xe_uptake/0.2)/(gcmc_data_10bar$Kr_uptake/0.8)

# create save path
save_path <- "Results/XeKr_1_10Bar_Selectivity_Distribution/"
if (!dir.exists(save_path)){
  dir.create(save_path)
}
#gcmc_data_1bar$Pres <- gcmc_data_1bar$Pres/1e5 # convert to bar
gcmc_data_1bar$Pres <- rep("1 Bar", times = nrow(gcmc_data_1bar))
#gcmc_data_10bar$Pres <- gcmc_data_10bar$Pres/1e5
gcmc_data_10bar$Pres <- rep("10 Bar", times = nrow(gcmc_data_10bar))
gcmc_data_1bar <- gcmc_data_1bar %>% select(Selectivity, Pres)
gcmc_data_10bar <- gcmc_data_10bar %>% select(Selectivity, Pres)
# merge selectivities 
gcmc_select <- rbind(gcmc_data_1bar, gcmc_data_10bar)
gcmc_select$Pres <- as.factor(gcmc_select$Pres)
names(gcmc_select)[names(gcmc_select) == "Pres"] <- "Pressure"
p <- ggplot(gcmc_select, aes(x=Pressure, y = Selectivity, fill = Pressure)) + 
  xlab(paste0("Pressure")) + 
  ylab(paste0("Selectivity")) + 
  ylim(0, 400) +
  theme_classic() + 
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  scale_y_continuous(breaks = c(10, 50, 100, 200, 300, 350)) + 
  geom_violin()

#p + stat_summary(fun.y=median, geom="point", shape=23, size=2, color = 'red')
fig <- p + geom_boxplot(width=0.1)
save_plot(paste0(save_path, "violin_dist.png"), fig, base_width = 10, base_height = 10, dpi = 600)

gcmc_select$Selectivity[gcmc_select$Selectivity < 0.01] <- 0.01
p_log <- ggplot(gcmc_select, aes(x=Pressure, y = Selectivity, fill = Pressure)) + 
  xlab(paste0("Pressure [Bar]")) + 
  ylab(paste0("Selectivity")) + 
  theme_classic() + 
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  geom_violin() + geom_boxplot(width=0.1) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 1, 4, 10, 20, 50, 100, 200, 300, 400))

save_plot(paste0(save_path, "violin_dist_log.png"), p_log, base_width = 10, base_height = 10, dpi = 600)
