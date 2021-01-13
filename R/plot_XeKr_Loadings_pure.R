rm(list=ls())
source("R/refined_bins_calc.R")
source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/save_train_test_data.R")
source("R/read_textural_data.R")
source("R/package_verification.R")

unit_for_plot <<- "" # Selectivity is unitless
poster <<- TRUE # save this as a global variable: no need to pass it around
# poster plots should have big dots
if (poster)
{
  dot_size <<- 3
}else{
  dot_size <<- 1 
}
Xe_1bar <- TRUE
Xe_10bar <- TRUE
Kr_1bar <- TRUE
Kr_10bar <- TRUE

save_path <- "Results/Xe_Kr_SingleComponent_Analysis/"

if (!dir.exists(save_path)){
  dir.create(save_path)
}

if(Xe_1bar){
  gcmc_Xe1bar_train <- read.csv("Results/Xe_273K_1Bar_Kr_1A/cm3overcm3/Xe_273K_1Bar_Kr_1A_RF_train.csv")
  gcmc_Xe1bar_test <- read.csv("Results/Xe_273K_1Bar_Kr_1A/cm3overcm3/Xe_273K_1Bar_Kr_1A_RF_test.csv")
  gcmc_Xe1bar <- rbind(gcmc_Xe1bar_train, gcmc_Xe1bar_test)
  #names(gcmc_Xe1bar)[names(gcmc_Xe1bar) == "y_act"] <- "Xe_1Bar"
  gcmc_data <- read_excel("All_data/SI.xlsx", sheet = "Xe_273K_1Bar")
  gcmc_data$ID <- as.character(gcmc_data$ID)
  if(sum(grepl("Error", names(gcmc_data))) == 1){
    gcmc_data <- gcmc_data %>% select(-`GCMC_Uptake_Error_Bar [cm3/cm3]`)
  }
  names(gcmc_data) <- c("ID", "Temp", "Pres", "Uptake", "Heat_of_Ads")
  heat_Xe1bar <- gcmc_data %>% filter(ID %in% gcmc_Xe1bar$id) %>% select(Heat_of_Ads)
  heat_Xe1bar <- heat_Xe1bar %>% na.omit()
  heat_Xe1bar$Condition <- strrep("Xe_1bar",times = rep(1,nrow(heat_Xe1bar)))
}
if(Xe_10bar){
  gcmc_Xe10bar_train <- read.csv("Results/Xe_273K_10Bar_Kr_1A/cm3overcm3/Xe_273K_10Bar_Kr_1A_RF_train.csv")
  gcmc_Xe10bar_test <- read.csv("Results/Xe_273K_10Bar_Kr_1A/cm3overcm3/Xe_273K_10Bar_Kr_1A_RF_test.csv")
  gcmc_Xe10bar <- rbind(gcmc_Xe10bar_train, gcmc_Xe10bar_test)
  #names(gcmc_Xe10bar)[names(gcmc_Xe10bar) == "y_act"] <- "Xe_10Bar"
  gcmc_data <- read_excel("All_data/SI.xlsx", sheet = "Xe_273K_10Bar")
  gcmc_data$ID <- as.character(gcmc_data$ID)
  if(sum(grepl("Error", names(gcmc_data))) == 1){
    gcmc_data <- gcmc_data %>% select(-`GCMC_Uptake_Error_Bar [cm3/cm3]`)
  }
  names(gcmc_data) <- c("ID", "Temp", "Pres", "Uptake", "Heat_of_Ads")
  heat_Xe10bar <- gcmc_data %>% filter(ID %in% gcmc_Xe10bar$id) %>% select(Heat_of_Ads)
  heat_Xe10bar <- heat_Xe10bar %>% na.omit()
  heat_Xe10bar$Condition <- strrep("Xe_10bar",times = rep(1,nrow(heat_Xe10bar)))
}
if(Kr_1bar){
  gcmc_Kr1bar_train <- read.csv("Results/Kr_273K_1Bar_Kr_1A/cm3overcm3/Kr_273K_1Bar_Kr_1A_RF_train.csv")
  gcmc_Kr1bar_test <- read.csv("Results/Kr_273K_1Bar_Kr_1A/cm3overcm3/Kr_273K_1Bar_Kr_1A_RF_test.csv")
  gcmc_Kr1bar <- rbind(gcmc_Kr1bar_train, gcmc_Kr1bar_test)
  #names(gcmc_Kr1bar)[names(gcmc_Kr1bar) == "y_act"] <- "Kr_1Bar"
  gcmc_data <- read_excel("All_data/SI.xlsx", sheet = "Kr_273K_1Bar")
  gcmc_data$ID <- as.character(gcmc_data$ID)
  if(sum(grepl("Error", names(gcmc_data))) == 1){
    gcmc_data <- gcmc_data %>% select(-`GCMC_Uptake_Error_Bar [cm3/cm3]`)
  }
  names(gcmc_data) <- c("ID", "Temp", "Pres", "Uptake", "Heat_of_Ads")
  heat_Kr1bar <- gcmc_data %>% filter(ID %in% gcmc_Kr1bar$id) %>% select(Heat_of_Ads)
  heat_Kr1bar <- heat_Kr1bar %>% na.omit()
  heat_Kr1bar$Condition <- strrep("Kr_1bar",times = rep(1,nrow(heat_Kr1bar)))
}
if(Kr_10bar){
  gcmc_Kr10bar_train <- read.csv("Results/Kr_273K_10Bar_Kr_1A/cm3overcm3/Kr_273K_10Bar_Kr_1A_RF_train.csv")
  gcmc_Kr10bar_test <- read.csv("Results/Kr_273K_10Bar_Kr_1A/cm3overcm3/Kr_273K_10Bar_Kr_1A_RF_test.csv")
  gcmc_Kr10bar <- rbind(gcmc_Kr10bar_train, gcmc_Kr10bar_test)
  #names(gcmc_Kr10bar)[names(gcmc_Kr10bar) == "y_act"] <- "Kr_10Bar"
  gcmc_data <- read_excel("All_data/SI.xlsx", sheet = "Kr_273K_10Bar")
  gcmc_data$ID <- as.character(gcmc_data$ID)
  if(sum(grepl("Error", names(gcmc_data))) == 1){
    gcmc_data <- gcmc_data %>% select(-`GCMC_Uptake_Error_Bar [cm3/cm3]`)
  }
  names(gcmc_data) <- c("ID", "Temp", "Pres", "Uptake", "Heat_of_Ads")
  heat_Kr10bar <- gcmc_data %>% filter(ID %in% gcmc_Kr10bar$id) %>% select(Heat_of_Ads)
  heat_Kr10bar <- heat_Kr10bar %>% na.omit()
  heat_Kr10bar$Condition <- strrep("Kr_10bar",times = rep(1,nrow(heat_Kr10bar)))
}
# extract the probe name, molecule name, pressure and temperature for picture naming
# the naming convention of rds file should be: "probe"_"size of grid"_whatever.rds
# have an option here to make the 4 distributions into 1 plot
if(Xe_1bar & Xe_10bar & Kr_1bar & Kr_10bar){
  gcmc_c <- gcmc_Xe1bar
  gcmc_d <- gcmc_Xe10bar
  gcmc_c$Condition <- strrep("Xe_1bar",times = rep(1,nrow(gcmc_c)))
  gcmc_c <- gcmc_c %>% select(y_act, Condition)
  gcmc_d$Condition <- strrep("Xe_10bar",times = rep(1,nrow(gcmc_d)))
  gcmc_d <- gcmc_d %>% select(y_act, Condition)
  gcmc_a <- rbind(gcmc_c, gcmc_d)
  gcmc_e <- gcmc_Kr1bar
  gcmc_f <- gcmc_Kr10bar
  gcmc_e$Condition <- strrep("Kr_1bar",times = rep(1,nrow(gcmc_e)))
  gcmc_e <- gcmc_e %>% select(y_act, Condition)
  gcmc_f$Condition <- strrep("Kr_10bar",times = rep(1,nrow(gcmc_f)))
  gcmc_f <- gcmc_f %>% select(y_act, Condition)
  gcmc_b <- rbind(gcmc_e, gcmc_f)
  stringtopaste <- "XeKr_1_10bar_dist"
  qst_a <- rbind(heat_Xe1bar, heat_Xe10bar)
  qst_b <- rbind(heat_Kr1bar, heat_Kr10bar)
}else if(Xe_1bar & Xe_10bar){
  gcmc_a <- gcmc_Xe1bar
  gcmc_b <- gcmc_Xe10bar
  gcmc_a$Condition <- strrep("Xe_1bar",times = rep(1,nrow(gcmc_a)))
  gcmc_b$Condition <- strrep("Xe_10bar",times = rep(1,nrow(gcmc_b)))
  stringtopaste <- "Xe_1_10bar_dist"
  qst_a <- heat_Xe1bar
  qst_b <- heat_Xe10bar
}else if(Xe_1bar & Kr_1bar){
  gcmc_a <- gcmc_Xe1bar
  gcmc_b <- gcmc_Kr1bar
  gcmc_a$Condition <- strrep("Xe_1bar",times = rep(1,nrow(gcmc_a)))
  gcmc_b$Condition <- strrep("Kr_1bar",times = rep(1,nrow(gcmc_b)))
  stringtopaste <- "XeKr_1bar_dist"
  qst_a <- heat_Xe1bar
  qst_b <- heat_Kr1bar
}else if(Kr_1bar & Kr_10bar){
  gcmc_a <- gcmc_Kr1bar
  gcmc_b <- gcmc_Kr10bar
  gcmc_a$Condition <- strrep("Kr_1bar",times = rep(1,nrow(gcmc_a)))
  gcmc_b$Condition <- strrep("Kr_10bar",times = rep(1,nrow(gcmc_b)))
  stringtopaste <- "Kr_1_10bar_dist"
  qst_a <- heat_Kr1bar
  qst_b <- heat_Kr10bar
}else if(Xe_10bar & Kr_10bar){
  gcmc_a <- gcmc_Xe10bar
  gcmc_b <- gcmc_Kr10bar
  gcmc_a$Condition <- strrep("Xe_10bar",times = rep(1,nrow(gcmc_a)))
  gcmc_b$Condition <- strrep("Kr_10bar",times = rep(1,nrow(gcmc_b)))
  stringtopaste <- "XeKr_10bar_dist"
  qst_a <- heat_Xe10bar
  qst_b <- heat_Kr10bar
}
gcmc_a <- gcmc_a %>% select(y_act, Condition)
gcmc_b <- gcmc_b %>% select(y_act, Condition)

# merge selectivities 
gcmc_select <- rbind(gcmc_a, gcmc_b)
#gcmc_select$sim <- as.factor(gcmc_select$sim)
p <- ggplot(gcmc_select, aes(x=Condition, y = y_act, fill = Condition)) + 
  xlab(paste0("Condition")) + 
  ylab(paste0("GCMC Loading (cm\u00B3/cm\u00B3)")) + 
  ylim(0,300) + 
  theme_classic() + 
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  geom_violin(scale = "width")
#p + stat_summary(fun.y=median, geom="point", shape=23, size=2, color = 'red')
fig <- p + geom_boxplot(width=0.1) + coord_flip()
save_plot(paste0(save_path, stringtopaste, "_Loading.png"), fig, base_width = 10, base_height = 10, dpi = 600)

# make plot for qst
qst <- rbind(qst_a, qst_b)
p <- ggplot(qst, aes(x=Condition, y = -1*Heat_of_Ads, fill = Condition)) + 
  xlab(paste0("Condition")) + 
  ylab(paste0("Heat of Adsorption (KJ/mol)")) + 
  ylim(0,40) + 
  theme_classic() + 
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20)) + 
  geom_violin(scale = "width")
#p + stat_summary(fun.y=median, geom="point", shape=23, size=2, color = 'red')
fig <- p + geom_boxplot(width=0.1) + coord_flip()
save_plot(paste0(save_path, stringtopaste,"_Qst", ".png"), fig, base_width = 10, base_height = 10, dpi = 600)
# also make plots like the ones in Ben Sikora's paper
orig_tobacco_data <-  na.omit(read_textual_data())
colnames(orig_tobacco_data)[colnames(orig_tobacco_data)=="MOF.ID"] <- "id"
Xe_diameter <- 4.1
Kr_diameter <- 3.636
if(Xe_1bar){
Xe_1bar_text <- merge(orig_tobacco_data, gcmc_Xe1bar)
p <- ggplot(Xe_1bar_text, aes(x=LCD, y = y_act)) + 
  xlab(paste0("Largest Cavity Diameter (\305)")) + 
  ylab(paste0("GCMC Loading (cm\u00B3/cm\u00B3)")) + 
  geom_point(color=I("#0070C0"))+ theme_classic() + 
  theme(axis.text = element_text(size = 30),
        axis.title=element_text(size=30,face="bold")) + 
  geom_vline(xintercept = 2*Xe_diameter, linetype="dashed", 
             color = "red", size=1.5) + 
  xlim(0, 80) +
  ylim(0, 300)
save_plot(paste0(save_path, "Xe1bar_Loading_vs_LCD", ".png"), p, base_width = 10, base_height = 10, dpi = 600)
}
if(Xe_10bar){
  Xe_10bar_text <- merge(orig_tobacco_data, gcmc_Xe10bar)
  p <- ggplot(Xe_10bar_text, aes(x=LCD, y = y_act)) + 
    xlab(paste0("Largest Cavity Diameter (\305)")) + 
    ylab(paste0("GCMC Loading (cm\u00B3/cm\u00B3)")) + 
    geom_point(color=I("#0070C0")) + theme_classic() + 
    theme(axis.text = element_text(size = 30),
          axis.title=element_text(size=30,face="bold")) + 
    geom_vline(xintercept = 2*Xe_diameter, linetype="dashed", 
               color = "red", size=1.5) + 
    xlim(0, 80) +
    ylim(0, 300)
  save_plot(paste0(save_path, "Xe10bar_Loading_vs_LCD", ".png"), p, base_width = 10, base_height = 10, dpi = 600)
}
if(Kr_1bar){
  Kr_1bar_text <- merge(orig_tobacco_data, gcmc_Kr1bar)
  p <- ggplot(Kr_1bar_text, aes(x=LCD, y = y_act)) + 
    xlab(paste0("Largest Cavity Diameter (\305)")) + 
    ylab(paste0("GCMC Loading (cm\u00B3/cm\u00B3)")) + 
    geom_point(color=I("#0070C0"))+ theme_classic() + 
    theme(axis.text = element_text(size = 30),
          axis.title=element_text(size=30,face="bold")) + 
    geom_vline(xintercept = 2*Kr_diameter, linetype="dashed", 
               color = "red", size=1.5) + 
    xlim(0, 80) +
    ylim(0, 300)
  save_plot(paste0(save_path, "Kr1bar_Loading_vs_LCD", ".png"), p, base_width = 10, base_height = 10, dpi = 600)
}
if(Kr_10bar){
  Kr_10bar_text <- merge(orig_tobacco_data, gcmc_Kr10bar)
  p <- ggplot(Kr_10bar_text, aes(x=LCD, y = y_act)) + 
    xlab(paste0("Largest Cavity Diameter (\305)")) + 
    ylab(paste0("GCMC Loading (cm\u00B3/cm\u00B3)")) + 
    geom_point(color=I("#0070C0"))+ theme_classic() + 
    theme(axis.text = element_text(size = 30),
          axis.title=element_text(size=30,face="bold")) + 
    geom_vline(xintercept = 2*Kr_diameter, linetype="dashed", 
               color = "red", size=1.5) + 
    xlim(0, 80) +
    ylim(0, 300)
  save_plot(paste0(save_path, "Kr10bar_Loading_vs_LCD", ".png"), p, base_width = 10, base_height = 10, dpi = 600)
}