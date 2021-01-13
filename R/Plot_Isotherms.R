library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(plotly)
library(cowplot)
library(readxl)
#gcmc_data <- read_table2("file.txt")
#
molecule <- "Ethane"
if(molecule == "Ethane"){
  gcmc_data <- read_excel("All_data/SI.xlsx", sheet = "Ethane_Isotherm")
  VP=40e5 # vapor pressure
  maxy = 300 # maximum loading range
}else if(molecule == "Propane"){
  gcmc_data <- read_excel("All_data/SI.xlsx", sheet = "Propane_Isotherm")
  VP=10e5
  maxy = 300
}
if(sum(grepl("Error", names(gcmc_data))) == 1){
  gcmc_data <- gcmc_data %>% select(-`GCMC_Uptake_Error_Bar [cm3/cm3]`)
}
names(gcmc_data) <- c("ID", "Temp", "Pres", "Uptake")
save_path <- "Results/Ethane_Propane_Isotherms/"
if (!dir.exists(save_path)){
  dir.create(save_path)
}
gcmc_data <- na.omit(gcmc_data)

gcmc_data$ID <- as.character(gcmc_data$ID)
gg <- gcmc_data %>% ggplot(aes(x = Pres/1e5, y = Uptake, color = ID)) + 
  geom_vline(xintercept=0.1*VP/1e5, linetype = "dashed", size = 1.5) +
  geom_vline(xintercept=0.5*VP/1e5, linetype = "dashed", size = 1.5) +
  geom_vline(xintercept=VP/1e5, linetype = "dashed", size = 1.5) +
  annotate(geom = 'text', label = molecule, x = 0.002*log10(VP), y = 0.8*maxy, 
           size = 20, fontface = 'bold') +
  geom_line(size = 1) + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                                      labels = scales::trans_format("log10", 
                                                                    scales::math_format(10^.x))) + 
                                        ylab("GCMC loading (cm\u00B3/cm\u00B3)") + xlab("Pressure (Bar)")
gg + theme_set(theme_cowplot(10))
gg <- gg +theme(axis.text=element_text(size=30),
          axis.title=element_text(size=30,face="bold"), legend.position = "none")

save_plot(paste0(save_path, molecule, "_isotherm.png"),gg, base_width = 10, base_height = 10, dpi = 600)