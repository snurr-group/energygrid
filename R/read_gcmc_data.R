library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
#gcmc_data <- read_table2("new_propane.txt")
#gcmc_data <- read_table("new_propane.txt")
#gcmc_data <- read_table("propane_60bar.txt")
#gcmc_data <- read_table("propane_1bar_60K_Cycle.txt")
read_data <- function(filename){ 
gcmc_data <- read_table(filename)
#gcmc_data <- read_table("All_data/propane_1000pa_140K.txt")
#gcmc_data <- read_table("tobacco_propane_5k.txt")
#gcmc_data <- read_table("propane_9000_cycle.txt")
#gcmc_data <- read_table("propane_1000pa.txt")
gcmc_data <- separate(gcmc_data, col = `ID  Temp   Pres`, into = c("ID Temp", "Pres"), sep = "  ")
gcmc_data <- separate(gcmc_data, col = `ID Temp`, into = c("ID",  "Temp"), sep = " ")


gcmc_data <- na.omit(gcmc_data)
gcmc_data
}
# low_loading <- read_table2("All_data/propane_1000pa_140K.txt")
# low_loading_heat  <- low_loading %>% select(ID, Heat_of_Ads)
# names(low_loading_heat)[2] <- "Qst_low_loading"
# low_loading_heat <- na.omit(low_loading_heat)
# hggg <- read_table2("All_data/HGGG_data_per_molec_1bar_28K.txt")
# #hggg <- read_table2("HGGG_data_1bar_28K_float.txt")
# hggg <- na.omit(hggg)
# #write.xlsx(gcmc_data, "propane_1000pa.xlsx")
# names(hggg)[1] <- "ID"
# gcmc_data <- merge(gcmc_data, hggg, by="ID")
# gcmc_data <- merge(gcmc_data, low_loading_heat, by="ID")

