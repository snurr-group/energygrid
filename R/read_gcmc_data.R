library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
read_data <- function(filename){ 
gcmc_data <- read_table(filename)
gcmc_data <- na.omit(gcmc_data)
gcmc_data <- separate(gcmc_data, col = `ID  Temp   Pres`, into = c("ID Temp", "Pres"), sep = "  ")
gcmc_data <- separate(gcmc_data, col = `ID Temp`, into = c("ID",  "Temp"), sep = " ")
mc_0_ids <- read_table("All_data/IDs_with_mc_0.txt")
#mc_0_ids <- as.character(mc_0_ids)
# for tobacco mofs, we exclude the dubious mc_0 node, all mofs with that ID will be filtered

gcmc_data <- gcmc_data %>% filter(!ID %in% mc_0_ids$IDs_with_mc_0)

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

