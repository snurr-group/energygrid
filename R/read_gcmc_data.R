library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
read_data <- function(filename, unit_of_ads="cm3overcm3"){ 
gcmc_data <- read_table2(filename)
gcmc_data <- na.omit(gcmc_data)
gcmc_data$ID <- as.character(gcmc_data$ID)
#gcmc_data <- separate(gcmc_data, col = `ID  Temp   Pres`, into = c("ID Temp", "Pres"), sep = "  ")
#gcmc_data <- separate(gcmc_data, col = `ID Temp`, into = c("ID",  "Temp"), sep = " ")
mc_0_ids <- read_table("All_data/IDs_with_mc_0.txt")
#mc_0_ids <- as.character(mc_0_ids)
# for tobacco mofs, we exclude the dubious mc_0 node, all mofs with that ID will be filtered

if (unit_of_ads == "Molec_cm3overcm3"){
  gcmc_data <- gcmc_data %>% select(ID, Temp, Pres, Molec_cm3overcm3, Heat_of_Ads, Heat_fluc)
  # rename the column using a general name
  colnames(gcmc_data)[colnames(gcmc_data) == "Molec_cm3overcm3"] <- "Uptake"
  unit_for_plot <<- "cm\u00B3/cm\u00B3"
  
} else if (unit_of_ads == "Mill_gram"){
  gcmc_data <- gcmc_data %>% select(ID, Temp, Pres, Mill_gram, Heat_of_Ads, Heat_fluc)
  unit_for_plot <<- "mg/gram framework"
  colnames(gcmc_data)[colnames(gcmc_data) == "Mill_gram"] <- "Uptake"

} else if (unit_of_ads == "Cm3_gram"){
  gcmc_data <- gcmc_data %>% select(ID, Temp, Pres, Cm3_gram, Heat_of_Ads, Heat_fluc)
  unit_for_plot <<- "cm\u00B3/gram framework"
  colnames(gcmc_data)[colnames(gcmc_data) == "Cm3_gram"] <- "Uptake"
  
} else if (unit_of_ads == "Mol_kg"){
  gcmc_data <- gcmc_data %>% select(ID, Temp, Pres, Mol_kg, Heat_of_Ads, Heat_fluc)
  unit_for_plot <<- "mol/kg framework"
  colnames(gcmc_data)[colnames(gcmc_data) == "Mol_kg"] <- "Uptake"
  
} else{
  unit_for_plot <<- " "
}
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

