library(readr)
library(readxl)
library(dplyr)
library(stringr)

### From setup_data.R
read_textual_data <- function(option = "ToBaCCo", keep = FALSE, keepVF = FALSE){
  if (option == "ToBaCCo"){
    h2_types <- paste0("y.h2.", c("g.L", "mol.kg", "wtp"))  # Prefix data with its source (Yamil, Scotty, etc.)
    orig_tobacco_data <- read_xlsx(
      "Data/CrystGrowthDesign_SI_without_tob_cleaner_py.xlsx",
      sheet = "data",
      skip = 3, na = "inf",
      col_names = c(
        "MOF.ID",
        "VF", "VSA", "GSA", "PLD", "LCD",
        paste0(h2_types, ".100.77"),
        paste0(h2_types, ".100.130"),
        paste0(h2_types, ".100.200"),
        paste0(h2_types, ".100.243"),
        paste0("h2.qst.6.", c(77, 130, 200, 243)),
        paste0("y.ch4.", rep(c("v.v.", "mg.g."), 2), c("100.298", "100.298", "65.298", "65.298")),
        "ch4.qst.6.298",
        paste0("y.xe.kr.1.", c("xe", "kr", "select")),
        paste0("y.xe.kr.5.", c("xe", "kr", "select")),
        "topology",
        paste0("n1.", c("sym", "character", "ID")),
        paste0("n2.", c("sym", "character", "ID")),
        "cbb.ID"
      )
    )
    tobacco_codes <- read_table2("All_Data/mofs_map_without_tobacco_cleaner_py.dat", col_names = c("MOF.ID", "python.id"), col_types="ic")
    
    orig_tobacco_data <- orig_tobacco_data %>% left_join(tobacco_codes, by="MOF.ID")
    
    orig_tobacco_data <- orig_tobacco_data %>% 
      mutate("n1.name" = str_c("sym", n1.sym, ifelse(n1.character=="organic", "on", "mc"), n1.ID, sep="_")) %>% 
      mutate("n2.name" = str_c("sym", n2.sym, ifelse(n2.character=="organic", "on", "mc"), n2.ID, sep="_"))

  
  if (keep == FALSE){
    # just select the structrual data
    original_data <- orig_tobacco_data %>% select(MOF.ID, VF, VSA, GSA, PLD, LCD)
  } else{
    original_data <- orig_tobacco_data
  }
  if (keepVF == FALSE){
    original_data <- original_data %>% select(-VF)
  }
  }
  if (option == "CoRE"){
    # read textural data from CoRE_MOF database
    original_CoRE_data <- read_xlsx("All_Data/COREMOF_Textual_data.xlsx")
    # just select the data set same as the tobaccos
    # vf, vsa, gsa, pld, lcd
    original_data <- original_CoRE_data %>% select(MOF.ID, LCD, PLD, ASA_m2_cm3, ASA_m2_g, AV_VF)
    names(original_data) <- c("MOF.ID", "LCD", "PLD", "VSA", "GSA", "VF")
  }
  
  original_data
}
