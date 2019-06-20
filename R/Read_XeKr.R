library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
library(readxl)
source("R/refined_bins_calc.R")
source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")
source("R/save_train_test_data.R")
library(randomForest)
gcmc_data <- read_table2("All_data/KrXe_mix_273_1bar.txt")

gcmc_data <- na.omit(gcmc_data)
gcmc_data$Selectivity <- (gcmc_data$Xe_uptake/0.2)/(gcmc_data$Kr_uptake/0.8)
gcmc_data <- gcmc_data %>% filter(Selectivity < 50)
# convert IDs to characters
gcmc_data$ID <- as.character(gcmc_data$ID)

# 
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


hmof_h2_grid <- read_rds("Kr_greater_range.Rds")
gcmc_data <- mutate(gcmc_data, id=ID)
new_grid <- hmof_h2_grid[hmof_h2_grid[, "id"] %in% gcmc_data$id,] # select out the grids with id same as gcmc_data id
binbounds <- automatic_bins(new_grid,1000000) # for non-uniform bins
ch4_binspec <- list("from" = -40, to = 40, bounds = binbounds)
DATA_SPLIT <- 1000 # Number of data points used in training data, from setup_data.R
p_ch4_sets <- partition_data_subsets(hmof_h2_grid, gcmc_data, DATA_SPLIT)

p_ch4_vol <- gcmc_data %>%
  mutate(g.L = Selectivity) %>%
  run_model_on_partitions(p_ch4_sets, ., ch4_binspec, plot_units="unitness", db_name = "tobacco")

p_ch4_vol$plots$parity_testing + scale_x_continuous(limits = c(0,30)) + scale_y_continuous(limits = c(0,30))+ xlab("GCMC selectivity (unitless)") + ylab("Predicted selectivity (unitless)")


export_data(p_ch4_sets$training, rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), "tobacco_XeKr_train.csv", ch4_binspec)
export_data(p_ch4_sets$testing, rename(gcmc_data %>% mutate(g.L = Selectivity), y_act=g.L), "tobacco_XeKr_test.csv", ch4_binspec)


train_data <- read_csv("tobacco_XeKr_train.csv")
test_data <- read_csv("tobacco_XeKr_test.csv")

model <- randomForest(x = train_data %>% select(-y_act), y = train_data$y_act, ntree = 500)

tested <- predict(model, test_data)

tested_rmse <- postResample(pred = tested, obs = test_data$y_act) # original model
gg<- qplot(x = model$y, y = model$predicted) + geom_point() + geom_abline(slope = 1, intercept = 0) + scale_x_continuous(limits = c(0,50)) + scale_y_continuous(limits = c(0,50)) + coord_equal()
gg
