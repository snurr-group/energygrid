# From tobacco_exploration.Rmd, modified Dec 6, 2017
# Note: the working directory needs to be the root of the Git repository
# (e.g. the parent directory of this script)

### From setup block

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)

library(purrr)
library(readxl)
library(magrittr)

library(testthat)

library(caret)
library(glmnet)

#theme_set(theme_bw(base_size=16) + theme(aspect.ratio=1))  # Default ggplot2 params, now overridden with cowplot


### From load.data code block

h2_types <- paste0("y.h2.", c("g.L", "mol.kg", "wtp"))  # Prefix data with its source (Yamil, Scotty, etc.)
tobacco_data <- read_xlsx(
  "BigData/CrystGrowthDesign_SI.xlsx",
  sheet = "data",
  skip = 3, na = "inf",
  col_names = c(
    "MOF.ID",
    "vf", "vsa", "gsa", "pld", "lcd",
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
tobacco_codes <- read_table2("BigData/mofs_map.dat", col_names = c("MOF.ID", "python.id"), col_types="ic")

tobacco_data <- tobacco_data %>% left_join(tobacco_codes, by="MOF.ID")

tobacco_data <- tobacco_data %>% 
  mutate("n1.name" = str_c("sym", n1.sym, ifelse(n1.character=="organic", "on", "mc"), n1.ID, sep="_")) %>% 
  mutate("n2.name" = str_c("sym", n2.sym, ifelse(n2.character=="organic", "on", "mc"), n2.ID, sep="_"))

# Translate between Scotty's ID and the traditional ToBaCCo filename
scotty_codes <- read_tsv("BigData/tobacco-20171114/key.tsv", col_names = c("filename", "id"), col_types="cc")
scotty_codes <- scotty_codes %>% mutate(python.id = str_sub(filename, 1, -5))  # trim .cif suffix
tobacco_data <- tobacco_data %>% inner_join(scotty_codes, by="python.id")

# Add in Scotty's LJ-only GCMC
tobacco_data <-
  read_tsv(
    "BigData/Emails/tobacco-gcmc-20171214/converted_gcmc_5bar_160k.tsv",
    col_names = c("tob.num", "id", "lj.h2.v.v.5.160", "lj.h2.err.v.v.5.160"),
    col_types = "icnn"
  ) %>% 
  select(-tob.num) %>% 
  mutate(`lj.h2.g.L.5.160` = `lj.h2.v.v.5.160` * 2.0 / 22.4) %>% 
  mutate(`lj.h2.err.g.L.5.160` = `lj.h2.err.v.v.5.160` * 2.0 / 22.4) %>% 
  left_join(tobacco_data, ., by="id")
tobacco_data <-
  read_tsv(
    "BigData/Emails/tobacco-gcmc-20171214/converted_gcmc_100bar_77k.tsv",
    col_names = c("tob.num", "id", "lj.h2.v.v.100.77", "lj.h2.err.v.v.100.77"),
    col_types = "icnn"
  ) %>% 
  select(-tob.num) %>% 
  mutate(`lj.h2.g.L.100.77` = `lj.h2.v.v.100.77` * 2.0 / 22.4) %>% 
  mutate(`lj.h2.err.g.L.100.77` = `lj.h2.err.v.v.100.77` * 2.0 / 22.4) %>% 
  left_join(tobacco_data, ., by="id")

# Now with Feynman-Hibbs
tobacco_data <-
  read_table2(
    "BigData/Emails/tobaccoFH-20180103/tobacco-2bar-FH-comb-volume.txt",
    col_names = c("tob.num", "id", "fh.h2.v.v.2.77", "fh.h2.err.v.v.2.77"),
    col_types = "icnn"
  ) %>% 
  select(-tob.num) %>% 
  mutate(`fh.h2.g.L.2.77` = `fh.h2.v.v.2.77` * 2.0 / 22.4) %>% 
  mutate(`fh.h2.err.g.L.2.77` = `fh.h2.err.v.v.2.77` * 2.0 / 22.4) %>% 
  left_join(tobacco_data, ., by="id")
tobacco_data <- 
  read_table2(
    "BigData/Emails/tobaccoFH-20180103/tobacco-100bar-FH-comb-volume.txt",
    col_names = c("tob.num", "id", "fh.h2.v.v.100.77", "fh.h2.err.v.v.100.77"),
    col_types = "icnn"
  ) %>% 
  select(-tob.num) %>% 
  mutate(`fh.h2.g.L.100.77` = `fh.h2.v.v.100.77` * 2.0 / 22.4) %>% 
  mutate(`fh.h2.err.g.L.100.77` = `fh.h2.err.v.v.100.77` * 2.0 / 22.4) %>% 
  left_join(tobacco_data, ., by="id")



#DATA_SPLIT <- c(0.4, 0.4, 0.2)  # Fractional split between hyperparams, training, and test data
DATA_SPLIT <- 1000  # Number of data points used in training data
R_GAS_KJ <- 8.314 / 1000
source("R/regression_and_features.R")  # Get regression utilities


### From load_low_p code block
# Also import the new low pressure data

# Update 12/14: import Scotty's (revised, mostly complete) GCMC results for H2 ToBaCCo, low P
# Also a few columns at higher P for verification

low_p_h2_data <- read_tsv(
  "BigData/Emails/tobacco-gcmc-20171214/converted_gcmc_2bar_77k.tsv",
  col_names = c("tob.num", "id", "nofh.h2.v.v.2.77", "nofh.h2.err.v.v.2.77"),
  col_types = "icnn"
)

low_p_h2_data <- low_p_h2_data %>% 
  mutate(`nofh.h2.g.L.2.77` = `nofh.h2.v.v.2.77` * 2.0 / 22.4) %>% 
  mutate(`nofh.h2.err.g.L.2.77` = `nofh.h2.err.v.v.2.77` * 2.0 / 22.4)
tobacco_data <- tobacco_data %>% left_join(low_p_h2_data, by="id")
# Warning: ToBaCCo data will contain rows with missing deliverable capcity.
# Training and testing code will need to handle these NAs.
tobacco_data <- tobacco_data %>% mutate(redo.h2.deliv.77 = y.h2.g.L.100.77 - nofh.h2.g.L.2.77)


# 77 K / 6 bar from EES paper, Yamil, 2017-11-15
# Similar format to the CrystGrowthDesign spreadsheet above, but with fewer columns to import
low_p_h2_data_raw <- read_xlsx(
  "BigData/Emails/tobacco-yamil-20171115/EES-SI-07-18-2016.xlsx",
  sheet = "data",
  skip = 3, na = "inf",
  col_names = c(
    "MOF.ID",
    "vf", "vsa", "gsa", "pld", "lcd",
    paste0(h2_types, ".100.77"),
    paste0(h2_types, ".5.160"),
    paste0(h2_types, ".6.77"),
    paste0("h2.qst.", c("6.77", "5.160")),
    "topology",
    paste0("n1.", c("sym", "character", "ID")),
    paste0("n2.", c("sym", "character", "ID")),
    "cbb.ID"
  )
)
# tobacco_data[1,1] = 234052409235  # Tested the next piece of test code by messing up one of the ID's
# Ensure data consistency: that compositions for the first MOF.ID column are the same
# anti_join will look for rows that mismatch within these columns
temp_mof_id_mismatches <- low_p_h2_data_raw %>%
  anti_join(tobacco_data, by=c(
    "MOF.ID", "topology",
    paste0("n1.", c("sym", "character", "ID")),
    paste0("n2.", c("sym", "character", "ID"))
  )) %>% 
  nrow
expect_equal(temp_mof_id_mismatches, 0)  # Now we can safely use MOF.ID to join
# Alternatively, we could have just performed an inner join on all of these columns
low_p_h2_data <- low_p_h2_data_raw %>% 
  select(c(1, 10:15, 17))  # MOF.ID and the new H2 data (not 100 bar cryo)
tobacco_data <- tobacco_data %>% inner_join(low_p_h2_data, by="MOF.ID")
#tobacco_data <- tobacco_data %>% mutate(h2.deliv.77 = h2.g.L.100.77 - h2.g.L.6.77)

# Let's also grab CH4 data at 6 bar (noting that this file also provides H2 at 6 bar, 243 K)
raw_ch4_6_bar <- read_xlsx(
  "BigData/Emails/tobacco-yamil-20171115/6bar.xlsx",
  sheet = "CH4",
  skip = 1, na = "inf",
  col_names = c("MOF.ID", paste0("y.ch4.", c("v.v", "mg.g"), ".6.298"))
)
tobacco_data <- tobacco_data %>% inner_join(raw_ch4_6_bar, by="MOF.ID")


### From load_grids code block
if (!exists("raw_grids_h2")) {
  raw_grids_h2 <- read_rds("BigData/Robj/tobacco_h2.Rds")
}


### From partition_data code block

tob_y_to_join <- tobacco_data %>% 
  mutate(g.L = fh.h2.g.L.100.77 - fh.h2.g.L.2.77) %>% 
  select(id, g.L) %>% 
  filter(!(is.na(g.L) | g.L < 0))  # Clean up unphysical MOFs
complete_ids <- raw_grids_h2 %>% 
  select(id) %>% 
  unique() %>% 
  inner_join(tob_y_to_join, by="id") %>% 
  select(id) %>% 
  unlist

tob_y_to_join <- tob_y_to_join %>% filter(id %in% complete_ids)
grids_h2 <- raw_grids_h2 %>% filter(id %in% complete_ids)

# Remove empty files when the RASPA simulation crashed
crashed_grid_ids <-
  grids_h2 %>% 
  group_by(id) %>% 
  summarize(total_pts = sum(counts)) %>% 
  ungroup %>% 
  filter(total_pts == 0) %>% 
  .$id
grids_h2 <- grids_h2 %>% filter(!(id %in% crashed_grid_ids))

tob_hist_sets <- partition_data_subsets(grids_h2, tob_y_to_join, DATA_SPLIT)

### From model_proof_of_concept code block (partial: removed test fits, but kept importing R functions)

source("R/plot_hists.R")
source("R/regression_and_features.R")
source("R/get_energy_stats.R")
source("R/plot_diagnostics.R")

source("R/load_data.R")
# Also load CH4 data for the hMOFs from Scotty.  We don't have this data from Wilmer's work, since they were only concerned with absolute uptake at that point.
raw_2500_cols <- c("kJ.mol", "cm3.cm3", "g.kg")
raw_ch4_hmof <- read_tsv(
  "BigData/2500hmof-data/2500hmofdata.tsv",
  skip = 2,
  na = c("", "#NAME?", "-"),
  col_names = c(
    "id",
    paste("h2", raw_2500_cols, "2.160", sep="."),
    paste("ch4", raw_2500_cols, "5_8.298", sep="."),
    paste("ch4", raw_2500_cols, "65.298", sep=".")
  ),
  col_types = "innnnnnnnn" # TODO update simulations or grids for the full database
)
# Note: some columns only have a "#NAME?" flag in the kJ/mol column, but zero in the other two, so clean them up here:
# deletes rows without any GCMC and erroneous zero entries associated with a corresponding NA
raw_ch4_hmof <- raw_ch4_hmof %>% 
  select(-starts_with("h2")) %>% 
  na.omit() %>% 
  full_join(na.omit(select(raw_ch4_hmof, id, starts_with("h2"))), by="id")
raw_ch4_hmof <- raw_ch4_hmof %>% 
  mutate(h2.g.L.2.160 = h2.cm3.cm3.2.160 * 2.0 / 22.4)
gcmc_data <- gcmc_data %>% 
  left_join(raw_ch4_hmof, by="id")

hmof_h2_grid <- read_rds("BigData/Robj/hmof_h2.Rds")
hmof_h2_grid <- hmof_h2_grid %>%
  mutate(str_id = str_sub(id, 6, -1)) %>% # strip the hMOF- prefix
  mutate(id = as.integer(str_id)) %>% 
  select(-str_id)
hmof_y_to_join <- gcmc_data %>% 
  select(id, g.L) %>% 
  filter(!(is.na(g.L) | g.L < 0))
hmof_hist_sets <- partition_data_subsets(hmof_h2_grid, hmof_y_to_join, DATA_SPLIT)
# Note: since we're starting with 2500 hMOFs, 40% for training is still 1000, so that works out great.


### From hmof_betas code block (partial: mostly wanted a color.R import)

source("R/color.R")
base_hist <- hmof_h2_grid %>% 
  filter(id == 71) %>% 
  plot_hist_bins(default_binspec)


### From h2_hmof_2bar code block (partial for hmof_methane block below)

p_2bar_data <- na.omit(gcmc_data)  # probably overzealous, but it'll remove all of the junk


### From hmof_methane code block (mostly wanted the CH4 setup)

# Let's use the data sets from the previous block.  There are some cases where we have methane data but not H2, and vice versa, but let's go with the MOFs in common for simplicitly.
raw_hmof_grids_ch4 <- read_rds("BigData/Robj/hmof_ch4.Rds") %>% rename(dirname = id)
p_ch4_grids <- raw_hmof_grids_ch4 %>% 
  mutate(id = as.integer(str_sub(dirname, 6))) %>% 
  filter(id %in% p_2bar_data$id)

ch4_binspec <- c(from=-26, to=0.0, step=2.0, width=2.0)

p_ch4_sets <- partition_data_subsets(p_ch4_grids, p_2bar_data, DATA_SPLIT)


### From ch4_tobacco code block (partial, again for loading data)

# Repeat the CH4 hMOF analysis for ToBaCCo.  Does it behave like H2 deliverable capacity?
# First, attempt training/testing only on ToBaCCo data
grids_ch4_tobacco <- read_rds("BigData/Robj/tobacco_ch4.Rds")
tob_ch4_sets <- partition_data_subsets(grids_ch4_tobacco, tobacco_data, DATA_SPLIT)


### From ccdc_applicability code block (partial for CCDC data)

# Derived from notebook 20171106_immediate_followup.Rmd
ccdc_h2_grids <- read_rds("BigData/Robj/ccdc_hist_vals.Rds")
ccdc_gcmc <- 
  read_tsv(
    "BigData/Emails/ccdc-random-20180111/completed_20180518_ccdc_vol.tsv",
    skip = 3,
    na = c("", "null", "#NAME?", "#VALUE!", "-"),
    col_names = c("id", paste("fh.h2.g.L", c(2, 5, 5, 100), c(77, 77, 160, 77), sep=".")),
    col_types = "cnnnn"
    ) %>% 
  mutate_at(vars(starts_with("fh.h2")), funs(. * 2.0 / 22.4)) %>% 
  # See https://stackoverflow.com/questions/39279724/use-mutate-at-to-change-multiple-column-types
  mutate(g.L = fh.h2.g.L.100.77 - fh.h2.g.L.2.77) %>% 
  na.omit()
ccdc_h2_160k_grids <- read_rds("BigData/Robj/ccdc_hist_160k.Rds")


# Import 160 K, 5 bar data for the 2500 hMOFs.  Use it to train (and validate) a model for CCDC screening.
p_160k_5bar_data <- 
  read_table2(
    "BigData/Mateo/EnergyGrid/GCMC-160k-2500hMOF/volume-uptake_combined.txt",
    col_names = c("id", "filebasename", "fh.h2.v.v.5.160", "fh.h2.err.v.v.5.160"),
    col_types = "icnn"
  ) %>% 
  select(-filebasename) %>% 
  mutate(`fh.h2.g.L.5.160` = `fh.h2.v.v.5.160` * 2.0 / 22.4) %>% 
  mutate(`fh.h2.err.g.L.5.160` = `fh.h2.err.v.v.5.160` * 2.0 / 22.4) %>% 
  left_join(gcmc_data, by="id")

# Set up datasets that combine hMOFs and ToBaCCo in a 50:50 proportion
# Can't do this merely by sampling 1000 from the combined database, because then the proportion will be different
mixed_h2_hist_sets <- list(
  hmof = partition_data_subsets(
    mutate(hmof_h2_grid, id=paste0("h", id)),
    mutate(hmof_y_to_join, id=paste0("h", id)),
    DATA_SPLIT/2
    ),
  tob = partition_data_subsets(grids_h2, tob_y_to_join, DATA_SPLIT/2)
  ) %>% 
  transpose %>% 
  map(~bind_rows(.x))
expect_equal(mixed_h2_hist_sets$training$id %>% unique %>% length, DATA_SPLIT)
expect_equal(mixed_h2_hist_sets$training %>% select(id) %>% unique %>% filter(str_detect(id, "^h")) %>% nrow, DATA_SPLIT/2)
mixed_h2_y_to_join <- bind_rows(
  mutate(hmof_y_to_join, id=paste0("h", id)),
  tob_y_to_join
  )


