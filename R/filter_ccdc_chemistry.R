# Parse CCDC MOF subset for elemental composition
# Implements a less stringent version of CoRE MOF's organic-organic bond requirement.
# Command-line tool to be used in the Make build system so the (somewhat expensive)
# CHNOSZ operations are cached.

# hMOF model doesn't work well on some CSD MOFs because they're metal atoms floating in space after solvent removal
# PREREQ: make the target BigData/csd_formula.txt to extract
# chemical information from [Open Babel](https://openbabel.org/docs/dev/Command-line_tools/babel.html#append-option)
# Note, the output has multiple warnings about atom labels like '0', '1', and '1.'.

library(dplyr)
library(readr)
library(purrr)
library(CHNOSZ)

raw_csd_formulas <- read_tsv(
  "BigData/csd_formula.txt",
  col_names = c("cif", "str_formula"),
  col_types = "cc"
)

# We can convert a purrr list to data.frame using dplyr::bind_rows.
# See also https://github.com/tidyverse/dplyr/issues/1676 for the as.list invocation,
# which converts each named vector into its own list
# TODO: consider using `safely` to parse or hide errors
csd_formulas <- 
  raw_csd_formulas$str_formula %>% 
  as.list %>% 
  map(makeup) %>% 
  map(as.list) %>% 
  bind_rows %>% 
  bind_cols(raw_csd_formulas, .)

csd_formulas <- 
  csd_formulas %>% 
  mutate(missing_cnps = (is.na(C) & is.na(N) & is.na(P) & is.na(S)))
no_cnps <- csd_formulas %>%
  filter(missing_cnps) %>% 
  select(cif, str_formula)

no_cnps %>% 
  select(cif) %>% 
  write_tsv("BigData/Robj/csd_mofs_without_cnps.txt", col_names = FALSE)
csd_formulas %>% 
  write_rds("BigData/Robj/csd_formulas.Rds")

