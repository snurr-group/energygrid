# get the heat map (E vs grad.E)
library(tidyverse)
library(stringr)
library(R.utils)
library(manipulate)
library(slam) # for calculating row norm
library(gplots)
source("R/hist2d_breaks.R")

heat_maps <- function(MOF_ID, binbounds) {
  energy <- read_table(paste0("GenericMOF_Energies/", MOF_ID, ".grid"), col_names = "V1", col_types = "d")$V1  # faster than read_tsv for similar purposes
  gradients <- read_table2(paste0("GenericMOF_Gradients/", MOF_ID, ".grid"), col_names = FALSE)
  
  eng_grad <- gradients
  eng_grad$energy <- energy
  new_grad <- eng_grad[eng_grad$X1 != 1e+23,]
  new_grad <- new_grad/120.28 # converting factor from Kelvin to KJ/mol
  # calculate the norms of the gradients
  new_grad$grad_size <- new_grad %>% select(-energy) %>% row_norms()
  new_grad <- new_grad %>% select(-X1,-X2,-X3)
  
  # latest_grad <- new_grad %>% filter(energy < 0)
  # ggplot(latest_grad, aes(x = energy, y = grad_size)) + stat_bin_2d()
  a <- append(binbounds$upper, -Inf, 0)
  a <- append(a, Inf, length(a))
  energy_breaks <- a
  
  #grad_breaks <- c(0, 7, 350, 11000, Inf) # taken from the quantile values of 1.grid from tobacco
  #grad_breaks <- c(0, 1.5, 2, 3, 50, 350, Inf) # finer values
  
  grad_breaks <- c(0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 20, 30, 40, 50, 350, Inf)
  abc <- hist2d_break(x = new_grad$energy, y = new_grad$grad_size, energy_breaks, grad_breaks, show = FALSE)
  # convert matrix to a 1D dataframe with MOFID
  abc$counts <- abc$counts/sum(abc$counts)
  saveRDS(abc$counts, paste0("HeatMaps/", MOF_ID, ".rds"))
  ab <- abc$counts #get the counts from the 2d histogram
  colnames(ab) <- seq(from = 1, to = 15, by = 1)
  rownames(ab) <- seq(from = 1, to = nrow(ab), by = 1)
  abv <- as.vector(ab)
  
  names(abv) <- seq(from = 1, to = length(abv), by = 1)
  abvt <- t(abv)
  abvtf <- as.data.frame(abvt)
  # finally, append MOF_ID to that vector
  abvtf$ID <- MOF_ID
  abvtf
}

namess <- list.files(path="Energies/") %>% str_remove(., ".grid")

count = 0
for (MOF_ID in namess){
  abvtf <- heat_maps(MOF_ID, binbounds)
  count = count + 1
  if (count == 1) {
    abvtfs <- abvtf
  } else{
    abvtfs <- rbind(abvtfs, abvtf)
  }
}
saveRDS(abvtfs, "heat_map_more_grads.rds")

