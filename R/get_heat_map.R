# get the heat map (E vs grad.E)
library(tidyverse)
library(stringr)
library(R.utils)
library(manipulate)
library(slam) # for calculating row norm

source("R/hist2d_breaks.R")
QUICK_TEST <- TRUE  # Set to true to do a "practice run" instead of all of the files

energy <- read_table("Energies/1.grid", col_names = "V1", col_types = "d")$V1  # faster than read_tsv for similar purposes
gradients <- read_table2("Gradients/1.grid", col_names = FALSE)

eng_grad <- gradients
eng_grad$energy <- energy
new_grad <- eng_grad[eng_grad$X1 != 1e+23,]
new_grad <- new_grad/120.28
# calculate the norms of the gradients
new_grad$grad_size <- new_grad %>% select(-energy) %>% row_norms()
new_grad <- new_grad %>% select(-X1,-X2,-X3)
library(gplots)
# latest_grad <- new_grad %>% filter(energy < 0)
# ggplot(latest_grad, aes(x = energy, y = grad_size)) + stat_bin_2d()
binbounds$upper <- append(binbounds$upper, -Inf, 0)
binbounds$upper <- append(binbounds$upper, Inf, length(binbounds$upper))
energy_breaks <- binbounds$upper

grad_breaks <- c(0, 7, 350, 11000, Inf)

abc <- hist2d_break(x = new_grad$energy, y = new_grad$grad_size, energy_breaks, grad_breaks, show = FALSE)
