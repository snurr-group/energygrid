# Functions to read relevant data
# See also data loading procedures from the old RF notebooks: C:\Users\Benjamin\Git\ManyNets\Data\RandomForest

library(dplyr)
library(stringr)

base_dir <- "Data/NNTraining/"


id_from_cif <- function(filename) {
  as.integer(str_sub(filename, 2, -5))
}

load_gcmc <- function(filename) {
  gcmc_data <- read.delim(filename,
                          header = TRUE,
                          skip = 1,  # Give a generic header row at the top
                          colClasses = c(rep("integer", 7), rep("numeric", 14)),
                          na.strings=c("", "null", "#NAME?", "#VALUE!")
                          )
  rename(gcmc_data, id = ID)
}

load_timing <- function(filename) {
  timing_data <- read.delim(filename,
                            header = FALSE,
                            colClasses = c("character", "numeric"),
                            col.names = c("filename", "run.time")
                            )
  timing_data$id <- id_from_cif(timing_data$filename)
  timing_data$filename <- NULL
  timing_data
}

load_metric <- function(filename) {
  metric_data <- read.delim(filename,
                            header = FALSE,
                            colClasses = c("character", "numeric"),
                            col.names = c("filename", "metric")
                            )
  metric_data$filename <- str_sub(metric_data$filename, 25, -2)
  metric_data$id <- id_from_cif(metric_data$filename)
  metric_data$filename <- NULL
  metric_data
}

load_agnostic <- function(filename) {
  agnostic_data <- read.delim(filename,
                              header = TRUE,
                              colClasses = c("integer", "numeric", "numeric", "numeric"),
                              col.names = c("id", "ag.geo.VF", "ag.bind.frac", "ag.sparse.frac")
                              )
  agnostic_data
}

# TODO: clean up the bottom of this, probably putting it in another script, notebook, etc.
gcmc_data <- load_gcmc(paste0(base_dir, "condensed-master-list.tsv.gz"))
timing_data <- load_timing(paste0(base_dir, "timer.txt"))
metric_data <- load_metric(paste0(base_dir, "Metric.txt"))
agnostic_data <- load_agnostic(paste0(base_dir, "10k-hMOFs-agnosticscreening.txt"))
# <- right_join(gcmc_data, metric_data, by = c("ID"="id")) %>%
scotty_data <- right_join(gcmc_data, metric_data) %>% right_join(timing_data) %>% left_join(agnostic_data)
# ggplot(scotty_data, aes(hv, metric)) + geom_point() + theme_bw(16)
scotty_no_na_data <- na.omit(scotty_data)
scotty_na_data <- scotty_data[attr(scotty_no_na_data, "na.action"),]
write.table(scotty_data,
            paste0(base_dir, "scotty_test.tsv"),
            sep = "\t",
            row.names = FALSE,
            na = "#VALUE!"
)

