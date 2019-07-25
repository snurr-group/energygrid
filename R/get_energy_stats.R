# Gets statistics on calculated energies from the scotty_map.py code

QUICK_TEST <- FALSE  # Set to true to do a "practice run" instead of all of the files


energy_stats <- function(data_dir, bin_width, min_max) {
  bins <- seq(from = bin_width * min_max[1],
              to   = bin_width * (min_max[2] + 1),
              by   = bin_width
  )
  stats_fcn <- function(energy) {
    # First cap the maximum energy so that the last bin (energy > max) includes the high end
    energy[energy >= (bin_width * (min_max[2]) + 1)] <- bin_width * (min_max[2] + 0.5)
    # Do the same for the low end.  If we care about the strongly attractive region, check
    # that both of the lower two bins are unoccupied.  Otherwise they're lumped in with the smallest.
    energy[energy <= (bin_width * min_max[1])] <- bin_width * (min_max[1] + 0.5)
    # Now that preprocessing is done, run the histogram binning
    raw_hist <- hist(energy, breaks = bins, plot = FALSE)
    cbind(lower = raw_hist$breaks[1:(length(raw_hist$breaks)-1)],
          upper = raw_hist$breaks[2:length(raw_hist$breaks)],
          counts = raw_hist$counts
    )
  }
  num_rows <- length(bins)-1
  # Run a stats_fcn on all raspa grid files
  # Returns a data.frame, where each grid's stats has num_rows by df_prototype entries
  files <- list.files(data_dir)
  if (QUICK_TEST) {
    files <- files[1:500]  # Debugging trick to only run the first 100 folders
  }
  num_cifs <- length(files)
  if (all(str_detect(files, "\\.grid$"))) {
    files <- str_sub(files, 0, -6)
  }
  
  # Preallocate a df using an R.utils function: http://r.789695.n4.nabble.com/idiom-for-constructing-data-frame-td4705353.html
  # Names for the dataframe will automatically be copied from the prototype
  df_prototype <- rep("integer", 3)  # Lower bound, upper bound, count
  names(df_prototype) <- c("lower", "upper", "counts")
  all_stats <- dataFrame(df_prototype, nrow = num_cifs*num_rows)
  
  write(paste0("Compiling statistics from the energy files in ", data_dir, "..."), "")
  pb <- txtProgressBar(min = 1, max = num_cifs, style = 3)  # Add a progress bar, courtesy of https://www.r-bloggers.com/r-monitoring-the-function-progress-with-a-progress-bar/
  current_row = 1
  for (raspa_grid in files) {
      energy_file <- file.path(data_dir, paste0(raspa_grid, ".grid"))
    
    # File I/O (next line) is the bottleneck.  read_table is faster than read_tsv
    energy <- read_table(energy_file, col_names = "V1", col_types = "d")$V1  # faster than read_tsv for similar purposes
    end_row <- current_row + num_rows - 1
    all_stats[current_row:end_row, ] <- stats_fcn(energy)
    current_row = current_row + num_rows
    setTxtProgressBar(pb, (current_row-1)/num_rows)
  }
  close(pb)
  
  # Assign a numeric ID column to hMOFs, character otherwise
  if (all(str_detect(files, "^h[^M]"))) {  # regular expressions(RegEx), prevent false detection of hMOF-number
    ids <- as.integer(str_sub(files, 2, -1))  # strip off leading "h" for hMOF designation
  } else {
    ids <- as.character(files)
  }
  ids <- rep(ids, each=num_rows)
  all_stats$id <- ids
  all_stats
}

# converting hists to metric stuff
metric_from_hists <- function(hist_df, lower_bound = -200, upper_bound = 0, warn = TRUE) {
  # Compute the "LJ metric" based on given cutoffs
  # Set a bound to NA for open intervals (e.g., energy > -200, but no upper bound)
  # TODO: NA code
  
  if (warn & !(lower_bound %in% hist_df$lower & upper_bound %in% hist_df$upper)) {
    warning("Metric bounds do not exactly line up with a histogram bin")
  }
  
  # be safer with float comparisons, and avoid double counting
  good_counts <- hist_df %>%
    filter(near(lower, lower_bound) | lower > lower_bound) %>%
    filter(near(upper, upper_bound) | upper < upper_bound) %>%
    group_by(id) %>% summarize(good = sum(counts))
  #print("Yes we got here") # for debugging purposes, if you wanna use, just uncomment this line
  total_counts <- hist_df %>% group_by(id) %>% summarize(total = sum(counts))

  lj_metric <- total_counts %>%
    inner_join(good_counts, by="id") %>% 
    mutate(metric = good / total) %>% 
    select(id, metric)
  lj_metric
}

# convert tibble into a matrix
# then calculate the correlation coefficient R
# Note: Value given by Post_resample is R^2 
# this function needs tibble with mof_id in it
# then the id column is converted to the row names and then get transposed
# finally it calls cor() function and returns a R-value-matrix
correlation_of_energy_histograms <- function(hists_with_id, write_to_csv = FALSE){
  a <- hists_with_id %>% select(-y_act)
  ab <- a %>% select(-MOF.ID)
  abc <- as.matrix(ab)
  namesss <- t(as.vector(a %>% select(MOF.ID)))
  rownames(abc) <- namesss
  correlation <- cor(t(abc))
  if (write_to_csv == TRUE){
    write.csv(correlation, file = "Results/correlation_between_histograms.csv")
  }
  correlation
}