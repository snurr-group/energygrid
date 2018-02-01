# Generates the hyper_tuned_hist analysis to test binning parameters
# Saves as an Rds instead of running the hypertuning process upon every notebook compile.
# Must be called from the root of the Git project (parent directory of this script).

source("Notebooks/setup_data.R")

hyper_tuned_hist <-
  seq(0.25, 2.5, 0.25) %>%
  map_dfr(
    function(x) run_bin_model(
      e_data = hmof_hist_sets$training,  # we no longer need hyperparameter data for other purposes
      y_with_id = hmof_y_to_join,
      step = x, width = x,
      bin_lims = c(default_binspec["from"], default_binspec["to"]),
      lambda = NULL, alpha = DEFAULT_ALPHA,
      align_bins = "downward"
    )
  ) %>% 
  # Warning: the above output cannot be viewed without deleting the fitted_model and id_list list columns
  (function(x) {
    y <- tibble()
    for (rownum in 1:nrow(x)) {
      curr_row <- x[rownum,]
      y <- coef_tbl(curr_row[[1,"fitted_model"]]$mod) %>% 
        mutate(q2 = curr_row$q2, binwidth = curr_row$width) %>% 
        inner_join(
          bin_loc_from_spec(
            c(
              from=curr_row$bin_lo,
              to=curr_row$bin_hi,
              step=curr_row$step,
              width=curr_row$width
              ),
            align_bins = "downward"
            ),
          by="bin"
        ) %>% 
        bind_rows(y, .)
    }
    y
  })(.)  # run the anonymous function

write_rds(hyper_tuned_hist, "BigData/Robj/hyper_tuned.Rds")
