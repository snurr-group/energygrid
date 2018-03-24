# Generates the hyper_tuned_hist analysis to test binning parameters
# Saves as an Rds instead of running the hypertuning process upon every notebook compile.
# Must be called from the root of the Git project (parent directory of this script).

source("Notebooks/setup_data.R")

hyper_tuned_hist <-
  seq(0.5, 2.5, 0.25) %>%
  as.list() %>% 
  map(.,
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
  map_dfr(
    .,
    function(x) {
      coef_tbl(x$fitted_model[[1]]$mod) %>% 
        mutate(q2 = x$q2, binwidth = x$width) %>% 
        inner_join(
          bin_loc_from_spec(
            c(
              from=x$bin_lo,
              to=x$bin_hi,
              step=x$step,
              width=x$width
              ),
            align_bins = "downward"
            ),
          by="bin"
          )
    })

write_rds(hyper_tuned_hist, "BigData/Robj/hyper_tuned.Rds")
