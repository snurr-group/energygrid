# Generates the hyper_tuned_hist analysis to test binning parameters
# Saves as an Rds instead of running the hypertuning process upon every notebook compile.
# Must be called from the root of the Git project (parent directory of this script).

source("Notebooks/setup_data.R")

# Ideally, we'd just throw seq(0.5, 2.5, 0.25) into as.list() then plain map to get our step/width, but
# unfortunately the way that bin_below/"downward" currently works is to leave a gap between the original
# "from" variable and the bottom bin, since the bin_below is from c(-Inf, from)
# To fix this, manually calculate what the bottom bins should be, then feed it in using map2

widths <- seq(0.25, 2.5, 0.25) %>% as.list()
hyper_tuned_alpha <- DEFAULT_ALPHA

bottoms <- map(
  widths,
  function(x) {
    bounds_from_params(
      step = x, width = x,
      from = default_binspec["from"], default_binspec["to"],
      align_bins = "downward"
      )$lower[1]  # we only care about the limit of the bottommost bin
    }
  )

get_training_fit_stats <- function(package_mod) {
  postResample(pred=predict(package_mod$mod, as.matrix(package_mod$x)), obs=package_mod$y)
}

hyper_tuned_hist <-
  map2(widths, bottoms,
    function(x, y) run_bin_model(
      e_data = hmof_hist_sets$training,  # we no longer need hyperparameter data for other purposes
      y_with_id = hmof_y_to_join,
      step = x, width = x,
      bin_lims = c(y, default_binspec["to"]),
      lambda = NULL, alpha = hyper_tuned_alpha,
      align_bins = "strict"  # used to be "downward", but now that we're specifying the bottom bin, strict will be a more rigorous check
    )
  ) %>% 
  # Warning: the above output cannot be viewed without deleting the fitted_model and id_list list columns
  map_dfr(
    .,
    function(x) {
      coef_tbl(x$fitted_model[[1]]$mod) %>% 
        mutate(q2 = x$q2, binwidth = x$width, MAE=get_training_fit_stats(x$fitted_model[[1]])["MAE"]) %>% 
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
