# Functions to calculate a (stepped) energy histogram and run ridge regression models

library(dplyr)
library(purrr)
library(glmnet)
library(caret)

# TODO: long-term, consider writing a data frame with beta, bin number, bin extremes, etc., for easy transfer of models

### FEATURE GENERATION ###

DEFAULT_ALPHA <- 1  # 0 for ridge regression, 1 for LASSO, otherwise elastic net.

bounds_from_params <- function(step, width, from, to, align_bins="strict") {
  # Consolidates calculations of bins from stepped_hist() and bin_loc_from_spec()
  # align_bins designates how to handle cases when from, to-width, and step do not align
  # "strict" means they must align exactly, going either way giving the same results.
  # This will still work in cases like bounds_from_params(2, 6, -40, 0).
  lower_bounds <- seq(from, to - width, step)
  upper_bounds <- lower_bounds + width
  
  upper_downward <- seq(to, from + width, -1*step) %>% rev
  lower_downward <- upper_downward - width
  
  if (align_bins == "strict") {
    stopifnot(lower_bounds == lower_downward)
    stopifnot(upper_bounds == upper_downward)
  } else if (align_bins == "downward") {  # direction of "to" to "from"
    lower_bounds <- lower_downward
    upper_bounds <- upper_downward
  } else if (align_bins != "upward") {  # if it's not a recognized keyword
    stop(paste("Unrecognized align_bins flag:", align_bins))
  }
  
  list(lower=lower_bounds, upper=upper_bounds)
}

stepped_hist <- function(df, step, width, from, to, bin_above=TRUE, align_bins="strict") {
  # Calculates a stepped histogram, where bins contain the overlap from the parent distribution
  # Usage: stepped_hist(only_one, 0.5, 0.5, -8, 0) %>% View
  # bin_above designates a bin to infinity above the "to" argument.
  
  lower_bounds <- bounds_from_params(step, width, from, to, align_bins)$lower
  upper_bounds <- bounds_from_params(step, width, from, to, align_bins)$upper

  #.id grabs the bin number, which should be sufficient as an identifier.
  # This approach will also likely be easier for "spread"-ing back into columns for plsr
  result_hist <- map2_dfr(lower_bounds, upper_bounds, metric_from_hists, hist_df=df, warn=FALSE, .id="bin")
  
  if (bin_above) {
    # If specified, add a bin above containing the remainder.
    # For whatever reason, the "bin" column from purrr is of type "character"
    result_hist <- result_hist %>% 
      bind_rows(mutate(metric_from_hists(df, to, Inf, warn=FALSE), bin="Inf"))
  }
  
  # < 2 sec per bin calculation  Just run metric_from_hists (an integral-type equation)
  result_hist
}

stepped_hist_spec <- function(df, binspec, ...) {
  # Runs stepped_hist given a binspec instead of manually specifying the arguments
  stepped_hist(df, binspec["step"], binspec["width"], binspec["from"], binspec["to"], ...)
}

bin_loc_from_spec <- function(binspec, bin_above=TRUE, align_bins="strict") {
  # Retrieves the locations of bins specified from stepped_hists
  both_bounds <- bounds_from_params(
    binspec["step"], binspec["width"],
    binspec["from"], binspec["to"],
    align_bins
    )
  lower_bounds <- both_bounds$lower
  upper_bounds <- both_bounds$upper
  
  bins <- tibble(
    lower = lower_bounds,
    upper = upper_bounds,
    bin = as.character(1:length(upper_bounds))
    )
  bins <- bins %>% mutate(loc = (lower + upper)/2.0)
  if (bin_above) {
    bins <- bins %>%
      add_row(
        lower = binspec["to"],
        upper = Inf,
        bin = "Inf",
        loc = binspec["to"] + 0.5*binspec["width"]
        )
  }
  bins
}

### DATA PROCESSING ###

partition_data_subsets <- function(unprocessed_x_with_id, y_with_id, data_split) {
  # Splits the raw histogram data into hyperparameter tuning, training, and test subsets.
  # That ensures the results are actually performed on independent sets of data (so the results
  # don't merely state that the model can fit itself)
  # As a side benefit, this speeds up the model processing since each section contains fewer rows.
  set.seed(20171017)
  
  # Get a list of unique ID's for splitting
  filtered_hist <- unprocessed_x_with_id %>% 
    filter(id %in% y_with_id$id)
  filtered_ids <- filtered_hist$id %>% unique
  num_ids <- length(filtered_ids)
  
  n_split <- integer(3)
  
  # Parse data_split into proportions for hyperparameter tuning, testing, and training
  if (length(data_split) == 1) {
    # A single number: the number of samples used for training.  Everything else is tested
    n_split <- c(0, data_split, num_ids - data_split)
  } else if (length(data_split) == 3) {
    # A vector of proportions for the three subsets
    n_split[1] <- num_ids*data_split[1]
    n_split[2] <- num_ids*data_split[2]
    n_split[3] <- NA  # not actually used
  } else {
    stop(paste("Improper data_split arg.  Must be a single number (training MOFs) or the fractional split.  Provided", data_split))
  }
  
  # For hyperparameter tuning
  hyperparam_rows <- sample(num_ids, n_split[1])
  hyperparam_ids <- filtered_ids[hyperparam_rows]
  hyperparam_data <- filtered_hist %>% filter(id %in% hyperparam_ids)
  # For training
  if (length(hyperparam_rows) > 0) {
    remaining_ids <- filtered_ids[-hyperparam_rows]
  } else {
    remaining_ids <- filtered_ids  # -numeric(0) returns numeric(0) or something similar
  }
  training_rows <- sample(length(remaining_ids), n_split[2])
  training_ids <- remaining_ids[training_rows]
  training_data <- filtered_hist %>% filter(id %in% training_ids)
  # For testing
  testing_ids <- remaining_ids[-training_rows]
  testing_data <- filtered_hist %>% filter(id %in% testing_ids)
  #y_to_join is already set in the calling function
  
  list(
    hyperparam = hyperparam_data,
    training = training_data,
    testing = testing_data
    )
}

standardize_from_props <- function(x, remove_cols, meanz, stdz) {
  # Apply standardization procedures
  x <- x[,!(names(x) %in% remove_cols)]
  # Thanks to this gem on the mailing list: https://stat.ethz.ch/pipermail/r-help/2006-May/104690.html
  x <- (x - meanz[rep(1,times=nrow(x)),]) / (stdz[rep(1,times=nrow(x)),])
  x
}

standardize <- function(glm_mod, x) {
  standardize_from_props(x, glm_mod$removed_cols, glm_mod$meanz, glm_mod$stdz)
}

standardization_from_x <- function(x, zscore = TRUE) {
  # Gets the standardization properties for an x matrix (removed_cols, meanz, and stdz)
  # By default, this transforms data with mean-center and unit variance
  
  # Remove columns with zero std dev
  removed_cols <- x %>%
    summarize_all(funs(sd)) %>% 
    gather("bin", "sd") %>% 
    filter(sd == 0) %>% # Could also replace with a tolerance
    .$bin
  x <- x[,!(names(x) %in% removed_cols)]
  
  # Mean-center
  if (zscore) {
    meanz <- x %>% summarize_all(funs(mean))
  } else {
    meanz <- x %>% summarize_all(function(a) {0})
  }
  x <- x - meanz[rep(1, times=nrow(x)),]
  
  # Unit variance
  if (zscore) {
    varz <- x %>% summarize_all(funs(sd))
  } else {
    varz <- x %>% summarize_all(function(a) {1})
  }
  x <- x / varz[rep(1, times=nrow(x)),]  # unused, but kept for completeness
  
  list(
    removed_cols = removed_cols,
    meanz = meanz,
    stdz = varz
    )
}

shuffle <- function(vec) {
  # Shuffles the values around in a vector
  sample(vec, size=length(vec))
}

### MODEL FITTING ###
# First, some potentially useful documentation links:
# * [caret's implementation](https://topepo.github.io/caret/available-models.html) of ridge regression is the [elasticnet package](https://cran.r-project.org/web/packages/elasticnet/elasticnet.pdf)
# * My current model uses [glmnet](https://cran.r-project.org/web/packages/glmnet/glmnet.pdf)
# * See also a helpful [vignette](https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet_beta.pdf) and [blog post](https://www.r-bloggers.com/ridge-regression-and-the-lasso/)

fit_glmnet <- function(x, y, lambda = NULL, alpha = DEFAULT_ALPHA, fit_intercept = TRUE, zscore = FALSE, ...) {
  # Fits a ridge regression model to a dataframe of predictors
  # By default, only removes zero variance columns without mean-centering or unit variance.
  # The alpha parameter specifies the type of model (ridge=0, LASSO=1, others=elastic net)
  
  orig_x <- x
  
  # Standardize x by removing zero variance columns.
  # Optionally with zscore flag do mean-centering and unit variance
  x_standardization <- standardization_from_x(x, zscore = zscore)  # By default, no longer zscore
  removed_cols <- x_standardization$removed_cols
  meanz <- x_standardization$meanz
  stdz <- x_standardization$stdz
  x <- standardize_from_props(x, removed_cols, meanz, stdz)
  
  # If lambda is not defined, let's calculate the largest (most regularized) value within one SE of the min. cross-validated error
  cvfit <- NULL
  if (is.null(lambda)) {
    trial_lambdas <- 10^seq(10, -6, length = 81)  # Idea from https://www.r-bloggers.com/ridge-regression-and-the-lasso/
    cvfit <- cv.glmnet(as.matrix(x), y, alpha=alpha, nfolds=10, type.measure="mse", lambda=trial_lambdas, intercept=fit_intercept, ...)
    lambda <- cvfit$lambda.min  # This is also where you could alternatively set lambda.1se, depending on which resource
  }
  
  # Actually run the model
  # alpha=0 is ridge regression,
  # We do not need to restandardize x (or y), which glmnet is likely doing.
  mod <- glmnet(as.matrix(x), y, alpha=alpha, lambda=lambda, intercept=fit_intercept, ...)
  
  # Return the relevant model details
  list(
    mod = mod,
    x = x,
    y = y,
    orig_x = orig_x,
    nfit = nrow(orig_x),
    removed_cols = removed_cols,
    meanz = meanz,
    stdz = stdz,
    lambda = lambda,
    alpha = alpha,
    coefs = NULL,
    cv_for_lambda = cvfit
  )
}

pred_glmnet <- function(glm_mod, x_tbl) {
  # Standardize data in x_tbl and apply the glmnet ridge regression model
  # glm_mod must be the form as returned by fit_glmnet, including the mod and x/y/stdz data
  #x <- x_tbl[,!(names(x_tbl) %in% glm_mod$removed_cols)]
  #x <- (x - glm_mod$meanz[rep(1,times=nrow(x)),]) / (glm_mod$stdz[rep(1,times=nrow(x)),])
  x <- standardize(glm_mod, x_tbl)
  predict(glm_mod$mod, as.matrix(x)) %>% as.numeric
}

run_bin_model <- function(e_data, y_with_id, step, width, bin_lims=c(-8, 0.5), lambda=NULL, alpha=DEFAULT_ALPHA, align_bins="strict", ...) {
  # Runs the ridge regression model, transforming e_data into a stepped histogram with appropriate y columns.
  # Also runs cross-validation and returns the model for later consideration
  x <- e_data %>% 
    stepped_hist(step, width, from=bin_lims[1], to=bin_lims[2], align_bins=align_bins) %>% 
    spread(key=bin, value=metric)
  y <- x %>% 
    left_join(y_with_id, by="id") %>%
    rename(y = g.L) %>% 
    .$y
  x_id <- x %>% select(id)
  x <- x %>% select(-id)
  
  q2 <- calc_q2(x, y, alpha=alpha)
  fitted_mod <- fit_glmnet(x, y, lambda=lambda, alpha=alpha, ...)
  tibble(
    fitted_model = list(fitted_mod),
    id_list = list(x_id),
    q2=q2,
    lambda = fitted_mod$lambda, 
    step = step,
    width = width,
    bin_lo = bin_lims[1],
    bin_hi = bin_lims[2]
  )
}

### MODEL ASSESSMENT ###

# First, some generally helpful plots:
# How the lambda parameter is calculated: plot(first_mod$cv_for_lambda)
# Basic parity plot (fixme): qplot(first_mod$y, predict(first_mod$mod, as.matrix(first_mod$x)), alpha=I(0.15))
# Beta's: coef(first_mod$mod)

# Q2, aka predicted R2
# Use caret's `createFolds` helper function to generate the CV fold indicies
calc_q2 <- function(x, y, lambda = NULL, alpha = DEFAULT_ALPHA, ...) {
  # Calculate Q2 for a given model fit
  # Alternatively, we could use the results from glmnet, but the code is messy and hard to read
  q2 <- NULL
  
  # First, obtain the global model parameters, such as lambda and mean-centering
  global_mod <- fit_glmnet(x, y, lambda, alpha, ...)
  global_lambda <- global_mod$lambda
  
  folds <- createFolds(y, k=10)
  
  press <- 0
  tss <- 0
  
  for (fold in folds) {
    # Can't use my fit and prediction functions, since we don't want to standardize inconsistently
    x_fold <- standardize(global_mod, x[fold,])
    y_fold <- y[fold]
    x_others <- standardize(global_mod, x[-fold,])
    y_others <- y[-fold]  # We can use the "-" notation in this case, since it's integers instead of names
    
    # Make glm model from others
    fold_model <- glmnet(as.matrix(x_others), y_others, alpha=alpha, lambda=global_lambda, ...)
    
    # Predict on current fold
    y_jack <- predict(fold_model, as.matrix(x_fold))
    #y_bar <- mean(y_fold)
    y_bar <- mean(y_others)
    
    # Add to press and tss accordingly
    press <- press + sum((y_fold - y_jack)^2)
    tss <- tss + sum((y_fold - y_bar)^2)
  }
  
  1 - press/tss  # Definition of Q2
}

coef_tbl <- function(mod_ridge) {
  # Extract beta coefficients from a trained (ridge/LASSO) regression model
  intermediate_coef <- coef(mod_ridge) %>% as.matrix
  tibble(bin = row.names(intermediate_coef), beta = as.numeric(intermediate_coef))
}

