# Makes diagnostic plots, such as model fit or distributions of bins

library(ggplot2)
library(dplyr)

default_binspec <- c(from=-20, to=1.0, step=1.0, width=1.0)

parity_line <- geom_abline(slope=1, intercept=0, linetype="dashed")

parity_plot <- function(act, pred, color=1, alpha=0.10) {
  # Parity plot between actual and predicted data, on square axes for g/L H2
  qplot(act, pred, alpha=I(alpha), color=I(color)) +
    xlab("'Actual' uptake (GCMC simulations)") +
    ylab("Predicted uptake (ridge regression)") +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(limits = c(0,50)) +
    scale_y_continuous(limits = c(0,50)) +
    parity_line
}


plot_bin_z_vs_y <- function(zs, y, betas) {
  # Plots the 2D distribution of z-score and y within each histogram bin, colored by model beta coef
  mod_betas <- betas %>% mutate(rbeta = round(beta, 2))
  
  bind_cols(z=zs, y=y) %>% 
    gather(key="bin", value="qty", -y) %>% 
    left_join(mod_betas, by="bin") %>%
    ggplot(aes(qty, y, col=rbeta)) +
    geom_point() +
    facet_wrap(~bin, scales="free_x") +
    xlab("z-score in bin") +
    theme(
      text = element_text(size=8),
      aspect.ratio = 0.75
    ) +
    scale_color_gradientn(colors = c("red", "darkgray", "blue"),
                          #guide=guide_colorbar(title=expression(beta)))
                          guide=guide_colorbar(title="beta"))
}
# Example: # plot_bin_z_vs_y(p_test_mod$x, p_test_mod$y, coef_tbl(p_test_mod$mod))

pred_grid <- function(glmnet_mod, test_grid, binspec) {
  # Runs predictions on grids by using binspec to calculate the energy histogram
  # and pred_glmnet(glmnet_mod, x) for z-scoring and evaluation.
  grid_desc <- test_grid %>%
    stepped_hist_spec(binspec) %>% 
    spread(key=bin, value=metric)
  
  grid_ids <- grid_desc$id
  grid_desc <- grid_desc %>% select(-id)
  
  grid_y_pred <- pred_glmnet(glmnet_mod, grid_desc)
  
  tibble(id = grid_ids, y_pred = grid_y_pred)
}

eval_test_grid <- function(glmnet_mod, test_grid, binspec, df_with_y_act) {
  # Runs pred_grid and calculates statistics/plots on model predictivity.
  # Requires a base grid and data.frame including a (renamed) y_act column.
  # TODO: what do do about color?
  # TODO: implement, then refactor run_model_on_partitions
  NULL
}

run_model_on_partitions <- function(partitioned_hists, y_with_id, binspec, alpha) {
  # Calculates the ridge or LASSO model on training data, then evaluates the performance on test data
  # Returns a large object with several plots
  
  results <- list(
    binspec = binspec, alpha = alpha,
    data = partitioned_hists, y = y_with_id,
    plots=list()
    )
  class(results) <- "partitioned_glmnet"
  
  if (alpha %in% c(1, 0)) {
    perf_models <- c("ridge regression", "LASSO")
    perf_model <- perf_models[alpha+1]  # R is one-indexed
  } else {
    perf_model <- "elasticnet"
  }
  results$model_name <- perf_model

  # Training a model, given our histogram parameters of 0.5 and 0.5
  trained_model <- run_bin_model(
    partitioned_hists$training, y_with_id,
    binspec["step"], binspec["width"],
    binspec[c("from", "to")],
    alpha=alpha
    )
  trained_mod <- trained_model$fitted_model[[1]]
  results$trained_model <- trained_model
  results$trained_mod <- trained_mod
  
  # Performance on the test data, which hasn't yet been used for anything
  predictions <- pred_grid(trained_mod, partitioned_hists$testing, binspec)
  testing_ids <- predictions$id
  y_act <- tibble(id = testing_ids) %>% 
    left_join(y_with_id, by="id") %>% 
    rename(y = g.L) %>% 
    .$y
  y_pred <- predictions$y_pred

  results$predictions <- list(
    y_act = y_act,
    y_pred = y_pred,
    test_ids = testing_ids
  )
  
  # What are the stats for R2, MAE, RMSE?
  results$training_fit <- postResample(pred=predict(trained_mod$mod, as.matrix(trained_mod$x)), obs=trained_mod$y)
  results$testing_fit <- postResample(pred=y_pred, obs=y_act)

  # Begin the plots
  results$plots$parity_bw <- parity_plot(y_act, y_pred) + ylab(paste0("Predicted uptake (", perf_model,")"))
  
  # Let's add two more plots, one with just the training data, and another with both training and test (color-coded)
  results$plots$parity_training <-
    parity_plot(
      trained_mod$y,
      pred_glmnet(trained_mod, trained_mod$orig_x),
      "#CA7C1B"
    ) +
    ylab(paste0("Predicted uptake (", perf_model,")"))
  # Add test data
  results$plots$parity_full <- results$plots$parity_training +
    geom_point(
      data=tibble(act=y_act, pred=y_pred),
      aes(act, pred),
      alpha=I(0.10),  # Make slightly darker
      color=I("#0070C0")
    )
  
  results  # return the partitioned_glmnet object
}

print.partitioned_glmnet <- function(x) {
  # Pretty-print results from run_model_on_partitions, in case we just want to view instead of saving
  #attach(x)
  cat(paste("Analysis for alpha =", x$alpha, x$model_name, "\n\n"))
  
  cat(paste("Q2 for the trained model is", x$trained_model$q2), fill=TRUE)
  cat("Fit of the training data\n")
  print(x$training_fit)
  
  cat("\nModel coefficients\n")
  coef(x$trained_mod$mod) %>% print
  
  cat("\nAccuracy of the test data\n")
  print(x$testing_fit)
  
  print(x$plots$parity_training)
  print(x$plots$parity_full)
  
  cat("\n")
  
  #detach("x")
}

# Other plots in old notebooks:
# Q2 vs. bin settings (short code and doesn't need to be generalized)
