# Makes diagnostic plots, such as model fit or distributions of bins

library(ggplot2)
library(dplyr)
library(grid)  # textGrob, since geom_text/annotate won't let you adjust an x=Inf, etc.
library(stringr)

parity_line <- geom_abline(slope=1, intercept=0, linetype="dashed")

parity_plot <- function(act, pred, color=1, alpha=0.50) {
  # Parity plot between actual and predicted data, on square axes for g/L H2
  qplot(act, pred, alpha=I(alpha), color=I(color)) +
    xlab("'Actual' capacity (GCMC)") +
    ylab("Predicted capacity (LASSO)") +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(limits = c(0,60)) +
    scale_y_continuous(limits = c(0,60)) +
    coord_fixed() +
    parity_line
}


pred_grid <- function(glmnet_mod, test_grid, binspec, align_bins = "strict") {
  # Runs predictions on grids by using binspec to calculate the energy histogram
  # and pred_glmnet(glmnet_mod, x) for z-scoring and evaluation.
  # TODO: consider incorporating binspec as part of glmnet_mod
  grid_desc <- test_grid %>%
    stepped_hist_spec(binspec, align_bins = align_bins) %>% 
    spread(key=bin, value=metric)
  
  grid_ids <- grid_desc$id
  grid_desc <- grid_desc %>% select(-id)
  
  grid_y_pred <- pred_glmnet(glmnet_mod, grid_desc)
  
  tibble(id = grid_ids, y_pred = grid_y_pred)
}

eval_test_grid <- function(glmnet_mod, test_grid, binspec, df_with_y_act, db_name = "MOFs", plot_units = "g/L", q2 = NULL, align_bins = "strict") {
  # Runs pred_grid and calculates statistics/plots on model predictivity.
  # Requires a base grid and data.frame including a (renamed) y_act column.
  # TODO: what do do about color?  Custom color for the dataframe?  Maybe just a manual parity plot.
  
  alpha_glm <- glmnet_mod$alpha
  trained_mod <- glmnet_mod
  
  results <- list(
    binspec = binspec,
    align_bins = align_bins,
    alpha = alpha_glm,
    trained_mod = glmnet_mod,
    plots=list()
  )
  class(results) <- "partitioned_glmnet"
  
  if (alpha_glm %in% c(1, 0)) {
    perf_models <- c("ridge", "LASSO")
    perf_model <- perf_models[alpha_glm+1]  # R is one-indexed
  } else {
    perf_model <- "elasticnet"
  }
  perf_label <- paste0(perf_model, ", ", plot_units)  # Label axes with corresponding units
  results$model_name <- perf_model

  predictions <- pred_grid(glmnet_mod, test_grid, binspec, align_bins = align_bins)
  df_with_ys <- predictions %>% 
    inner_join(df_with_y_act, by="id") %>% 
    mutate(y_err = abs(y_act - y_pred))
  results$pred_df <- df_with_ys
  
  # Convenience variables for plotting, etc.
  y_act <- df_with_ys$y_act
  y_pred <- df_with_ys$y_pred
  
  # What are the stats for R2, MAE, RMSE?
  results$training_fit <- postResample(pred=predict(trained_mod$mod, as.matrix(trained_mod$x)), obs=trained_mod$y)
  results$testing_fit <- postResample(pred=y_pred, obs=y_act)
  
  # Begin the plots
  results$plots$parity_bw <- parity_plot(y_act, y_pred) +
    xlab(paste0("'Actual' capacity (GCMC, ", plot_units, ")")) +
    ylab(paste0("Predicted capacity (", perf_label,")"))
  
  # Let's add two more plots, one with just the training data, and another with both training and test (color-coded)
  results$plots$parity_training <-
    parity_plot(
      trained_mod$y,
      pred_glmnet(trained_mod, trained_mod$orig_x),
      "#CA7C1B"
    ) +
    xlab(paste0("'Actual' capacity (GCMC, ", plot_units, ")")) +
    ylab(paste0("Fitted capacity (", perf_label,")"))
  # Add test data
  results$plots$parity_full <- results$plots$parity_training +
    geom_point(
      data=tibble(act=y_act, pred=y_pred),
      aes(act, pred),
      alpha=I(0.10),  # Make slightly darker
      color=I("#0070C0")
    )
  
  # Label stats on the training and testing plots
  label_stats <- function(p, postresample_results, label_q2=NULL, do_label_r2=FALSE) {
    training_stats <- ""
    if (!is.null(label_q2)) {
      training_stats <- paste0(
        training_stats, "\n",
        "Q\u00B2 = ", format(label_q2, digits=2, nsmall=2)
        )
    }
    if (do_label_r2) {
      training_stats <- paste0(
        training_stats, "\n",
        "R\u00B2 = ", format(postresample_results["Rsquared"], digits=2, nsmall=2)
        )
    }
    training_stats <- paste0(
      training_stats, "\n",
      "MAE = ", format(postresample_results["MAE"], digits=2, nsmall=1), " ", plot_units, "\n",
      "RMSE = ", format(postresample_results["RMSE"], digits=2, nsmall=1), " ", plot_units
       )  # https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts
    p %>% annotate_plot(training_stats, "bottom.right")
  }

  results$plots$parity_training <- 
    results$plots$parity_training %>% 
    annotate_plot(paste0("Training data\n", glmnet_mod$nfit," ", db_name), "top.left", "#CA7C1B") %>% 
    #label_stats(results$training_fit, label_q2=NULL, do_label_r2=TRUE)  # skip Q2 stats
    label_stats(results$training_fit, label_q2=NULL)  # decided against labeling R2 on any of the figures
  
  results$plots$parity_testing <- 
    parity_plot(y_act, y_pred, "#0070C0") %>% 
    # We could use a thousands separator within the figures, but it looks weird and out of place
    #annotate_plot(paste0("Testing data\n", format(nrow(df_with_ys), , big.mark=",")," ", db_name), "top.left", "#0070C0") %>%
    annotate_plot(paste0("Testing data\n", nrow(df_with_ys)," ", db_name), "top.left", "#0070C0") %>% 
    label_stats(results$testing_fit) +
    xlab(paste0("'Actual' capacity (GCMC, ", plot_units, ")")) +
    ylab(paste0("Predicted capacity (", perf_label,")"))
  
  # Check the normality of the residuals.  Though recall that ridge/LASSO are biased estimators
  results$plots$resid_normality <- df_with_ys %>% 
    ggplot(aes(y_act - y_pred)) +
    geom_histogram(bins = 30) +
    xlab(expression('Residuals'~y-hat(y)))
  
  # How well do the models rank the test data
  results$test_spearman <- results$pred_df %>% 
    select(y_act, y_pred) %>% 
    cor(method = "spearman")
  results$test_kendall <- results$pred_df %>% 
    select(y_act, y_pred) %>% 
    cor(method = "kendall")
  results$plots$test_ranking <- results$pred_df %>% 
    # rank() by default is in ascending order, but for some reason dplyr will accept desc
    # https://stackoverflow.com/questions/26106408/create-a-ranking-variable-with-dplyr
    mutate(r_act = rank(desc(y_act)), r_pred = rank(desc(y_pred))) %>% 
    ggplot(aes(r_act, r_pred)) +
    geom_point() +
    xlab("'Actual' ranking (GCMC): 1 is the best") +
    ylab(paste0("Predicted ranking (", perf_model,")")) +
    scale_x_reverse() + scale_y_reverse()  # Consistent interpretation as top right for best MOFs
  results$plots$test_ranking <- 
    results$plots$test_ranking %>% 
    annotate_plot(
      as_plotmath(paste0(
        "atop(",  # plotmath can't handle \n, so use a fraction per their recommendation
        "r[s] == ", format(results$test_spearman["y_pred", "y_act"], digits=2, nsmall=2),
        ",tau == ", format(results$test_kendall["y_pred", "y_act"], digits=2, nsmall=2),
        ")"
        )),
      "bottom.right"
    )
  
  results  # return the partitioned_glmnet object
}

run_model_on_partitions <- function(partitioned_hists, y_with_id, binspec, alpha=DEFAULT_ALPHA, lambda=NULL, db_name = "MOFs", plot_units="g/L", align_bins="strict", ...) {
  # Calculates the ridge or LASSO model on training data, then evaluates the performance on test data
  # Returns a large object with several plots
  
  # Training a model, given our histogram parameters
  trained_model <- run_bin_model(
    partitioned_hists$training, y_with_id,
    binspec["step"], binspec["width"],
    binspec[c("from", "to")],
    alpha=alpha,
    lambda=lambda,
    align_bins = align_bins,
    ...
    )
  trained_mod <- trained_model$fitted_model[[1]]
  
  # Performance on the test data, which hasn't yet been used for anything
  predictions <- y_with_id %>% 
    mutate(y_act = g.L) %>%
    eval_test_grid(trained_mod, partitioned_hists$testing, binspec, ., db_name, plot_units, trained_model$q2, align_bins=align_bins)
  predictions$trained_model <- trained_model
  predictions  # return the partitioned_glmnet object
}

print.partitioned_glmnet <- function(x) {
  # Pretty-print results from run_model_on_partitions, in case we just want to view instead of saving
  #attach(x)
  cat(paste("Analysis for alpha =", x$alpha, x$model_name, "\n\n"))
  
  cat(paste("Q2 for the trained model is", x$trained_model$q2), fill=TRUE)
  cat(paste("Lambda for the trained model is", x$trained_mod$lambda), fill=TRUE)
  cat("Fit of the training data\n")
  print(x$training_fit)
  
  cat("\nModel coefficients\n")
  coef(x$trained_mod$mod) %>% print
  
  cat("\nAccuracy of the test data\n")
  print(x$testing_fit)
  
  print(x$plots$parity_training)
  print(x$plots$parity_testing)
  print(x$plots$parity_full)
  print(x$plots$resid_normality)
  
  cat(paste("Spearman correlation for ranking of the test data is", x$test_spearman["y_pred", "y_act"]), fill=TRUE)
  cat(paste("Kendall tau is", x$test_kendall["y_pred", "y_act"]), fill=TRUE)
  cat("\n(On the ranking plot, a higher number indicates higher capacity)\n")
  print(x$plots$test_ranking)
  
  cat("\n")
  
  #detach("x")
}

annotate_plot <- function(p, label, pos="top.left", col="black", fontsize = 10, ...) {
  # Annotates a plot, including margins
  # A great resource is this cheatsheet: http://zevross.com/blog/2014/08/04/beautiful-plotting-in-r-a-ggplot2-cheatsheet-3/
  
  # Set default justification and positioning to center
  xpos <- 0.5
  ypos <- 0.5
  hj <- 0.5
  vj <- 0.5
  # Margins to end of document
  plot_margin <- 0.05  # fraction of the plot area (0.02 is nearly flush)
  
  if (str_detect(pos, "left")) {
    xpos <- plot_margin
    hj <- 0
  }
  if (str_detect(pos, "right")) {
    xpos <- 1.0 - plot_margin
    hj <- 1
  }
  if (str_detect(pos, "top")) {
    ypos <- 1.0 - plot_margin
    vj <- 1
  }
  if (str_detect(pos, "bot")) {
    ypos <- plot_margin
    vj <- 0
  }
  
  p + annotation_custom(textGrob(
    label,
    x=xpos, y=ypos, hjust=hj, vjust=vj,
    gp=gpar(fontsize=fontsize, col=col, ...)
    ))
}

as_plotmath <- function(x) {
  # Converts a string x to the necessary `expression` format for plotmath annotation, etc.
  # Can use in `annotate_plot(as_plotmath("Q^2 == 0.95))``
  scales::parse_format()(x)[[1]]
  # Plotmath is such a pain.  At least scales helps take care of it: https://jangorecki.gitlab.io/data.table/library/scales/html/parse_format.html
}

rescale_ch4_parity <- function(p, lim=250) {
  p +
    scale_x_continuous(limits = c(0,lim)) +
    scale_y_continuous(limits = c(0,lim))
}
