# Makes diagnostic plots, such as model fit or distributions of bins

parity_line <- geom_abline(slope=1, intercept=0, linetype="dashed")

parity_plot <- function(act, pred, color=1, alpha=0.50) {
  # Parity plot between actual and predicted data, on square axes for g/L H2
  qplot(act, pred, alpha=I(alpha), color=I(color), size = I(dot_size)) +
    xlab("'Actual' capacity (GCMC)") +
    ylab("Predicted capacity (LASSO)") +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(limits = c(0,60)) +
    scale_y_continuous(limits = c(0,60)) +
    coord_fixed() +
    parity_line + 
    theme_classic()
}


pred_grid <- function(glmnet_mod, test_grid, binspec, align_bins = "strict") {
  # Runs predictions on grids by using binspec to calculate the energy histogram
  # and pred_glmnet(glmnet_mod, x) for z-scoring and evaluation.
  # TODO: consider incorporating binspec as part of glmnet_mod
  grid_desc <- test_grid %>%
    stepped_hist(binspec, align_bins = align_bins) %>% 
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
  print("Plotting and applying model to test data")
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
  
  do_label_low_loading <- FALSE
  low_loading_postresample_train <- NULL
  low_loading_postresample_test <- NULL
  low_loading_threshold <- 50
  LASSO_train_pred <- predict(results$trained_mod$mod, as.matrix(results$trained_mod$x))
  # result from predict function is a list, convert it to a vector
  LASSO_train_pred <- as.vector(LASSO_train_pred)
  # if the population within the threshold is large enough
    do_label_low_loading <-  TRUE
    low_loading_postresample_train <- low_loading_stat(pred = LASSO_train_pred, 
                                                       obs = results$trained_mod$y)
    low_loading_postresample_test <- low_loading_stat(pred = results$pred_df$y_pred, 
                                                       obs = results$pred_df$y_act)
    
  # calculate the boundaries of a parity plot: pick the lowest
  # put all together: training/testing, prediction/actual, get the lowest and highest
  pred_list <- c(LASSO_train_pred, results$trained_mod$y, results$pred_df$y_pred, results$pred_df$y_act)
  plot_limit <- c(min(pred_list), max(pred_list))
  #plot_limit[2] <- plot_limit[2]-mod(plot_limit[2],100)+100
  plot_limit[2] <- plot_limit[2]-mod(plot_limit[2],50)+50
  if (previous_plot_lim < plot_limit[2])
  {
    previous_plot_lim <<- plot_limit[2] # use global variable
  }else{
    plot_limit[2] <- previous_plot_lim
  }
  print(plot_limit[2])
  # get to see if we are doing XeKr selectivity
  if(XeKr){
    if (poster){
      results$plots$parity_training <- 
        results$plots$parity_training %>% rescale_ch4_parity(lims = plot_limit) %>% 
        annotate_plot(paste0("Training data\n", glmnet_mod$nfit," ", db_name), "top.left", "#CA7C1B") %>% 
        #label_stats(results$training_fit, label_q2=NULL, do_label_r2=TRUE)  # skip Q2 stats
        label_stats(results$training_fit, label_q2=NULL, plot_units = plot_units, do_label_r2=TRUE, label_position = "forXeKr")
    }else{
      results$plots$parity_training <- 
        results$plots$parity_training %>% rescale_ch4_parity(lims = plot_limit) %>% 
        annotate_plot(paste0("Training data\n", glmnet_mod$nfit," ", db_name), "top.left", "#CA7C1B") %>% 
        #label_stats(results$training_fit, label_q2=NULL, do_label_r2=TRUE)  # skip Q2 stats
        label_stats(results$training_fit, label_q2=NULL, plot_units = plot_units, do_label_r2=TRUE,
                    low_loading_label = do_label_low_loading, 
                    low_loading_stat = low_loading_postresample_train, label_position = "forXeKr")  # decided against labeling R2 on any of the figures
    }
  }else{
    if (poster){
      results$plots$parity_training <- 
        results$plots$parity_training %>% rescale_ch4_parity(lims = plot_limit) %>% 
        annotate_plot(paste0("Training data\n", glmnet_mod$nfit," ", db_name), "top.left", "#CA7C1B") %>% 
        #label_stats(results$training_fit, label_q2=NULL, do_label_r2=TRUE)  # skip Q2 stats
        label_stats(results$training_fit, label_q2=NULL, plot_units = plot_units, do_label_r2=TRUE)
    }else{
      results$plots$parity_training <- 
      results$plots$parity_training %>% rescale_ch4_parity(lims = plot_limit) %>% 
      annotate_plot(paste0("Training data\n", glmnet_mod$nfit," ", db_name), "top.left", "#CA7C1B") %>% 
      #label_stats(results$training_fit, label_q2=NULL, do_label_r2=TRUE)  # skip Q2 stats
      label_stats(results$training_fit, label_q2=NULL, plot_units = plot_units, do_label_r2=TRUE,
                  low_loading_label = do_label_low_loading, 
                  low_loading_stat = low_loading_postresample_train)  # decided against labeling R2 on any of the figures
    }
  }
  if(XeKr){
    if (poster){
      results$plots$parity_testing <- 
        parity_plot(y_act, y_pred, "#0070C0") %>% rescale_ch4_parity(lims = plot_limit) %>% 
        # We could use a thousands separator within the figures, but it looks weird and out of place
        #annotate_plot(paste0("Testing data\n", format(nrow(df_with_ys), , big.mark=",")," ", db_name), "top.left", "#0070C0") %>%
        annotate_plot(paste0("Testing data\n", nrow(df_with_ys)," ", db_name), "top.left", "#0070C0") %>% 
        label_stats(results$testing_fit, plot_units = plot_units, do_label_r2=TRUE, label_position = "forXeKr")
    }else{
      results$plots$parity_testing <- 
        parity_plot(y_act, y_pred, "#0070C0") %>% rescale_ch4_parity(lims = plot_limit) %>% 
        # We could use a thousands separator within the figures, but it looks weird and out of place
        #annotate_plot(paste0("Testing data\n", format(nrow(df_with_ys), , big.mark=",")," ", db_name), "top.left", "#0070C0") %>%
        annotate_plot(paste0("Testing data\n", nrow(df_with_ys)," ", db_name), "top.left", "#0070C0") %>% 
        label_stats(results$testing_fit, plot_units = plot_units, do_label_r2=TRUE,
                    low_loading_label = do_label_low_loading, 
                    low_loading_stat = low_loading_postresample_test, label_position = "forXeKr")
    }
  }else{
    # Zhao's idea: maybe labeling R2 is worthy...
    if (poster){
      results$plots$parity_testing <- 
        parity_plot(y_act, y_pred, "#0070C0") %>% rescale_ch4_parity(lims = plot_limit) %>% 
        # We could use a thousands separator within the figures, but it looks weird and out of place
        #annotate_plot(paste0("Testing data\n", format(nrow(df_with_ys), , big.mark=",")," ", db_name), "top.left", "#0070C0") %>%
        annotate_plot(paste0("Testing data\n", nrow(df_with_ys)," ", db_name), "top.left", "#0070C0") %>% 
        label_stats(results$testing_fit, plot_units = plot_units, do_label_r2=TRUE)
   }else{
    results$plots$parity_testing <- 
      parity_plot(y_act, y_pred, "#0070C0") %>% rescale_ch4_parity(lims = plot_limit) %>% 
      # We could use a thousands separator within the figures, but it looks weird and out of place
      #annotate_plot(paste0("Testing data\n", format(nrow(df_with_ys), , big.mark=",")," ", db_name), "top.left", "#0070C0") %>%
      annotate_plot(paste0("Testing data\n", nrow(df_with_ys)," ", db_name), "top.left", "#0070C0") %>% 
      label_stats(results$testing_fit, plot_units = plot_units, do_label_r2=TRUE,
                  low_loading_label = do_label_low_loading, 
                  low_loading_stat = low_loading_postresample_test)
   }
  }
  
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
  print("Partitioning data")
  trained_model <- run_bin_model(
    partitioned_hists$training, y_with_id,
    binspec,
    alpha=alpha,
    lambda=lambda,
    align_bins = align_bins,
    ...
    )
  trained_mod <- trained_model$fitted_model[[1]]
  print("training model complete")
  # Performance on the test data, which hasn't yet been used for anything
  predictions <- y_with_id %>% 
    mutate(y_act = g.L) %>%
    eval_test_grid(trained_mod, partitioned_hists$testing, binspec, ., db_name, plot_units, trained_model$q2, align_bins=align_bins)
  predictions$trained_model <- trained_model
  predictions  # return the partitioned_glmnet object
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
  if (poster){
    fontsize= 30
  }
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
  if (str_detect(pos, "forXeKr")) {
    xpos <- 1.0 - plot_margin 
    ypos <- plot_margin + 0.05
    vj <- 0
    hj <- 1
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

rescale_ch4_parity <- function(p, lims=c(0,250)) {
  if (lims[2] < 100)
  {
    by_value <- round((lims[2] - lims[1])/5) 
  }else
  {
    by_value <- 50
  }
  p +
    scale_x_continuous(limits = c(lims[1],lims[2]), breaks = seq(from = 0, to = lims[2], by = by_value)) +
    scale_y_continuous(limits = c(lims[1],lims[2]), breaks = seq(from = 0, to = lims[2], by = by_value))
}

rescale_ch4_parity_XeKr <- function(p, lims=c(0,250)) {
  if (lims[2] < 100)
  {
    by_value <- round((lims[2] - lims[1])/5) 
  }else
  {
    by_value <- 50
  }
  p +
    scale_x_continuous(limits = c(lims[1],lims[2]), breaks = seq(from = 0, to = lims[2], by = by_value)) +
    scale_y_continuous(limits = c(lims[1],lims[2]), breaks = seq(from = 0, to = lims[2], by = by_value))
}
# compute statistics for low_loading region
low_loading_stat <- function(pred, obs){
  # compute statistics based on GCMC loading
  threshold = quantile(obs)[[2]] # use the first quartile of the vector as threshold
  combined_df <- as.data.frame(row.names = NULL, pred) %>% cbind(obs)

  combined_df <- combined_df %>% filter(obs <= threshold)
  low_loading_statistics <- postResample(pred = combined_df$pred, obs = combined_df$obs)
  # add low loading threshold to the named numbers
  low_loading_statistics <- c(low_loading_statistics, low_loading_Threshold=threshold)
}

# Label stats on the training and testing plots
label_stats <- function(p, postresample_results, label_q2=NULL, do_label_r2=FALSE, plot_units = "g/L", low_loading_label=FALSE, low_loading_statistics = NULL, label_position = "bottom.right") {
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
  if (low_loading_label) {
    training_stats <- paste0(
      training_stats, "\n",
      "Low Loading R\u00B2 = ", format(low_loading_statistics["Rsquared"], digits=2, nsmall=1), "\n", 
      "Low Loading MAE = ", format(low_loading_statistics["MAE"], digits=2, nsmall=1), " ", plot_units, "\n", 
      "Low Loading RMSE = ", format(low_loading_statistics["RMSE"], digits=2, nsmall=1), " ", plot_units, "\n",
      "Note: Low Loading Statistics generated by First Quartile of GCMC Loading Data", "\n", 
      "Low Loading threshold = ", format(low_loading_statistics["low_loading_Threshold"], digits=2, nsmall=1), " ", plot_units
    )
  }
  training_stats <- paste0(
    training_stats, "\n",
    "MAE = ", format(postresample_results["MAE"], digits=2, nsmall=1), " ", plot_units, "\n",
    "RMSE = ", format(postresample_results["RMSE"], digits=2, nsmall=1), " ", plot_units
  )  # https://en.wikipedia.org/wiki/Unicode_subscripts_and_superscripts
  p %>% annotate_plot(training_stats, label_position)
}

# generate rf predictions and plots
make_rf_prediction_plots <- function(condition_name, plot_name, rf_model, test_data, lims = c(0,250), axis_label = paste0("capacity", "(", unit_for_plot, ")")){
  # condition names should go like: molecule_temperature_pressure
  new_string <- paste0(condition_name, "_")
  if(poster){
    new_string <- paste0(new_string, "_poster")
  }
  # make prediction
  tested <- predict(rf_model, test_data)
  train_rmse <- postResample(pred = rf_model$predicted, obs = rf_model$y)
  test_rmse <- postResample(pred = tested, obs = test_data$y_act)
  # consider doing low loading stats
  do_label_low_loading <- FALSE
  low_loading_postresample_train <- NULL
  low_loading_postresample_test <- NULL
  # if the population within the threshold is large enough
    do_label_low_loading <-  TRUE
    low_loading_postresample_train <- low_loading_stat(pred = rf_model$predicted, 
                                                       obs = rf_model$y)
    low_loading_postresample_test <- low_loading_stat(pred = tested, 
                                                      obs = test_data$y_act)
  # get the lowest and highest 
  pred_list <- c(rf_model$predicted, rf_model$y, tested, test_data$y_act)
  # for here, we also need to scale
  plot_limit <- c(min(pred_list), max(pred_list))
  
  plot_limit[2] <- plot_limit[2]-mod(plot_limit[2],100)+100
  if (plot_limit[2] < previous_plot_lim)
  {
    plot_limit[2] <- previous_plot_lim
  }
  print(plot_limit)
  
  gg_rf_train <- qplot(x = rf_model$y, y = rf_model$predicted) %>% rescale_ch4_parity(., lims=plot_limit) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    xlab(paste0("GCMC ", axis_label)) + 
    ylab(paste0("Predicted ", axis_label)) + 
    geom_point(color="#CA7C1B", size = dot_size) + 
    theme_classic()
  
  if(poster){
    gg_rf_train <- gg_rf_train %>% annotate_plot(paste0("Training data\n", length(rf_model$y)," ", "points"), "top.left", "#CA7C1B") %>% 
      label_stats(., train_rmse, plot_units = unit_for_plot, do_label_r2=TRUE)
    gg_rf_train <- gg_rf_train + theme(axis.text=element_text(size=30),
                                       axis.title=element_text(size=30,face="bold"))
  }else{
    gg_rf_train <- gg_rf_train %>% annotate_plot(paste0("Training data\n", length(rf_model$y)," ", "points"), "top.left", "#CA7C1B") %>% 
      label_stats(., train_rmse, plot_units = unit_for_plot, do_label_r2=TRUE,  
                  low_loading_label = do_label_low_loading, 
                  low_loading_stat = low_loading_postresample_train)
  }
  save_plot(paste(save_path, paste0(new_string, paste0(plot_name, "_train.png")), sep = ""), 
            gg_rf_train, base_width = 10, base_height = 10, dpi = 600)
  
  number_of_data_smaller_than_zero <- sum(rf_model$predicted < 0)
  if (number_of_data_smaller_than_zero > 0){
    file_name <- file(paste(save_path,"warning.txt"), 'a')
    writeLines(paste0(plot_name, " has ", number_of_data_smaller_than_zero, " smaller than zero"), file_name)
  }
  
  
  gg_rf_test <- qplot(x = test_data$y_act, y = tested) %>% rescale_ch4_parity(., lims=plot_limit) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    xlab(paste0("GCMC ", axis_label)) + 
    ylab(paste0("Predicted ", axis_label)) + 
    geom_point(color="#0070C0", size = dot_size) + 
    theme_classic()
 
  if(poster){
    gg_rf_test <- gg_rf_test %>% annotate_plot(paste0("Testing data\n", nrow(test_data)," ", "points"), "top.left", "#0070C0") %>% 
      label_stats(., test_rmse, plot_units = unit_for_plot, do_label_r2=TRUE)
    gg_rf_test <- gg_rf_test + theme(axis.text=element_text(size=30),
                                     axis.title=element_text(size=30,face="bold"))
  }else{
    gg_rf_test <- gg_rf_test %>% annotate_plot(paste0("Testing data\n", nrow(test_data)," ", "points"), "top.left", "#0070C0") %>% 
      label_stats(., test_rmse, plot_units = unit_for_plot, do_label_r2=TRUE,
                  low_loading_label = do_label_low_loading, 
                  low_loading_stat = low_loading_postresample_test)
  }
  save_plot(paste(save_path, paste0(new_string, paste0(plot_name, "_test.png")), sep = ""), 
            gg_rf_test, base_width = 10, base_height = 10, dpi = 600)
  ## add another section that checks if predicted values has negative elements
  # if so, write a warning text 
  # now, check on test data
  number_of_data_smaller_than_zero <- sum(tested < 0)
  if (number_of_data_smaller_than_zero > 0){
    file_name <- file(paste(save_path,"warning.txt"), 'a')
    writeLines(paste0(plot_name, " has ", number_of_data_smaller_than_zero, " smaller than zero"), file_name)
  }
}

# get RandomForest Variable Importance
getVarImp <- function(RFmodel, modelshort = "RF", howmany = 10, condition, modelname, rename_bins = TRUE){
  rf_varimp <- importance(RFmodel)
  rf_varimp <- as.data.frame(rf_varimp)
  rf_varimp$varnames <- rownames(rf_varimp)
  rownames(rf_varimp) <- NULL
  
  # use mid-points for energy bins
  rf_varimp$midpoints <- rf_varimp$varnames
  # rename the energy bins to ranges of energies?
  if (rename_bins){
    for (i in 1:length(rf_varimp$varnames)){
      if (rf_varimp$varnames[i] %>% grepl("^[0-9]+$", ., perl = T)){ # if it is number-only
        upper <- round(ch4_binspec$bounds$upper[strtoi(rf_varimp$varnames[i])], 1)
        lower <- round(ch4_binspec$bounds$lower[strtoi(rf_varimp$varnames[i])], 1)
        if (lower == -50){
          lower = -Inf
        }
        varstr <- paste0(toString(lower), " to\n", toString(upper))
        rf_varimp$varnames[i] <- varstr
        rf_varimp$midpoints[i] <- (lower + upper)/2
      } else if (rf_varimp$varnames[i] == "Inf") {
        rf_varimp$varnames[i] <- "Inf"
        rf_varimp$midpoints[i] <- "2"
      }
    }
  }
  rf_varimp %>% write_csv(., paste(save_path, paste0(modelshort, "_VarImp.csv"), sep = ""))
  #rf_varimp %>% write_csv(paste(save_path, paste0(modelshort, "_VarImp.png"), sep = ""))
  rf_varimp_purity <- rf_varimp %>% top_n(., IncNodePurity, n = howmany)
  plt_rfimp <- ggplot(rf_varimp_purity, aes(x = reorder(varnames, IncNodePurity), 
                                     weight = IncNodePurity))+
    geom_bar() +
    scale_fill_discrete(name="Feature") +
    ylab("Feature Importance") +
    xlab("Feature Name")
  texttopaste <- paste0(condition, "\n", modelname)
  plt_rfimp <- plt_rfimp %>% annotate_plot(., texttopaste, "top.left", 
                                           col="black", fontsize = 10) + 
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25,face="bold"))
  # convert variable name to string
  # deparse(substitute(RFmodel))
  
  save_plot(paste(save_path, paste0(modelshort, "_VarImp.png"), sep = ""), 
            plt_rfimp, base_width = 15, base_height = 10, dpi = 600)
  # make a similar figure, but with %IncMSE
  rf_varimp_IncMSE <- rf_varimp %>% top_n(., `%IncMSE`, n = howmany)
  plt_rfimp <- ggplot(rf_varimp_IncMSE, aes(x = reorder(varnames, `%IncMSE`), 
                                     weight = `%IncMSE`))+
    geom_bar() +
    scale_fill_discrete(name="Feature") +
    ylab("Feature %IncMSE") +
    xlab("Feature Name")
  texttopaste <- paste0(condition, "\n", modelname)
  
  plt_rfimp <- plt_rfimp %>% annotate_plot(., texttopaste, "top.left", 
                                           col="black", fontsize = 10) + 
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25,face="bold"))
  # convert variable name to string
  # deparse(substitute(RFmodel))
  
  save_plot(paste(save_path, paste0(modelshort, "_Var_PIncMSE.png"), sep = ""), 
            plt_rfimp, base_width = 15, base_height = 10, dpi = 600)
}

make_rf_prediction_plots_XeKr <- function(condition_name, plot_name, rf_model, test_data, lims = c(0,250), axis_label = paste0("capacity", "(", unit_for_plot, ")"), Errorlabel = "", interval = 50){
  # condition names should go like: molecule_temperature_pressure
  new_string <- paste0(condition_name, "_")
  if(poster)
  {
    new_string <- paste0(new_string, "_poster")
  }
  # make prediction
  tested <- predict(rf_model, test_data)
  train_rmse <- postResample(pred = rf_model$predicted, obs = rf_model$y)
  test_rmse <- postResample(pred = tested, obs = test_data$y_act)
  # consider doing low loading stats
  do_label_low_loading <- FALSE
  low_loading_postresample_train <- NULL
  low_loading_postresample_test <- NULL
  # if the population within the threshold is large enough
  do_label_low_loading <-  TRUE
  low_loading_postresample_train <- low_loading_stat(pred = rf_model$predicted, 
                                                     obs = rf_model$y)
  low_loading_postresample_test <- low_loading_stat(pred = tested, 
                                                    obs = test_data$y_act)
  # get the lowest and highest 
  pred_list <- c(rf_model$predicted, rf_model$y, tested, test_data$y_act)
  # for here, we also need to scale
  plot_limit <- c(min(pred_list), max(pred_list))
  
  plot_limit[2] <- plot_limit[2]-mod(plot_limit[2],interval)+interval
  if (plot_limit[2] < previous_plot_lim)
  {
    plot_limit[2] <- previous_plot_lim
  }
  print(plot_limit)
  
  gg_rf_train <- qplot(x = rf_model$y, y = rf_model$predicted) %>% rescale_ch4_parity_XeKr(., lims=plot_limit) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    xlab(paste0("GCMC ", axis_label)) + 
    ylab(paste0("Predicted ", axis_label)) + 
    geom_point(color="#CA7C1B", size = dot_size) + 
    theme_classic()
  
  if(poster){
    gg_rf_train <- gg_rf_train %>% annotate_plot(paste0("Training data\n", length(rf_model$y)," ", "points"), "top.left", "#CA7C1B") %>% 
      label_stats(., train_rmse, plot_units = unit_for_plot, do_label_r2=TRUE, label_position = "forXeKr")
    gg_rf_train <- gg_rf_train + theme(axis.text=element_text(size=30),
                                       axis.title=element_text(size=30,face="bold"))
  }else{
    gg_rf_train <- gg_rf_train %>% annotate_plot(paste0("Training data\n", length(rf_model$y)," ", "points"), "top.left", "#CA7C1B") %>% 
      label_stats(., train_rmse, plot_units = unit_for_plot, do_label_r2=TRUE,  
                  low_loading_label = do_label_low_loading, 
                  low_loading_stat = low_loading_postresample_train, label_position = "forXeKr")
  }
  save_plot(paste(save_path, paste0(new_string, paste0(plot_name, "_train.png")), sep = ""), 
            gg_rf_train, base_width = 10, base_height = 10, dpi = 600)
  
  number_of_data_smaller_than_zero <- sum(rf_model$predicted < 0)
  if (number_of_data_smaller_than_zero > 0){
    file_name <- file(paste(save_path,"warning.txt"), 'a')
    writeLines(paste0(plot_name, " has ", number_of_data_smaller_than_zero, " smaller than zero"), file_name)
  }
  
  
  gg_rf_test <- qplot(x = test_data$y_act, y = tested) %>% rescale_ch4_parity_XeKr(., lims=plot_limit) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + 
    xlab(paste0("GCMC ", axis_label)) + 
    ylab(paste0("Predicted ", axis_label)) + 
    geom_point(color="#0070C0", size = dot_size) + 
    theme_classic()
  
  if(poster){
    gg_rf_test <- gg_rf_test %>% annotate_plot(paste0("Testing data\n", nrow(test_data)," ", "points"), "top.left", "#0070C0") %>% 
      label_stats(., test_rmse, plot_units = unit_for_plot, do_label_r2=TRUE, label_position = "forXeKr")
    gg_rf_test <- gg_rf_test + theme(axis.text=element_text(size=30),
                                     axis.title=element_text(size=30,face="bold"))
  }else{
    gg_rf_test <- gg_rf_test %>% annotate_plot(paste0("Testing data\n", nrow(test_data)," ", "points"), "top.left", "#0070C0") %>% 
      label_stats(., test_rmse, plot_units = unit_for_plot, do_label_r2=TRUE,
                  low_loading_label = do_label_low_loading, 
                  low_loading_stat = low_loading_postresample_test, label_position = "forXeKr")
  }
  save_plot(paste(save_path, paste0(new_string, paste0(plot_name, "_test.png")), sep = ""), 
            gg_rf_test, base_width = 10, base_height = 10, dpi = 600)
  ## add another section that checks if predicted values has negative elements
  # if so, write a warning text 
  # now, check on test data
  number_of_data_smaller_than_zero <- sum(tested < 0)
  if (number_of_data_smaller_than_zero > 0){
    file_name <- file(paste(save_path,"warning.txt"), 'a')
    writeLines(paste0(plot_name, " has ", number_of_data_smaller_than_zero, " smaller than zero"), file_name)
  }
}

make_topology_histograms <- function(condition_name, topo_data){
  # set a folder for these histograms
  new_string <- paste0(condition_name, "_")
  save_path_2 <- paste0(save_path, "topology_histograms/")
  if (!dir.exists(save_path_2)){
    dir.create(save_path_2)
  }
  # vf
  #gg_histo_vf <- qplot(topo_data$VF, geom = "histogram", xlim = c(0,1)) + 
   # xlab("Void Fraction") + ylab("Counts") + theme_classic() + 
   # theme(axis.text=element_text(size=30),
   #       axis.title=element_text(size=30,face="bold"))
  # vsa
  gg_histo_vsa <- qplot(topo_data$VSA, geom = "histogram", xlim = c(0,3000)) + 
    xlab("Volumetric Surface Area (m\u00B2/cm\u00B3)") + ylab("Counts") + theme_classic() + 
    theme(axis.text=element_text(size=30),
          axis.title=element_text(size=30,face="bold"))
  # gsa
  gg_histo_gsa <- qplot(topo_data$GSA, geom = "histogram", xlim = c(0,10000)) + 
    xlab("Gravitational Surface Area (m\u00B2/g)") + ylab("Counts") + theme_classic() + 
    theme(axis.text=element_text(size=30),
          axis.title=element_text(size=30,face="bold"))
  # pld
  gg_histo_pld <- qplot(topo_data$PLD, geom = "histogram", xlim = c(0,70)) + 
    xlab("Pore Limiting Diameter (\305)") + ylab("Counts") + theme_classic() + 
    theme(axis.text=element_text(size=30),
          axis.title=element_text(size=30,face="bold"))
  # lcd
  gg_histo_lcd <- qplot(topo_data$LCD, geom = "histogram", xlim = c(0,80)) + 
    xlab("Largest Cavity Diameter (\305)") + ylab("Counts") + theme_classic() + 
    theme(axis.text=element_text(size=30),
          axis.title=element_text(size=30,face="bold"))
  if(poster){
    
  }
  # save the plots
  save_plot(paste(save_path_2, paste0(new_string, "_vf.png"), sep = ""), 
            gg_histo_vf, base_width = 10, base_height = 10, dpi = 600)
  save_plot(paste(save_path_2, paste0(new_string, "_vsa.png"), sep = ""), 
            gg_histo_vsa, base_width = 10, base_height = 10, dpi = 600)
  save_plot(paste(save_path_2, paste0(new_string, "_gsa.png"), sep = ""), 
            gg_histo_gsa, base_width = 10, base_height = 10, dpi = 600)
  save_plot(paste(save_path_2, paste0(new_string, "_pld.png"), sep = ""), 
            gg_histo_pld, base_width = 10, base_height = 10, dpi = 600)
  save_plot(paste(save_path_2, paste0(new_string, "_lcd.png"), sep = ""), 
            gg_histo_lcd, base_width = 10, base_height = 10, dpi = 600)
}