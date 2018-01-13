# Figures for Scotty 1/12/18

# First load everything by running Notebooks/figures.Rmd
# Might consider moving this analysis to the main Figure 1, too.

ad_hoc_grids <- read_rds("BigData/Robj/ad_hoc.Rds")

list(
  `A` = 1898,
  `B` = 5072244
  ) %>% 
  map2(
    ., names(.),
    ~ filter(ad_hoc_grids, id == paste0("hMOF-", .x)) %>% 
    plot_hist_bins(default_binspec, "Fraction of unit cell") %>% 
    overlay_cat_betas(coef_tbl(trained_mod$mod), default_binspec) %>% 
    save_ben_fig(paste0("scotty_", .y, ".png"))
    )

annotate_scotty_parity <- function(p, top_left_label, data_color, mae) {
  (p + ylab("Predicted uptake (ridge regression)")) %>% 
    annotate_plot(top_left_label, "top.left", data_color) %>% 
    annotate_plot(
      paste("MAE =", format(mae, digits=2), "g/L"),
      "bottom.right"
    )
}

parity_plot(trained_mod$y, pred_glmnet(trained_mod, trained_mod$orig_x), "#CA7C1B") %>% 
  annotate_scotty_parity(
    paste0("Training data\n", trained_mod$nfit, " hMOFs"),
    "#CA7C1B",
    hmof_partitioned_mod$training_fit["MAE"]
    ) %>% 
  save_ben_fig("scotty_parity_training.png")

hmof_partitioned_mod$pred_df %$%
  parity_plot(y_act, y_pred, "#0070C0") %>% 
  annotate_scotty_parity(
    paste0("Testing data\n", nrow(hmof_partitioned_mod$pred_df), " hMOFs"),
    "#0070C0",
    hmof_partitioned_mod$testing_fit["MAE"]
    ) %>% 
  save_ben_fig("scotty_parity_testing.png")

