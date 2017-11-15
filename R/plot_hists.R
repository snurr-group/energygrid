# Plots either of the energy histogram or containing the energy histogram as a base layer

library(ggplot2)
library(dplyr)

# e.g. plot_hist_bins(filter(hist_vals, id==55), default_binspec)
plot_hist_bins <- function(one_grid, binspec) {
  # Returns an energy histogram plot, possibly as a base layer to beta coefficients
  # TODO: add in bin colors, which will require calculating them first
  one_grid %>% 
    stepped_hist_spec(binspec) %>% 
    mutate(height = metric) %>%
    inner_join(bin_loc_from_spec(binspec), by="bin") %>% 
    # also join on a color matrix, or at least cbind it
    # left_join(trained_betas, by="bin") %>% 
    ggplot(aes(loc, height)) +
    geom_col() +
    #geom_col(fill=rgb(HIST_COLORS_JUNK)) +  # THIS IS THE BIGGEST ONE
    labs(
      x = "Energy (kJ/mol)",
      y = NULL
    ) +
    theme(aspect.ratio = 0.6
    ) +
    scale_y_continuous(
      labels = NULL,
      breaks = NULL
    )
}
# plot_hist_bins(filter(hist_sets$training, id==1), default_binspec) + coord_cartesian(ylim=c(-0.6, 0.6))

# Add in betas as an additive plot layer?  TODO
overlay_color_betas <- function(hist_plot, betas, binspec) {
  # Function that takes a histogram plot and overlays the model betas in a single color
  # TODO: stub
  hist_plot
}

overlay_cat_betas <- function(hist_plot, betas, binspec, scaling = 10.0, hist_max = 0.5) {
  # Overlay betas from one or more categories, saved in the column `cat`, on top of a histogram plot
  beta_data <- betas %>%
    inner_join(bin_loc_from_spec(binspec), by="bin") %>% 
    mutate(beta = beta / scaling)
  if (!("cat" %in% colnames(beta_data))) {
    beta_data <- beta_data %>% mutate(`cat` = "1")
  }
  
  hist_plot +
    geom_point(
      data = beta_data,
      aes(y = beta, col = cat),
      size = 3
      ) +
    coord_cartesian(ylim = c(-1.1*hist_max, 1.1*hist_max)) +
    scale_y_continuous(
      labels = NULL,
      breaks = NULL,
      minor_breaks = c(-1*hist_max, hist_max),
      sec.axis = sec_axis(
        ~.*scaling,
        name = expression(beta),
        breaks = c(-1*hist_max*scaling, hist_max*scaling)
        )
      ) +
    theme(axis.title.y.right = element_text(angle=0, vjust = 0.5)) +
    geom_hline(yintercept=0)
}

# More functions:
# histogram density: merge content from partial_plot_histogram_for_hists.R (figures in model explanation)

# Other plots in old notebooks:
# Heatmap randomization cartoon
