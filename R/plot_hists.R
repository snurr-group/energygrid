# Plots either of the energy histogram or containing the energy histogram as a base layer

library(ggplot2)
library(dplyr)

source("R/regression_and_features.R")

# e.g. plot_hist_bins(filter(hist_vals, id==55), default_binspec)
plot_hist_bins <- function(one_grid, binspec) {
  # Returns an energy histogram plot, possibly as a base layer to beta coefficients
  one_grid %>% 
    stepped_hist_spec(binspec) %>% 
    mutate(height = metric) %>%
    inner_join(color_from_binloc(bin_loc_from_spec(binspec)), by="bin") %>% 
    ggplot(aes(loc, height)) +
    geom_col(aes(fill = I(color))) +  # need I() so ggplot doesn't think the rgb strings are a factor
    guides(fill = FALSE) +
    labs(
      x = "Energy (kJ/mol)",
      y = NULL
    ) +
    theme(aspect.ratio = 0.6
    ) +
    scale_y_continuous(
      labels = NULL,
      breaks = NULL
    ) +
    scale_x_continuous(minor_breaks = NULL)
}
# plot_hist_bins(filter(hist_sets$training, id==1), default_binspec) + coord_cartesian(ylim=c(-0.6, 0.6))

overlay_cat_betas <- function(hist_plot, betas, binspec, scaling = 10.0, hist_max = 0.5) {
  # Overlay betas from one or more categories, saved in the column `cat`, on top of a histogram plot
  # If `cat` is not defined, use a colorbar based on the magnitude of beta.
  
  cat_missing <- FALSE  # Category not specified.  What a sad variable name  :(
  
  beta_data <- betas %>%
    inner_join(bin_loc_from_spec(binspec), by="bin") %>% 
    mutate(beta = beta / scaling)
  if ("cat" %in% colnames(beta_data)) {
    beta_data <- beta_data %>% mutate(color = `cat`)
  } else {
    cat_missing <- TRUE
    beta_data <- beta_data %>% mutate(color = beta)
  }
  
  p <- hist_plot +
    geom_point(
      data = beta_data,
      aes(y = beta, col = color),
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
        breaks = c(-1*hist_max*scaling, 0, hist_max*scaling)
        )
      ) +
    theme(axis.title.y.right = element_text(angle=0, vjust = 0.5)) +
    geom_hline(yintercept=0)
  
  if (cat_missing) {  # use colorbar
    p <- p +
      scale_color_gradientn(colors = c("red", "darkgray", "blue"), guide = "none") +
      theme(axis.text.y.right = element_text(color = c("red", "darkgray", "blue")))
  } else {
    p <- p + theme(legend.title = element_blank())
  }
  p
}

theme_diagram_min <- 
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(size=2)
  )

make_minimal <- function(p, ymax = 0.25) {
  p +
    theme_diagram_min +
    scale_y_continuous() +
    coord_cartesian(ylim = c(0, ymax))
}

plot_hist_curve <- function(one_grid) {
  # Plots an abstract depiction of the potential energy distribution in a single MOF.
  # Note: there used to be a related geom_histogram, but that's been deprecated for the main histogram function above
  truncated_grid <- one_grid %>% 
    mutate(dens = counts / sum(counts))
  top_dens <- truncated_grid[nrow(truncated_grid), "dens"]
  truncated_grid <- truncated_grid[-nrow(truncated_grid),]
  truncated_grid %>% 
    ggplot(aes(lower, weight = dens)) +
    stat_density(
      aes(col=I(rgb(78, 42, 132, maxColorValue=255))),
      adjust = 1/5,
      size = 3,
      geom = "line"  # Can't use geom_density directly without an ugly line along the x axis
      )
}

plot_mof_minimal <- function(all_grid, mof_id, binspec, ymax = 0.25) {
  # Plots the energy distribution and histogram cartoons for a given id
  mof_data <- all_grid %>% filter(id == mof_id)
  mof_data %>% 
    plot_hist_curve %>% 
    make_minimal(ymax) %>% 
    print  # this will raise a warning about densities!=1.
  # Can hide the warnings/messages by saving the plot and wrapping print with suppressWarnings, etc.
  mof_data %>% 
    plot_hist_bins(binspec) %>% 
    make_minimal(ymax) %>% 
    print
  invisible(NULL)  # return something so the second plot isn't shown twice
}
# Example: plot_mof_minimal(hmof_h2_grid, 1, default_binspec)

# Other plots in old notebooks:
# Heatmap randomization cartoon
