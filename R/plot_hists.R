# Plots either of the energy histogram or containing the energy histogram as a base layer

library(ggplot2)
library(dplyr)

source("R/regression_and_features.R")

BETA_H2_SCALING <- 1000.0
BETA_CH4_SCALING <- 10 * BETA_H2_SCALING
BETA_LEGEND_MARGIN <- 0.75  # pt.

# e.g. plot_hist_bins(filter(hist_vals, id==55), default_binspec)
plot_hist_bins <- function(one_grid, binspec, y_title = NULL, extend_top=FALSE) {
  # Returns an energy histogram plot, possibly as a base layer to beta coefficients
  # If y_title is NULL, skip plotting the primary y axis.
  # If extend_top is TRUE, make the top bin wider to denote its width
  result_data <- 
    one_grid %>% 
    stepped_hist_spec(binspec) %>% 
    mutate(height = metric) %>%
    filter(!near(height, 0)) %>%  # Remove columns where height == 0, floating point
    inner_join(color_from_binloc(bin_loc_from_spec(binspec)), by="bin") %>% 
    mutate(bar_width=1)
  if (extend_top) {
    result_data <- 
      result_data %>% 
      mutate(bar_width = ifelse((is.infinite(upper) & upper > 0), 1.5, bar_width)) %>% 
      # mutate(height = ifelse((is.infinite(upper) & upper > 0), 0.1, height)) %>%  # useful for debugging widths in this section
      mutate(loc = ifelse((is.infinite(upper) & upper > 0), loc+0.5/2.0, loc))
  }
  result <- 
    result_data %>% 
    ggplot(aes(loc, height)) +
    geom_col(aes(fill = I(color), width=bar_width)) +  # need I() so ggplot doesn't think the rgb strings are a factor
    guides(fill = FALSE) +
    labs(
      x = "Energy (kJ/mol)",
      y = y_title
    ) +
    theme(aspect.ratio = 0.6
    ) +
    scale_x_continuous(minor_breaks = NULL)
  if (is.null(y_title)) {
    result <- result +
      scale_y_continuous(
        labels = NULL,
        breaks = NULL
      )
  }
  result
}
# plot_hist_bins(filter(hist_sets$training, id==1), default_binspec) + coord_cartesian(ylim=c(-0.6, 0.6))

replot_hist_and_beta <- function(hist_beta_plot, hist_max = 0.5, undo_scaling=1.0) {
  # Conceptually similar to `overlay_cat_betas`, but resets the histogram heights to the secondary axis
  # from the combined plot.
  # Digs into code from `plot_hist_bins` and `overlay_cat_betas` (not DRY!).
  # Undoes the scaling applied to betas in the hist_beta_plot if specified.  Otherwise leaves them along.
  hist_data <- hist_beta_plot$data
  beta_data <- hist_beta_plot$layers[[2]]$data %>% mutate(beta = beta * undo_scaling)
  auto_scaling <- max(beta_data$beta) / hist_max
  
  # Scale everything in terms of the primary plot axis (beta)
  scaled_hist_data <- hist_data %>% mutate(height = height * auto_scaling)
  scaled_beta_data <- beta_data
  
  result <- ggplot(mapping = aes(x=loc)) +
    geom_col(data=scaled_hist_data, aes(fill=I(color), width=bar_width, y=height)) +
    geom_point(data=scaled_beta_data, aes(shape = shp, col = color, y = beta), size = 2) +
    guides(fill=FALSE) +
    scale_x_continuous(minor_breaks = NULL) +
    scale_y_continuous(
      name = expression(beta),
      sec.axis = sec_axis(
        ~ . / auto_scaling,
        name = "Volume fraction",
        breaks = c(0, hist_max)
       )
      ) +
    geom_hline(yintercept=0) +
    theme(axis.text.y.right = element_text(color="#8A226A"), axis.title.y.right = element_text(color="#8A226A")) +
    theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
    theme(
      legend.box.background = element_rect(color = "black"),
      legend.margin = margin(BETA_LEGEND_MARGIN, unit(5, "pt"), BETA_LEGEND_MARGIN, 0, "pt"),
      legend.title = element_blank()
    ) +
    labs(x="Energy (kJ/mol)", shape="", color="")
  result
}

overlay_violin_distr <- function(hist_plot, all_grids, binspec, color="gray", ...) {
  # Overlay a violin plot of the bin distributions on top of the bins.
  # Optionally pass other options such as color
  hist_plot + geom_violin(
    data = all_grids %>%
      stepped_hist_spec(binspec) %>% mutate(height = metric) %>%
      left_join(bin_loc_from_spec(binspec), by="bin"),
    aes(group = loc, y = height),
    scale="width",  # full width of bars, http://ggplot2.tidyverse.org/reference/geom_violin.html
    alpha=0, color=color,
    ...  # size=1, etc.
    )
}

plot_avg_with_distr <- function(all_grids, binspec, print_violin = FALSE, ...) {
  # Combines plot_hist_bins with overlay_violin_distr, and simplifies summary statistics calculations
  mean_grid <- 
    all_grids %>% 
    select(-id) %>% 
    group_by(lower, upper) %>% 
    summarize(counts = mean(counts)) %>% 
    ungroup %>% 
    mutate(id = "avg_grid")  # field is later used internally
  result <- mean_grid %>% plot_hist_bins(binspec, "Fraction of unit cell", extend_top=FALSE)
  if (print_violin) {
    result <- result %>% overlay_violin_distr(all_grids, binspec, ...)
  }
  result
}

plot_cat_betas <- function(betas, binspec) {
  # Successor to overlay_cat_betas (below) which only plots the points/symbols, not underlying histogram bars
  # This will simplify plots and explanations in the text
  # On second thought, maybe it would be clearer to keep the bars and just explain them better.
  beta_data <- betas %>%
    inner_join(bin_loc_from_spec(binspec), by="bin")
  if ("cat" %in% colnames(beta_data)) {
    beta_data <- beta_data %>% mutate(color = `cat`, shp = `cat`)
  } else {
    beta_data <- beta_data %>% mutate(color = I("darkgray"), shp = I(16))  # See `overlay_cat_betas` for info on shape 16
  }
  
  p <- beta_data %>% 
    ggplot(aes(x = loc)) +
    geom_point(aes(shape = shp, col = color, y = beta), size = 2) +
    scale_y_continuous(name = expression(beta)) +
    theme(axis.title.y = element_text(angle=0, vjust = 0.5)) +
    geom_hline(yintercept=0) +
    labs(shape="", color="", x="Energy (kJ/mol)")
  
  #p <- p + theme(legend.title = element_blank())
  p
}

overlay_cat_betas <- function(hist_plot, betas, binspec, scaling = BETA_H2_SCALING, hist_max = 0.5) {
  # Overlay betas from one or more categories, saved in the column `cat`, on top of a histogram plot
  # If `cat` is not defined, use a colorbar based on the magnitude of beta.
  
  cat_missing <- FALSE  # Category not specified.  What a sad variable name  :(
  color_if_cat_missing <- FALSE  # If true, color the beta points by their magnitude.  Else, darkgray.
  
  beta_data <- betas %>%
    inner_join(bin_loc_from_spec(binspec), by="bin") %>% 
    mutate(beta = beta / scaling)
  if ("cat" %in% colnames(beta_data)) {
    beta_data <- beta_data %>% mutate(color = `cat`, shp = `cat`)
  } else {
    cat_missing <- TRUE
    if (color_if_cat_missing) {
      beta_data <- beta_data %>% mutate(color = beta, shp = I(16))
    } else {
      beta_data <- beta_data %>% mutate(color = I("darkgray"), shp = I(16))
    }
    # Consider adding shapes based on color cat (or I(16)).
    # Why shape 16?  It's the default and drawn in the outline color (not fill)
    # See also http://sape.inf.usi.ch/quick-reference/ggplot2/shape
    # But also ggplot doesn't play nicely with combining the aesthetics, so maybe there's another way, like using a column as-is with I.
    # It looks like you can handle the shape below automatically by disabling the relevant guide titles or setting them equal.
    # https://stackoverflow.com/questions/37140266/how-to-merge-color-line-style-and-shape-legends-in-ggplot
    # http://environmentalcomputing.net/plotting-with-ggplot-colours-and-symbols/
  }
  
  p <- hist_plot +
    geom_point(
      data = beta_data,
      aes(y = beta, col = color, shape = shp),
      size = 2
      ) +
    coord_cartesian(ylim = c(-1.1*hist_max, 1.1*hist_max)) +
    scale_y_continuous(
      name = "Volume fraction",
      breaks = c(0, hist_max),
      sec.axis = sec_axis(
        ~.*scaling,
        name = expression(beta),
        breaks = c(-1*hist_max*scaling, 0, hist_max*scaling)
        )
      ) +
    theme(axis.title.y.right = element_text(angle=0, vjust = 0.5)) +
    geom_hline(yintercept=0) +
    labs(shape="", color="")
  
  if (cat_missing & color_if_cat_missing) {  # use colorbar
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
