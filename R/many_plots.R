plot_lj <- function(lower, upper) {
  upper_lab <- ifelse(upper == 25, Inf, upper)
  hist_vals %>%
    metric_from_hists(lower, upper) %>%
    left_join(gcmc_data, by="id") %>%
    ggplot(aes(metric, g.L, col=void.frac)) +
    geom_point() +
    theme_bw(16) +
    annotate("text", Inf, Inf, label = paste0("LJ cutoff: (", lower, ", ", upper_lab, ")"), hjust=1, vjust=1)  # https://github.com/tidyverse/ggplot2/issues/1244
}

export_many <- function(many_lower, many_upper) {
  dir.create("Scotty")
  for (lower in many_lower) {
    for (upper in many_upper) {
        ggsave(paste0("Scotty/plot_", lower, "_", upper, ".png"), plot_lj(lower, upper))
    }
  }
}

# For earlier exporting, I used a "Lower" scan of (-30:10:-500, -20)
# and an upper of (-1000, 0:5:-300)

# use seq(bot, top, step) to feed into the export cmd

### COLORATION WORK

# What columns do we have available?
# hist_vals %>% metric_from_hists(-1000, -10) %>% left_join(gcmc_data, by="id") %>% colnames
# To add a color to the plot, just change the ggplot command:
# ggplot(aes(metric, g.L, col=void.frac))



