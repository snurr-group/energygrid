# Plots either of the energy histogram or containing the energy histogram as a base layer

library(ggplot2)
library(dplyr)

parity_line <- geom_abline(slope=1, intercept=0, linetype="dashed")

parity_plot <- function(act, pred, color=1, alpha=0.05) {
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


# More functions:
# beta plot (with adjustable (automatic?) histogram bins)
# histogram density
