# Plots either of the energy histogram or containing the energy histogram as a base layer

library(ggplot2)
library(dplyr)

parity_line <- geom_abline(slope=1, intercept=0, linetype="dashed")

parity_plot <- function(act, pred, color=1, alpha=0.05) {
  qplot(act, pred, alpha=I(alpha), color=I(color)) +
    xlab("'Actual' uptake (GCMC simulations)") +
    ylab("Predicted uptake (ridge regression)") +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(limits = c(0,50)) +
    scale_y_continuous(limits = c(0,50)) +
    parity_line
}


# More functions:
# beta plot (with adjustable (automatic?) histogram bins)
# histogram density
