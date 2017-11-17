# Sets the colormap for energy histogram bins

library(viridis)  # Use the "inferno" colormap, which has nice printing properties, etc.
library(dplyr)
library(purrr)

full_e_colormap <- viridis.map %>% 
  filter(opt == "B") %>%   #inferno
  select(-opt)

interpolate_row <- function(df, double_row, allow_overflow = TRUE) {
  # Interpolates the values between two rows of a dataframe
  # If allow_overflow is set to true, then double_row that exceeds dimensions of the data frame will be set to the extremes
  if (allow_overflow) {
    if (double_row < 1) {
      return(df[1,])
    } else if (double_row >= nrow(df)) {
      return(df[nrow(df),])
    }
  }  # otherwise, the code will just naturally error out
  mixing <- double_row - floor(double_row)  # fraction of second in the linear combination
  sm_df <- df[c(floor(double_row), ceiling(double_row)), ]
    #summarize_all(funs(mean))  # mean of the two columns
  sm_df[1,] * (1-mixing) + sm_df[2,] * mixing
}
# Example:
# full_e_colormap[2:3,]
# full_e_colormap %>% interpolate_row(2.99)

color_from_binloc <- function(binloc, threshold = 10, top_col = 200) {
  # Assigns RGB colors to binloc, likely generated from regression_and_features.R:bin_loc_from_spec.
  # threshold specifies the lower energy bound (kJ/mol) for maxing out the color intensity.
  # top_col (if not false) specifies which color (1-256) should be used for an Inf bin.
  
  # 256 rows of colors, so midpoint (0 kJ/mol) is 128
  binloc %>% 
    mutate(loc = ifelse(top_col & bin=="Inf", Inf, loc)) %>%  # flag top bin if applicable
    .$loc %>% as.list %>%
    # Get fraction of threshold (+/-), move to positive frac, scale by the half colorbar, and make it one-indexed
    map_dbl(~ifelse(is.infinite(.), top_col, (./threshold + 1)*128 + 1)) %>%  # If top bin, manually color
    map_dfr(interpolate_row, df=full_e_colormap) %>% 
    mutate(color = rgb(R, G, B)) %>% 
    bind_cols(binloc, .)
}
# Example:
# default_binspec %>% bin_loc_from_spec %>% color_from_binloc
