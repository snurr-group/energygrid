packages = c("plotly","randomForest","cowplot", "ggplot2", "dplyr", "stringr", "magrittr", "readr", "tidyr", "openxlsx", "tidyverse", "R.utils", "manipulate", "slam", "gplots", "purrr", "glmnet", "caret", "testthat", "readxl","grid")
package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
    }
})
search()