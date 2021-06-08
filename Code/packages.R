silent.require <- function(x) suppressMessages(require(package=x, character.only=TRUE, quietly=TRUE))

# Load packages that are already installed
packages <- c("minpack.lm", "xtable", "tidyr", "dplyr", "ggplot2", "nls.multstart"
              )

ready <- sapply(packages, silent.require)
