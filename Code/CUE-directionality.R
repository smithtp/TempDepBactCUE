###############################
# CUE-directionality.R
# Tom Smith
# 2021-05-19
#
# Determine directionality of CUE response based on bootstrapped fits
#

source("Code/packages.R")

bootstrapped_fits <- read.csv("Results/bootstrap_results_BA.csv")

names(bootstrapped_fits)[2] <- "E_growth"
names(bootstrapped_fits)[8] <- "E_resp"

# function to check the overlaps
check_overlaps <- function(dataset){
  
  overlaps <- dataset
  
  # this is simple, does the CI overlap zero?
  
  overlaps$result <- "Invariant"
  overlaps[overlaps$E_CUE_min > 0,]$result <- "Increasing"
  overlaps[overlaps$E_CUE_max < 0,]$result <- "Decreasing"
  
  return(overlaps)
  
}

overlap_results <- check_overlaps(bootstrapped_fits)

# final results table is then:
results_table <- overlap_results[,c("strain", "phylum", "class", "order", "family", "genus", "species", "E_growth", "E_resp", "E_CUE_median", "E_CUE_min", "E_CUE_max", "result")]

# order it from incresing to decreasing
results_table$result <- as.factor(results_table$result)

results_table$result <- factor(results_table$result, levels = c("Increasing", "Invariant", "Decreasing"))

results_table <- results_table[order(results_table$result, -results_table$E_CUE_median),]

# write it to file
write.csv(results_table, "Results/CUE_directionality.csv", row.names = FALSE)

# print a nice latex table version for supplement
print(xtable(results_table[,c("strain", "E_growth", "E_resp", "E_CUE_median", "E_CUE_min", "E_CUE_max", "result")]))
