###############################
# Filtering.R
# Tom Smith
# 2021-05-19
#
# Filter the raw data to remove poor quality rate estimates
# and keep only the aerobic strains for CUE analyses
#

source("Code/packages.R")

# read the raw data
raw_data <- read.csv("Data/dataset_rates_taxonomy.csv")

# add CUE column
raw_data$CUE <- raw_data$biomass_increase / (raw_data$biomass_increase + raw_data$C_respired)

# and get rid of the crappy data
good_data <- raw_data[raw_data$cells_per_ml_start > 2000000,] # too few cells measured = unreliable rate estimates
good_data <- good_data[good_data$growth_rate > 0,] # keep only positive rates
good_data <- good_data[good_data$growth_rate_biomass > 0,]
good_data <- good_data[good_data$Resp_rate > 0,]
good_data <- good_data[!(good_data$Strain %in% c("30_30_01", "30_30_06", "30_30_07")),] # cell counts from these streptomyces unreliable
good_data <- good_data[!is.na(good_data$Strain),]

# merge with functional data
functional_traits <- read.csv("Data/isolates_growth_temperatures.csv")

combined_data <- merge(x = good_data, y = functional_traits[,c("ID_code", "Gram_type", "Motility", "Spore_forming", "Shape", "Aero_anaero")], 
                       by.x = "Strain", by.y = "ID_code", all.x = TRUE)

# cut to the aerobes only
aerobic_data <- combined_data[combined_data$Aero_anaero == "Aerobic",]

# write this out
write.csv(aerobic_data, "Data/filtered_data.csv", row.names = FALSE)
