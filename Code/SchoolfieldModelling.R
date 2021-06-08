###############################
# SchoolfieldModelling.R
# Tom Smith
# 2021-05-19
#
# Fit Sharpe-Schoolfield TPC model to filtered rate data
#

library(minpack.lm) # for NLLS
library(ggplot2)

# assign Boltzmann constant (units of eV * K^-1) as a global:
assign("k", 8.617 * 10^-5, envir = .GlobalEnv)  

#### Functions #####

SchoolFTpk <- function(B0, E, E_D, T_pk, T_ref, temp)
{ # Sharpe-Schoolfield model with explicit T_pk parameter
  
  # PARAMETERS/INPUTS (all temperatures in Kelvin) -
  # temp   : temperature values to evaluate function at (single, scalar or vector of values)
  # B0     : Normalisation constant (log transformed)
  # E      : Activation energy (> 0)
  # E_D    : High temperature de-activation energy (> 0) 
  # T_ref  : Standardization (reference) temperature; set to 0 if not wanted   
  # T_pk   : Temperature at which trait reaches peak value
  
  return(B0 + log(exp((-E/k) * ((1/temp)-(1/T_ref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_pk - 1/temp)))))
}

#### Estimate STARTING VALUES for the nls

GetE <- function(tmp, rate, T_p, k)
{
  # Estimate starting value for E, taking linear regression using the rise part
  # of the curve only.
  # ~~~ Parameters ~~~
  # tmp  : temperature data (in K).
  # rate : rate data corresponding to temperature above.
  # T_p  : temperature at which rate peaks, used as a cutoff point.
  # k    : Boltzmann constant.
  
  tmp.w <- which(tmp <= T_p)
  if (length(tmp.w) > 1)
  {
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w])))) # include ludicrous error catching device
    msg <- tryCatch(return(abs(summary(m)$coefficients[2, 1])), error = function(e)e)
    if(length(msg) < 2){
      return(abs(summary(m)$coefficients[2, 1]))
    } else 
    {
      return(0.6)
    }
  } else
  {
    return(0.6)
  }
}

GetB0 <- function(tmp, rate)
{
  # Estimate starting value for the normalising constant.
  # ~~~ Parameters ~~~
  # tmp   : temperature data (in K).
  # rate  : rate data corresponding to temperature above.
  # T_ref : estimate normalising constant at this temperature (in K).
  
  if (min(tmp,na.rm=TRUE) > T_ref)
  {
    return(log(min(rate[1],na.rm=TRUE)))
  } else
  {
    return(log(max(rate[which(tmp <= T_ref)],na.rm=TRUE)))
  }
}


GetTpk <- function(tmp, rate)
{
  # Temperature at which the rate is maximised (estimate of T.peak).
  # ~~~ Parameters ~~~
  # tmp  : Temperature data (in K).
  # rate : Rate data corresponding to temperature above.
  
  return(max(tmp[which.max(rate)]))
}


##-------------------------------------------------------------------------------------##

# load in the data

data <- read.csv("Data/filtered_data.csv")

# set up for fitting
data$K <- data$Temperature+273.15 # better make a temperature column in kelvin

# set a reference temperature (note this is in Kelvin, set here at 0C)
T_ref <- 273.15

# pull out the IDs
IDs <- as.character(unique(data$Strain))

# Initialize empty vectors to store the parameter estimates
# that will be obtained.

volume <- c()
strain <- c()
trait <- c()
B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_pk_sch <- c()
P_pk_sch <- c()
r_sq_sch <- c()
n_temps <- c()
temps_before_peak <- c()
temps_after_peak <- c()

# now loop through the IDs and fit...

for(i in 1:length(IDs)){
  # subset the data
  subs <- data[data$Strain == IDs[i],]
  
  # only inlcude positive rates!
  subs <- subs[!is.na(subs$growth_rate_biomass),]
  subs <- subs[subs$growth_rate_biomass > 0,]
  subs <- subs[is.finite(subs$growth_rate_biomass),]
  
  strain <- c(strain, as.character(subs$Strain[[1]]))
  volume <- c(volume, (mean(subs$cell_volume_start, na.rm = TRUE) + mean(subs$cell_volume_end, na.rm = TRUE) )/2 )
  n_temps <- c(n_temps, length(unique(subs$Temperature)))
  trait <- "Growth Rate"
  
  # generate starting values for the model
  T_pk_st  <- GetTpk(tmp=subs$K, rate=subs$growth_rate_biomass)
  E_st    <- GetE(tmp=subs$K, rate=subs$growth_rate_biomass, T_p=T_pk_st, k = k)
  B_st <- GetB0(tmp=subs$K, rate=subs$growth_rate_biomass)
  
  # try to fit Schoolfield
  # first reset the loop
  schoolfield_nls <- NA
  try(schoolfield_nls <- nlsLM(
    log(growth_rate_biomass) ~ SchoolFTpk(B0, E, E_D, T_pk, T_ref, temp = K), data= subs, 
    start=list(B0 = B_st, E = E_st, E_D = 2*E_st, T_pk=T_pk_st)))
  
  # if fit didn't work, try again with different E_D estimate (I think this is causing problems)
  
  if(is.na(schoolfield_nls[1])){
    try(schoolfield_nls <- nlsLM(
      log(growth_rate_biomass) ~ SchoolFTpk(B0, E, E_D, T_pk, T_ref, temp = K), data= subs,
      start=list(B0 = B_st, E = E_st, E_D = 4*E_st, T_pk=T_pk_st)))
  }
  
  # If fitting worked ...
  if(!is.na(schoolfield_nls[1])){ 
    
    # Collect the parameter estimates...
    B0_sch <- c(B0_sch, exp(coef(schoolfield_nls)["B0"]))
    E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
    E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
    T_pk_sch <- c(T_pk_sch, coef(schoolfield_nls)["T_pk"])
    
    # Calculate the R squared value as: 1 - (rss/tss)
    rss <- sum((exp(predict(schoolfield_nls)) - 
                  subs$growth_rate_biomass)^2, 
               na.rm = TRUE)
    tss <- sum(
      (subs$growth_rate_biomass - 
         mean(subs$growth_rate_biomass, na.rm = TRUE))^2, 
      na.rm = TRUE)
    
    if ( tss != 0 )
    {
      r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
    } else
    {
      r_sq_sch <- c(r_sq_sch, 1)
    }
  
  
  # Calculate the peak of the curve
  max_rate <- NA
  try(max_rate <- exp(SchoolFTpk(B0 = coef(schoolfield_nls)["B0"],
                                 E = coef(schoolfield_nls)["E"],
                                 E_D = coef(schoolfield_nls)["E_D"],
                                 T_pk = coef(schoolfield_nls)["T_pk"],
                                 T_ref = 273.15,
                                 temp = coef(schoolfield_nls)["T_pk"])))
    
  P_pk_sch <- c(P_pk_sch, max_rate)
  
  # also check points before/after the peak
  
  temps_before_peak <- c(temps_before_peak, 
                         length(unique(subs[subs$Temperature <=  (coef(schoolfield_nls)[["T_pk"]] - 273.15),]$Temperature)))
  temps_after_peak <- c(temps_after_peak,
                        length(unique(subs[subs$Temperature >  (coef(schoolfield_nls)[["T_pk"]] - 273.15),]$Temperature)))
  
  
  ##############################
  # Plotting Schoolfield's fit #
  ########################
  ######
  
  # Create a name for the output file using:
  #	- the original id number
  #   - the species name
  #   - the model
  output_name <- paste(
    i, 
    subs$Strain[1],
    subs$Genus[1],
    subs$Species[1],
    "Growth_Rate_Biomass",
    #'Schoolfield',
    sep = "_"
  )
  
  # Remove any characters that won't look good in a file name,
  # using a regular expression.
  output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
  
  # Convert spaces to underscores.
  output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
  
  # CHANGE THIS to set an alternative output directory.
  outdir <- "Results/Schoolfield_fits/"
  
  # Generate predictions from the model fit...
  tmp_temps <- seq(min(
    floor(subs$K), na.rm = TRUE), 
    ceiling(max(subs$K, na.rm = TRUE)
    ), length = 200)
  
  tmp_model <- exp(SchoolFTpk(
    coef(schoolfield_nls)["B0"],
    coef(schoolfield_nls)["E"],
    coef(schoolfield_nls)["E_D"],
    coef(schoolfield_nls)["T_pk"],
    T_ref,
    tmp_temps
  ))
  
  ModelToPlot <- data.frame(
    Temperature = tmp_temps - 273.15, 
    TraitValue = tmp_model
  )
  
  # Prepare the data points of the original values.
  DataToPlot <- data.frame(
    Temperature = subs$K - 273.15, 
    TraitValue = subs$growth_rate_biomass
  )
  DataToPlot <- na.omit(DataToPlot)
  
  # Plot!
  p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                    y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                             alpha = 0.7, pch = 21) + 
    geom_line(data = ModelToPlot, 
              aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
              lwd = 1.3) + 
    annotate("text", x = min(DataToPlot$Temperature)+5, y = max(DataToPlot$TraitValue), label = paste("E = ", round(coef(schoolfield_nls)[["E"]], 3), sep = "")) +
    # add an annotation with E!
    ggtitle(paste(subs$Strain[1], subs$Genus[1], subs$Species[1], "Growth", sep = " ")) +
    xlab(expression(paste("Temperature (", degree, C, ")"))) + 
    ylab("Growth Rate (Biomass)") +
    theme_bw() + theme(plot.title = element_text(size = 12), 
                       axis.title = element_text(size = 10))
  
  # Save it as an svg file.
  svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
  ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)
  
} else # If fitting failed ...
{
  
  # Populate the vectors with missing values.
  B0_sch <- c(B0_sch, NA)
  E_sch <- c(E_sch, NA)
  E_D_sch <- c(E_D_sch, NA)	
  T_pk_sch <- c(T_pk_sch, NA)
  P_pk_sch <- c(P_pk_sch, NA)
  r_sq_sch <- c(r_sq_sch, NA)
  temps_before_peak <- c(temps_before_peak, NA)
  temps_after_peak <- c(temps_after_peak, NA)
}
}
  
  
# Compile all data into a data frame.
results <- data.frame(
    strain, volume, trait, E_sch, B0_sch, E_D_sch, 
    T_pk_sch, P_pk_sch, r_sq_sch, n_temps,
    temps_before_peak, temps_after_peak
)

# subsequently add the taxonomic and incubation information back onto the dataframe
strain_info <- read.csv("Data/isolates_growth_temperatures.csv")
growth_results <- merge(x = results, y = strain_info[,c("ID_code", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Iso_temp", "Incu_temp",
                                                "Gram_type", "Motility", "Spore_forming", "Shape", "Aero_anaero")],
                        by.x = "strain", by.y = "ID_code", all.x = TRUE)
  
# write this out
write.csv(growth_results, file = "Results/Schoolfield_growth.csv", row.names = FALSE)


# ----------------------------------------------------- #
# Repeat process for respiration rate fits
# ----------------------------------------------------- #

# Initialize empty vectors to store the parameter estimates
# that will be obtained.

volume <- c()
strain <- c()
trait <- c()
B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_pk_sch <- c()
P_pk_sch <- c()
r_sq_sch <- c()
n_temps <- c()
temps_before_peak <- c()
temps_after_peak <- c()

# now loop through the IDs and fit...

for(i in 1:length(IDs)){
  # subset the data
  subs <- data[data$Strain == IDs[i],]
  
  # only inlcude positive rates!
  subs <- subs[!is.na(subs$Resp_rate_per_biomass),]
  subs <- subs[subs$Resp_rate_per_biomass > 0,]
  subs <- subs[is.finite(subs$Resp_rate_per_biomass),]
  
  strain <- c(strain, as.character(subs$Strain[[1]]))
  volume <- c(volume, (mean(subs$cell_volume_start, na.rm = TRUE) + mean(subs$cell_volume_end, na.rm = TRUE) )/2 )
  n_temps <- c(n_temps, length(unique(subs$Temperature)))
  trait <- "Respiration Rate"
  
  # generate starting values for the model
  T_pk_st  <- GetTpk(tmp=subs$K, rate=subs$Resp_rate_per_biomass)
  E_st    <- GetE(tmp=subs$K, rate=subs$Resp_rate_per_biomass, T_p=T_pk_st, k = k)
  B_st <- GetB0(tmp=subs$K, rate=subs$Resp_rate_per_biomass)
  
  # try to fit Schoolfield
  # first reset the loop
  schoolfield_nls <- NA
  try(schoolfield_nls <- nlsLM(
    log(Resp_rate_per_biomass) ~ SchoolFTpk(B0, E, E_D, T_pk, T_ref, temp = K), data= subs, 
    start=list(B0 = B_st, E = E_st, E_D = 2*E_st, T_pk=T_pk_st)))
  
  # if fit didn't work, try again with different E_D estimate (I think this is causing problems)
  
  if(is.na(schoolfield_nls[1])){
    try(schoolfield_nls <- nlsLM(
      log(Resp_rate_per_biomass) ~ SchoolFTpk(B0, E, E_D, T_pk, T_ref, temp = K), data= subs,
      start=list(B0 = B_st, E = E_st, E_D = 4*E_st, T_pk=T_pk_st)))
  }
  
  # If fitting worked ...
  if(!is.na(schoolfield_nls[1])){ 
    
    # Collect the parameter estimates...
    B0_sch <- c(B0_sch, exp(coef(schoolfield_nls)["B0"]))
    E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
    E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
    T_pk_sch <- c(T_pk_sch, coef(schoolfield_nls)["T_pk"])
    
    # Calculate the R squared value as: 1 - (rss/tss)
    rss <- sum((exp(predict(schoolfield_nls)) - 
                  subs$Resp_rate_per_biomass)^2, 
               na.rm = TRUE)
    tss <- sum(
      (subs$Resp_rate_per_biomass - 
         mean(subs$Resp_rate_per_biomass, na.rm = TRUE))^2, 
      na.rm = TRUE)
    
    if ( tss != 0 )
    {
      r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
    } else
    {
      r_sq_sch <- c(r_sq_sch, 1)
    }
    
    
    # Calculate the peak of the curve
    max_rate <- NA
    try(max_rate <- exp(SchoolFTpk(B0 = coef(schoolfield_nls)["B0"],
                           E = coef(schoolfield_nls)["E"],
                           E_D = coef(schoolfield_nls)["E_D"],
                           T_pk = coef(schoolfield_nls)["T_pk"],
                           T_ref = 273.15,
                           temp = coef(schoolfield_nls)["T_pk"])))
    
    P_pk_sch <- c(P_pk_sch, max_rate)
    
    # also check points before/after the peak
    
    temps_before_peak <- c(temps_before_peak, 
                           length(unique(subs[subs$Temperature <=  (coef(schoolfield_nls)[["T_pk"]] - 273.15),]$Temperature)))
    temps_after_peak <- c(temps_after_peak,
                          length(unique(subs[subs$Temperature >  (coef(schoolfield_nls)[["T_pk"]] - 273.15),]$Temperature)))
    
    
    ##############################
    # Plotting Schoolfield's fit #
    ##############################
    
    # Create a name for the output file using:
    #	- the original id number
    #   - the species name
    #   - the model
    output_name <- paste(
      i, 
      subs$Strain[1],
      subs$Genus[1],
      subs$Species[1],
      "Resp_Rate",
      #'Schoolfield',
      sep = "_"
    )
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "Results/Schoolfield_fits/"
    
    # Generate predictions from the model fit...
    tmp_temps <- seq(min(
      floor(subs$K), na.rm = TRUE), 
      ceiling(max(subs$K, na.rm = TRUE)
      ), length = 200)
    
    tmp_model <- exp(SchoolFTpk(
      coef(schoolfield_nls)["B0"],
      coef(schoolfield_nls)["E"],
      coef(schoolfield_nls)["E_D"],
      coef(schoolfield_nls)["T_pk"],
      T_ref,
      tmp_temps
    ))
    
    ModelToPlot <- data.frame(
      Temperature = tmp_temps - 273.15, 
      TraitValue = tmp_model
    )
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
      Temperature = subs$K - 273.15, 
      TraitValue = subs$Resp_rate_per_biomass
    )
    DataToPlot <- na.omit(DataToPlot)
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                      y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot, 
                aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      annotate("text", x = min(DataToPlot$Temperature)+5, y = max(DataToPlot$TraitValue), label = paste("E = ", round(coef(schoolfield_nls)[["E"]], 3), sep = "")) +
      # add an annotation with E!
      ggtitle(paste(subs$Strain[1], subs$Genus[1], subs$Species[1], "Resp", sep = " ")) +
      xlab(expression(paste("Temperature (", degree, C, ")"))) + 
      ylab("Respiration Rate (per biomass)") +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10))
    
    # Save it as an svg file.
    svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
    ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)
    
  } else # If fitting failed ...
  {
    
    # Populate the vectors with missing values.
    B0_sch <- c(B0_sch, NA)
    E_sch <- c(E_sch, NA)
    E_D_sch <- c(E_D_sch, NA)	
    T_pk_sch <- c(T_pk_sch, NA)
    P_pk_sch <- c(P_pk_sch, NA)
    r_sq_sch <- c(r_sq_sch, NA)
    temps_before_peak <- c(temps_before_peak, NA)
    temps_after_peak <- c(temps_after_peak, NA)
  }
}


# Compile all data into a data frame.
results <- data.frame(
  strain, volume, trait, E_sch, B0_sch, E_D_sch, 
  T_pk_sch, P_pk_sch, r_sq_sch, n_temps,
  temps_before_peak, temps_after_peak
)

# subsequently add the taxonomic and incubation information back onto the dataframe
resp_results <- merge(x = results, y = strain_info[,c("ID_code", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Iso_temp", "Incu_temp",
                                                        "Gram_type", "Motility", "Spore_forming", "Shape", "Aero_anaero")],
                        by.x = "strain", by.y = "ID_code", all.x = TRUE)

# write this out
write.csv(resp_results, file = "Results/Schoolfield_resp.csv", row.names = FALSE)


# ------------------------------------------------------ #
# Merge resp and growth results, clean up and write out
# ------------------------------------------------------ #

# a bit of column renaming
names(resp_results) <- c("strain", "volume", "trait", "E_resp", "B0_resp", "E_D_resp", "Tpk_resp", "Ppk_resp", "r_sq_resp", "n_temps",
                         "temps_before_peak_resp", "temps_after_peak_resp", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                         "iso_temp", "incu_temp", "gram_type", "motility", "spore_forming", "shape", "aero_anaero")

names(growth_results) <- c("strain", "volume", "trait", "E_growth", "B0_growth", "E_D_growth", "Tpk_growth", "Ppk_growth", "r_sq_growth", "n_temps",
                         "temps_before_peak_growth", "temps_after_peak_growth", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                         "iso_temp", "incu_temp", "gram_type", "motility", "spore_forming", "shape", "aero_anaero")

combined_dataset <- merge(x = resp_results[, !names(resp_results) == "trait"], y = growth_results[c("strain", "E_growth", "B0_growth", "E_D_growth", "Tpk_growth",
                                                                                                    "Ppk_growth", "temps_before_peak_growth", "temps_after_peak_growth",
                                                                                                    "r_sq_growth")],
                          by = "strain")

# now cut it down to the strains which had a sensible Schoolfield fit
# and enough data points for estimates to be realistic
good_fits <- combined_dataset[combined_dataset$E_growth > 0 & combined_dataset$E_resp > 0 & combined_dataset$n_temps > 3 &
                                combined_dataset$temps_before_peak_growth > 2,]
good_fits <- good_fits[!is.na(good_fits$strain),]

# save combined data for further analysis

write.csv(good_fits, "Results/Schoolfield_results.csv", row.names = FALSE)
