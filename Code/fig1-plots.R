###############################
# fig1-plots.R
# Tom Smith
# 2021-05-19
#
# plot out all of the individual CUE responses for fig 1
#

library("ggplot2")

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

BoltzmannTref <- function(B0, E, T_ref, temp) {
  # B0 is the normalization constant (logged).  
  # E is the activation energy.
  calc <- B0 - E/k * (1/temp - 1/T_ref)
  return(calc)
}


# read the data
# read in the Sch fitting results for Tpk
good_fits <- read.csv("Results/Schoolfield_results.csv")

# also read the raw data
filtered_data <- read.csv("Data/filtered_data.csv")
# add CUE column
filtered_data$CUE <- filtered_data$biomass_increase / (filtered_data$biomass_increase + filtered_data$C_respired)

# also read the bootstrapped BA fits
bootstrapped_fits <- read.csv("Results/bootstrapped_BA_fits.csv")
# and the raw BA fits
raw_BA_fits <- read.csv("Results/BA_fits_phase1.csv")


# the for loop isn't working properly for geom_rect as
# R apparently is lazy and doesn't evaluate i until later
# so we have to make the plotting a function and run a
# for loop over the function later

plot_fun <- function(i) {
  force(i)
  subs_data <- filtered_data[filtered_data$Strain == strains[i],]
  boot_results <- bootstrapped_fits[bootstrapped_fits$strain == strains[i],]
  
  # get Sch so we can cut it to OTR
  sch_subs <- good_fits[good_fits$strain == strains[i],]

  # had to add in a weird catch here for the strain where Tpk doesnt occur within the data range
  temperatures <- seq(min(subs_data$Temperature), min(c(max(subs_data$Temperature), sch_subs$Tpk_growth-273.15)), 0.1)
  
  full_temperatures <- c()
  full_responses <- c()
  response_group <- c()
  
  for(ix in 1:length(boot_results$rep_num)){
    
    boot_subs_row <- boot_results[ix,]
    growth_curve <- exp(BoltzmannTref(B0 = log(boot_subs_row$B0_growth), E = boot_subs_row$E_growth, T_ref = 273.15,
                                      temp = temperatures+273.15))
    resp_curve <- exp(BoltzmannTref(B0 = log(boot_subs_row$B0_resp), E = boot_subs_row$E_resp, T_ref = 273.15,
                                    temp = temperatures+273.15))
    cue_curve <- growth_curve / (growth_curve + resp_curve)

    # collect results
    full_temperatures <- c(full_temperatures, temperatures)
    full_responses <- c(full_responses, cue_curve)
    response_group <- c(response_group, rep(i, length(temperatures)))
  }
  
  bootstrapped_cue <- data.frame(full_temperatures, full_responses, response_group)
  
  # iteratively take the 95% CIs at each temperature point
  # and median for line
  Temperature <- unique(bootstrapped_cue$full_temperatures)
  
  Upper_CUE <- c()
  Lower_CUE <- c()
  Median_CUE <- c()
  
  BA_results <- raw_BA_fits[raw_BA_fits$strain == strains[i],]
  
  for(ixx in 1:length(Temperature)){
    
    Median_CUE <- c(Median_CUE, median(bootstrapped_cue[bootstrapped_cue$full_temperatures == Temperature[ixx],]$full_responses,
                                       na.rm = TRUE))
    Upper_CUE <- c(Upper_CUE, quantile(bootstrapped_cue[bootstrapped_cue$full_temperatures == Temperature[ixx],]$full_responses,
                                       probs = c(0.975), na.rm = TRUE))
    Lower_CUE <- c(Lower_CUE, quantile(bootstrapped_cue[bootstrapped_cue$full_temperatures == Temperature[ixx],]$full_responses,
                                       probs = c(0.025), na.rm = TRUE))
  }
  
  new_CIs_frame <- data.frame(Temperature, Upper_CUE, Lower_CUE, Median_CUE)
  names(new_CIs_frame) <- c("Temperature", "Upper", "Lower", "Median")
  
  # define the plotting axes 
  data_min <- min(subs_data$CUE, na.rm = TRUE)
  bootstrap_min <- min(new_CIs_frame$Lower, na.rm = TRUE)
  data_max <- max(subs_data$CUE, na.rm = TRUE)
  bootstrap_max <- max(new_CIs_frame$Upper, na.rm = TRUE)
  
  plot_min <- min(c(data_min, bootstrap_min))
  plot_max <- max(c(data_max, bootstrap_max))
  
  p <- ggplot(subs_data, aes(x = Temperature))
  p <- p + geom_point(aes(y = CUE), size = 5)
  # only add the line and the rectangle if it actually fits in the plotting area
  if(sch_subs$Tpk_growth-273.15 < max(subs_data$Temperature)){
    #p <- p + geom_rect(aes(xmin = growth_sch$T_pk_est_sch-273.15, xmax = Inf, ymin = -Inf, ymax = Inf), col = "grey", alpha = 0.1)
    p <- p + annotate("rect", xmin = sch_subs$Tpk_growth-273.15, xmax = Inf, ymin = -Inf, ymax = Inf, col = "grey", alpha = 0.5)
    p <- p + geom_vline(xintercept = sch_subs$Tpk_growth-273.15, linetype = "dashed")
  }
  p <- p + geom_line(data = new_CIs_frame, aes(x = Temperature, y = Median), col = "red", size = 3, linetype = "dashed")
  #p <- p + geom_line(data = new_CIs_frame, aes(x = Temperature, y = Upper), col = "red", linetype = "dashed", size = 1)
  #p <- p + geom_line(data = new_CIs_frame, aes(x = Temperature, y = Lower), col = "red", linetype = "dashed", size = 1)
  p <- p + geom_ribbon(data = new_CIs_frame, aes(x = Temperature, ymin = Lower, ymax = Upper), alpha = 0.3, fill = "red")
  p <- p + coord_cartesian(ylim = c(plot_min, plot_max),
                           xlim = c(min(subs_data$Temperature, na.rm = TRUE), max(subs_data$Temperature, na.rm = TRUE)))
  p <- p + ggtitle(paste(strains[i], as.character(unique(subs_data$Genus)), sep = " "))
  p <- p + theme_bw() 
  p <- p + theme(aspect.ratio=1,
                 axis.title = element_blank(),
                 axis.text.x = element_text(size = 28),
                 axis.text.y = element_text(size = 28),
                 plot.title = element_text(size=28, vjust=1))
  
  # print out a figure for each strain individually
  # in order to have more control over arranging them in a figure later
  
  # Create a name for the output file using:
  #	- the original id number
  #   - the species name
  #   - the model
  output_name <- paste(strains[i],
                      as.character(unique(subs_data$Genus)),
                      sep = "_")
  
  # Remove any characters that won't look good in a file name,
  # using a regular expression.
  output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
  
  # Convert spaces to underscores.
  output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
  
  # CHANGE THIS to set an alternative output directory.
  outdir <- "Results/CUE_response_figures/"
  
  # Save it as an svg file.
  svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
  ggsave(filename = svg_file, plot = p, height = 6, width = 6)
}

# sneakily here we are going to use a list of strains in a
# different order, so they're plotted out the way we want
ordered_data <- read.csv("Results/CUE_directionality.csv")

strains <- as.character(unique(ordered_data$strain))

for (i in 1:length(strains)) {
  plot_fun(i)
}
