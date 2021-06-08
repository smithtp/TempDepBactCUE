###############################
# BAModelling.R
# Tom Smith
# 2021-05-19
#
# Fit Boltzmann-Arrhenius model to OTR region of TPCs
#

library(minpack.lm) # for NLLS
library(nls.multstart)

# assign Boltzmann constant (units of eV * K^-1) as a global:
assign("k", 8.617 * 10^-5, envir = .GlobalEnv)  

BoltzmannTref <- function(B0, E, T_ref, temp) {
  # B0 is the normalization constant (logged).  
  # E is the activation energy.
  calc <- B0 - E/k * (1/temp - 1/T_ref)
  return(calc)
}

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


# load in the data

data <- read.csv("Data/filtered_data.csv")

# set up for fitting
data$K <- data$Temperature+273.15 # better make a temperature column in kelvin

# also read Schoolfield fits
good_fits <- read.csv("Results/Schoolfield_results.csv")

# bring data down to only those that we already determined Schoolfield parameters for
data <- data[data$Strain %in% good_fits$strain,]

# set a reference temperature (note this is in Kelvin, set here at 0C)
T_ref <- 273.15

# pull out the IDs
IDs <- as.character(unique(good_fits$strain))

strain <- c()
B0_growth <- c()
E_growth <- c()
E_growth_min <- c()
E_growth_max <- c()
B0_resp <- c()
E_resp <- c()
E_resp_min <- c()
E_resp_max <- c()
phylum <- c()
class <- c()
order <- c()
family <- c()
genus <- c()
species <- c()

for(i in 1: length(IDs)){
  
  subs <- data[data$Strain == IDs[i],] # subset the data
  # only inlcude positive rates!
  subs <- subs[!is.na(subs$growth_rate_biomass),]
  subs <- subs[subs$growth_rate_biomass > 0,]
  subs <- subs[is.finite(subs$growth_rate_biomass),]
  
  phylum <- c(phylum, as.character(unique(subs$Phylum)))
  class <- c(class, as.character(unique(subs$Class)))
  order <- c(order, as.character(unique(subs$Order)))
  family <- c(family, as.character(unique(subs$Family)))
  genus <- c(genus, as.character(unique(subs$Genus)))
  species <- c(species, as.character(unique(subs$Species)))
  strain <- c(strain, IDs[i])
  
  growth_subs <- subs[subs$Temperature <= good_fits[good_fits$strain == IDs[i],]$Tpk_growth-273.15,]
  resp_subs <- subs[subs$Temperature <= good_fits[good_fits$strain == IDs[i],]$Tpk_resp-273.15,]
  
  # now fit BA for growth 
  
  # first reset the loop
  boltzmann_nls <- NA
  # try to use loads of start parameter combinations
  # from dan padfields package
  boltzmann_nls <- nls_multstart(
    log(growth_rate_biomass) ~ BoltzmannTref(B0, E, T_ref = 273.15, temp = K), data= growth_subs,
    iter = 100,
    start_lower = c(B0 = -10, E = 0.1),
    start_upper = c(B0 = 10, E = 2),
    supp_errors = 'Y',
    convergence_count = 50,
    na.action = na.omit,
    lower = c(B0 = -20, E = 0))
  
  # If fitting worked ...
  if(length(boltzmann_nls) > 0){ 
    
    # Collect the parameter estimates...
    B0_growth <- c(B0_growth, exp(coef(boltzmann_nls)["B0"]))
    E_growth <- c(E_growth, coef(boltzmann_nls)["E"])
    
    growth_confints <- confint(boltzmann_nls)
    E_growth_min <- c(E_growth_min, growth_confints[2,1])
    E_growth_max <- c(E_growth_max, growth_confints[2,2])
    
  } else # If fitting failed ...
  {
    
    # Populate the vectors with missing values.
    B0_growth <- c(B0_growth, NA)
    E_growth <- c(E_growth, NA)
    E_growth_min <- c(E_growth_min, NA)
    E_growth_max <- c(E_growth_max, NA)
    
  }
  
  # repeat the whole thing for respiration rate
  
  boltzmann_nls <- NA
  # try to use loads of start parameter combinations
  # from dan padfields package
  boltzmann_nls <- nls_multstart(
    log(Resp_rate_per_biomass) ~ BoltzmannTref(B0, E, T_ref = 273.15, temp = K), data= resp_subs,
    iter = 100,
    start_lower = c(B0 = -10, E = 0.1),
    start_upper = c(B0 = 10, E = 2),
    supp_errors = 'Y',
    convergence_count = 50,
    na.action = na.omit,
    lower = c(B0 = -10, E = 0))
  
  # If fitting worked ...
  if(length(boltzmann_nls) > 0){ 
    
    # Collect the parameter estimates...
    B0_resp <- c(B0_resp, exp(coef(boltzmann_nls)["B0"]))
    E_resp <- c(E_resp, coef(boltzmann_nls)["E"])
    
    resp_confints <- confint(boltzmann_nls)
    E_resp_min <- c(E_resp_min, resp_confints[2,1])
    E_resp_max <- c(E_resp_max, resp_confints[2,2])
    
    
  } else # If fitting failed ...
  {
    
    # Populate the vectors with missing values.
    B0_resp <- c(B0_resp, NA)
    E_resp <- c(E_resp, NA)
    E_resp_min <- c(E_resp_min, NA)
    E_resp_max <- c(E_resp_max, NA)
  }
}

BA_results <- data.frame(
  strain, phylum, class, order, family, genus, species,
  E_growth, E_growth_min, E_growth_max, B0_growth, E_resp, E_resp_min, E_resp_max, B0_resp)

write.csv(BA_results, "Results/BA_fits_phase1.csv", row.names = FALSE)


# --------------------------------------------------------- #
# Now bootstrap these fits to generate confidence intervals
# --------------------------------------------------------- #

strain <- c()
rep_num <- c()
B0_growth <- c()
E_growth <- c()
B0_resp <- c()
E_resp <- c()
phylum <- c()
class <- c()
order <- c()
family <- c()
genus <- c()
species <- c()

for(i in 1: length(IDs)){
  
  subs <- data[data$Strain == IDs[i],] # subset the data
  # only inlcude positive rates!
  subs <- subs[!is.na(subs$growth_rate_biomass),]
  subs <- subs[subs$growth_rate_biomass > 0,]
  subs <- subs[is.finite(subs$growth_rate_biomass),]
  
  
  for(ix in 1:1000){ # number of bootstraps
    
    rep_num <- c(rep_num, ix)
    strain <- c(strain, IDs[i])
    phylum <- c(phylum, as.character(unique(subs$Phylum)))
    class <- c(class, as.character(unique(subs$Class)))
    order <- c(order, as.character(unique(subs$Order)))
    family <- c(family, as.character(unique(subs$Family)))
    genus <- c(genus, as.character(unique(subs$Genus)))
    species <- c(species, as.character(unique(subs$Species)))
    
    # now resample
    # doing it with rows means that the pair-wise measurements are
    # kept together.. though I'm not sure how much this matters
    sample_data <- subs[sample(nrow(subs), replace = TRUE),]
    
    # cut it to the OTR using previous Schoolfield fits
    growth_sample <- sample_data[sample_data$Temperature <= good_fits[good_fits$strain == IDs[i],]$Tpk_growth-273.15,]
    resp_sample <- sample_data[sample_data$Temperature <= good_fits[good_fits$strain == IDs[i],]$Tpk_resp-273.15,]
    
    # now fit BA for growth 
    
    # first reset the loop
    boltzmann_nls <- NA

    # do it the old-fashioned way rather than multstart
    # generate starting values for the model
    E_st    <- GetE(tmp=growth_sample$K, rate=growth_sample$growth_rate_biomass, T_p = max(growth_sample$K), k = k)
    B_st <- GetB0(tmp=growth_sample$K, rate=growth_sample$growth_rate_biomass)
    # try to fit
    try(boltzmann_nls <- nlsLM(
      log(growth_rate_biomass) ~ BoltzmannTref(B0, E, T_ref = 273.15, temp = K), data= growth_sample,
      start = list(B0 = B_st, E = E_st)), silent = TRUE)
    
    # If fitting worked ...
    if(!is.na(boltzmann_nls[1])){ 
      
      # Collect the parameter estimates...
      B0_growth <- c(B0_growth, exp(coef(boltzmann_nls)["B0"]))
      E_growth <- c(E_growth, coef(boltzmann_nls)["E"])
      
    } else # If fitting failed ...
    {
      
      # Populate the vectors with missing values.
      B0_growth <- c(B0_growth, NA)
      E_growth <- c(E_growth, NA)
    }
    
    #### repeat the whole thing for respiration rate ####
    
    boltzmann_nls <- NA

    E_st    <- GetE(tmp=resp_sample$K, rate=resp_sample$Resp_rate_per_biomass, T_p = max(resp_sample$K), k = k)
    B_st <- GetB0(tmp=resp_sample$K, rate=resp_sample$Resp_rate_per_biomass)
    # try to use loads of start parameter combinations
    # from dan padfields package
    try(boltzmann_nls <- nlsLM(
      log(Resp_rate_per_biomass) ~ BoltzmannTref(B0, E, T_ref = 273.15, temp = K), data= resp_sample,
      start = list(B0 = B_st, E = E_st)), silent = TRUE)
    
    
    # If fitting worked ...
    if(!is.na(boltzmann_nls[1])){  
      
      # Collect the parameter estimates...
      B0_resp <- c(B0_resp, exp(coef(boltzmann_nls)["B0"]))
      E_resp <- c(E_resp, coef(boltzmann_nls)["E"])
      
    } else # If fitting failed ...
    {
      
      # Populate the vectors with missing values.
      B0_resp <- c(B0_resp, NA)
      E_resp <- c(E_resp, NA)
    }
    
    # lets record how many times we've iterated, for the sake of sanity
    print(paste("Strain = ", i, ", Rep = ", ix, sep = ""))
  }
  
}

BA_results <- data.frame(
  strain, rep_num, phylum, class, order, family, genus, species,
  E_growth, B0_growth, E_resp, B0_resp)

# add CUE_E and CUE_B0 values here
BA_results$E_CUE <- BA_results$E_growth - ((BA_results$E_growth * BA_results$B0_growth +
                                              BA_results$E_resp * BA_results$B0_resp) / 
                                             (BA_results$B0_growth + BA_results$B0_resp))

BA_results$B0_CUE <- BA_results$B0_growth / (BA_results$B0_growth + BA_results$B0_resp)

write.csv(BA_results, "Results/bootstrapped_BA_fits.csv", row.names = FALSE)

# now collect the medians and CIs from the distributions for each parameter
BA_results <- read.csv("Results/bootstrapped_BA_fits.csv")

IDs <- as.character(unique(BA_results$strain))

# initialise a lot of parameters
strain <- c()
E_growth_median <- c()
E_growth_min <- c()
E_growth_max <- c()
B0_growth_median <- c()
B0_growth_min <- c()
B0_growth_max <- c()
E_resp_median <- c()
E_resp_min <- c()
E_resp_max <- c()
B0_resp_median <- c()
B0_resp_min <- c()
B0_resp_max <- c()
E_CUE_median <- c()
E_CUE_min <- c()
E_CUE_max <- c()
B0_CUE_median <- c()
B0_CUE_min <- c()
B0_CUE_max <- c()
phylum <- c()
class <- c()
order <- c()
family <- c()
genus <- c()
species <- c()

for(i in 1:length(IDs)){
  
  subs_data <- BA_results[BA_results$strain == IDs[i],]
  
  strain <- c(strain, IDs[i])
  phylum <- c(phylum, as.character(unique(subs_data$phylum)))
  class <- c(class, as.character(unique(subs_data$class)))
  order <- c(order, as.character(unique(subs_data$order)))
  family <- c(family, as.character(unique(subs_data$family)))
  genus <- c(genus, as.character(unique(subs_data$genus)))
  species <- c(species, as.character(unique(subs_data$species)))
  
  E_growth_median <- c(E_growth_median, median(subs_data$E_growth, na.rm = TRUE))
  E_growth_min <- c(E_growth_min, quantile(subs_data$E_growth, probs = c(0.025), na.rm = TRUE))
  E_growth_max <- c(E_growth_max, quantile(subs_data$E_growth, probs = c(0.975), na.rm = TRUE))
  B0_growth_median <- c(B0_growth_median, median(subs_data$B0_growth, na.rm = TRUE))
  B0_growth_min <- c(B0_growth_min, quantile(subs_data$B0_growth, probs = c(0.025), na.rm = TRUE))
  B0_growth_max <- c(B0_growth_max, quantile(subs_data$B0_growth, probs = c(0.975), na.rm = TRUE))
  E_resp_median <- c(E_resp_median, median(subs_data$E_resp, na.rm = TRUE))
  E_resp_min <- c(E_resp_min, quantile(subs_data$E_resp, probs = c(0.025), na.rm = TRUE))
  E_resp_max <- c(E_resp_max, quantile(subs_data$E_resp, probs = c(0.975), na.rm = TRUE))
  B0_resp_median <- c(B0_resp_median, median(subs_data$B0_resp, na.rm = TRUE))
  B0_resp_min <- c(B0_resp_min, quantile(subs_data$B0_resp, probs = c(0.025), na.rm = TRUE))
  B0_resp_max <- c(B0_resp_max, quantile(subs_data$B0_resp, probs = c(0.975), na.rm = TRUE))
  E_CUE_median <- c(E_CUE_median, median(subs_data$E_CUE, na.rm = TRUE))
  E_CUE_min <- c(E_CUE_min, quantile(subs_data$E_CUE, probs = c(0.025), na.rm = TRUE))
  E_CUE_max <- c(E_CUE_max, quantile(subs_data$E_CUE, probs = c(0.975), na.rm = TRUE))
  B0_CUE_median <- c(B0_CUE_median, median(subs_data$B0_CUE, na.rm = TRUE))
  B0_CUE_min <- c(B0_CUE_min, quantile(subs_data$B0_CUE, probs = c(0.025), na.rm = TRUE))
  B0_CUE_max <- c(B0_CUE_max, quantile(subs_data$B0_CUE, probs = c(0.975), na.rm = TRUE))
}

BA_boot_results <- data.frame(strain, E_growth_median, E_growth_min, E_growth_max, B0_growth_median, B0_growth_min,
                              B0_growth_max, E_resp_median, E_resp_min, E_resp_max, B0_resp_median,
                              B0_resp_min, B0_resp_max, E_CUE_median, E_CUE_min, E_CUE_max, B0_CUE_median,
                              B0_CUE_min, B0_CUE_max, phylum, class, order, family, genus, species)

write.csv(BA_boot_results, "Results/bootstrap_results_BA.csv", row.names = FALSE)

