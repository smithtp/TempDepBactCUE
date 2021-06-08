###############################
# make-figures.R
# Tom Smith
# 2021-05-19
#
# Create all the plots for figs 2 and 3 in the manuscript
#

source("Code/packages.R")

# read the necessary data/results files
results_table <- read.csv("Results/CUE_directionality.csv")


#############################################
# ------ CUE directionality Figures --------#
#############################################

main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16))

# reverse the factor to plot from highest at the top to lowest at the bottom
results_table$strain <- factor(results_table$strain,
                               levels = as.character(results_table$strain))
results_table$strain <- factor(results_table$strain, levels=rev(results_table$strain))
results_table$result <- factor(results_table$result, levels = c("Increasing", "Invariant", "Decreasing"))


CUE_E_plot <- ggplot(results_table, aes(x = E_CUE_median, y = strain)) + 
  geom_errorbarh(aes(xmax=E_CUE_max, xmin=E_CUE_min, col = phylum), size = 1, show.legend = FALSE, height = 0, alpha = 0.8) +
  geom_point(aes(col = phylum), size = 4, show.legend = TRUE) +
  geom_vline(aes(xintercept=0), linetype = 'dashed', size = 1) +
  scale_colour_manual(values=c("red", "blue", "orange")) +
  scale_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
  main_theme +
  # ylab("Strain") +
  # xlab(expression(italic(E)[CUE])) +
  labs(x = "", y = "") +
  theme(axis.title = element_blank(),
        legend.title=element_blank(), 
        legend.position=c(0.7,0.5),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 12))
CUE_E_plot

ggsave("Results/Figures/Fig2A.png", CUE_E_plot, height = 6, width = 6)


# need to join max CUE with the E_CUE results for fig 2B
good_data <- read.csv("Data/filtered_data.csv")

max_CUEs <- good_data %>%
  group_by(Strain) %>%
  dplyr::summarise(max_CUE = max(CUE, na.rm = TRUE))


# read the phase-1 boltzmann results
BA_fits <- read.csv("Results/BA_fits_phase1.csv")

# E of the combined BA model, according to Tom C's maths
# can be described by this:

CUE_model_E <- function(B0_mu, E_mu, B0_R, E_R){
  
  return(E_mu - (E_mu*B0_mu + E_R*B0_R)/(B0_mu + B0_R))
  
}

# B0 of the combined BA model
# is described as:

CUE_model_B0 <- function(B0_mu, B0_R){
  
  return(B0_mu / (B0_mu + B0_R))
  
}

# find E for the combined model
BA_fits$CUE_E <- CUE_model_E(B0_mu = BA_fits$B0_growth, E_mu = BA_fits$E_growth, 
                             B0_R = BA_fits$B0_resp, E_R = BA_fits$E_resp)

# find B0 for the combined model

BA_fits$CUE_B0 <- CUE_model_B0(B0_mu = BA_fits$B0_growth, B0_R = BA_fits$B0_resp)
aerobic_results <- merge(BA_fits, max_CUEs, by.x = "strain", by.y = "Strain")

# linear or polynomial regression
lin_mod <- lm(CUE_E ~ max_CUE, data = aerobic_results)
poly_mod <- lm(CUE_E ~ poly(max_CUE, 2), data = aerobic_results)

summary(lin_mod)
summary(poly_mod)

AIC(lin_mod)
AIC(poly_mod)
# linear is preferred by AIC

CUE_E_vs_maxCUE <- ggplot(aerobic_results, aes(x = max_CUE, y = CUE_E)) + 
  geom_point(aes(fill = phylum), size = 3, shape = 21) +
  scale_fill_manual(values=c("red", "blue", "orange")) +
  geom_smooth(method = lm,  col = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # labs(x = "Maximum CUE" , y = expression(italic(E)[CUE])) +
  labs(x = "", y = "") +
  main_theme +
  theme(legend.position = c(0.75, 0.85),
        legend.title = element_blank(),
        aspect.ratio = 1)
CUE_E_vs_maxCUE

ggsave("Results/Figures/Fig2B.png", CUE_E_vs_maxCUE, height = 6, width = 6)


################################################
# ------ Growth vs Respiration Figures --------#
################################################

######################
####    Tpk?    ######
######################
# 

# merge schoolfield data in for peaks
sch_fits <- read.csv("Results/Schoolfield_results.csv")
aerobic_results <- merge(aerobic_results, sch_fits[,c("strain", "temps_before_peak_growth", "temps_after_peak_growth", "temps_before_peak_resp",
                                                      "temps_after_peak_resp", "Tpk_resp", "Tpk_growth")], by = "strain")

# # subset so peaks make sense:
Tpk_data_biomass <- aerobic_results[aerobic_results$temps_after_peak_resp > 0 &  aerobic_results$temps_after_peak_growth > 0,]

# check distribution
ggplot(Tpk_data_biomass, aes(x = Tpk_resp)) + geom_histogram()
ggplot(Tpk_data_biomass, aes(x = Tpk_growth)) + geom_histogram()
# a bit right skewed - maybe not t-test?

Tpk_plot <- ggplot(Tpk_data_biomass, aes(x = Tpk_resp-273.15, y = Tpk_growth-273.15)) +
  geom_point(size = 3) +
  geom_abline(slope = 1) +
  xlim(20, 60) +
  ylim(20, 60) +
  xlab("") +
  ylab("") +
  main_theme +
  theme(aspect.ratio = 1)
Tpk_plot

# Tpk stats
linear_model <- lm(Tpk_data_biomass$Tpk_growth ~ Tpk_data_biomass$Tpk_resp)
summary(linear_model) 
confint(linear_model)

# Pearson's correlation test
cor.test(Tpk_data_biomass$Tpk_resp, Tpk_data_biomass$Tpk_growth)
# paired t test
t.test(Tpk_data_biomass$Tpk_resp, Tpk_data_biomass$Tpk_growth, paired = TRUE)
# two-tailed test is significant that means arent different
t.test(Tpk_data_biomass$Tpk_resp, Tpk_data_biomass$Tpk_growth, paired = TRUE, alternative = "greater")
# however lower p-value when we ask whether resp is lower than growth

# ######################
# ####     E?     ######
# ######################
# 
# # subset differently, but so there's enough points to estimate E from
E_data_biomass <- aerobic_results[aerobic_results$temps_before_peak_growth > 2 & aerobic_results$temps_before_peak_resp > 2,]

E_plot <- ggplot(E_data_biomass, aes(x = E_resp, y = E_growth)) +
  geom_point(size = 3) +
  geom_abline(slope = 1) +
  xlim(0.25, 1.4) +
  ylim(0.25, 1.75) +
  labs(x = "",
       y = "") +
  main_theme +
  theme(aspect.ratio = 1)
E_plot


# E stats
linear_model <- lm(E_data_biomass$E_growth ~ E_data_biomass$E_resp)
summary(linear_model) # is significant
confint(linear_model) # confidence interval is somewhat wide, but does include 1

# Pearson's correlation test
cor.test(E_data_biomass$E_resp, E_data_biomass$E_growth)
# paired t test
t.test(E_data_biomass$E_resp, E_data_biomass$E_growth, paired = TRUE)
# two-tailed test is significant that means arent different
t.test(E_data_biomass$E_growth, E_data_biomass$E_resp, paired = TRUE, alternative = "greater")
# however lower p-value when we ask whether resp is lower than growth E

########################
# Distribution
#########################

# plot the distribution of E for growth and resp
# first need to make the E data "long" rather than "wide"

E_data_growth <- E_data_biomass[,c("strain", "E_growth", "phylum")]
names(E_data_growth) <- c("strain", "E", "Phylum")
E_data_growth$Trait <- "Growth Rate"

E_data_resp <- E_data_biomass[,c("strain", "E_resp", "phylum")]
names(E_data_resp) <- c("strain", "E", "Phylum")
E_data_resp$Trait <- "Respiration Rate"

E_data <- rbind(E_data_growth, E_data_resp)

E_distribution_label <- ggplot(E_data, aes(x = E)) +
  geom_density(aes(fill = Trait), alpha = 0.6) +
  labs(x = "", y = "") +
  xlim(0, 1.6) +
  main_theme +
  #geom_vline(xintercept = mean(E_data[E_data$Trait == "Growth Rate",]$E), col = "red") +
  geom_vline(xintercept = median(E_data[E_data$Trait == "Growth Rate",]$E), col = "red", linetype = "dashed") +
  #geom_vline(xintercept = mean(E_data[E_data$Trait == "Respiration Rate",]$E), col = "blue") +
  geom_vline(xintercept = median(E_data[E_data$Trait == "Respiration Rate",]$E), col = "blue", linetype = "dashed") +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.8),
        aspect.ratio = 1) 
E_distribution_label

# save plots seperately

ggsave("Results/Figures/Fig3A.png", Tpk_plot, width = 6, height = 6)
ggsave("Results/Figures/Fig3B.png", E_plot, width = 6, height = 6)
ggsave("Results/Figures/Fig3C.png", E_distribution_label, width = 6, height = 6)


###############################
# Data-synthesis Distribution
###############################

# this is the flux data only
flux_means <- read.csv("Data/group_means_non_aggregated_fluxes.csv")
flux_means <- flux_means[flux_means$ConKingdom == "Bacteria",]

# this is the growth data
growth_means <- read.csv("Data/all_means_non_aggregated.csv")
growth_means <- growth_means[growth_means$ConKingdom == "Bacteria",]
growth_means$Trait <- "Growth Rate"

# These data are from our data synthesis paper:
# "Community-level respiration of prokaryotic microbes may rise with global warming"
# these can be obtained from: https://github.com/smithtp/hotterbetterprokaryotes

# make the fluxes all one "trait"
flux_means$Trait <- "Metabolic Flux"

# subset down to only the few columns we need
flux_data <- flux_means[,c("Species", "Trait", "E")]
growth_data <- growth_means[,c("Species", "Trait", "E")]

hib_data <- rbind(growth_data, flux_data)

E_distribution <- ggplot(hib_data, aes(x = E)) +
  geom_density(aes(fill = Trait), alpha = 0.6) +
  labs(x = "", y = "") +
  xlim(0, 1.6) + 
  #geom_vline(xintercept = mean(E_data[E_data$Trait == "Growth Rate",]$E), col = "red") +
  geom_vline(xintercept = median(hib_data[hib_data$Trait == "Growth Rate",]$E), col = "red", linetype = "dashed") +
  #geom_vline(xintercept = mean(E_data[E_data$Trait == "Respiration Rate",]$E), col = "blue") +
  geom_vline(xintercept = median(hib_data[hib_data$Trait == "Metabolic Flux",]$E), col = "blue", linetype = "dashed") +
  main_theme +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.8),
        aspect.ratio = 1) 
E_distribution

ggsave("Results/Figures/Fig3D.png", E_distribution, height = 6, width = 6)

# t-test
t.test(E ~ Trait, data = hib_data)
