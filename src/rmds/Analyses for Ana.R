

###########################
## SET WORKING DIRECTORY ##
###########################

setwd("~/Ana SDM")

######################
## LIBRARY PACKAGES ##
######################

library(plot3D)
library(mapproj)
library(brms)
library(bayesplot)
library(ggplot2)
library(dplyr)

#########################################################
## FUNCTION TO STANDARDIZE VARIABLES TO MEAN 0, SD 0.5 ##
#########################################################

z_score_2sd <- function(x){ (x - mean(x,na.rm=T)) / (2*sd(x,na.rm=T))}

###################
## LOAD THE DATA ##
###################

shark_data <- read.table("shark_data.txt")

###############################################
## DO ANY VARIABLES NEED TO BE TRANSFORMED ? ##
###############################################

num_vars <- shark_data %>%
  select(dist_mangrove:pop_20km)

graphics.off()
par(mfrow=c(3,3))
for(i in 1:ncol(num_vars)){
  
  num_var_min <- min(num_vars[,i], na.rm=TRUE)
  
  num_var_log <- if(num_var_min > 0) { log(num_vars[,i])
  } else {
    log(num_vars[,i]+ceiling(abs(num_var_min))+1)
  }
  
  hist(num_vars[,i], main=NA)
  title(colnames(num_vars)[i])
  
  hist(num_var_log, main=NA)
  title(paste("log", colnames(num_vars)[i]))
  
  hist((num_vars[,i]^2), main=NA)
  title(paste("square", colnames(num_vars)[i]))
  
}

# GRAVITY METRICS SHOULD BE LOGGED, POP_20, MAYBE DIST TO MANGROVE, NEARPOP AND NEARTT IF WE USED THEM

summarise(shark_data,
          across(c('dist_mangrove',
                   'neartt',
                   'nearpop',
                   'neargrav1',
                   'neargrav2',
                   'neargrav3',
                   'gravtot1',
                   'gravtot2',
                   'gravtot3',
                   'pop_20km'),
                 min, na.rm=TRUE))

shark_data <- mutate(shark_data,
                     across(c('dist_mangrove',
                              'neartt',
                              'nearpop',
                              'neargrav1',
                              'neargrav2',
                              'neargrav3',
                              'gravtot1',
                              'gravtot2',
                              'gravtot3'),
                            log))

shark_data$pop_20km <- log(shark_data$pop_20km + 1)

#################################
##  LOOK AT ONLY SHARK PRESENT ##
#################################

shark_present <- shark_data %>%
  filter(shark_present==1)

#################################
## IDENTIFY CONSERVATION AREAS ##
#################################

# AREAS NEAR MANGROVES WITH HIGH HUMAN IMPACTS - GOOD POTENTIAL FOR THE SHARK, BUT IN DANGER ??
# TOP 25% GRAVITY, TOP 25% CLOSEST TO MANGROVE

cons_areas <- shark_present %>%
  filter(poptot >= quantile(poptot, probs=0.75) &
           mangrove_area <= quantile(mangrove_area, probs=0.25))

graphics.off()
plot(shark_present$Lon, shark_present$Lat, cex=0, ylab=NA, xlab=NA)
map(database = "world",fill=TRUE,col="beige", add=T)
corners <- par("usr")
polygon(c(corners[1], corners[1], corners[2], corners[2]),
        c(corners[3], corners[4], corners[4], corners[3]),lwd=1.35)
points(shark_present$Lon, shark_present$Lat, cex=2, pch=21, col=1, bg="grey")
points(cons_areas$Lon, cons_areas$Lat, cex=2, pch=21, col=1, bg="red")
legend("topright",cex=1,pch=19, pt.cex=1.25, col=c("grey","red"),
       legend = c("Shark Present", "Priority Area"))


#########################################
## TRY TO DO A MODEL AT THE SITE LEVEL ##-------------------------------------------------------------------
#########################################
###########################################################
## SIMPLE LINEAR MODEL (NO RANDOM EFFECTS) AT SITE LEVEL ##
###########################################################

sapply(shark_data, function(x) paste(round(sum(is.na(x))/length(x),2)*100,"%",sep=""))

shark_data <- mutate(shark_data,
                       across(c('dist_mangrove':'pop_20km'),
                              z_score_2sd)) 

site_model <- brm(shark_present ~ 
                    dist_mangrove +
                    neargrav2 + 
                    gravtot2 +
                    nearpop +
                    poptot,
                  #c(set_prior(class="Intercept", "normal(0,1)"),
                  #  set_prior(class="b", "normal(0,1)")),
                  family=bernoulli, data=shark_data, chains=3, iter=2000)

summary(site_model, prob=0.90)
site_drivers <- data.frame(as.matrix(site_model))
site_drivers <- site_drivers %>%
  select('b_dist_mangrove','b_gravtot2', 'b_neargrav2', 'b_nearpop', 'b_poptot')
mcmc_intervals(site_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

conditional_effects(site_model)

##########################################
## IDENTIFY COUNTRY FOR EACH COORDINATE ##
##########################################

# GOING TO NEED TO GET AN EEZ SHAPEFILE FOR THIS

##################################
## LOAD COUNTRY LEVEL VARIABLES ##
##################################

country_vars <- read.csv("country_variables.csv")
names(country_vars)
country_vars <- country_vars %>%
  rename(country = ï..country)

##########################################
## SIMPLE LINEAR MODEL AT COUNTRY LEVEL ##
##########################################

country_vars <- mutate(country_vars,
                     across(c('hdi':'hd_intertidal'),
                            z_score_2sd)) 

country_model <- brm(occurrence ~ 
                       hdi +
                       voice +
                       strength_of_governance +
                       habitat_condition_of_mangroves,
                     family=bernoulli, data=country_vars, chains=3, iter=2000)

summary(country_model, prob=0.90)
country_drivers <- data.frame(as.matrix(country_model))
country_drivers <- country_drivers %>%
  select('b_hdi', 'b_voice', 'b_strength_of_governance',
         'b_habitat_condition_of_mangroves')
mcmc_intervals(country_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0)

conditional_effects(country_model)


#######################################
## MERGE COUNTRY AND SITE LEVEL DATA ##
#######################################

point_country <- read.table("point_country.txt")
point_country <- point_country %>%
  select(cell, TERRITORY1) %>%
  rename(country = TERRITORY1)

shark_data <- merge(shark_data, point_country, by="cell")

shark_data <- merge(shark_data, country_vars, by="country", all=TRUE)

##############################
## TRY A HIERARCHICAL MODEL ##
##############################

hier_model <- site_model <- brm(shark_present ~ 
                                  dist_mangrove +
                                  #dist_port +
                                  poptot + 
                                  hdi +
                                  #gdp +
                                  voice +
                                  habitat_condition_of_mangroves +
                                  (1 | country),
                                c(set_prior(class="Intercept", "normal(0,1)"),
                                  set_prior(class="b", "normal(0,1)")),
                                family=bernoulli, data=shark_data, chains=2, iter=4000,
                                control = list(adapt_delta=0.90))

color_scheme_set("darkgray")
summary(hier_model, prob=0.90)
hier_drivers <- data.frame(as.matrix(hier_model))
hier_drivers <- hier_drivers %>%
  select('b_dist_mangrove', 'b_poptot',
         'b_hdi', 'b_voice', 'b_habitat_condition_of_mangroves')
hier_drivers <- hier_drivers %>%
  rename('Distance to Mangroves' = b_dist_mangrove,
         'Human Population Density' = b_poptot,
         'Human Development Index' = b_hdi,
         'Voice and Accountability' = b_voice,
         'Mangroves Habitat Condition' = b_habitat_condition_of_mangroves)
mcmc_intervals(hier_drivers, point_est = "median", prob = 0.5, prob_outer = 0.90,
               outer_size = 1,
               inner_size = 4,
               point_size = 6) + geom_vline(xintercept = 0) +
  labs(x = "Standardized Effect Size") +
  theme(axis.text.y = element_text(size = 13, family="Sans",face="plain"),
        axis.text.x = element_text(size = 12, family="Sans",face="plain"),
        axis.title.x = element_text(size = 13, family="Sans",face="plain")) 
  
conditional_effects(hier_model)
performance::r2_bayes(hier_model)



###################################
## DIGGING INTO DISTANCE TO PORT ##
###################################

graphics.off()
par(mfrow=c(1,2))
scatter2D(shark_data$Lon, shark_data$Lat, colvar=shark_data$dist_mangrove, cex=0, ylab=NA, xlab=NA)
map(database = "world",fill=TRUE,col="beige", add=T)
corners <- par("usr")
polygon(c(corners[1], corners[1], corners[2], corners[2]),
        c(corners[3], corners[4], corners[4], corners[3]),lwd=1.35)
scatter2D(shark_data$Lon, shark_data$Lat, colvar=shark_data$dist_mangrove, cex=1, ylab=NA, xlab=NA, add=T)
title("Distance from Port")


plot(shark_data$Lon, shark_data$Lat, cex=0, ylab=NA, xlab=NA)
map(database = "world",fill=TRUE,col="beige", add=T)
corners <- par("usr")
polygon(c(corners[1], corners[1], corners[2], corners[2]),
        c(corners[3], corners[4], corners[4], corners[3]),lwd=1.35)
points(shark_data$Lon, shark_data$Lat, cex=1, pch=21, col=1, bg="grey")
points(shark_present$Lon, shark_present$Lat, cex=1, pch=21, col=1, bg="red")
legend("topright",cex=1,pch=19, pt.cex=1.25, col=c("grey","red"),
       legend = c("Possible Habitat", "Shark Occurence"))


