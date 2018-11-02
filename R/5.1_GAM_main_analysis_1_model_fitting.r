################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Oct 23 2018
################################################################################

#
# Description: This is where the two models are fit to the data, both in
# using maximum likelihood (package 'mgcv') and Hamiltonian Monte Carlo 
# (package 'brms'). The fitted objects are exported to files, so that they
# can be further used for predictions and inference.

# Note: All predictors are scaled to 0 mean and variance of 1. Also, some
# predictors (Area_km, Tree_dens) are log transformed.

################################################################################ 
# LOAD THE DATA AND THE PACKAGES
################################################################################

source("0_libraries_functions_settings.r")
source("4.1_Data_loading_standardization_and_centering.r")


################################################################################ 
# DEFINING THE R MODEL FORMULAS
################################################################################

REALM.formula <- S ~ 
  REALM + poly(Area_km,3):REALM + # biogeographic realm
  Tree_dens + Tree_dens:Area_km + # tree density
  min_DBH + min_DBH:Area_km +     # minimum DBH
  GPP + GPP:Area_km +             # gross primary productivity
  ANN_T + ANN_T:Area_km +         # mean annual temperature
  ISO_T + ISO_T:Area_km +         # mean isothermality 
  MIN_P + MIN_P:Area_km +         # minimum precipitation in the driest quarter of the year
  P_SEAS + P_SEAS:Area_km +       # mean precipitation seasonality
  ALT_DIF + ALT_DIF:Area_km +     # altitude difference
  ISLAND + ISLAND:Area_km         # insularity


SMOOTH.formula <- S ~ 
  s(Lat, Lon, by=DAT_TYPE, bs="sos", k=14) + # smooth region effect
  poly(Area_km, 3) +              # species-area relationship
  Tree_dens + Tree_dens:Area_km + # tree density
  min_DBH + min_DBH:Area_km +     # minimum DBH
  GPP + GPP:Area_km +             # gross primary productivity
  ANN_T + ANN_T:Area_km +         # mean annual temperature
  ISO_T + ISO_T:Area_km +         # mean isothermality 
  MIN_P + MIN_P:Area_km +         # minimum precipitation in the driest quarter of the year
  P_SEAS + P_SEAS:Area_km +       # mean precipitation seasonality
  ALT_DIF + ALT_DIF:Area_km +     # altitude difference
  ISLAND + ISLAND:Area_km         # insularity



################################################################################ 
# FIT THE MODELS using 'mgcv' function 'gam'
################################################################################

gam.REALM <- gam(REALM.formula, data=DAT, family="nb")
    summary(gam.REALM)
    
    # save to a file
    save(gam.REALM, file="../Models/gam_REALM.Rdata")


gam.SMOOTH <- gam(SMOOTH.formula, data = DAT, family="nb")
    summary(gam.SMOOTH)
    
    # save to a file
    save(gam.SMOOTH, file="../Models/gam_SMOOTH.Rdata")


# additionally, fit model with only the intercept (a 'null' model)
gam.NULL <- gam(S~1, data=DAT, family="nb")


    # compare the models using AIC and BIC
    AIC(gam.NULL, gam.REALM, gam.SMOOTH)
    BIC(gam.NULL, gam.REALM, gam.SMOOTH)



################################################################################
# FIT THE MODELS IN STAN, using 'brms' function 'brm'
################################################################################

# Beware, this will take > 1 hour to run on a 2GHz 2-core laptop.

brm.REALM <- brm(REALM.formula, family="negbinomial", data=DAT,
                     cores=2,
                     seed=12355,
                     chains=3, iter=3000, warmup=1000, thin=10)
    
# save to a file
save(brm.REALM, file="../Models/brms_REALM.RData")
    
    
brm.SMOOTH <- brm(SMOOTH.formula, family="negbinomial", data=DAT,
                  cores=3,
                  seed=12355,
                  chains=3, iter=3000, warmup=1000, thin=10)
    
    # save to a file
    save(brm.SMOOTH , file="../Models/brms_SMOOTH.RData")
    

# ------------------------------------------------------------------------------
# FIRST DIAGNOSTICS    

# are the predictions from the brms and mgcv models the same?
prd.gam <- predict(gam.SMOOTH, type="response")
prd <- predict(brm.SMOOTH, newdata=DAT)[,'Estimate']
plot(prd, prd.gam)
abline(a=0, b=1)

prd.gam <- predict(gam.REALM, type="response")
prd <- predict(brm.REALM, newdata=DAT)[,'Estimate']
plot(prd, prd.gam)
abline(a=0, b=1)
    
        
    
# ------------------------------------------------------------------------------
# TRACEPLOTS and CATERPILLAR PLOTS of model parameters

pars.REALM <- rownames(data.frame(summary(brm.REALM$fit)[1]))[1:46]
pars.SMOOTH <- rownames(data.frame(summary(brm.SMOOTH$fit)[1]))[1:51]

# traceplots - this indicate if the convergence is good (which it is)

# save to a file
png("../Figures/traceplot_REALM.png", width = 2000, height=2000, res=150)
rstan::traceplot(brm.REALM$fit, pars=pars.REALM)
dev.off()

# save to a file
png("../Figures/traceplot_SMOOTH.png", width = 2000, height=2000, res=150)
rstan::traceplot(brm.SMOOTH$fit, pars=pars.SMOOTH)
dev.off()

# caterpillar plots

# save to a file
png("../Figures/caterpillar_REALM.png", width = 2000, height=2000, res=150)
rstan::plot(brm.REALM$fit, pars=pars.REALM)
dev.off()

png("../Figures/caterpillar_SMOOTH.png", width = 2000, height=2000, res=150)
rstan::plot(brm.SMOOTH$fit, pars=pars.SMOOTH)
dev.off()




