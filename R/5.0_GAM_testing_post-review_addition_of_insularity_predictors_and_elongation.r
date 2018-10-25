################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Oct 23 2018
################################################################################

# Description: Here we test performance of predictors that were added to the
# dataset after the first round of review in Nature Ecology & Evolution.
# Specifically, these are the 3 different definition of island/mainland, and
# the ELONG variable that measures elongation of the spatial units.

# The way these scripts were used is an interactive one, i.e. there is no linear
# workflow, but we were simply manually fiddling with the code, noting perfomance
# metrics of the models with different predictors.

################################################################################

# clean the workspace and load the libraries and the data
source("0_libraries_functions_settings.r")
source("4.1_Data_loading_standardization_and_centering.r")



################################################################################ 
# MODEL FORMULAS
################################################################################

NULL.formula <- S ~ 1


REALM.formula <- S ~ REALM + poly(Area_km,3):REALM +  
  Tree_dens + Tree_dens:Area_km + 
  min_DBH + min_DBH:Area_km +
  GPP + GPP:Area_km + 
  ANN_T + ANN_T:Area_km +
  ISO_T + ISO_T:Area_km + 
  MIN_P + MIN_P:Area_km +  
  P_SEAS + P_SEAS:Area_km + 
  ALT_DIF + ALT_DIF:Area_km +
  ISL_LS + ISL_LS:Area_km +
  # ISL_ST + ISL_ST:Area_km #+
  # ISL_DIS + ISL_DIS:Area_km #+
  # ELONG + ELONG:Area_km


SMOOTH.formula <- S ~ s(Lat, Lon, by=DAT_TYPE, bs="sos", k=14) +
  poly(Area_km, 3) +  
  Tree_dens + Tree_dens:Area_km + 
  min_DBH + min_DBH:Area_km +
  GPP + GPP:Area_km + 
  ANN_T + ANN_T:Area_km +
  ISO_T + ISO_T:Area_km + 
  MIN_P + MIN_P:Area_km +  
  P_SEAS + P_SEAS:Area_km + 
  ALT_DIF + ALT_DIF:Area_km +
  ISL_LS + ISL_LS:Area_km #+
  # ISL_ST + ISL_ST:Area_km# +
  # ISL_DIS + ISL_DIS:Area_km# +
  # ELONG + ELONG:Area_km



################################################################################

summary(gam(S ~ Tree_dens, data=DAT, family="nb"))

# manually assessing the models

model.formula = SMOOTH.formula
m1 <- gam(model.formula, data=DAT, family="nb")
summary(m1)
AIC(m1)
BIC(m1)
















