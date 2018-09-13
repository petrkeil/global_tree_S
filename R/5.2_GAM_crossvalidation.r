################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Sep 11 2018
################################################################################

#
# Description: 
#

################################################################################ 
# 0. LOAD THE DATA AND THE PACKAGES
################################################################################


# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")

# ------------------------------------------------------------------------------

# load the data
PLT <- read.csv("../Data/Main_dataset_full_detail.csv")


################################################################################ 
# 1. PREPARE THE DATA FOR THE ANALYSES
################################################################################

# calculate tree density (note the x+1 step!!)
PLT$Tree_dens <- (PLT$N + 1) / PLT$Area_km 

# select only the variables of interest from the larger data.frame
DAT <- dplyr::select(PLT, S, Area_km, Tree_dens, min_DBH=min_DBH_cm, 
                     GPP, ET, ANN_T, WARM_T, ISO_T, MIN_P, P_SEAS, ALT_DIF,
                     ISLAND, REALM=REALM_PK, Lat, Lon, DAT_TYPE, Loc_ID) 

# order the data frame by regions and area
DAT <- DAT[order(DAT$REALM),]

# transform area and tree density
DAT$Area_km <- log(DAT$Area_km)
DAT$Tree_dens <- log(DAT$Tree_dens)

# remove NAs
DAT <- DAT[rowSums(is.na(DAT)) == 0,]

# mean and sd of Area (will be used later to bring these to their origina scale)
A.mean <- mean(DAT$Area_km)
A.sd <- sd(DAT$Area_km)

# means and sd of the rest of the variables
centr <- attributes(scale(DAT[,2:12]))$'scaled:center'
scale <- attributes(scale(DAT[,2:12]))$'scaled:scale'
scale.tab <- data.frame(var=names(centr), centr, scale)
write.csv(scale.tab, file="scale_tab.csv", row.names = FALSE)

# do the actual scaling
DAT[,2:12] <- scale(DAT[,2:12])



################################################################################ 
# 2. THE MODEL FORMULAS
################################################################################

REALM.formula <- S ~ REALM + poly(Area_km,3):REALM +  
  Tree_dens + Tree_dens:Area_km + 
  min_DBH + min_DBH:Area_km +
  GPP + GPP:Area_km + 
  ANN_T + ANN_T:Area_km +
  ISO_T + ISO_T:Area_km + 
  MIN_P + MIN_P:Area_km +  
  P_SEAS + P_SEAS:Area_km + 
  ALT_DIF + ALT_DIF:Area_km +
  ISLAND + ISLAND:Area_km


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
  ISLAND + ISLAND:Area_km

# ------------------------------------------------------------------------------

nfolds <- 4
case.folds <- rep(1:nfolds,length.out=nrow(DAT))
# divide the cases as evenly as possible
case.folds <- sample(case.folds) # randomly permute the order

par(mfrow=c(2,2))
for (fold in 1:nfolds) {
  # What are the training cases and what are the test cases?
  train <- DAT[case.folds!=fold,]
  test <- DAT[case.folds==fold,]
  
  gam.train <- gam(SMOOTH.formula, data=train, family="nb")
  preds <- predict.gam(gam.train, newdata = test, type="response")
  plot( log(test$S),log(preds), col=test$DAT_TYPE); abline(a=0, b=1)
}



################################################################################ 
# 3. FIT THE TRAINING MODELS
################################################################################

gam.REALM <- gam(REALM.formula, data=DAT, family="nb")
summary(gam.REALM)
save(gam.REALM, file="../STAN_models/gam_REALM.Rdata")


gam.SMOOTH <- gam(SMOOTH.formula, data = DAT, family="nb")
summary(gam.SMOOTH)
save(gam.SMOOTH, file="../STAN_models/gam_SMOOTH.Rdata")

