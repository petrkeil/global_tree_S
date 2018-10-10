################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Sep 11 2018
################################################################################

#
# Description: 
#

################################################################################ 
# LOAD THE DATA AND THE PACKAGES
################################################################################

source("0_libraries_functions_settings.r")
source("4_Data_loading_standardization_and_centering.r")



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


# CROSSVALIDATING MODEL SMOOTH

preds <- list()
par(mfrow=c(2,2))
for (fold in 1:nfolds) 
{
  print(fold)
  
  train <- DAT[case.folds!=fold,]
  test <- DAT[case.folds==fold,]
  
  gam.train <- gam(SMOOTH.formula, data=train, family="nb")
  res <- predict.gam(gam.train, newdata = test, type="response")
  res <- data.frame(predicted = res,
                    observed = test$S,
                    grain = test$DAT_TYPE, 
                    fold = paste("Fold", fold))
  preds[[fold]] <- res
}

preds <- ldply(preds)

SMOOTH.cross <- ggplot(data = preds, aes(x=observed, y=predicted)) +
                geom_point(aes(colour = grain), shape = 1) +
                scale_x_continuous(trans = "log10") + 
                scale_y_continuous(trans = "log10") +
                facet_grid(.~fold) +
                geom_abline(intercept = 0, slope = 1) +
                labs(x = "Observed test S", 
                     y = "Predicted test S",
                     title = "Model SMOOTH") +
                theme_bw()
SMOOTH.cross

# CROSSVALIDATING MODEL REALM
preds <- list()
par(mfrow=c(2,2))
for (fold in 1:nfolds) 
{
  print(fold)
  
  train <- DAT[case.folds!=fold,]
  test <- DAT[case.folds==fold,]
  
  gam.train <- gam(REALM.formula, data=train, family="nb")
  res <- predict.gam(gam.train, newdata = test, type="response")
  res <- data.frame(predicted = res,
                    observed = test$S,
                    grain = test$DAT_TYPE, 
                    fold = paste("Fold", fold))
  preds[[fold]] <- res
}

preds <- ldply(preds)

REALM.cross <- ggplot(data = preds, aes(x=observed, y=predicted)) +
               geom_point(aes(colour = grain), shape = 1) +
               scale_x_continuous(trans = "log10") + 
               scale_y_continuous(trans = "log10") +
               facet_grid(.~fold) +
               geom_abline(intercept = 0, slope = 1) +
               labs(x = "Observed test S", y = "Predicted test S",
                    title = "Model REALM") +
               theme_bw()
REALM.cross


png("../Figures/crossvalidation.png", width=2500, height=1600, res=250)
  grid.arrange(REALM.cross, SMOOTH.cross, ncol=1, nrow=2)
dev.off()
