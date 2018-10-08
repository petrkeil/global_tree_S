

# clean the workspace and load the libraries and the data
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
  ISLAND + ISLAND:Area_km +
  #INSULARITY + INSULARITY:Area_km +
  ELONGATION + ELONGATION:Area_km


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
  ISLAND + ISLAND:Area_km +
  INSULARITY + INSULARITY:Area_km +
  #ELONGATION + ELONGATION:Area_km


# ------------------------------------------------------------------------------


# crossvalidation function
crossvalidate <- function(model.formula, DAT, nfolds = 10)
{

  case.folds <- rep(1:nfolds,length.out=nrow(DAT))
  
  # divide the cases as evenly as possible
  case.folds <- sample(case.folds) # randomly permute the order
  
  preds <- list()
  for (fold in 1:nfolds) 
  {
   # print(paste("Fold", fold))
    
    train <- DAT[case.folds!=fold,]
    test <- DAT[case.folds==fold,]
    
    gam.train <- gam(model.formula, data=train, family="nb")
    res <- predict.gam(gam.train, newdata = test , type="response")
    res <- data.frame(predicted = log10(res),
                      observed = log10(test$S),
                      grain = test$DAT_TYPE, 
                      fold = paste("Fold", fold))
    preds[[fold]] <- res
  }
  
  preds <- ldply(preds)
 # plot(preds$observed, preds$predicted)
  
  null.sum.sq <- sum((res$observed - mean(res$observed))^2)
  
  resid.sum.sq <- sum((res$observed - res$predicted)^2)
  
  r2s <- 1 - resid.sum.sq/null.sum.sq
  return(r2s)
}

# performing the crossvalidation

model.formula = REALM.formula
model.formula
m1 <- gam(model.formula, data=DAT, family="nb")
summary(m1)
AIC(m1)
BIC(m1)
crossvalidate(model.formula, DAT)

