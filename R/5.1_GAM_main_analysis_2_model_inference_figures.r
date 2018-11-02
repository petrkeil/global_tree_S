################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Oct 25 2018
################################################################################

#
# Description: Here we produce figures related to model inference (parameters,
# region effects), observed-vs-predicted figures, 
# and also the correlation pairplot for all the predictors.
#

################################################################################ 
# LOAD THE PACKAGES, THE DATA, AND THE MODELS
################################################################################

source("0_libraries_functions_settings.r")

source("4.1_Data_loading_standardization_and_centering.r")

# load the models
load("../Models/gam_REALM.Rdata")
load("../Models/gam_SMOOTH.Rdata")
load("../Models/brms_REALM.RData")
load("../Models/brms_SMOOTH.RData")


################################################################################
# DISSECTING THE STAN MODEL RESULTS
################################################################################

# ------------------------------------------------------------------------------
# Observed vs predicted values with the full Bayesian uncertainty

pred.REALM.brm <- predict(brm.REALM, type="response", 
                          newdata = DAT, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pred.REALM.brm <- data.frame(pred.REALM.brm, 
                             S = DAT$S, 
                             grain = DAT$DAT_TYPE,
                             model = "Model REALM")

pred.SMOOTH.brm <- predict(brm.SMOOTH, type="response", 
                          newdata = DAT, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
pred.SMOOTH.brm <- data.frame(pred.SMOOTH.brm,  
                              S = DAT$S, 
                              grain = DAT$DAT_TYPE,
                              model = "Model SMOOTH")

pred.brm <- rbind(pred.REALM.brm, pred.SMOOTH.brm)


# observed vs predicted plots
obs.pred.brm <- ggplot(pred.brm, aes(x = S, y = Q50)) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5, colour = grain), alpha = 0.2) +
  geom_linerange(aes(ymin = Q25, ymax = Q75, colour = grain), size=1, alpha = 0.4) +
  geom_point(aes(colour=grain), shape = 1) +
  xlab("Observed S") + 
  ylab("Predicted S") +
  geom_abline(intercept = 0, slope = 1, colour="black") + theme_bw() +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_y_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000)) +
  facet_grid(.~model) +
  labs(colour = "grain") +
  guides(colour = guide_legend(override.aes = list(size=4, shape=19), title="grain"))

obs.pred.brm

# export the figure
png("../Figures/observed_vs_predicted.png", width=2000, height=950, res=250)
obs.pred.brm
dev.off()


# this is how to get the data out:
standata(brm.REALM)




# ------------------------------------------------------------------------------
# brm.SMOOTH results
# ------------------------------------------------------------------------------

all.varnames.SMOOTH <- rownames(summary(brm.SMOOTH)$fixed)
data.frame(all.varnames.SMOOTH)
data.frame(1:52, summary(brm.SMOOTH$fit)$summary[,'50%'])

all.vars.SMOOTH <- standata(brm.SMOOTH)$X

all.splines.SMOOTH <- standata(brm.SMOOTH)$Zs_2_1
pars.country.splines.SMOOTH <- summary(brm.SMOOTH$fit)$summary[,'50%'][26:38]
pars.plot.splines.SMOOTH <- summary(brm.SMOOTH$fit)$summary[,'50%'][39:51]
pars.envivars.SMOOTH <- summary(brm.SMOOTH$fit)$summary[,'50%'][5:22]

SMOOTH.area.id <- c(2:4)
area.varnames.SMOOTH <- all.varnames.SMOOTH[SMOOTH.area.id]
area.parnames.SMOOTH <- paste("b", area.varnames.SMOOTH, sep="_")
area.vars.SMOOTH <- model.matrix.gam(gam.SMOOTH)[,SMOOTH.area.id]


# ------------------------------------------------------------------------------
# brm.REALM results
# ------------------------------------------------------------------------------

all.varnames.REALM <- rownames(summary(brm.REALM)$fixed)
data.frame(all.varnames.REALM)
all.vars.REALM <- model.matrix.gam(gam.REALM)
all.pars.REALM <- summary(brm.REALM)$fixed[,'Estimate']

REALM.hist.id <- c(2:7, 17:37)
hist.varnames.REALM <- all.varnames.REALM[REALM.hist.id]
hist.parnames.REALM <- paste("b", hist.varnames.REALM, sep="_")
hist.vars.REALM <- model.matrix.gam(gam.REALM)[,REALM.hist.id]



# ------------------------------------------------------------------------------
# EXTRACT ONE PARAMETER FROM THE BRMS FITTED OBJECT (BY NAME)
extract.par <- function(par.name, brm.fit)
{
  n.iter <- summary(brm.fit)$iter
  n.chains <- summary(brm.fit)$chains
  n.thin <- summary(brm.fit)$thin
  n.warm <- summary(brm.fit)$warmup
  
  res <- list()
  
  for(chain in 1:n.chains){
    one.chain <- brm.fit$fit@sim$samples[[chain]][par.name][[1]]
    res[[chain]] <- one.chain[((n.warm/n.thin)+1):length(one.chain)]
  }
  
  res <- unlist(res)
  return(res)
}
# b <- extract.par("b_Intercept", brm.REALM)

# EXTRACT MULTIPLE PARAMETERS FROM THE BRMS FITTED OBJECT (BY NAMES)
extract.pars <- function(par.names, brm.fit)
{
  
  res <- list()
  for(par in par.names){
      res[[par]] <- extract.par(par, brm.fit)
  }
  return(res)
}
#bb <- extract.pars(hist.parnames.REALM, brm.REALM)

# CHOOSE ONE SAMPLING ITERATION AND SAMPLE PARAMETER VALUES FROM THAT
draw.one <- function(par.list, i){
  myfun <- function(x, i) x[i]
  unlist(lapply(par.list, FUN=myfun, i=i))
}
#bbb <- draw.one(bb, 1)

# CALCULATE PREDICTION FOR EACH PARAMETER DRAW IN THE CHAIN
draw.all <- function(par.names, vars, brm.fit, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
{
  B <- extract.pars(par.names, brm.fit)
  N <- length(B[[1]])
  
  res <- matrix(nrow=nrow(vars), ncol=N)
  
  for(i in 1:N){
    b <- draw.one(B, i=i)
    res[,i] <- vars %*% b
  }
  
  require(matrixStats)
  res <- matrixStats::rowQuantiles(res, probs=probs)
  return(res)
}

################################################################################
# SUMMARIZE THE REALM PREDICTIONS
################################################################################


prd.brm.REALM <- all.vars.REALM %*% all.pars.REALM
#prd.brm.REALM <- predict(brm.REALM, newdata=DAT, type="link")[,'Estimate']
res.brm.REALM <- prd.brm.REALM - log(DAT$S)

prd.history.REALM <- draw.all(par.names = hist.parnames.REALM, 
                        vars = hist.vars.REALM, 
                        brm.fit = brm.REALM)

resid.REALM <- prd.history.REALM[,'50%'] + res.brm.REALM
  
  

HIST <- data.frame(prd.history.REALM, DAT, resid.REALM)


p.hist.multi <- ggplot(HIST, aes(exp(Area_km*A.sd + A.mean), X50.)) +
  geom_point(aes(exp(Area_km*A.sd + A.mean), resid.REALM), 
             colour="grey", shape=1) +
  geom_linerange(aes(ymin=X2.5., ymax=X97.5., colour=REALM), alpha=0.5) +
  geom_linerange(aes(ymin=X25., ymax=X75., colour=REALM), size=1) +
  geom_line(aes(colour=REALM)) + 
  scale_y_continuous(minor_breaks = NULL) +
  scale_x_continuous(trans = "log10",
                     minor_breaks = NULL,
                     breaks = c(0.01, 1, 100, 10000, 1000000),
                     labels = c("0.01", "1", "100", expression(10^4), expression(10^6))) +
  xlab(expression("Area" ~ (km^2))) +
  ylab("Region effect") + 
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  theme_bw()  +  
  theme(legend.position="none") +
  facet_grid(. ~ REALM) 
p.hist.multi


theme.legend <- theme(legend.position="right",#c(0.25,0.8), 
                      legend.title = element_blank(),
                      legend.background = element_rect(fill="white",
                                                       size=0.2, linetype="solid", 
                                                       colour ="black"))


p.hist.single <- ggplot(HIST, aes(exp(Area_km*A.sd + A.mean), X50.)) +
  geom_linerange(aes(ymin=X2.5., ymax=X97.5., colour=REALM), alpha=0.4) +
  geom_linerange(aes(ymin=X25., ymax=X75., colour=REALM), size=1) +
  geom_line( aes(colour=REALM)) +
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  scale_y_continuous(minor_breaks = NULL,
                     breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6 )) +
  scale_x_continuous(trans = "log10",
                     minor_breaks = NULL,
                     breaks = c(0.01, 1, 100, 10000, 1000000),
                     labels = c("0.01", "1", "100", expression(10^4), expression(10^6))) +
  xlab(expression("Area" ~ (km^2))) +
  ylab("Region effect") + 
  theme_bw()+
  theme.legend +
  theme(plot.title = element_text(hjust = -0.125))
p.hist.single


# Figure for the main text - save to a file
pdf("../Figures/region_effects_ovelaid_REALM.pdf", width=6, height=4)
p.hist.single
dev.off()

  
# Figure for supplementary material - save to a file
tiff("../Figures/region_effects_individually_REALM.tif", 
     width=2500, height=600, res=250,
     compression = "lzw")
p.hist.multi
dev.off()


################################################################################
#  TRIPHASIC SAR CURVE FROM THE SMOOTH PREDICTIONS
################################################################################

prd.area.SMOOTH <- draw.all(par.names = area.parnames.SMOOTH, 
                            vars = area.vars.SMOOTH, 
                            brm.fit = brm.SMOOTH)


partial.res.area <- prd.area.SMOOTH + residuals(gam.SMOOTH) 


AREA <- data.frame(prd.area.SMOOTH, DAT, partial.resid = partial.res.area)

# save to a file
png("../Figures/triphasic_SAR_SMOOTH.png", width=1500, height=1500, res=500)
ggplot(AREA, aes(x=exp(Area_km*A.sd + A.mean), y=X50.)) + 
       geom_point(aes(x=exp(Area_km*A.sd + A.mean), y=partial.resid.50.), 
                  color="grey", shape=1) +
       geom_line(size=0.8) +
  scale_x_continuous(trans = "log10",
                     minor_breaks = NULL,
                     breaks = c(0.01, 1, 100, 10000, 1000000),
                     labels = c("0.01", "1", expression(10^2), expression(10^4), expression(10^6))) +
       ylab("Partial effect of Area") +
       xlab(expression("Area" ~ (km^2))) +
       geom_ribbon(aes(ymin=X2.5., ymax=X97.5.), alpha=0.3) +
       theme_bw()
dev.off()




################################################################################
#  GRAIN-DEPENDENT ENVIRONMENTAL PREDICTORS
################################################################################

# indices of variable coefficients - model REALM
data.frame(all.varnames.REALM)
pure.id.REALM <- 8:16
area.id.REALM <- 38:46

# coefficient names
pure.varnames.REALM <- all.varnames.REALM[pure.id.REALM]
pure.coefnames.REALM <- paste("b", pure.varnames.REALM, sep="_")
area.varnames.REALM <- all.varnames.REALM[area.id.REALM]
area.coefnames.REALM <- paste("b", area.varnames.REALM, sep="_")

pure.coefs.REALM <- extract.pars(pure.coefnames.REALM, brm.REALM)
area.coefs.REALM <- extract.pars(area.coefnames.REALM, brm.REALM)

# -----------------
# indices of variable coefficients - model SMOOTH
data.frame(all.varnames.SMOOTH)
pure.id.SMOOTH <- 5:13
area.id.SMOOTH <- 14:22

# coefficient names
pure.varnames.SMOOTH <- all.varnames.SMOOTH[pure.id.SMOOTH]
pure.coefnames.SMOOTH <- paste("b", pure.varnames.SMOOTH, sep="_")
area.varnames.SMOOTH <- all.varnames.SMOOTH[area.id.SMOOTH]
area.coefnames.SMOOTH <- paste("b", area.varnames.SMOOTH, sep="_")

pure.coefs.SMOOTH <- extract.pars(pure.coefnames.SMOOTH, brm.SMOOTH)
area.coefs.SMOOTH <- extract.pars(area.coefnames.SMOOTH, brm.SMOOTH)

# -----------------
# function that extracts posterior quantiles of a parameter from a brms model object
get.one.var <- function(i, pure.coefs, area.coefs, varnames, model,
                         probs=c(0.025, 0.25, 0.5,0.75, 0.975))
{
  pure.coef <- pure.coefs[[i]]
  area.coef <- area.coefs[[i]]
  A <- seq(from=-0.94, to=2.83, by=0.1)
  N <- length(pure.coef)
  
  res <- matrix(nrow=length(A), ncol=N)
  
  for(j in 1:N)
  {
    b.pure <- pure.coef[j]
    b.area <- area.coef[j]
    res[,j] <- b.pure + b.area*A
  }
  
  require(matrixStats)
  res <- matrixStats::rowQuantiles(res, probs=probs)
  res <- data.frame(Model=model,
                    variable=rep(varnames[i], times=length(A)),
                    A, res)
  return(res)
}

# ------------------
# function that extracts posterior quantiles of 
# multiple parameters from a brms model object
get.vars <- function(pure.coefs, area.coefs, varnames, model,
                     probs=c(0.025, 0.25, 0.5,0.75, 0.975))
{
  res <- list()
  for(i in 1:length(varnames))
  {
    res[[i]] <- get.one.var(i=i, 
                            pure.coefs=pure.coefs, 
                            area.coefs=area.coefs, 
                            varnames= varnames,
                            model=model)
  }

  require(plyr)
  res <- plyr::ldply(res)
  return(res)
}

# ------------------------------------------------------------------------------
# get the quantiles of the posterior distributions using the functions above

REALM.scale.coefs <- get.vars(pure.coefs=pure.coefs.REALM, 
                              area.coefs=area.coefs.REALM, 
                              varnames= pure.varnames.REALM,
                              model="REALM")

SMOOTH.scale.coefs <- get.vars(pure.coefs=pure.coefs.SMOOTH, 
                               area.coefs=area.coefs.SMOOTH, 
                               varnames= pure.varnames.SMOOTH,
                               model="SMOOTH")

scale.coefs <- rbind(REALM.scale.coefs, SMOOTH.scale.coefs)
scale.coefs$variable <- as.character(scale.coefs$variable)

# ------------------------------------------------------------------------------
# use better names for scale coefficients

scale.coefs[scale.coefs$variable=="Tree_dens",'variable'] <- "Tree density"
scale.coefs[scale.coefs$variable=="min_DBH",'variable'] <- "Minimum DBH"
scale.coefs[scale.coefs$variable=="ANN_T",'variable'] <- "Annual T"
scale.coefs[scale.coefs$variable=="ISO_T",'variable'] <- "Isothermality"
scale.coefs[scale.coefs$variable=="MIN_P",'variable'] <- "Minimum P"
scale.coefs[scale.coefs$variable=="P_SEAS",'variable'] <- "P Seasonality"
scale.coefs[scale.coefs$variable=="ALT_DIF",'variable'] <- "Elevation span"
scale.coefs[scale.coefs$variable=="ISLANDmainland",'variable'] <- "Mainland"
# scale.coefs[scale.coefs$variable=="ELONG",'variable'] <- "Elongation"

scale.coefs$variable <- factor(scale.coefs$variable, 
                               levels=c("Mainland",
                                        "Elevation span", "GPP", 
                                        "P Seasonality", "Minimum P",
                                        "Annual T", "Isothermality", 
                                        "Tree density", "Minimum DBH") )

# ------------------------------------------------------------------------------
# plot the coefficients


coef.plot <- ggplot(scale.coefs, aes(x=exp(A*A.sd + A.mean), y=X50.)) + 
             geom_ribbon(aes(ymin=X2.5., ymax=X97.5., fill=Model), alpha=0.3) +
             #geom_ribbon(aes(ymin=X25., ymax=X75., fill=NA), alpha=0.3) +
             geom_line(aes(colour=Model), size=1) +
             geom_line(aes(x=exp(A*A.sd + A.mean), y=X25., colour=Model), 
                       linetype="dashed") +
             geom_line(aes(x=exp(A*A.sd + A.mean), y=X75., colour=Model), 
                       linetype="dashed") +
             facet_grid(.~variable) +
             geom_hline(yintercept=0, colour="darkgrey") +
             xlab(expression("Area" ~ (km^2))) + 
             ylab("Environment effect") +
             scale_colour_brewer(palette = "Set1") +
             scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(trans = "log10",
                     minor_breaks = NULL,
                     breaks = c(0.001, 1, 1000, 1000000),
                     labels = c(expression(10^-3), "1", expression(10^3), expression(10^6))) +
             theme_bw() +
             theme(legend.position="right")
coef.plot

# save to a file
png("../Figures/environment_effects.png", width=2200, height=600, res=200)
coef.plot
dev.off()

# save to a file
pdf("../Figures/environment_effects.pdf", width=11, height=3)
coef.plot
dev.off()




################################################################################
#  PAIRPLOTS OF THE PREDICTORS
################################################################################


for.pairs <- dplyr::select(DAT, S, Tree_dens, min_DBH, GPP, ANN_T, ISO_T, 
                           MIN_P, P_SEAS, ALT_DIF, ISLAND, REALM, DAT_TYPE)

p <- ggpairs(for.pairs, aes(colour=DAT_TYPE)) + theme_bw() +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))

# save to a file
png("../Figures/predictors_pairplot.png", width=2000, height=2300, res=100)
ggpairs(for.pairs, aes(colour=DAT_TYPE)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



################################################################################
# CORRELATION MATRIX OF MODEL PARAMETERS
################################################################################

library(corrplot)

one.chain = data.frame(as.mcmc(brm.REALM)[1][[1]])
data.frame(names(one.chain))
one.chain = one.chain[, c(1, 2:7, 17:37, 8:16, 38:46)]

par.cors <- cor(one.chain)

png("../Figures/parameter_correlation_matrix.png", width=2000, height=2000, res = 150)
corrplot(par.cors, method = "square", tl.col = "black")
dev.off()


