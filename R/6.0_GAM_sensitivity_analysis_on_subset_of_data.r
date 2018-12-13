################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: April 26 2018
################################################################################


# Description: THIS IS WHERE WE DO THE ANALYSIS WITH JUST A SUBSET OF THE DATA, 
# OR THE SO CALLED 'SENSITIVITY ANALYSIS' (SEE THE METHODS).
# THE GOAL IS TO MAKE SURE THAT THE RESULTS ARE ROBUST TO DATA SOURCES
# AND TO DIFFERENT DEFINITIONS OF TREES.
# THE PROTOCOL IS THE SAME AS IN THE FULL DATASET, ONLY THE BASE DATA ARE SMALLER.



################################################################################ 
# LOAD THE DATA AND THE PACKAGES
################################################################################


# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")

# ------------------------------------------------------------------------------

# load the data
PLT <- read.csv("../Data/Main_dataset_subset.csv")


################################################################################ 
# PREPARE THE DATA FOR THE ANALYSES
################################################################################

# cakculate tree density (note the x+1 step!!)
PLT$Tree_dens <- (PLT$N + 1) / PLT$Area_km 

# select only the variables of interest from the larger data.frame
DAT <- dplyr::select(PLT, S, Area_km, Tree_dens, min_DBH=min_DBH_cm, 
                     GPP, ET, ANN_T, WARM_T, ISO_T, MIN_P, P_SEAS, ALT_DIF,
                     ISLAND = ISL_LS, REALM=REALM_PK, Lat, Lon, DAT_TYPE, Loc_ID) 

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
# THE MODEL FORMULAS
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


################################################################################ 
# FIT THE MODELS IN mgcv
################################################################################

gam.REALM <- gam(REALM.formula, data=DAT, family="nb")
summary(gam.REALM)

gam.SMOOTH <- gam(SMOOTH.formula, data = DAT, family="nb")
summary(gam.SMOOTH)

# additionally, fit model with only the intercept (a 'null' model)
gam.NULL <- gam(S~1, data=DAT, family="nb")

# compare the models using AIC
AIC(gam.NULL, gam.REALM, gam.SMOOTH)



################################################################################
# FIT THE MODELS IN STAN, using 'brms' function 'brm'
################################################################################

# Beware, this will take several hours to run.
# YOU CAN SKIP THIS AND LAOD THE FITTED MODELS IN THE NEXT STEP (DEFAULT BEHAVIOR)

fit.brm = FALSE # Should we skip this step? TRUE/FALSE
if(fit.brm){
  brm.SMOOTH <- brm(SMOOTH.formula, family="negbinomial", data=DAT,
                    cores=3,
                    seed=12355,
                    chains=3, iter=3000, warmup=1000, thin=10)
  
  save(brm.SMOOTH , file="../Models/subset_brms_SMOOTH.RData")
  
  brm.REALM <- brm(REALM.formula, family="negbinomial", data=DAT,
                   cores=3,
                   seed=12355,
                   chains=3, iter=3000, warmup=1000, thin=10)
  
  save(brm.REALM, file="../Models/subset_brms_REALM.RData")
}


################################################################################
# LOAD THE FITTED STAN OBJECTS
################################################################################

load("../Models/subset_brms_SMOOTH.RData")
load("../Models/subset_brms_REALM.RData")


################################################################################ 
# PLOT OBSERVED VS. PREDICTED VALUES
################################################################################

pred.REALM <- data.frame( DAT, pred=predict(gam.REALM, type="response"),
                          model="Model REALM")
pred.SMOOTH <- data.frame( DAT, pred=predict(gam.SMOOTH, type="response"),
                           model="Model SMOOTH")
pred <- rbind(pred.REALM, pred.SMOOTH)

obs.glm <- ggplot(pred, aes(log10(S), log10(pred))) +
  geom_point(aes(colour=REALM), shape=1) +
  xlab("log10 Observed S") + ylab("log10 Predicted S") +
  geom_abline(intercept = 0, slope = 1, colour="black") + theme_bw() +
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  scale_x_continuous(limits = c(0, 4)) + scale_y_continuous(limits = c(0, 4)) +
  facet_grid(.~model) +
  labs(size = "log10 Area [km]") + labs(colour = "Realm") +
  guides(colour = guide_legend(override.aes = list(size=4, shape=19)))



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
# export the figure
png("../Figures/Subset_data_sensitivity_analysis/subset_observed_vs_predicted.png", 
    width=2000, height=950, res=250)
print(obs.pred.brm)
dev.off()





################################################################################ 
# PLOT MAPS
################################################################################

# extract the smooth GAM surfaces representing history
prd <- predict.gam(gam.SMOOTH, type="terms")
smooth.country <- prd[,"s(Lat,Lon):DAT_TYPECountry"]
smooth.plot <- prd[,"s(Lat,Lon):DAT_TYPEPlot"]
PRED <- data.frame(DAT, smooth.country, smooth.plot)


# extract predicted species richness
S.pred <- predict.gam(gam.SMOOTH, type="response")
PRED <- data.frame(PRED, S.pred)


# reproject geogrpahic coordinates of centroids to a Mollweide projection
latlon <- SpatialPoints(data.frame(PRED$Lon, PRED$Lat), CRS(WGS84))
xy <- data.frame(coordinates(spTransform(latlon, CRSobj = CRS(MOLLWEIDE))))
names(xy) <- c("X", "Y")
PRED <- data.frame(PRED, xy)


# select prediction data for countries
PRED.CNTRS <- PRED[PRED$DAT_TYPE=="Country",]


# global mainlands (not divided by country boundaries)
MAINL <- readOGR(dsn = "../Data/COUNTRIES", layer = "GSHHS_i_L1_simple")
MAINL <- spTransform(MAINL, CRSobj = CRS(MOLLWEIDE))
MAINL <- tidy(MAINL, region="id")

# equator, tropics, and polar circles
LINES <- readOGR(dsn = "../Data/COUNTRIES", layer = "ne_110m_geographic_lines")
LINES <- spTransform(LINES, CRSobj = CRS(MOLLWEIDE))
LINES <- tidy(LINES, region="name")

# load the country shapefile
CNTRS <- readOGR(dsn = "../Data/COUNTRIES", layer = "COUNTRIES_with_environment")
CNTRS <- CNTRS[CNTRS$Aggregated == "Original country",]
CNTRS@data <- dplyr::select(CNTRS@data, NAME)


CNTRS <- spTransform(CNTRS, CRSobj = CRS(MOLLWEIDE))

PRED.CNTRS <- merge(CNTRS, PRED.CNTRS, by.x="NAME", by.y="Loc_ID")
#PRED.CNTRS <- PRED.CNTRS[is.na(PRED.CNTRS$DAT_TYPE) == FALSE, ]



# DATA FOR PLOTS ---------------------------------------------------------------

PRED.PLOTS <- PRED[PRED$DAT_TYPE=="Plot",]

# prepare the polygons for GGPLOT2 ---------------
C.fort <- tidy(PRED.CNTRS, region="NAME")
C.fort <- merge(C.fort, PRED.CNTRS@data, by.x="id", by.y="NAME")


# BLANK THEME FOR GGPLOT2

blank.theme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())

blank.theme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.position=c(0.6, 0.05),
                     legend.direction = "horizontal",
                     legend.title = element_blank(),
                     #legend.title.align = 0,
                     panel.background=element_blank(),
                     panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())




# ------------------------------------------------------------------------------
#### SMOOTHER MAPS


g.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=smooth.country), colour="black", size=.2) + 
  scale_fill_viridis(option = "viridis", 
                       limits=c(-1, 0.8),
                       name="Region\neffect") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") + 
  ggtitle("Smooth region effects at country grain") +  
  theme_minimal() + blank.theme


g.plot <- ggplot(MAINL, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="white", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=smooth.plot)) +
  scale_colour_viridis(option = "viridis", 
                         limits=c(-1, 0.8),
                         name="Region\neffect") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") +
  ggtitle("Smooth region effects at plot grain") +
  theme_minimal() + blank.theme

png("../Figures/Subset_data_sensitivity_analysis/subset_region_effect_map_SMOOTH.png", width=2000, height=2000, res=250)
grid.arrange(g.cntr, g.plot, ncol=1, nrow=2)
dev.off()


# ------------------------------------------------------------------------------
# PREDICTED RICHNESS MAPS

s.pred.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=S.pred), colour="black", size=.2) + 
  scale_fill_viridis(option = "magma", name="Predicted S", 
                       trans="log10") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  ggtitle("Predicted richness at country grain (model SMOOTH)") + theme_minimal() +
  xlab("") + ylab("") + blank.theme

s.pred.plot <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="white", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=S.pred)) +
  scale_colour_viridis(option = "magma", name="Predicted S", 
                         trans="log10") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  ggtitle("Predicted richness at plot grain (model SMOOTH)") + theme_minimal() +
  xlab("") + ylab("") + blank.theme

png("../Figures/Subset_data_sensitivity_analysis/subset_predicted_richness_SMOOTH.png", 
    width=2000, height=2000, res=250)
grid.arrange(s.pred.cntr, s.pred.plot, ncol=1, nrow=2)
dev.off()

# ------------------------------------------------------------------------------
# RAW RICHNESS MAPS

s.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=S), colour="black", size=.2) + 
  scale_fill_viridis(option = "magma", name="S", 
                       trans="log10", limits=c(1,10000)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  ggtitle("Observed species richness - country grain") + theme_minimal() +
  # labs(subtitle = expression(S[country] ~ "(richness at the country grain)")) +
  xlab("") + ylab("") + blank.theme

s.plot <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="lightgrey", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=S), size=1) +
  #geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), size=1, colour="black", 
  #           shape =1 , size=.2) +
  scale_colour_viridis(option = "magma", name="S", 
                         trans="log10", limits=c(1,10000)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  ggtitle("Observed species richness - plot grain") + theme_minimal() +
  #labs(subtitle = expression(S[plot] ~ "(richness at the plot grain)")) +
  xlab("") + ylab("") + blank.theme

png("../Figures/Subset_data_sensitivity_analysis/subset_observed_richness_map.png", width=2000, height=2000, res=250)
grid.arrange(s.cntr, s.plot, ncol=1, nrow=2)
dev.off()


# ------------------------------------------------------------------------------
# BIOGEOGRAPHIC REALMS

blank.theme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())


png("../Figures/Subset_data_sensitivity_analysis/subset_realms_map.png", width=1800, height=1000, res=200)
realm.plot <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(aes(fill=REALM), colour="darkgray", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), shape=3, colour="black") +
  scale_fill_brewer(palette = "Dark2", name="Realm") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  ggtitle("Biogeographic realms") +
  theme(legend.position="none") 
realm.plot
dev.off()





################################################################################
# DISSECTING THE STAN MODEL RESULTS
################################################################################

# are the predictions from the brms and mgcv models the same?
prd.gam <- predict(gam.SMOOTH, type="response")
prd <- predict(brm.SMOOTH, newdata=DAT)[,'Estimate']
plot(prd, prd.gam)
abline(a=0, b=1)


# ------------------------------------------------------------------------------

# TRACEPLOTS and CATERPILLAR PLOTS of model parameters

pars.REALM <- rownames(data.frame(summary(brm.REALM$fit)[1]))[1:46]
pars.SMOOTH <- rownames(data.frame(summary(brm.SMOOTH$fit)[1]))[1:51]

# traceplots - this indicate if the convergence is good (which it is)

rstan::traceplot(brm.REALM$fit, pars=pars.REALM)
rstan::traceplot(brm.SMOOTH$fit, pars=pars.SMOOTH)


# caterpillar plots
rstan::plot(brm.REALM$fit, pars=pars.REALM)


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
png("../Figures/Subset_data_sensitivity_analysis/subset_region_effects_ovelaid_REALM.png", 
     width=1200, height=800, res=200)
p.hist.single
dev.off()



################################################################################
# TRIPHASIC SAR CURVE FROM THE SMOOTH PREDICTIONS
################################################################################

prd.area.SMOOTH <- draw.all(par.names = area.parnames.SMOOTH, 
                            vars = area.vars.SMOOTH, 
                            brm.fit = brm.SMOOTH)


partial.res.area <- prd.area.SMOOTH + residuals(gam.SMOOTH) 


AREA <- data.frame(prd.area.SMOOTH, DAT, partial.resid = partial.res.area)

png("../Figures/Subset_data_sensitivity_analysis/subset_partial_area_effect_SMOOTH.png", width=1500, height=1500, res=500)
ggplot(AREA, aes(x=log10(exp(Area_km*A.sd + A.mean)), y=X50.)) + 
  geom_point(aes(x=log10(exp(Area_km*A.sd + A.mean)), y=partial.resid.50.), 
             color="grey", shape=1) +
  geom_line(size=0.8) +
  ylab("Partial effect of area") +
  xlab("log10 Area [km^2]") +
  geom_ribbon(aes(ymin=X2.5., ymax=X97.5.), alpha=0.3) +
  theme_bw()
dev.off()


################################################################################
# GRAIN-DEPENDENT ENVIRONMENTAL PREDICTORS
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
png("../Figures/Subset_data_sensitivity_analysis/subset_environment_effects.png", width=2200, height=600, res=200)
coef.plot
dev.off()


