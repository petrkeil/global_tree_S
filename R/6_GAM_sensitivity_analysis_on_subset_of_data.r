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
# 0. LOAD THE DATA AND THE PACKAGES
################################################################################


# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")

# ------------------------------------------------------------------------------

# load the data
PLT <- read.csv("../Data/Main_dataset_subset.csv")


################################################################################ 
# 1. PREPARE THE DATA FOR THE ANALYSES
################################################################################

# cakculate tree density (note the x+1 step!!)
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



################################################################################ 
# 3. FIT THE MODELS
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
# 4. PLOT OBSERVED VS. PREDICTED VALUES
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

# export the figure
png("../Figures/Subset_data_sensitivity_analysis/subset_obs_vs_pred.png", 
    width=2000, height=950, res=250)
    print(obs.glm)
dev.off()



################################################################################ 
# 5. MORAN'S I CORRELOGRAMS OF THE RESIDUALS for both models at two scales
################################################################################

# width of the correlogram increment [km]
increment = 200

# Residual correlograms for Model REALM -----------------
res.REALM <- data.frame(Lat=DAT$Lat, Lon=DAT$Lon, S=DAT$S, 
                        DAT_TYPE=DAT$DAT_TYPE, 
                        resid=residuals(gam.REALM))

res.REALM.cntr <- res.REALM[res.REALM$DAT_TYPE=="Country",]
res.REALM.plot <- res.REALM[res.REALM$DAT_TYPE=="Plot",]

cor.REALM.cntr <- ncf::correlog(x=res.REALM.cntr$Lat, y=res.REALM.cntr$Lon, 
                                z=res.REALM.cntr$resid, 
                                latlon=TRUE, resamp=1, increment=increment)
cor.REALM.plot <- ncf::correlog(x=res.REALM.plot$Lat, y=res.REALM.plot$Lon, 
                                z=res.REALM.plot$resid, 
                                latlon=TRUE, resamp=1, increment=increment)


# Residual correlograms for Model SMOOTH -----------------
res.SMOOTH <- data.frame(Lat=DAT$Lat, Lon=DAT$Lon, S=DAT$S, 
                         DAT_TYPE=DAT$DAT_TYPE, 
                         resid=residuals(gam.SMOOTH))

res.SMOOTH.cntr <- res.SMOOTH[res.SMOOTH$DAT_TYPE=="Country",]
res.SMOOTH.plot <- res.SMOOTH[res.SMOOTH$DAT_TYPE=="Plot",]

cor.SMOOTH.cntr <- ncf::correlog(x=res.SMOOTH.cntr$Lat, y=res.SMOOTH.cntr$Lon, 
                                 z=res.SMOOTH.cntr$resid, 
                                 latlon=TRUE, resamp=1, increment=increment)
cor.SMOOTH.plot <- ncf::correlog(x=res.SMOOTH.plot$Lat, y=res.SMOOTH.plot$Lon, 
                                 z=res.SMOOTH.plot$resid, 
                                 latlon=TRUE, resamp=1, increment=increment)

# Correlograms for S
cor.S.cntr <- ncf::correlog(x=res.SMOOTH.cntr$Lat, y=res.SMOOTH.cntr$Lon, 
                            z=res.SMOOTH.cntr$S, 
                            latlon=TRUE, resamp=1, increment=increment)
cor.S.plot <- ncf::correlog(x=res.SMOOTH.plot$Lat, y=res.SMOOTH.plot$Lon, 
                            z=res.SMOOTH.plot$S, 
                            latlon=TRUE, resamp=1, increment=increment)


# extract the correlograms for further plotting
N <- length(cor.REALM.plot$mean.of.class)

RP <- data.frame(Dist=cor.REALM.plot$mean.of.class, 
                 Cor=cor.REALM.plot$correlation,
                 Scale=rep("Plot"), 
                 Variable=rep("REALM residuals"))
RC <- data.frame(Dist=cor.REALM.cntr$mean.of.class, 
                 Cor=cor.REALM.cntr$correlation,
                 Scale=rep("Country"), 
                 Variable=rep("REALM residuals"))
SP <- data.frame(Dist=cor.SMOOTH.plot$mean.of.class, 
                 Cor=cor.SMOOTH.plot$correlation,
                 Scale=rep("Plot"), 
                 Variable=rep("SMOOTH residuals"))
SC <- data.frame(Dist=cor.SMOOTH.cntr$mean.of.class, 
                 Cor=cor.SMOOTH.cntr$correlation,
                 Scale=rep("Country"), 
                 Variable=rep("SMOOTH residuals"))
RIC <- data.frame(Dist=cor.S.cntr$mean.of.class, 
                  Cor=cor.S.cntr$correlation,
                  Scale=rep("Country"), 
                  Variable=rep("Species richness S"))
RIP <- data.frame(Dist=cor.S.plot$mean.of.class, 
                  Cor=cor.S.plot$correlation,
                  Scale=rep("Plot"), 
                  Variable=rep("Species richness S"))

cor.data <- rbind(RP, RC, SP, SC, RIP, RIC)


# plot the correlograms
png("../Figures/Subset_data_sensitivity_analysis/subset_correlograms.png", width=2000, height=900, res=250)
ggplot(cor.data, aes(x=Dist, y=Cor)) + 
  geom_point(aes(colour=Variable, shape=Variable)) + 
  geom_line(aes(colour=Variable)) + 
  geom_hline(yintercept = 0, colour="darkgrey") +
  xlim(c(0, 3000)) + ylim(c(-0.2,0.5)) +
  scale_colour_brewer(palette="Set1") +
  facet_grid(.~Scale) +
  xlab("Distance [km]") + ylab("Moran's I") +
  theme_bw()
dev.off()



################################################################################ 
# 6. PLOT MAPS
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
  scale_fill_distiller(palette = "Spectral", 
                       limits=c(-1, 0.8),
                       name="Region\neffect") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") + 
  ggtitle("A") +  
  theme_minimal() + blank.theme


g.plot <- ggplot(MAINL, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="white", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=smooth.plot)) +
  scale_colour_distiller(palette = "Spectral", 
                         limits=c(-1, 0.8),
                         name="Region\neffect") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") +
  ggtitle("B") +
  theme_minimal() + blank.theme

png("../Figures/Subset_data_sensitivity_analysis/subset_smooth_map.png", width=2000, height=2000, res=250)
grid.arrange(g.cntr, g.plot, ncol=1, nrow=2)
dev.off()


# ------------------------------------------------------------------------------
# PREDICTED RICHNESS MAPS

s.pred.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=S.pred), colour="black", size=.2) + 
  scale_fill_distiller(palette = "Spectral", name="Predicted S", 
                       trans="log10") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  ggtitle("A") + theme_minimal() +
  xlab("") + ylab("") + blank.theme

s.pred.plot <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="white", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=S.pred)) +
  scale_colour_distiller(palette = "Spectral", name="Predicted S", 
                         trans="log10") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  ggtitle("B") + theme_minimal() +
  xlab("") + ylab("") + blank.theme

png("../Figures/Subset_data_sensitivity_analysis/subset_predictions.png", 
    width=2000, height=2000, res=250)
grid.arrange(s.pred.cntr, s.pred.plot, ncol=1, nrow=2)
dev.off()

# ------------------------------------------------------------------------------
# RAW RICHNESS MAPS

s.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=S), colour="black", size=.2) + 
  scale_fill_distiller(palette = "Spectral", name="S", 
                       trans="log10", limits=c(1,10000)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  ggtitle("A") + theme_minimal() +
  # labs(subtitle = expression(S[country] ~ "(richness at the country grain)")) +
  xlab("") + ylab("") + blank.theme

s.plot <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="lightgrey", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=S), size=1) +
  #geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), size=1, colour="black", 
  #           shape =1 , size=.2) +
  scale_colour_distiller(palette = "Spectral", name="S", 
                         trans="log10", limits=c(1,10000)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  ggtitle("B") + theme_minimal() +
  #labs(subtitle = expression(S[plot] ~ "(richness at the plot grain)")) +
  xlab("") + ylab("") + blank.theme

png("../Figures/Subset_data_sensitivity_analysis/subset_richness_map.png", width=2000, height=2000, res=250)
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


png("../Figures/Subset_data_sensitivity_analysis/realms_map.png", width=2300, height=1000, res=200)
realm.plot <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(aes(fill=REALM), colour="darkgray", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), shape=3, colour="black") +
  scale_fill_brewer(palette = "Dark2", name="Realm") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  theme(legend.position="none") 
realm.plot
dev.off()



################################################################################
# 7. FIT THE MODELS IN STAN, using 'brms' function 'brm'
################################################################################

# Beware, this will take several hours to run.
# YOU CAN SKIP THIS AND LAOD THE FITTED MODELS IN STEP 8!!!

brm.SMOOTH <- brm(SMOOTH.formula, family="negbinomial", data=DAT,
                  cores=3,
                  seed=12355,
                  chains=3, iter=3000, warmup=1000, thin=10)

save(brm.SMOOTH , file="../STAN_models/subset_brms_SMOOTH.RData")

brm.REALM <- brm(REALM.formula, family="negbinomial", data=DAT,
                 cores=3,
                 seed=12355,
                 chains=3, iter=3000, warmup=1000, thin=10)

save(brm.REALM, file="../STAN_models/subset_brms_REALM.RData")



################################################################################
# 8. LOAD THE FITTED STAN OBJECTS
################################################################################

load("../STAN_models/subset_brms_SMOOTH.RData")
load("../STAN_models/subset_brms_REALM.RData")



################################################################################
# 9. DISSECTING THE STAN MODEL RESULTS
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

png("../Figures/Subset_data_sensitivity_analysis/subset_traceplot_REALM.png", width = 2000, height=2000, res=150)
rstan::traceplot(brm.REALM$fit, pars=pars.REALM)
dev.off()

png("../Figures/Subset_data_sensitivity_analysis/subset_traceplot_SMOOTH.png", width = 2000, height=2000, res=150)
rstan::traceplot(brm.SMOOTH$fit, pars=pars.SMOOTH)
dev.off()

# caterpillar plots

png("../Figures/caterpillar_REALM.png", width = 2000, height=2000, res=150)
rstan::plot(brm.REALM$fit, pars=pars.REALM)
dev.off()

# this is how to get the data out:
standata(brm.REALM)

# ------------------------------------------------------------------------------
# brm.SMOOTH results
# ------------------------------------------------------------------------------

all.varnames.SMOOTH <- rownames(summary(brm.SMOOTH)$fixed)
data.frame(all.varnames.SMOOTH)
all.vars.SMOOTH <- standata(brm.SMOOTH)$X

all.splines.SMOOTH <- standata(brm.SMOOTH)$Zs_2_1
pars.country.splines.SMOOTH <- summary(brm.SMOOTH$fit)$summary[,'50%'][26:38]
pars.plot.splines.SMOOTH <- summary(brm.SMOOTH$fit)$summary[,'50%'][39:51]
pars.envivars.SMOOTH <- summary(brm.SMOOTH$fit)$summary[,'50%'][39:51]


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
# 10. SUMMARIZE THE REALM PREDICTIONS
################################################################################


prd.brm.REALM <- all.vars.REALM %*% all.pars.REALM
#prd.brm.REALM <- predict(brm.REALM, newdata=DAT, type="link")[,'Estimate']
res.brm.REALM <- prd.brm.REALM - log(DAT$S)

prd.history.REALM <- draw.all(par.names = hist.parnames.REALM, 
                              vars = hist.vars.REALM, 
                              brm.fit = brm.REALM)

resid.REALM <- prd.history.REALM[,'50%'] + res.brm.REALM



HIST <- data.frame(prd.history.REALM, DAT, resid.REALM)


p.hist.multi <- ggplot(HIST, aes(log10(exp(Area_km*A.sd + A.mean)), X50.)) +
  geom_point(aes(log10(exp(Area_km*A.sd + A.mean)), resid.REALM), 
             colour="grey", shape=1) +
  geom_linerange(aes(ymin=X2.5., ymax=X97.5., colour=REALM), alpha=0.5) +
  geom_linerange(aes(ymin=X25., ymax=X75., colour=REALM), size=1) +
  # geom_hline(yintercept = 0, linetype=2) +
  geom_line(aes(colour=REALM)) + 
  xlab(expression(log[10] ~ "Area" ~ (km^2))) +
  ylab("Region effect") + 
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  theme_bw()  +  
  theme(legend.position="none") +
  ggtitle("B") +
  facet_grid(. ~ REALM) 
p.hist.multi

p.hist.single <- ggplot(HIST, aes(log10(exp(Area_km*A.sd + A.mean)), X50.)) +
  geom_linerange(aes(ymin=X2.5., ymax=X97.5., colour=REALM), alpha=0.4) +
  geom_linerange(aes(ymin=X25., ymax=X75., colour=REALM), size=1) +
  #geom_hline(yintercept = 0, linetype=2) +
  #geom_point( aes(colour=REALM)) + 
  geom_line( aes(colour=REALM)) +
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  # xlab("log10 Area [km^2]") + 
  xlab(expression(log[10] ~ "Area" ~ (km^2))) +
  ylab("Region effect") + 
  ggtitle("C") +
  theme_bw() + 
  theme(legend.position="none") 
p.hist.single

lay.mat <- matrix(c(1, 1, 2, 
                    3, 3, 3), 2, 3, byrow = TRUE)

png("../Figures/Subset_data_sensitivity_analysis/subset_historical_effects_curves.png", width=3500, height=2200, res=350)
grid.arrange(realm.plot + ggtitle("A") , 
             p.hist.single, 
             p.hist.multi, 
             #p.hist.linear,
             heights=c(1.2,0.8),
             # widths=c(4,2),
             layout_matrix=lay.mat)
dev.off()




################################################################################
# 11. TRIPHASIC SAR CURVE FROM THE SMOOTH PREDICTIONS
################################################################################

prd.area.SMOOTH <- draw.all(par.names = area.parnames.SMOOTH, 
                            vars = area.vars.SMOOTH, 
                            brm.fit = brm.SMOOTH)


partial.res.area <- prd.area.SMOOTH + residuals(gam.SMOOTH) 


AREA <- data.frame(prd.area.SMOOTH, DAT, partial.resid = partial.res.area)

png("../Figures/Subset_data_sensitivity_analysis/subset_SAR_SMOOTH.png", width=1500, height=1500, res=500)
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
# 12. GRAIN-DEPENDENT ENVIRONMENTAL PREDICTORS
################################################################################

data.frame(all.varnames.REALM)
data.frame(all.varnames.SMOOTH)


# indices of variable coefficients - model REALM
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
  A <- seq(from=-0.94, to=2.82, by=0.1)
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
scale.coefs[scale.coefs$variable=="min_DBH_cm",'variable'] <- "Minimum DBH"
scale.coefs[scale.coefs$variable=="ANN_T",'variable'] <- "Annual T"
scale.coefs[scale.coefs$variable=="ISO_T",'variable'] <- "Isothermality"
scale.coefs[scale.coefs$variable=="MIN_P",'variable'] <- "Minimum P"
scale.coefs[scale.coefs$variable=="P_SEAS",'variable'] <- "P Seasonality"
scale.coefs[scale.coefs$variable=="ALT_DIF",'variable'] <- "Elevation span"
scale.coefs[scale.coefs$variable=="ISLAND",'variable'] <- "Island"
scale.coefs$variable <- factor(scale.coefs$variable, 
                               levels=c("Island", "Elevation span", "GPP", 
                                        "P Seasonality", "Minimum P",
                                        "Annual T", "Isothermality", 
                                        "Tree density", "Minimum DBH") )

# ------------------------------------------------------------------------------
# plot the coefficients

coef.plot <- ggplot(scale.coefs, aes(x=log10(exp(A*A.sd + A.mean)), y=X50.)) + 
  geom_ribbon(aes(ymin=X2.5., ymax=X97.5., fill=Model), alpha=0.3) +
  #geom_ribbon(aes(ymin=X25., ymax=X75., fill=NA), alpha=0.3) +
  geom_line(aes(colour=Model), size=1) +
  geom_line(aes(x=log10(exp(A*A.sd + A.mean)), y=X25., colour=Model), 
            linetype="dashed") +
  geom_line(aes(x=log10(exp(A*A.sd + A.mean)), y=X75., colour=Model), 
            linetype="dashed") +
  facet_grid(.~variable) +
  geom_hline(yintercept=0, colour="darkgrey") +
  xlab("log10 Area [km^2]") + 
  ylab("Standardized coefficient") +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position="top")
coef.plot

png("../Figures/Subset_data_sensitivity_analysis/subset_envir_effects.png", width=3000, height=800, res=250)
coef.plot
dev.off()


################################################################################
# 13. PAIRPLOTS SHOWING COLLINEARITY BETWEEN THE PREDICTORS
################################################################################


for.pairs <- dplyr::select(DAT, Tree_dens, min_DBH, GPP, ANN_T, ISO_T, 
                           MIN_P, P_SEAS, ALT_DIF, ISLAND, REALM, DAT_TYPE)
for.pairs$ISLAND <- as.factor(for.pairs$ISLAND)
p <- ggpairs(for.pairs, aes(colour=DAT_TYPE)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

png("../Figures/Subset_data_sensitivity_analysis/subset_pairplot.png", width=1500, height=1700, res=100)
p
dev.off()




################################################################################
# 14. ADDITIONAL RANDOM FOREST ANALYSIS
################################################################################

# fit the random forest models
rf.plot <- randomForest(S ~ REALM + Area_km + Tree_dens + min_DBH + 
                          GPP + ANN_T + ISO_T +
                          MIN_P + P_SEAS + ALT_DIF + ISLAND, 
                        data=DAT[DAT$DAT_TYPE=="Plot",])
rf.country <- randomForest(S ~ REALM + Area_km + Tree_dens + min_DBH + 
                             GPP + ANN_T + ISO_T +
                             MIN_P + P_SEAS + ALT_DIF + ISLAND, 
                           data=DAT[DAT$DAT_TYPE=="Country",])

# plot variable importance
par(mfrow=c(1,2))
varImpPlot(rf.plot)
varImpPlot(rf.country)

# prepare data for ggplot variable importance plots
plot.imp <- data.frame(variable=rownames(importance(rf.plot)), 
                       importance=importance(rf.plot),
                       type = as.factor(c(1,0,0,0,0,0,0,0,0,0,0)))
cntr.imp <- data.frame(variable=rownames(importance(rf.country)), 
                       importance=importance(rf.country),
                       type = as.factor(c(1,0,0,0,0,0,0,0,0,0,0)))


# variable importance plots
cntr.p <- ggplot(cntr.imp, aes(variable, IncNodePurity)) + 
  geom_col(aes(fill=type)) + 
  theme_bw() + xlab("Predictor") + ylab("Importance") + coord_flip() +
  labs(title = "A", subtitle = "Country")  +
  theme(legend.position="none")


plot.p <- ggplot(plot.imp, aes(variable, IncNodePurity)) + 
  geom_col(aes(fill=type)) + 
  theme_bw() + xlab("Predictor") + ylab("Importance") + coord_flip() +
  labs(title = "B", subtitle = "Plot") +
  theme(legend.position="none")


png("../Figures/Subset_data_sensitivity_analysis/subset_random_forest_variable_importance.png", width=1000, height=500, res=150)
grid.arrange(cntr.p, plot.p, ncol=2)
dev.off()







