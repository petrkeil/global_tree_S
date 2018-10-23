################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Dec 7 2017
################################################################################

# Description: Here is where I perform the deviance partitioning. In short, model
# REALM is fitted to different subset of the data, and the deviance explained by
# environmental variables or biogeographic realms is calculated.

# Also, here is where the PCA plots of the environmental conditions within realms
# are calculated.

################################################################################

# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")

# load the data
source("4_Data_loading_standardization_and_centering.r")




################################################################################
# THE DEVIANCE PARTITIONING
################################################################################
  
# define the family that will be used in the exercise (so that we can easily
# also check the effect of this on the obtained partitions)
fam = "nb"

  
# function that extracts R2 and deviance explained from a model object
perf <- function(obj)
{
 res <- c(dev.exp = summary(obj)$dev.exp, 
          r.sq = summary(obj)$r.sq)
 return(res)
}


################################################################
# DEVIANCE PARTITIONING - MODEL REALM, PLOT GRAIN, GLOBAL EXTENT
################################################################

plot.hist <- gam(S ~ REALM + poly(Area_km,3):REALM, 
                 data=DAT[DAT$DAT_TYPE=="Plot",],
                 family=fam)


plot.env <- gam(S ~   Tree_dens + Tree_dens:Area_km + 
                  min_DBH + min_DBH:Area_km +
                  GPP + GPP:Area_km + 
                  ANN_T + ANN_T:Area_km +
                  ISO_T + ISO_T:Area_km + 
                  MIN_P + MIN_P:Area_km +  
                  P_SEAS + P_SEAS:Area_km + 
                  ALT_DIF + ALT_DIF:Area_km +
                  ISLAND + ISLAND:Area_km, 
                data=DAT[DAT$DAT_TYPE=="Plot",], 
                family=fam)

plot.full <- gam(S ~   REALM + poly(Area_km,3):REALM +
                   Tree_dens + Tree_dens:Area_km + 
                   min_DBH + min_DBH:Area_km +
                   GPP + GPP:Area_km + 
                   ANN_T + ANN_T:Area_km +
                   ISO_T + ISO_T:Area_km + 
                   MIN_P + MIN_P:Area_km +  
                   P_SEAS + P_SEAS:Area_km + 
                   ALT_DIF + ALT_DIF:Area_km +
                   ISLAND + ISLAND:Area_km, 
                 data=DAT[DAT$DAT_TYPE=="Plot",], 
                 family=fam)

env <- perf(plot.env)
env
hist <- perf(plot.hist) 
hist
full  <- perf(plot.full)
full

# independent - environment
indep.env = full - hist
indep.env
# independent - history
indep.hist <- full - env
indep.hist
# overlap
overlap = full - (indep.hist + indep.env)
overlap
# unexplained
1 - full

###################################################################
# DEVIANCE PARTITIONING - MODEL REALM, COUNTRY GRAIN, GLOBAL EXTENT
###################################################################

cntr.hist <- gam(S ~ REALM + poly(Area_km,3):REALM, 
                 data=DAT[DAT$DAT_TYPE=="Country",],
                 family=fam)
cntr.env <- gam(S ~   Tree_dens + Tree_dens:Area_km + 
                  min_DBH + min_DBH:Area_km +
                  GPP + GPP:Area_km + 
                  ANN_T + ANN_T:Area_km +
                  ISO_T + ISO_T:Area_km + 
                  MIN_P + MIN_P:Area_km +  
                  P_SEAS + P_SEAS:Area_km + 
                  ALT_DIF + ALT_DIF:Area_km +
                  ISLAND + ISLAND:Area_km, 
                data=DAT[DAT$DAT_TYPE=="Country",], 
                family=fam)

cntr.full <- gam(S ~   REALM + poly(Area_km,3):REALM +
                   Tree_dens + Tree_dens:Area_km + 
                   min_DBH + min_DBH:Area_km +
                   GPP + GPP:Area_km + 
                   ANN_T + ANN_T:Area_km +
                   ISO_T + ISO_T:Area_km + 
                   MIN_P + MIN_P:Area_km +  
                   P_SEAS + P_SEAS:Area_km + 
                   ALT_DIF + ALT_DIF:Area_km +
                   ISLAND + ISLAND:Area_km, 
                 data=DAT[DAT$DAT_TYPE=="Country",], 
                 family=fam)


env <- perf(cntr.env)
env
hist <- perf(cntr.hist) 
hist
full  <- perf(cntr.full)
full

# independent - environment
indep.env = full - hist
indep.env
# independent - history
indep.hist <- full - env
indep.hist
# overlap
overlap = full - (indep.hist + indep.env)
overlap
# unexplained
1 - full



################################################################################
# DEVIANCE PARTITIONING, NORTHERN HEMISPHERE
################################################################################

DAT.north <- DAT[DAT$REALM %in% c("Western Palearctic","Eastern Palearctic","Nearctic"),]
DAT.north <- dplyr::filter(.data=DAT,
                           DAT$REALM %in% c("Western Palearctic","Eastern Palearctic","Nearctic"))



######################################################################
# DEVIANCE PARTITIONING - MODEL REALM, PLOT GRAIN, NORTHERN HEMISPHERE
######################################################################

plot.hist <- gam(S ~ REALM + poly(Area_km,3):REALM, 
                 data=DAT.north[DAT.north$DAT_TYPE=="Plot",],
                 family="nb")
plot.env <- gam(S ~   Tree_dens + Tree_dens:Area_km + 
                  min_DBH + min_DBH:Area_km +
                  GPP + GPP:Area_km + 
                  ANN_T + ANN_T:Area_km +
                  ISO_T + ISO_T:Area_km + 
                  MIN_P + MIN_P:Area_km +  
                  P_SEAS + P_SEAS:Area_km + 
                  ALT_DIF + ALT_DIF:Area_km +
                  ISLAND + ISLAND:Area_km, 
                  data=DAT.north[DAT.north$DAT_TYPE=="Plot",], 
                  family="nb")

plot.full <- gam(S ~   REALM + poly(Area_km,3):REALM +
                   Tree_dens + Tree_dens:Area_km + 
                   min_DBH + min_DBH:Area_km +
                   GPP + GPP:Area_km + 
                   ANN_T + ANN_T:Area_km +
                   ISO_T + ISO_T:Area_km + 
                   MIN_P + MIN_P:Area_km +  
                   P_SEAS + P_SEAS:Area_km + 
                   ALT_DIF + ALT_DIF:Area_km +
                   ISLAND + ISLAND:Area_km, 
                   data=DAT.north[DAT.north$DAT_TYPE=="Plot",], 
                   family="nb")

env <- perf(plot.env)
env
hist <- perf(plot.hist) 
hist
full  <- perf(plot.full)
full

# independent - environment
indep.env = full - hist
indep.env
# independent - history
indep.hist <- full - env
indep.hist
# overlap
overlap = full - (indep.hist + indep.env)
overlap
# unexplained
1 - full

#########################################################################
# DEVIANCE PARTITIONING - MODEL REALM, COUNTRY GRAIN, NORTHERN HEMISPHERE
#########################################################################

cntr.hist <- gam(S ~ REALM + poly(Area_km,3):REALM, 
                 data=DAT.north[DAT.north$DAT_TYPE=="Country",],
                 family="nb")

cntr.env <- gam(S ~   Tree_dens + Tree_dens:Area_km + 
                  min_DBH + min_DBH:Area_km +
                  GPP + GPP:Area_km + 
                  ANN_T + ANN_T:Area_km +
                  ISO_T + ISO_T:Area_km + 
                  MIN_P + MIN_P:Area_km +  
                  P_SEAS + P_SEAS:Area_km + 
                  ALT_DIF + ALT_DIF:Area_km +
                  ISLAND + ISLAND:Area_km, 
                data=DAT.north[DAT.north$DAT_TYPE=="Country",], 
                family="nb")

cntr.full <- gam(S ~   REALM + poly(Area_km,3):REALM +
                   Tree_dens + Tree_dens:Area_km + 
                   min_DBH + min_DBH:Area_km +
                   GPP + GPP:Area_km + 
                   ANN_T + ANN_T:Area_km +
                   ISO_T + ISO_T:Area_km + 
                   MIN_P + MIN_P:Area_km +  
                   P_SEAS + P_SEAS:Area_km + 
                   ALT_DIF + ALT_DIF:Area_km +
                   ISLAND + ISLAND:Area_km, 
                 data=DAT.north[DAT.north$DAT_TYPE=="Country",], 
                 family="nb")


env <- perf(cntr.env)
env
hist <- perf(cntr.hist) 
hist
full  <- perf(cntr.full)
full
# independent - environment
indep.env = full - hist
indep.env
# independent - history
indep.hist <- full - env
indep.hist
# overlap
overlap = full - (indep.hist + indep.env)
overlap
# unexplained
1 - full


DAT.north <- DAT[DAT$REALM %in% c("Western Palearctic","Eastern Palearctic","Nearctic"),]
DAT.north <- dplyr::filter(.data=DAT,
                           DAT$REALM %in% c("Western Palearctic","Eastern Palearctic","Nearctic"))



################################################################################
# DEVIANCE PARTITIONING, TROPICS
################################################################################

DAT.south <- DAT[DAT$REALM %in% c("Western Palearctic","Eastern Palearctic","Nearctic") == FALSE,]

##########################################################
# DEVIANCE PARTITIONING - MODEL REALM, PLOT GRAIN, TROPICS
##########################################################

plot.hist <- gam(S ~ REALM + poly(Area_km,3):REALM, 
                 data=DAT.south[DAT.south$DAT_TYPE=="Plot",],
                 family="nb")
plot.env <- gam(S ~   Tree_dens + Tree_dens:Area_km + 
                  min_DBH + min_DBH:Area_km +
                  GPP + GPP:Area_km + 
                  ANN_T + ANN_T:Area_km +
                  ISO_T + ISO_T:Area_km + 
                  MIN_P + MIN_P:Area_km +  
                  P_SEAS + P_SEAS:Area_km + 
                  ALT_DIF + ALT_DIF:Area_km +
                  ISLAND + ISLAND:Area_km, 
                data=DAT.south[DAT.south$DAT_TYPE=="Plot",], 
                family="nb")

plot.full <- gam(S ~   REALM + poly(Area_km,3):REALM +
                   Tree_dens + Tree_dens:Area_km + 
                   min_DBH + min_DBH:Area_km +
                   GPP + GPP:Area_km + 
                   ANN_T + ANN_T:Area_km +
                   ISO_T + ISO_T:Area_km + 
                   MIN_P + MIN_P:Area_km +  
                   P_SEAS + P_SEAS:Area_km + 
                   ALT_DIF + ALT_DIF:Area_km +
                   ISLAND + ISLAND:Area_km, 
                 data=DAT.south[DAT.south$DAT_TYPE=="Plot",], 
                 family="nb")

env <- perf(plot.env)
env
hist <- perf(plot.hist) 
hist
full  <- perf(plot.full)
full

# independent - environment
indep.env = full - hist
indep.env
# independent - history
indep.hist <- full - env
indep.hist
# overlap
overlap = full - (indep.hist + indep.env)
overlap
# unexplained
1 - full


#############################################################
# DEVIANCE PARTITIONING - MODEL REALM, COUNTRY GRAIN, TROPICS
#############################################################

cntr.hist <- gam(S ~ REALM + poly(Area_km,3):REALM, 
                 data=DAT.south[DAT.south$DAT_TYPE=="Country",],
                 family="nb")
cntr.env <- gam(S ~   Tree_dens + Tree_dens:Area_km + 
                  min_DBH + min_DBH:Area_km +
                  GPP + GPP:Area_km + 
                  ANN_T + ANN_T:Area_km +
                  ISO_T + ISO_T:Area_km + 
                  MIN_P + MIN_P:Area_km +  
                  P_SEAS + P_SEAS:Area_km + 
                  ALT_DIF + ALT_DIF:Area_km +
                  ISLAND + ISLAND:Area_km, 
                data=DAT.south[DAT.south$DAT_TYPE=="Country",], 
                family="nb")

cntr.full <- gam(S ~   REALM + poly(Area_km,3):REALM +
                   Tree_dens + Tree_dens:Area_km + 
                   min_DBH + min_DBH:Area_km +
                   GPP + GPP:Area_km + 
                   ANN_T + ANN_T:Area_km +
                   ISO_T + ISO_T:Area_km + 
                   MIN_P + MIN_P:Area_km +  
                   P_SEAS + P_SEAS:Area_km + 
                   ALT_DIF + ALT_DIF:Area_km +
                   ISLAND + ISLAND:Area_km, 
                 data=DAT.south[DAT.south$DAT_TYPE=="Country",], 
                 family="nb")


env <- perf(cntr.env)
env
hist <- perf(cntr.hist) 
hist
full  <- perf(cntr.full)
full

# independent - environment
indep.env = full - hist
indep.env
# independent - history
indep.hist <- full - env
indep.hist
# overlap
overlap = full - (indep.hist + indep.env)
overlap
# unexplained
1 - full






################################################################################
# PCA of the data to show the environmental overlap
################################################################################

library(ade4)


# PLOT-LEVEL PCA ---------------------------------------------------------------

# convert factor ISLAND to 0/1
island = ((DAT$ISLAND == "island") *1)
DAT <- data.frame(DAT, island)

# select only the relevant environmental predictors
for.PCA.plot <- DAT[DAT$DAT_TYPE=="Plot",
                    c("REALM", "Tree_dens", "GPP", "ANN_T", "ISO_T", 
                      "MIN_P", "P_SEAS", "ALT_DIF", "island")]



# do the PCA:
PCA.plot <- dudi.pca(for.PCA.plot[,2:ncol(for.PCA.plot)], scannf=FALSE, nf=4,
                     scale=TRUE, center=TRUE)
summary(PCA.plot)

# export the plot coordinates
PCA.plot <- data.frame(PCA.plot$li, for.PCA.plot)

# calculate medians (these will be the large colourful symbols)
PCA.plot.sum <- ddply(.data=PCA.plot,
                      .variables = "REALM",
                      .fun=summarize,
                      Axis1 = median(Axis1),
                      Axis2 = median(Axis2))

# make the PCA biplot for the 1st and 2nd axis
PCA.A <- ggplot(data=PCA.plot, aes(x=Axis1, y=Axis2)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept=0) +
  geom_point(aes(colour=REALM)) + 
  geom_point(data=PCA.plot.sum, 
             aes(x=Axis1, y=Axis2, colour=REALM), size=4) +
  geom_point(data=PCA.plot.sum, 
             aes(x=Axis1, y=Axis2), size=4, shape=1) +
  geom_text(data=PCA.plot.sum, 
            aes(x=Axis1, y=Axis2, label=REALM), size=4, nudge_y = 0.25) +
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  labs(title = "a", subtitle = "Plot grain")  +
  theme_bw() + 
  theme(legend.position="none", 
        plot.title = element_text(face=quote(bold)),
        legend.background = element_rect(fill="lightgrey",
                                         size=0.5, linetype="solid", 
                                         colour ="darkgrey"))

# COUNTRY-LEVEL PCA ---------------------------------------------------------------

# select only the relevant environmental predictors
for.PCA.cntr <- DAT[DAT$DAT_TYPE=="Country",
                    c("REALM", "Tree_dens", "GPP", "ANN_T", "ISO_T", 
                      "MIN_P", "P_SEAS", "ALT_DIF", "island")]
# do the PCA:
PCA.cntr <- dudi.pca(for.PCA.cntr[,2:ncol(for.PCA.cntr)], scannf=FALSE, nf=4,
                     scale=TRUE, center=TRUE)

PCA.cntr <- data.frame(PCA.cntr$li, for.PCA.cntr)
PCA.cntr.sum <- ddply(.data=PCA.cntr,
                      .variables = "REALM",
                      .fun=summarize,
                      Axis1 = median(Axis1),
                      Axis2 = median(Axis2))

PCA.B <- ggplot(data=PCA.cntr, aes(x=Axis1, y=Axis2)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept=0) +
  geom_point(aes(colour=REALM)) + 
  geom_point(data=PCA.cntr.sum, 
             aes(x=Axis1, y=Axis2, colour=REALM), size=4) +
  geom_point(data=PCA.cntr.sum, 
             aes(x=Axis1, y=Axis2), size=4, shape=1) +
  geom_text(data=PCA.cntr.sum, 
            aes(x=Axis1, y=Axis2, label=REALM), size=4, nudge_y = 0.25) +
  scale_colour_brewer(palette = "Dark2", name="Realm") +
  labs(title = "b", subtitle = "Country grain")  +
  theme_bw() +  theme(legend.position="none",
                      plot.title = element_text(face=quote(bold)))

png(file="../Figures/PCA_of_environmental_variables.png", width=2000, height=1000, res=180)
grid.arrange(PCA.A, PCA.B, ncol=2)
dev.off()





