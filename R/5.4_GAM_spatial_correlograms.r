
source("0_libraries_functions_settings.r")

source("4_Data_loading_standardization_and_centering.r")

load("../Models/gam_REALM.Rdata")
load("../Models/gam_SMOOTH.Rdata")


################################################################################ 
# MORAN'S I CORRELOGRAMS OF THE RESIDUALS for both models at two scales
################################################################################

# width of the correlogram increment [km]
increment = 200

# how many times to resample for significance testing:
resample = 200

# ------------------------------------------------------------------------------

# Residual correlograms for Model REALM -----------------
res.REALM <- data.frame(Lat=DAT$Lat, Lon=DAT$Lon, S=DAT$S, 
                        DAT_TYPE=DAT$DAT_TYPE, 
                        resid=residuals(gam.REALM))

res.REALM.cntr <- res.REALM[res.REALM$DAT_TYPE=="Country",]
res.REALM.plot <- res.REALM[res.REALM$DAT_TYPE=="Plot",]

cor.REALM.cntr <- ncf::correlog(x=res.REALM.cntr$Lat, y=res.REALM.cntr$Lon, 
                                z=res.REALM.cntr$resid, 
                                latlon=TRUE, resamp=resample, increment=increment)
cor.REALM.plot <- ncf::correlog(x=res.REALM.plot$Lat, y=res.REALM.plot$Lon, 
                                z=res.REALM.plot$resid, 
                                latlon=TRUE, resamp=resample, increment=increment)


# Residual correlograms for Model SMOOTH -----------------
res.SMOOTH <- data.frame(Lat=DAT$Lat, Lon=DAT$Lon, S=DAT$S, 
                         DAT_TYPE=DAT$DAT_TYPE, 
                         resid=residuals(gam.SMOOTH))

res.SMOOTH.cntr <- res.SMOOTH[res.SMOOTH$DAT_TYPE=="Country",]
res.SMOOTH.plot <- res.SMOOTH[res.SMOOTH$DAT_TYPE=="Plot",]

cor.SMOOTH.cntr <- ncf::correlog(x=res.SMOOTH.cntr$Lat, y=res.SMOOTH.cntr$Lon, 
                                 z=res.SMOOTH.cntr$resid, 
                                 latlon=TRUE, resamp=resample, increment=increment)
cor.SMOOTH.plot <- ncf::correlog(x=res.SMOOTH.plot$Lat, y=res.SMOOTH.plot$Lon, 
                                 z=res.SMOOTH.plot$resid, 
                                 latlon=TRUE, resamp=resample, increment=increment)

# Correlograms for S
cor.S.cntr <- ncf::correlog(x=res.SMOOTH.cntr$Lat, y=res.SMOOTH.cntr$Lon, 
                            z=res.SMOOTH.cntr$S, 
                            latlon=TRUE, resamp=resample, increment=increment)
cor.S.plot <- ncf::correlog(x=res.SMOOTH.plot$Lat, y=res.SMOOTH.plot$Lon, 
                            z=res.SMOOTH.plot$S, 
                            latlon=TRUE, resamp=resample, increment=increment)


save(cor.REALM.cntr,
     cor.REALM.plot,
     cor.SMOOTH.cntr,
     cor.SMOOTH.plot,
     cor.S.cntr,
     cor.S.plot,
     file = "../Models/correlograms.Rdata")


# ------------------------------------------------------------------------------

load("../Models/correlograms.Rdata")

# extract the correlograms for further plotting
N <- length(cor.REALM.plot$mean.of.class)

RP <- data.frame(Dist=cor.REALM.plot$mean.of.class, 
                 Cor=cor.REALM.plot$correlation,
                 P = ifelse(cor.REALM.plot$p < 0.01, "P < 0.01", "not significant" ),
                 Scale=rep("Plot"), 
                 Variable=rep("REALM residuals"))
RC <- data.frame(Dist=cor.REALM.cntr$mean.of.class, 
                 Cor=cor.REALM.cntr$correlation,
                 P = ifelse(cor.REALM.cntr$p < 0.01, "P < 0.01", "not significant" ), 
                 Scale=rep("Country"), 
                 Variable=rep("REALM residuals"))
SP <- data.frame(Dist=cor.SMOOTH.plot$mean.of.class, 
                 Cor=cor.SMOOTH.plot$correlation,
                 P = ifelse(cor.SMOOTH.plot$p < 0.01, "P < 0.01", "not significant" ), 
                 Scale=rep("Plot"), 
                 Variable=rep("SMOOTH residuals"))
SC <- data.frame(Dist=cor.SMOOTH.cntr$mean.of.class, 
                 Cor=cor.SMOOTH.cntr$correlation,
                 P = ifelse(cor.SMOOTH.cntr$p < 0.01, "P < 0.01", "not significant" ), 
                 Scale=rep("Country"), 
                 Variable=rep("SMOOTH residuals"))
RIC <- data.frame(Dist=cor.S.cntr$mean.of.class, 
                  Cor=cor.S.cntr$correlation,
                  P = ifelse(cor.S.cntr$p < 0.01, "P < 0.01", "not significant" ), 
                  Scale=rep("Country"), 
                  Variable=rep("Species richness S"))
RIP <- data.frame(Dist=cor.S.plot$mean.of.class, 
                  Cor=cor.S.plot$correlation,
                  P = ifelse(cor.S.plot$p < 0.01, "P < 0.01", "not significant" ), 
                  Scale=rep("Plot"), 
                  Variable=rep("Species richness S"))

cor.data <- rbind(RP, RC, SP, SC, RIP, RIC)


# plot the correlograms
png("../Figures/Fig_S4_correlograms.png", width=2000, height=900, res=250)
ggplot(cor.data, aes(x=Dist, y=Cor)) + 
  geom_hline(yintercept = 0, colour="darkgrey") +
  geom_line(aes(colour=Variable)) + 
  geom_point(aes(colour=Variable, shape=P), size = 2) + 
  scale_shape_manual(values = c(1,19)) + 
  xlim(c(0, 3000)) + ylim(c(-0.2,0.5)) +
  scale_colour_brewer(palette="Set1") +
  facet_grid(.~Scale) +
  xlab("Distance [km]") + ylab("Moran's I") +
  theme_bw()
dev.off()
