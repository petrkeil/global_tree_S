################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Oct 25 2018
################################################################################

#
# Description: Here is where we produce the maps that use data observational units
# (the original plots and countries), both for mapping of the raw data and the
# model predictions.
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
CNTRS@data <- dplyr::select(CNTRS@data, NAME)

# remove the big countries
CNTRS <- CNTRS[(CNTRS$NAME %in% c("China", "United States", "Brazil")) == FALSE,]

# reproject the countries
CNTRS <- spTransform(CNTRS, CRSobj = CRS(MOLLWEIDE))

PRED.CNTRS <- merge(CNTRS, PRED.CNTRS, by.x="NAME", by.y="Loc_ID")

PRED.CNTRS[PRED.CNTRS$NAME == "Alaska","REALM"] <- "Nearctic"
PRED.CNTRS[PRED.CNTRS$NAME == "Alaska","ISLAND"] <- 0

# create the norhtern-southern hemisphere identifier
PRED.CNTRS@data <- data.frame(PRED.CNTRS@data, Hemisphere=rep("Palearctic & Nearctic",
                                                              times=nrow(PRED.CNTRS@data)))
is.Tropics <- (PRED.CNTRS@data$REALM %in% c("Western Palearctic","Eastern Palearctic","Nearctic")) == FALSE
a <- PRED.CNTRS@data
a$Hemisphere <- as.character(a$Hemisphere)
a$Hemisphere[is.Tropics] <- "Tropics"
PRED.CNTRS@data <- a



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
                     plot.title = element_text(face=quote(bold)),
                     #legend.title.align = 0,
                     panel.background=element_blank(),
                     panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())



# ------------------------------------------------------------------------------
# RAW RICHNESS MAPS

s.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=S), colour="black", size=.2) + 
  scale_fill_distiller(palette = "Spectral", name="S", 
                       trans="log10", limits=c(1,6500)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  ggtitle("a") + theme_minimal() +
  xlab("") + ylab("") + blank.theme

s.plot <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="white", size=0.2) +
  geom_polygon(fill="black", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=S), size=1, shape=1) +
  scale_colour_distiller(palette = "Spectral", name="S", 
                         trans="log10", limits=c(1,6500)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  ggtitle("b") + theme_minimal() +
  xlab("") + ylab("") + blank.theme +
  theme(panel.background = element_rect(fill = "grey",
                                        colour = "grey",
                                        size = 0.5, linetype = "solid"),)

tiff("../Figures/observed_richness_maps.tif", width=2000, height=2100, res=350, 
     compression = "lzw")
grid.arrange(s.cntr, s.plot, ncol=1, nrow=2)
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
  ggtitle("a") + theme_minimal() +
  xlab("") + ylab("") + blank.theme

s.pred.plot <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="white", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=S.pred)) +
  scale_colour_distiller(palette = "Spectral", name="Predicted S", 
                         trans="log10") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  ggtitle("b") + theme_minimal() +
  xlab("") + ylab("") + blank.theme



# ------------------------------------------------------------------------------
# SMOOTHER MAPS in observational units

g.cntr <- ggplot(C.fort, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=smooth.country), colour="black", size=.2) + 
  scale_fill_distiller(palette = "Spectral", 
                       limits=c(-2, 2),
                       name="Region\neffect") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") + 
  ggtitle("a") +  
  theme_minimal() + blank.theme


g.plot <- ggplot(MAINL, aes(long, lat, group=group)) + 
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(fill="white", colour="black", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour=smooth.plot)) +
  scale_colour_distiller(palette = "Spectral", 
                         limits=c(-2, 2),
                         name="Region\neffect") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") +
  ggtitle("b") +
  theme_minimal() + blank.theme

# save to a file
png("../Figures/smoothing_splines_at_data_points_map.png", width=2000, height=2000, res=250)
grid.arrange(g.cntr, g.plot, ncol=1, nrow=2)
dev.off()



# ------------------------------------------------------------------------------
# BIOGEOGRAPHIC REALMS

blank.theme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())

# save to a file
png("../Figures/realms_map.png", width=2300, height=1000, res=200)
realm.plot <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(aes(fill=REALM), colour="darkgray", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), shape=3, colour="black") +
  scale_fill_brewer(palette = "Dark2", name="Realm") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))
realm.plot
dev.off()

# save to a file
png("../Figures/realms_map.png", width=2300, height=1000, res=200)
realm.plot <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(aes(fill=REALM), colour="darkgray", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), shape=1, colour="black") +
  scale_fill_brewer(palette = "Dark2", name="Realm") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))
realm.plot
dev.off()

# save to a file
png("../Figures/realms_map_small.png", width=2400, height=1000, res=400)
realm.plot.small <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(aes(fill=REALM), colour="darkgray", size=.2) + 
  scale_fill_brewer(palette = "Dark2", name="Realm") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  theme(legend.position="right") +
  guides(fill=guide_legend(title=NULL))
realm.plot.small
dev.off()


# ------------------------------------------------------------------------------
# ISLANDS VS MAINLANDS

# The island/mainland/shelf classification from Weigelt et al.
ISL_MAINL <- raster("../Data/GRIDS/ALL_CLASSES_clean.tif")
ISL_MAINL[ISL_MAINL == 0] <- NA
ISL_MAINL_MOLL <- projectRaster(ISL_MAINL, crs = MOLLWEIDE, method = "ngb")
isl.plot.orig <- gplot(ISL_MAINL_MOLL, maxpixels = 1e6) + 
  geom_raster(aes(fill=factor(value))) + 
  #scale_fill_brewer(palette = "Set2", name="Island") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  scale_fill_manual(labels = c("mainland", "shelf island", "true island", ""), 
                    values = c("#377EB8","#4DAF4A","#E41A1C"))+
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  ggtitle("a")
isl.plot.orig


# The island/mainland classification that we end up using
isl.plot.countries <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(aes(fill=as.factor(ISLAND)), colour="darkgrey", size=.2) + 
  #geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour = as.factor(ISLAND)), shape=3) +
  #scale_fill_brewer(palette = "Set1", name="Island") +
  scale_fill_manual(labels = c("island", "mainland", ""), 
                    values = c("#E41A1C", "#377EB8"))+
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme +
  ggtitle("b")
isl.plot.countries

isl.plot.plots <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  geom_polygon(fill="white", colour="darkgrey", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL, colour = as.factor(ISLAND)), shape=3) +
  scale_colour_brewer(palette = "Set1", name="Island") +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") + theme_minimal() + blank.theme+
  ggtitle("c")
isl.plot.plots

# save to a file
png("../Figures/island_vs_mainland_map.png", width=2300, height=3450, res=250)
grid.arrange(isl.plot.orig, isl.plot.countries, isl.plot.plots, nrow=3, ncol=1)
dev.off()

# ------------------------------------------------------------------------------
# NORTHERN AND SOUTHERN HEMISPHERES

# save to a file
png("../Figures/north_south_hemispheres_map.png", width=1800, height=1000, res=200)
hemi.plot <- ggplot(C.fort[is.na(C.fort$REALM) == FALSE,], aes(long, lat, group=group)) +
  #geom_hline(yintercept = 0, colour="black", size=0.2) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(aes(fill=as.factor(Hemisphere)), colour="darkgrey", size=.2) + 
  geom_point(data=PRED.PLOTS, aes(x=X, y=Y, group=NULL), shape=3, colour="black") +
  scale_fill_brewer(palette = "Reds", name="", type=qual) +
  scale_x_continuous(limits = c(-13000000, 16000000)) +
  xlab("") + ylab("") + theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position="none") + blank.theme
hemi.plot
dev.off()
