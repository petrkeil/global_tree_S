################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Dec 13 2018
################################################################################

# Description: Here is where model SMOOTH is used to generate predictions to the
# regular global network of 1 ha plots, and to the grid of large hexagons.


################################################################################

# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")


################################################################################
### Read, transform and scale the data

# read the data
pts <- read.csv(file="../Data/GRIDS/Fine_points_with_environment.csv")
grid5 <- readOGR(dsn = "../Data/GRIDS", layer = "hex5_with_environment")
grid5 <- spTransform(x = grid5, CRSobj = WGS84)

# -----------------------------------------

pts$Tree_dens <- (pts$TREE_DENS + 1) / pts$A # calculate tree density (note the x+1 step!!)
pts <- data.frame(pts, 
                  Area_km = 0.01, 
                  min_DBH = 0, 
                  ELONG = 1,
                  DAT_TYPE = "Plot")

# tree density at the grid level
grid5$Tree_dens <- (grid5$TREE_DENS + 1) / grid5$LandArea
grid5@data <- data.frame(grid5@data, 
                         min_DBH = 0, 
                         ELONG = 1,
                         DAT_TYPE = "Country")

# -----------------------------------------

pts <- dplyr::select(pts, Area_km, Tree_dens, min_DBH, 
                     GPP, ANN_T, ISO_T, MIN_P, P_SEAS, ALT_DIF, ELONG,
                     ISLAND = ISL_LS, Lat, Lon, DAT_TYPE) %>%
              mutate(Area_km = log(Area_km), Tree_dens=log(Tree_dens))

grid5.dat <- dplyr::select(grid5@data, Area_km = LandArea, Tree_dens, min_DBH,
                           GPP, ANN_T, ISO_T, MIN_P, P_SEAS, ALT_DIF, ELONG,
                           ISLAND = ISL_LS, Lat, Lon, DAT_TYPE) %>%
                    mutate(Area_km = log(Area_km), Tree_dens=log(Tree_dens))

# get the scaling constants that were used to scale the raw plot and country data:
scal.tab <- read.csv("scale_tab.csv")
scal.tab <- scal.tab[scal.tab$var %in% c("ET","WARM_T") == FALSE,]

# scale the grid data in the same way as the original data
pts[,1:10] <- scale(pts[,1:10],
                   center = scal.tab$centr, 
                   scale = scal.tab$scale)

grid5.dat[,1:10] <- scale(grid5.dat[,1:10],
                   center = scal.tab$centr, 
                   scale = scal.tab$scale)

################################################################################
### Make the predictions

# load the saved SMOOTH model that will be used for the global predictions
library(mgcv)
load("../Models/gam_SMOOTH.Rdata")



################################################################################
### Predictions in hexagons

# predict S from the model SMOOTH
grid.pred.S <- predict(gam.SMOOTH, 
                       newdata = grid5.dat, 
                       type="response")
grid.preds.S <- round(grid.pred.S, 2)


# predict the regional effect from model SMOOTH
grid.pred.smth <- predict.gam(gam.SMOOTH, 
                              type = "terms",
                              newdata = grid5.dat)[,"s(Lat,Lon):DAT_TYPECountry"]

# merge with the original grid
grid5@data <- data.frame(grid5@data, S = grid.pred.S, smooth.country = grid.pred.smth)
grid5@data$id <- as.character(grid5@data$id)

# remove cells with little land area
good.cells <-  grid5@data$LandArea / grid5@data$CellArea > 0.5
good.cells[is.na(good.cells)] <- FALSE
grid5 <- grid5[good.cells,]

# remove cells with 0 or NA species predicted
good.cells <- grid5@data$S > 1
good.cells[is.na(good.cells)] <- FALSE
grid5 <- grid5[good.cells, ]


################################################################################
### Predictions in 1 ha plots

# predict S in the plots from the SMOOTH model
plot.pred.S <- predict(gam.SMOOTH, 
                       newdata = pts, 
                       type="response")
plot.pred.S <- round(plot.pred.S, 2)

# predict region effects in the plots from the SMOOTH model
plot.pred.smth <- predict.gam(gam.SMOOTH, 
                              type = "terms",
                              newdata = pts)[,"s(Lat,Lon):DAT_TYPEPlot"]

# put all together
plot.preds <- data.frame(pts, 
                         S = as.numeric(plot.pred.S), 
                         smooth.plot = plot.pred.smth)

# remove predictions of S < 0.8 (an arbitrarily selected threshold)
plot.preds$S[plot.preds$S < 0.8] <- NA

plot.preds <- plot.preds[rowSums(is.na(plot.preds)) == 0,]

# put predictions to a spatial object
plot.preds <- SpatialPointsDataFrame(coords = data.frame(plot.preds$Lon, plot.preds$Lat), 
                                     data = plot.preds, 
                                     proj4string = CRS(WGS84))

# ------------------------------------------------------------------------------
# calculate BETA DIVERSITY for PLOTS

# extract S values from the hexagonal grid to the points
gamma <- over(x=plot.preds, y=grid5)[,c("S", "ALT_DIF", "smooth.country")]
names(gamma) <- c("gamma", "ALT_DIF_grid", "smooth.country")
# calculate beta diversity per plot
plot.preds@data <- data.frame(plot.preds@data, gamma) %>%
  mutate(beta = gamma/S, reg.beta = exp(smooth.country)/exp(smooth.plot))

# ------------------------------------------------------------------------------
# write out data with no NA values
write.csv(plot.preds@data, 
          file="../Data/GRIDS/Fine_points_with_predictions.csv", row.names=FALSE)

# ------------------------------------------------------------------------------
# transform the data for fancy plotting
plot.preds.ml <- spTransform(plot.preds, CRSobj = MOLLWEIDE)
plot.preds.ml <- data.frame(plot.preds.ml@data, 
                            data.frame(X=coordinates(plot.preds.ml)[,1],
                                       Y=coordinates(plot.preds.ml)[,2]))

grid5.ml <- spTransform(grid5, CRSobj=MOLLWEIDE)
grid5.mlf <- tidy(grid5.ml, region="id")
grid5.mlf <- left_join(x=grid5.mlf, y=grid5.ml@data, by="id")



################################################################################
# PLOTTING THE MAPS 

# Read the shapefiles 

  # coutnry boundaries
  CNTR <- readOGR(dsn="../Data/COUNTRIES", layer="COUNTRIES")
  CNTRml <- spTransform(CNTR, CRSobj=MOLLWEIDE)
  CNTRml <- tidy(CNTRml, region="NAME")
  
  # global mainlands (not divided by country boundaries)
  MAINL <- readOGR(dsn = "../Data/COUNTRIES", layer = "GSHHS_i_L1_simple")
  MAINL <- spTransform(MAINL, CRSobj = CRS(MOLLWEIDE))
  MAINL <- tidy(MAINL, region="id")
  
  # equator, tropics, and polar circles
  LINES <- readOGR(dsn = "../Data/COUNTRIES", layer = "ne_110m_geographic_lines")
  LINES <- spTransform(LINES, CRSobj = CRS(MOLLWEIDE))
  LINES <- tidy(LINES, region="name")


blank.theme <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                     axis.text.y=element_blank(),axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.position=c(0.63, 0.09),
                     legend.direction = "horizontal",
                      legend.title = element_blank(),
                     legend.title.align = 0,
                     #plot.title = element_text(hjust = 0),
                     plot.subtitle = element_text(vjust=-3),
                     panel.background=element_blank(),
                     panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())


S.option = "magma"
R.option = "viridis"

# predicted S in hexagons
plot.gr.S <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=S)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_viridis(option = S.option, 
                       name=expression(S[hex]), 
                       #limits=c(1,5000),
                       trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  #ggtitle("A") +  
  labs(subtitle = expression(hat(S)[hex] ~ "- predicted richness in 209,903" ~ km^2 ~ "hexagons")) +
  theme_minimal() + blank.theme
#plot.gr.S


# predicted S in plots
plot.pl.S <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=S))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_viridis(option = S.option, 
                         name=expression(S[plot]),
                         #limits=c(1,5000),
                         trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  #ggtitle("B") +  
  labs(subtitle = expression(hat(S)[plot] ~ "- predicted richness in 1 ha plots")) +
  theme_minimal() + blank.theme


# predicted beta in plots
plot.beta.S <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=beta))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_viridis(option = S.option, 
                         name=expression(beta),
                         trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  #ggtitle("C") +  
  labs(subtitle = expression(beta ~ "=" ~ hat(S)[hex]/hat(S)[plot])) +
  theme_minimal() + blank.theme


# predicted region effects in the hexagons
plot.gr.smth <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=smooth.country)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_viridis(option = R.option, 
                       name="Region effect",
                       limits=c(-2, 2)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  #ggtitle("D") +  
  labs(subtitle = expression(s[2](Lat, Lon) ~ "- smooth region effects in 209,903" ~ km^2 ~ "hexagons")) +
  theme_minimal() + blank.theme


# predicted region effects in the plots
plot.pl.smth <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=smooth.plot))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_viridis(option = R.option, 
                       limits=c(-2, 2),
                       name= R.option) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(s[1](Lat, Lon) ~ "- smooth region effects in 1 ha plots")) +
  theme_minimal() + blank.theme 
#plot.pl.smth

# predicted ratios of region effects between local and hexagon grains
plot.beta.smth <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=reg.beta))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_viridis(option = R.option, 
                         name=expression("Region" ~ beta)) +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  #ggtitle("F") +  
  labs(subtitle = expression(Delta ~ "=" ~ e^s[2](Lat, Lon) / e^s[1](Lat, Lon) )) +
  theme_minimal() + blank.theme
#plot.beta.smth

# ------------------------------------------------------------------------------

# write to file

library(cowplot)
tiff("../Figures/SMOOTH_prediction_grids.tif", width=4000, height=3400, res=350,
     compression = "lzw")
  plot_grid(plot.gr.S, plot.gr.smth,
            plot.pl.S, plot.pl.smth, 
            plot.beta.S, plot.beta.smth,
            nrow=3, ncol=2,
            labels = c("a", "d", "b", "e", "c", "f"), vjust = 1.1, hjust = -0.1)
dev.off()

#pdf("../Figures/SMOOTH_prediction_grids.pdf", width=10, height=8.5)
#plot_grid(plot.gr.S, plot.gr.smth,
#          plot.pl.S, plot.pl.smth, 
#          plot.beta.S, plot.beta.smth,
#          nrow=3, ncol=2,
#          labels = c("a", "d", "b", "e", "c", "f"), vjust = 1.1, hjust = -0.1)
#dev.off()



################################################################################
# LATITUDINAL GRADIENTS OF THE PREDICTIONS
################################################################################

# data for latitudinal gradient plots
LG.data <-
  rbind(
    data.frame(Latitude = grid5.ml@data$Lat, S = grid5.ml@data$S, 
               Grain = "hexagons", Longitude = grid5.ml$Lon),
    data.frame(Latitude = plot.preds$Lat, S = plot.preds$S, 
               Grain = "plots",  Longitude = plot.preds$Lon)
  )

# plot the latitudinal gradients
LG.plot <- ggplot(LG.data, aes(x=Latitude, y=S)) +
            geom_vline(xintercept = 0, size=.2) +
            geom_vline(xintercept = 23.5, size=.2) +
            geom_vline(xintercept = -23.5, size=.2) +
            geom_point(aes(shape=Grain), alpha=0.3) +
            geom_smooth(colour="red", aes(linetype=Grain), method="loess", span=0.3) +
            scale_y_log10() +
            theme_bw() +
            scale_shape(solid = FALSE) +
            coord_flip() 

# write to file
png(file="../Figures/latitudinal_gradient.png", width=1500, height=1200, res=250)
LG.plot
dev.off()




################################################################################
# RELATIONSHIP BETWEEN BETA DIVERSITY AND ELEVATION SPAN
################################################################################

DAT <- plot.preds@data[is.na(plot.preds@data$beta) == FALSE, ]
DAT <- DAT[is.na(DAT$ALT_DIF_grid) == FALSE, ]
DAT <- DAT[is.na(DAT$ALT_DIF) == FALSE, ]
m1 <- lm(log10(beta)~ poly(ALT_DIF_grid,2) + poly(ALT_DIF, 2), 
          data=DAT,
          na.action=na.omit)
par(mfrow=c(1,2))
termplot(m1, se=T)

ggplot(data=plot.preds@data, aes(x=ALT_DIF_grid, y=beta)) + 
  geom_point(shape=1) + 
  geom_smooth(method="lm", formula= y ~ poly(x, 2)) +
  scale_y_continuous(trans="log10") +
  xlab("Elavation span within hexagon [m]") +
  ylab(expression(gamma / alpha)) +
  theme_bw()





















