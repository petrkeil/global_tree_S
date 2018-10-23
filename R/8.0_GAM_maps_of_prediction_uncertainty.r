################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: April 26 2018
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
load("../Models/brms_SMOOTH.RData")

################################################################################
### Predictions in hexagons

# predict S from the model SMOOTH
grid.pred.S.brm <- data.frame(predict(brm.SMOOTH, 
                              newdata = grid5.dat,
                              probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))

names(grid.pred.S.brm)[1] <- "S"
grid.pred.S.brm <- data.frame(grid.pred.S.brm, 
                              ratio.95 = grid.pred.S.brm$Q97.5 / grid.pred.S.brm$Q2.5,
                              std.diff = (grid.pred.S.brm$Q97.5 - grid.pred.S.brm$Q2.5) / grid.pred.S.brm$Q50)
  

# merge with the original grid
grid5@data <- data.frame(grid5@data, grid.pred.S.brm)
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
plot.pred.S.brm <- predict(brm.SMOOTH, 
                           newdata = pts, 
                           type="response",
                           probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
plot.pred.S.brm <- data.frame(plot.pred.S.brm)

names(plot.pred.S.brm)[1] <- "S"
plot.pred.S.brm <- data.frame(plot.pred.S.brm, 
                              ratio.95 = plot.pred.S.brm$Q97.5 / plot.pred.S.brm$Q2.5,
                              std.diff = (plot.pred.S.brm$Q97.5 - plot.pred.S.brm$Q2.5) / plot.pred.S.brm$Q50)

# put all together
plot.preds <- data.frame(pts, 
                         plot.pred.S.brm)

# remove predictions of S < 0.8 (an arbitrarily selected threshold)
plot.preds$S[plot.preds$S < 0.8] <- NA

plot.preds <- plot.preds[rowSums(is.na(plot.preds)) == 0,]

# put predictions to a spatial object
plot.preds <- SpatialPointsDataFrame(coords = data.frame(plot.preds$Lon, plot.preds$Lat), 
                                     data = plot.preds, 
                                     proj4string = CRS(WGS84))






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


# Set the minimalist theme for the plotting
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
                     panel.grid.minor=element_blank(),plot.background=element_blank(),
                     plot.title = element_text(face=quote(bold)))

# ------------------------------------------------------------------------------
# plot the standardized difference

plot.gr.diff <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=std.diff)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_distiller(palette = "Spectral") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(paste((q[97.5] - q[2.5]) / q[50])), title="k") +
    theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))
  
plot.gr.diff


plot.pl.diff <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=std.diff))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_distiller(palette = "Spectral") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(paste((q[97.5] - q[2.5]) / q[50])), title="l") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))

plot.pl.diff

# ------------------------------------------------------------------------------



plot.gr.high <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=Q97.5)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_distiller(palette = "Spectral", 
                       limits=c(1,max(grid5.mlf$X97.5.ile)),
                       trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(q[97.5]), title="a") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))
plot.gr.high

# predicted S in hexagons
plot.gr.S <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=Q50)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_distiller(palette = "Spectral",
                       limits=c(1,max(grid5.mlf$X97.5.ile)),
                       trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(q[50]), title="c") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))
plot.gr.S

plot.gr.low <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=Q2.5)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_distiller(palette = "Spectral", 
                       limits=c(1,max(grid5.mlf$X97.5.ile)),
                       trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(q[2.5]), title="e") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))
#plot.gr.low



# ------------------------------------------------------------------------------

plot.pl.high <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=X97.5.ile))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_distiller(palette = "Spectral", 
                         name=expression(S[plot]),
                         limits=c(1,max(plot.preds.ml$Q97.5)),
                         trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(q[97.5]), title="b") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))

#plot.pl.high

plot.pl.S <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=Q50))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_distiller(palette = "Spectral", 
                         name=expression(S[plot]),
                         limits=c(1,max(plot.preds.ml$X97.5.ile)),
                         trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(q[50]), title="d") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))

#plot.pl.S 

plot.pl.low <- ggplot(MAINL, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(colour=NA, fill="white", size=.2) + 
  geom_point(data=plot.preds.ml, size=0.01,
             aes(x=X, y=Y, group=NULL, colour=Q2.5))  +
  geom_polygon(colour="black", fill=NA, size=.2) + 
  scale_colour_distiller(palette = "Spectral", 
                         name=expression(S[plot]),
                         limits=c(1,max(plot.preds.ml$X97.5.ile)),
                         trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  labs(subtitle = expression(q[2.5]), title="f") +
  theme_minimal() + blank.theme + theme(plot.title = element_text(face=quote(bold)))

#plot.pl.low


# ------------------------------------------------------------------------------
# more uncertainty plots
grid.rank <- grid5@data[order(grid5$Q50), ]

rank.gr.lin <- ggplot(data = grid.rank, aes(x = order(Q50), y = Q50)) + 
  geom_linerange(aes(x = order(S), ymin = Q2.5, ymax = Q97.5), 
                 #size = 1, 
                 colour="lightgrey") +
  geom_linerange(aes(x = order(S), ymin = Q25, ymax = Q75), 
                 #size = 1, 
                 colour="darkgrey") +
  geom_point(size = 0.5) + theme_bw() +
  xlab("Rank") + ylab("Predicted S")  + ggtitle("g") + 
  theme(plot.title = element_text(face=quote(bold)))
rank.gr.lin


rank.gr.log <- ggplot(data = grid.rank, aes(x = order(Q50), y = Q50)) + 
               geom_linerange(aes(x = order(S), ymin = Q2.5, ymax = Q97.5), 
                              #size = 1, 
                              colour="lightgrey") +
               geom_linerange(aes(x = order(S), ymin = Q25, ymax = Q75), 
                              #size = 1, 
                              colour="darkgrey") +
               scale_y_continuous(trans = "log10", 
                                   # limits=c(0.1, max(grid5$X97.5.ile)),
                                  labels = c(1,10,100,1000, 10000),
                                  breaks = c(1,10,100,1000, 10000)) +
               xlab("Rank") + ylab("Predicted S") + ggtitle("h") +
               geom_point(size = 0.5) + theme_bw() + 
               theme(plot.title = element_text(face=quote(bold)))
rank.gr.log

               
plot.rank <- plot.preds.ml[order(plot.preds.ml$Q50), ]


rank.pl.lin <- ggplot(data = plot.rank, aes(x = order(Q50), y = Q50)) + 
  geom_linerange(aes(x = order(S), ymin = Q2.5, ymax = Q97.5), 
                 #size = 1, 
                 colour="lightgrey") +
  geom_linerange(aes(x = order(S), ymin = Q25, ymax = Q75), 
                 #size = 1, 
                 colour="darkgrey") +
  geom_point(size = 0.5) + theme_bw() +
  xlab("Rank") + ylab("Predicted S") + ggtitle("i") + 
  theme(plot.title = element_text(face=quote(bold)))
rank.pl.lin

rank.pl.log <- ggplot(data = plot.rank, aes(x = order(Q50), y = Q50)) + 
               geom_linerange(aes(x = order(S), ymin = Q2.5, ymax = Q97.5), 
                              #size = 1,
                              colour="lightgrey") +
               geom_linerange(aes(x = order(S), ymin = Q25, ymax = Q75), 
                              #size = 1, 
                              colour="darkgrey") +
               scale_y_continuous(trans = "log10", 
                                  # limits=c(0.1, max(grid5$X97.5.ile)),
                                  labels = c(1,10,100,1000, 10000),
                                  breaks = c(1,10,100,1000, 10000)) +
               xlab("Rank") + ylab("Predicted S") +
               geom_point(size = 0.5) + theme_bw() + ggtitle("j") + 
  theme(plot.title = element_text(face=quote(bold)))
rank.pl.log



# ------------------------------------------------------------------------------

# write to file

lay <- matrix(nrow=5, ncol=4, byrow = TRUE,
              c(1,1,2,2,
                3,3,4,4,
                5,5,6,6,
                7,8,9,10,
                11, 11, 12, 12))

tiff("../Figures/maps_of_prediction_uncertainty.tif", width=4000, height=5400, res=350,
     compression = "lzw")
  grid.arrange(plot.gr.high, plot.pl.high,
               plot.gr.S, plot.pl.S,
               plot.gr.low, plot.pl.low,
               rank.gr.lin, rank.gr.log, rank.pl.lin, rank.pl.log, 
               plot.gr.diff, plot.pl.diff,
               layout_matrix = lay,
               heights = c(1,1,1,0.8, 1)) 
dev.off()






################################################################################
# LATITUDINAL GRADIENTS OF THE PREDICTIONS - WITH UNCERTAINTY BARS
################################################################################

# data for latitudinal gradient plots
LG.data <-
  rbind(
    data.frame(Latitude = grid5.ml@data$Lat, 
               Q2.5 = grid5.ml$Q2.5,
               Q25 = grid5.ml$Q25,
               Q50 = grid5.ml$Q50, 
               Q75 = grid5.ml$Q75,
               Q97.5 = grid5.ml$Q97.5,
               Grain = "hexagons"),
    data.frame(Latitude = plot.preds$Lat, 
               Q2.5 = plot.preds$Q2.5,
               Q25 = plot.preds$Q25,
               Q50 = plot.preds$Q50, 
               Q75 = plot.preds$Q75,
               Q97.5 = plot.preds$Q97.5,
               Grain = "plots")
  )

# plot the latitudinal gradients
LG.plot <- ggplot(LG.data, aes(x=Latitude, y=Q50)) +
  geom_vline(xintercept = 0, size=.2) +
  geom_vline(xintercept = 23.5, size=.2) +
  geom_vline(xintercept = -23.5, size=.2) +
  #geom_linerange(aes(ymin = Q2.5, ymax = Q97.5, colour = Grain), alpha = 0.1) +
  #geom_linerange(aes(ymin = Q25, ymax = Q75, colour = Grain), size=1, alpha = 0.3) +
  geom_point(aes(colour=Grain), shape=1) +
  #geom_point(shape=1) +
  geom_smooth(colour = "black", method="loess", span=0.3, se=FALSE) +
  scale_y_log10() +
  theme_bw() +
  ylab("Predicted S") +
  facet_grid(.~Grain) +
  scale_shape(solid = FALSE)  
LG.plot



# write to file
png(file="../Figures/latitudinal_gradient.png", 
    width=2500, height=1200, res=240)
LG.plot
dev.off()

###########################################################################################
# CALCULATING CREDIBLE INTERVALS AROUND THE MEAN OF THE LATITUDINAL GRADIENT


# predict S from the model SMOOTH in the grid
grid.S <- data.frame(predict(brm.SMOOTH, 
                     newdata = grid5.dat,
                     summary = TRUE))

grid.mcmc <- data.frame(predict(brm.SMOOTH, 
                        newdata = grid5.dat,
                        summary = FALSE))

# non-NA cells with S > 1
good.cells <- grid.S$Estimate > 1
good.cells[is.na(good.cells)] <- FALSE

# cells with large enough land area
large.cells <-  grid5@data$LandArea / grid5@data$CellArea > 0.5
large.cells[is.na(large.cells)] <- FALSE

# subset the data
grid.mcmc.sub <- grid.mcmc[,(good.cells * large.cells) == 1]
grid.dat.sub <-  grid5.dat[(good.cells * large.cells) == 1,]

# fit a smoother to each iteration of the MCMC algorithm
res.grid <- matrix(nrow=nrow(grid.dat.sub), ncol=nrow(grid.mcmc.sub))
for(i in 1:nrow(grid.mcmc.sub))
{
  print(i)
  S.i <- unlist(c(grid.mcmc.sub[i,]))
  m.i <- gam(S.i ~ s(grid.dat.sub$Lat), family = "nb")
  #plot(grid.dat.sub$Lat, grid.mcmc.sub[i,])
  #points(grid.dat.sub$Lat, predict(m.i, type="response"), col = 'red')
  res.grid[,i] <- predict(m.i, type="response")
}

# summarize the fitted smoothers using quantiles
smoother.grid <- data.frame(Lat = grid.dat.sub$Lat,
                            Q2.5 = apply(X = res.grid, MARGIN = 1, FUN = quantile, probs = c(0.025)),
                            Q50  = apply(X = res.grid, MARGIN = 1, FUN = quantile, probs = c(0.5)),
                            Q97.5= apply(X = res.grid, MARGIN = 1, FUN = quantile, probs = c(0.975)),
                            Grain = "hexagons")

smoother.grid <- smoother.grid[order(smoother.grid$Lat),]

plot(Q97.5 ~ Lat, data = smoother.grid)
points(Q2.5 ~ Lat, data = smoother.grid)
points(Q50 ~ Lat, data = smoother.grid)


# ----------------------------------------------------------------------------
### Predictions in 1 ha plots

# predict S in the plots from the SMOOTH model
plot.S <- predict(brm.SMOOTH, 
                  newdata = pts, 
                  type="response",
                  summary = TRUE)
plot.S <- data.frame(plot.S)

plot.mcmc <- predict(brm.SMOOTH, 
                     newdata = pts, 
                     type="response",
                     summary = FALSE)

# remove NAs and predictions of S < 0.8 (an arbitrarily selected threshold)
good.cells <- plot.S$Estimate > 0.8
not.na.cells <- rowSums(is.na(plot.S)) == 0
good.cells <- (good.cells * not.na.cells) == 1
good.cells[is.na(good.cells)] <- FALSE

plot.S.sub <- plot.S[good.cells,]
plot.mcmc.sub <- plot.mcmc[,good.cells]
plot.dat.sub <- pts[good.cells,]

# the following parallel loop takes about 10 min to run on my Intel i5 4-core laptop
library(foreach)
library(doMC)
registerDoMC(cores=4)
res.plot <- foreach(i = 1:nrow(plot.mcmc.sub), 
                    .combine = cbind, 
                    .packages = "mgcv",
                    .verbose = TRUE) %dopar% 
{
  S.i <- unlist(c(plot.mcmc.sub[i,]))
  m.i <- gam(S.i ~ s(plot.dat.sub$Lat), family = "nb")
  predict(m.i, type="response")
}

# summarize the fitted smoothers using quantiles
smoother.plot <- data.frame(Lat = plot.dat.sub$Lat,
                            Q2.5 = apply(X = res.plot, MARGIN = 1, FUN = quantile, probs = c(0.025)),
                            Q50  = apply(X = res.plot, MARGIN = 1, FUN = quantile, probs = c(0.5)),
                            Q97.5= apply(X = res.plot, MARGIN = 1, FUN = quantile, probs = c(0.975)),
                            Grain = "plots")

plot(Q97.5 ~ Lat, data = smoother.plot)
points(Q2.5 ~ Lat, data = smoother.plot)
points(Q50 ~ Lat, data = smoother.plot)

smoother.plot <- smoother.plot[order(smoother.plot$Lat),]

# bring the smoothers from the plot and the grid together
smoothers.both <- rbind(smoother.grid, smoother.plot)


# ------------------------------------------------------------------------------

# plot IT
LGU.plot <- ggplot(LG.data, aes(x=Latitude, y=Q50)) +
  geom_vline(xintercept = 0, size=.2) +
  geom_vline(xintercept = 23.5, size=.2) +
  geom_vline(xintercept = -23.5, size=.2) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5, colour = Grain), alpha = 0.1) +
  #geom_linerange(aes(ymin = Q25, ymax = Q75, colour = Grain), size=1, alpha = 0.3) +
  geom_point(aes(colour=Grain), shape=1) +
  #geom_point(shape=1) +
  geom_line(data = smoothers.both, aes(x=Lat, y=Q2.5), linetype = 2) +
  geom_line(data = smoothers.both, aes(x=Lat, y=Q50)) +
  geom_line(data = smoothers.both, aes(x=Lat, y=Q97.5), linetype = 2) +
  scale_y_log10() +
  theme_bw() +
  ylab("Predicted S") +
  facet_grid(.~Grain) +
  scale_shape(solid = FALSE)  
LGU.plot


# write to file
png(file="../Figures/latitudinal_gradient_with_uncertainty.png", 
    width=2500, height=1200, res=240)
LGU.plot
dev.off()




