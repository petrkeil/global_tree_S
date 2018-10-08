# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")


################################################################################
### Read, transform and scale the data

# read the data
grid5 <- readOGR(dsn = "../Data/GRIDS", layer = "hex5_with_environment")
grid5 <- spTransform(x = grid5, CRSobj = WGS84)
grid5@data <- data.frame(grid5@data, ELONGATION = 1)

# -----------------------------------------

# tree density at the grid level
grid5$Tree_dens <- (grid5$TREE_DENS + 1) / grid5$LandArea
grid5@data <- data.frame(grid5@data, min_DBH = 0, DAT_TYPE = "Country")

# -----------------------------------------

grid5.dat <- dplyr::select(grid5@data, Area_km = LandArea, Tree_dens, min_DBH,
                           GPP, ANN_T, ISO_T, MIN_P, P_SEAS, ALT_DIF,
                           INSULARITY, ELONGATION, Lat, Lon, DAT_TYPE) %>%
  mutate(Area_km = log(Area_km), Tree_dens=log(Tree_dens))

# get the scaling constants that were used to scale the raw plot and country data:
scal.tab <- read.csv("scale_tab.csv")
scal.tab <- scal.tab[scal.tab$var %in% c("ET","WARM_T") == FALSE,]

# scale the grid data in the same way as the original data

grid5.dat[,1:9] <- scale(grid5.dat[,1:9],
                         center = scal.tab$centr, 
                         scale = scal.tab$scale)

################################################################################
# 1. The African data

RB <- read.csv("../Data/VALIDATION/RAINBIO_Africa_SN_grid5.csv")
EU <- read.csv("../Data/VALIDATION/EUForest_SN_grid5.csv")
BI <- read.csv("../Data/VALIDATION/BIEN_woody_grid5.csv")

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

# merge with the original grid

for.plot <- data.frame(grid5@data, grid.pred.S.brm)


################################################################################
### External data for the hexagon

RAINBIO <- left_join(x = for.plot, y = RB, by="id")
RAINBIO <- RAINBIO[is.na(RAINBIO$N) == FALSE, ]
RAINBIO <- RAINBIO[RAINBIO$N > 4000,]
RAINBIO <- data.frame(RAINBIO, Dataset = "Africa - RAINBIO")

EUForest <- left_join(x = for.plot, y = EU, by="id")
EUForest <- EUForest[is.na(EUForest$N) == FALSE, ]
EUForest <- EUForest[EUForest$N > 1000,]
EUForest <- data.frame(EUForest, Dataset = "Europe - EUForest")

BIEN <- left_join(x = for.plot, y = BI, by="id")
BIEN <- BIEN[is.na(BIEN$N) == FALSE, ]
BIEN <- BIEN[BIEN$N > 10000,]
BIEN <- data.frame(BIEN, Dataset = "New World - BIEN")


all.external <- rbind(RAINBIO, EUForest, BIEN)


p.valid <- ggplot(data = all.external, aes(x=S.y, y = X50.ile)) + 
            geom_linerange(aes(x=S.y, ymin = X2.5.ile, ymax = X97.5.ile, colour = Dataset), 
                           alpha = 0.5) +
            geom_linerange(aes(x=S.y, ymin = X25.ile, ymax = X75.ile, colour = Dataset), 
                           alpha = 0.5, size=1.5) +
            geom_point(aes(colour=Dataset)) +
              scale_y_log10(limits = c(10, max(all.external$X97.5.ile)), 
                            breaks = c(10,100, 1000, 10000))  +  
              scale_x_log10(limits = c(10, max(all.external$X97.5.ile)),
                            breaks = c(10,100, 1000, 10000)) +
            geom_abline(intercept = 0, slope=1) +
             # xlab("Observed S") + ylab("Predicted S") +
              labs(title = "b", x = "Observed S", y = "Predicted S (model SMOOTH)") +
              theme_bw() + theme(plot.title = element_text(face=quote(bold)),
                                 legend.position = c(0.24,0.82))
p.valid


# ------------------------------------------------------------------------------
# Plotting the external validation maps

grid5.ex <- grid5
grid5.ex@data <- left_join(x = grid5.ex@data, y = all.external, by = "id")
grid5.ex <- grid5.ex[is.na(grid5.ex@data$S.x) == FALSE, ]
grid5.ex$id <- as.character(grid5.ex$id)


grid5.ml <- spTransform(grid5.ex, CRSobj=MOLLWEIDE)
grid5.mlf <- tidy(grid5.ml, region="id")
grid5.mlf <- left_join(x=grid5.mlf, y=grid5.ex@data, by="id")

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
                     # legend.title = element_blank(),
                     legend.title.align = 0,
                     #plot.title = element_text(hjust = 0),
                     plot.subtitle = element_text(vjust=-3),
                     panel.background=element_blank(),
                     panel.border=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),plot.background=element_blank())



plot.gr.S <- ggplot(grid5.mlf, aes(long, lat, group=group)) +
  geom_polygon(data=LINES,  aes(long, lat, group=group), 
               colour="darkgrey", size=0.2) +
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill="white", colour=NA, size=.2) +
  geom_polygon(aes(fill=S.y)) + 
  geom_polygon(data=MAINL, aes(long, lat, group=group), 
               fill=NA, colour="black", size=.2) +
  scale_fill_distiller(palette = "Spectral", 
                       name=expression(S), 
                       #limits=c(1,5000),
                       trans="log10") +
  scale_x_continuous(limits = c(-12000000, 16000000)) +
  scale_y_continuous(limits = c(-6.4e+06, 8.8e+06)) +
  xlab("") + ylab("") +
  #ggtitle("A") +  
  labs(title = "a") +
  blank.theme + theme(plot.title = element_text(face=quote(bold)))
plot.gr.S

tiff("../Figures/Fig_external_validation.tif", width=4000, height=1400, res=350,
     compression = "lzw")
grid.arrange(plot.gr.S, p.valid, ncol=2, nrow = 1, widths = c(0.65, 0.35))
dev.off()

















