################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Dec 7 2017
################################################################################

# Description: Here the workspace is cleaned, all the necessary libraries are 
# loaded, and the geographic projections are defined.

# ------------------------------------------------------------------------------

# clean everything
rm(list=ls())

# ------------------------------------------------------------------------------
# LIBRARIES

# statistical modelling
library(dummies)
library(MASS)
library(rstan)
library(mcmcplots)
library(brms)
library(randomForest)
library(ncf)
library(mgcv)
library(broom)
library(matrixStats)

# tables
library(xtable)
library(R2HTML)

# GIS libraries
library(rgdal)
library(raster)
library(rgeos)
library(sp)
library(maps)

# data manipulation libraries
library(dplyr)
library(plyr)
library(reshape)

# plotting libraries
library(RColorBrewer)
library(latticeExtra)
library(gridExtra)
library(ggplot2)
library(colorRamps)
library(GGally)
library(ggthemes)


# ------------------------------------------------------------------------------
# PROJECTIONS for GIS operations

WGS84 <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
LAMBERT <- "+proj=cea +lon_0=Central Meridian +lat_ts=Standard Parallel +x_0=False Easting +y_0=False Northing"
MOLLWEIDE <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ROBINSON <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Equidistant Azimuthal Projection
AED <- "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=100000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


# ------------------------------------------------------------------------------
# Blank theme for ggplot2

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








































