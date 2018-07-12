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
