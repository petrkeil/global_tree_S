################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Oct 23 2018
################################################################################

# Description: These are scripts that generate the global layer of the 1-ha plots,
# and the hexagonal prediction grid.

################################################################################

# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")
library(dggridR)
library(dplyr)
library(mapproj)
library(rgdal)
library(sp)
library(raster)
library(gridExtra)
library(rgeos)



# read the COUNTRY data (as a shapefile with data in the attribute table)
COUNTR.shp <- readOGR(dsn="../Data/COUNTRIES", layer = "COUNTRIES")
COUNTR.shp <- spTransform(COUNTR.shp, CRSobj = WGS84)

LAND <- readOGR(dsn="../Data/COUNTRIES", layer = "GSHHS_i_L1")
LAND <- spTransform(LAND, CRSobj = WGS84)




# GENERATE THE POINTS ----------------------------------------------------------
dggs <- dgconstruct(res=8, metric=FALSE, resround='down')
#grid <- dgearthgrid(dggs, frame=TRUE, wrapcells=TRUE)
centers <- dgSEQNUM_to_GEO(dggs, in_seqnum=1:65612)

pts <- SpatialPoints(coords=data.frame(centers$lon_deg, 
                                       centers$lat_deg), 
                     proj4string = CRS(WGS84))

# subset only to data that lay on the mainland
points.land <- sp::over(pts, LAND)
points.land <- is.na(points.land$level) == FALSE
pts.land <- pts[points.land]

# plot the result
plot(pts.land, cex=0.1, col="red")
plot(LAND, add=TRUE)

# export the points to a .csv file
pts.latlon <- data.frame(coordinates(pts.land))
names(pts.latlon) <- c("Lon","Lat")
write.csv(pts.latlon, file="../Data/GRIDS/Fine_points.csv", row.names = FALSE)




# GENERATE THE HEXAGONAL COARSE GRID - RES 5 -----------------------------------
# area of a grid cell for res=5 is 209,903 km^2

dggs <- dgconstruct(res=5, metric=FALSE, resround='down')
grid <- dgearthgrid(dggs, frame=FALSE, wrapcells=TRUE)
grid <- spTransform(grid, WGS84)
grid <- SpatialPolygonsDataFrame(grid, data=data.frame(id=1:length(grid)))

# overlay over the global countries
x <- over(grid, COUNTR.shp)
good.cells <- is.na(x[,1]) == FALSE
grid.land <- grid[good.cells,]
plot(grid.land)


writeOGR(obj=grid.land, dsn = "../Data/GRIDS" , layer="hex5", driver="ESRI Shapefile")
# NOTE: This file will be manually adjusted in QGIS to delete the grid cells 
# with both negative and positive longitude values.





















