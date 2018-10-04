################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: April 26 2018
################################################################################

# Description: These scripts extract the environmental conditions for the set
# of regularly distributed 1 ha artificial plots, and for the large hexagons. 
# The raw environmental data are not provided in the GitHub repository, 
# since they are huge. Thus, if you want to replicate
# these specific parts of the analyses, you need to download the raw data
# from the repositories described in the Material and Methods section.


################################################################################

# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")

# ------------------------------------------------------------------------------

# load the mainland layer (i.e. everyting except for islands)
MAINL <- readOGR(dsn = "/media/pk33loci/Elements/GIS_data/Boundaries/GLOBAL_SHORELINE", 
                 layer = "main_landmasses")
proj4string(MAINL) <- WGS84

# ------------------------------------------------------------------------------
# LOAD THE HEXAGONAL GRID THAT WILL BE USED FOR PREDICTIONS

grid5 <- readOGR(dsn="../Data/GRIDS", layer = "hex5")
grid5 <- spTransform(grid5, CRSobj = WGS84)
gr.coords <- data.frame(coordinates(grid5))
names(gr.coords) <- c("Lon", "Lat")
grid5@data <- data.frame(grid5@data, gr.coords)


# read a global 10X10 km raster on annual temperature (will be used to create a grid)
# and convert it to a global raster of 10x10km cells with area of each cell
RAST <- raster("../Data/GRIDS/ANN_T10.tif")
RAST[is.na(RAST) == FALSE] <- 1
RAST.area <- raster::area(RAST)
RAST.area <- RAST.area * RAST
RAST.area[is.na(RAST.area)] <- 0
plot(RAST.area)

# the area of each hexagonal cell in km^2
grid5.totarea <- gArea(spTransform(grid5, CRSobj = LAMBERT), byid=TRUE) / 1000000

# calculate the area of land in km^2
grid5.area <- raster::extract(x=RAST.area, y=grid5, fun=sum)

# calculate area of mainland (non-island) in km
ML.rast <- rasterize(x = MAINL, y=RAST)
ML.rast[is.na(ML.rast) == FALSE] <- 1
ML.rast.area <- raster::area(ML.rast)
ML.rast.area <- ML.rast.area * ML.rast
ML.rast.area[is.na(ML.rast.area)] <- 0
plot(ML.rast.area)

grid5.ML.area <- raster::extract(x=ML.rast.area, y=grid5, fun=sum)

# put the data together
grid5@data <- data.frame(grid5@data, 
                         LandArea = grid5.area,
                         MainlArea = grid5.ML.area,
                         CellArea = grid5.totarea)


# ------------------------------------------------------------------------------
# LOAD THE POINTS THAT WILL BE USED FOR PREDICTIONS

pt.coords <- read.csv(file="../Data/GRIDS/Fine_points.csv")
pts <- SpatialPointsDataFrame(pt.coords, 
                              proj4string=CRS(WGS84),
                              data=data.frame(ptID = paste("pt", 1:nrow(pt.coords), sep="")))
pts@data <- data.frame(pts@data, pt.coords)


# ------------------------------------------------------------------------------
# FUNCTION THAT DOES THE ZONAL EXTRACTION
# It overlays the shapefile on a raster, and extracts raster values that lay 
# underneath each polygon.

# Arguments:
# SHAPEFILE - shapefile that is going to be extracted to
# RASTER - a raster object that is going to be extracted from
# fun - Function used to summarize the data.
# fact - How much do you want to aggregate the raster prior to the extraction?

env.extract <- function(SHAPEFILE, RASTER, fun, fact, na.threshold = FALSE)
{
  message("1 - Aggregating the raster.")
  if(fact > 1)
  {
    RASTER <- aggregate(RASTER, fact=fact, fun=fun)
  }
  
  if(na.threshold)
  {
    RASTER[RASTER > na.threshold] <- NA  
  }
  
  RASTER.px <- as(RASTER, "SpatialGridDataFrame")
  message("2 - Extracting the values")
  SUMMARY <- over(x=SHAPEFILE, y=RASTER.px, fn=fun, returnList=FALSE)
  SUMMARY <- round(SUMMARY, 2)
  return(SUMMARY[,1])
}


# ------------------------------------------------------------------------------

# Extracting the rasters, one by one. 

# IMPORTANT NOTE! The rasters do not come with these data (they are too large), 
# but their source is provided in the Methods section


###########################
## CROWTHER TREE DENSITY ##
###########################

TREE_DENS <- raster("/media/pk33loci/Elements/GIS_data/Crowther_TreeDensity/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif")
proj4string(TREE_DENS) <- WGS84

# caclulate pixel area
A <- raster::area(TREE_DENS)
A.pts <- raster::extract(x = A, y= SpatialPoints(coords=data.frame(pts$Lon, pts$Lat), 
                                                 proj4string=CRS(WGS84)))

# extract the hexagonal grid for predictions
grid5@data$TREE_DENS <- env.extract(grid5, TREE_DENS, fun=sum, fact=25)

# extract points for predictions
pts@data <- data.frame(pts@data, A = A.pts,
                       TREE_DENS=round(raster::extract(TREE_DENS, pts)))

# delete the raster to save memory
rm(TREE_DENS)


#########
## GPP ##
#########

GPP <- raster("/media/pk33loci/Elements/GIS_data/MOD17A3_1km/MOD17A3_Science_GPP_mean_00_15.tif")
proj4string(GPP) <- WGS84
#GPP[GPP > 50000] <- NA # this blows off the memory

# extract the hexagonal grid for predictions
grid5@data$GPP <- env.extract(grid5, GPP, fun=mean, fact=25, na.threshold = 50000)

# extract points for predictions
pts@data <- data.frame(pts@data, GPP=raster::extract(GPP, pts))
pts$GPP[pts$GPP > 50000] <- NA

rm(GPP)

########################
## EVAPOTRANSPIRATION ##
########################


ET <- raster("/media/pk33loci/Elements/GIS_data/Global_mean_evapotranspiration/w0010001.tif")
proj4string(ET) <- WGS84
# ET[ET > 50000] <- NA # this blows off the memory

# extract the hexagonal grid for predictions
grid5@data$ET <- env.extract(grid5, ET, fun=mean, fact=25, na.threshold=50000)

# extract points for predictions
pts@data <- data.frame(pts@data, ET=raster::extract(ET, pts))
pts$ET[pts$ET > 50000] <- NA

rm(ET)

############################
# MEAN ANNUAL TEMPERATURE - BIO 1
############################

ANN_T <- raster("/media/pk33loci/Elements/GIS_data/Bioclim_30s/bio_1.bil")
names(ANN_T) <- "ANN_T"

ANN_T10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio1.bil")

# extract the hexagonal grid for predictions
grid5@data$ANN_T <- env.extract(grid5, ANN_T10, fun=mean, fact=1)/10

# extract points for predictions
pts@data <- data.frame(pts@data, ANN_T=raster::extract(ANN_T, pts))
pts@data$ANN_T <- pts@data$ANN_T / 10

rm(ANN_T)

########################################
# MEAN TEMPERATURE OF WARMEST QUQRTER - BIO 10
########################################

WARM_T <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio10-19_30s_bil/bio_10.bil")
names(WARM_T) <- "WARM_T"

WARM_T10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio10.bil")/10

# extract the hexagonal grid for predictions
grid5@data$WARM_T <- env.extract(grid5, WARM_T10, fun=mean, fact=1)

# extract points for predictions
pts@data <- data.frame(pts@data, WARM_T=raster::extract(WARM_T, pts))
pts@data$WARM_T <- pts@data$WARM_T / 10

rm(WARM_T)

########################################
# ISOTHERMALITY - BIO 3
########################################

ISO_T <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio1-9_30s_bil/bio_3.bil")
names(ISO_T) <- "ISO_T"

ISO_T10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio3.bil")

# extract the hexagonal grid for predictions
grid5@data$ISO_T <- env.extract(grid5, ISO_T10, fun=mean, fact=1)

# extract points for predictions
pts@data <- data.frame(pts@data, ISO_T=raster::extract(ISO_T, pts))

rm(ISO_T)

########################################
# PRECIPITATION IN THE DRIEST QUARTER - BIO 17
########################################

MIN_P <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio10-19_30s_bil/bio_17.bil")
names(MIN_P) <- "MIN_P"

MIN_P10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio17.bil")

# extract the hexagonal grid for predictions
grid5@data$MIN_P <- env.extract(grid5, MIN_P10, fun=mean, fact=1)

# extract points for predictions
pts@data <- data.frame(pts@data, MIN_P=raster::extract(MIN_P, pts))

rm(MIN_P, MIN_P10)

########################################
# PRECIPITATION SEASONALITY - BIO 15
########################################

P_SEAS <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio10-19_30s_bil/bio_15.bil")
names(P_SEAS) <- "P_SEAS"

P_SEAS10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio15.bil")

# extract the hexagonal grid for predictions
grid5@data$P_SEAS <- env.extract(grid5, P_SEAS10, fun=mean, fact=1)

# extract points for predictions
pts@data <- data.frame(pts@data, P_SEAS=raster::extract(P_SEAS, pts))

rm(P_SEAS, P_SEAS10)

###############################
# TOPOGRAPHIC HETEROGENEITY  ##
###############################

MIN_ALT <- raster("/media/pk33loci/Elements/GIS_data/Elevation_GMTED2010_USGS/mi30_grd/mi30_grd")
MAX_ALT <- raster("/media/pk33loci/Elements/GIS_data/Elevation_GMTED2010_USGS/mx30_grd/mx30_grd")
# this takes about 4 minutes
ALT_DIF <- MAX_ALT - MIN_ALT

# extract the hexagonal grid for predictions
MIN_ALT_CN <- env.extract(grid5, MIN_ALT, fun=min, fact=25)
MAX_ALT_CN <- env.extract(grid5, MAX_ALT, fun=max, fact=25)
grid5@data$ALT_DIF <- MAX_ALT_CN - MIN_ALT_CN

# extract points for predictions
pts@data <- data.frame(pts@data, ALT_DIF=raster::extract(ALT_DIF, pts))

rm(MIN_ALT, MAX_ALT, ALT_DIF)

######################
# ISLAND/MAINLAND  ###
######################

## OLD APPROACH WHERE TRUE ISLANDS AND SHELF ILANDS WERE ALL LUMPED TOGETHER
 MAINL <- readOGR(dsn = "/media/pk33loci/Elements/GIS_data/Boundaries/GLOBAL_SHORELINE", 
                 layer = "main_landmasses")
 MAINL <- spTransform(MAINL,  CRSobj = WGS84)
# extract points for predictions
 CONTS <- over(x=pts, y=MAINL)
 is.island <- is.na(CONTS$continent == "<NA>")*1
 pts@data$ISLAND <- is.island
 plot(pts, col=pts@data$ISLAND+1); plot(MAINL, add=T)
# calculate the ISLAND status of the hexagonal cells
 hexISL <- 1 - grid5$MainlArea/grid5$LandArea
 hexISL <- ifelse(hexISL > 0.9, 1, 0)
 grid5@data <- data.frame(grid5@data, ISLAND=hexISL)


# UPDATED APPROACH IN WHICH SHELF ISLANDS ARE TREATED
# AS EFFECTIVELY MAINLANDS

ISLAND <- raster("/media/pk33loci/Elements/GIS_data/ISLANDNESS/rasters/ISLAND_clean.tif")
ALL.LAND <- raster("/media/pk33loci/Elements/GIS_data/ISLANDNESS/rasters/LAND_clean.tif")
MAINLAND <- ALL.LAND - ISLAND

is.island.plots <- raster::extract(x = ISLAND, y = pts)
is.island.plots <- ifelse(is.island.plots == 1, "island", "mainland")
pts@data$INSULARITY <- is.island.plots

is.mainl.hex <- raster::extract(x = MAINLAND, y = grid5, fun = max)
is.isl.hex <- as.vector((is.mainl.hex == 0) * 1)
grid5@data$INSULARITY <- is.isl.hex


################################################################################
# EXPORT THE DATA

write.csv(pts@data, file="../Data/GRIDS/Fine_points_with_environment.csv", row.names = FALSE)

writeOGR(obj = grid5, dsn = "../Data/GRIDS", 
         layer = "hex5_with_environment", 
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)






