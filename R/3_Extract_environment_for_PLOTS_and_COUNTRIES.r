################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Oct 23 2018
################################################################################

# Description: These scripts extract the environmental conditions for the 
# COUNTRY and PLOT spatial units. The raw environmental data are not provided 
# in the GitHub repository, since they are huge. Thus, if you want to replicate
# these specific parts of the analyses, you need to download the raw data
# from the repositories described in the Material and Methods section.

################################################################################


# clean the workspace and load the libraries
source("0_libraries_functions_settings.r")

# ------------------------------------------------------------------------------


# read the raw PLOT data
PLOTS <- read.csv("../Data/PLOTS/PLOTS_all_data_raw_tab_delimited.txt",
                  sep="\t", header=TRUE)

# convert the coordinates to 'sp' object 
plots <- SpatialPointsDataFrame(coords=data.frame(PLOTS$Lon, PLOTS$Lat),
                                data = PLOTS,
                                proj4string=CRS(WGS84))

# read the COUNTRY data (as a shapefile with data in the attribute table)
COUNTR.shp <- readOGR(dsn="../Data/COUNTRIES", layer = "COUNTRIES")
COUNTR.shp <- spTransform(x = COUNTR.shp, CRSobj = WGS84)


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

# extract plots
plots@data <- data.frame(plots@data, TREE_DENS=round(extract(TREE_DENS, plots)))

# extract countries
COUNTR.shp@data$TREE_DENS <- env.extract(COUNTR.shp, TREE_DENS, fun=sum, fact=25)

rm(TREE_DENS)

#########
## GPP ##
#########

GPP <- raster("/media/pk33loci/Elements/GIS_data/MOD17A3_1km/MOD17A3_Science_GPP_mean_00_15.tif")
proj4string(GPP) <- WGS84
#GPP[GPP > 50000] <- NA # this blows off the memory

# extract plots
plots@data <- data.frame(plots@data, GPP=raster::extract(GPP, plots))
plots$GPP[plots$GPP > 50000] <- NA

# extract countries
COUNTR.shp@data$GPP <- env.extract(COUNTR.shp, GPP, fun=mean, fact=25, na.threshold = 50000)

rm(GPP)

########################
## EVAPOTRANSPIRATION ##
########################


# extract plots
ET <- raster("/media/pk33loci/Elements/GIS_data/Global_mean_evapotranspiration/w0010001.tif")
proj4string(ET) <- WGS84
# ET[ET > 50000] <- NA # this blows off the memory

# extract plots
plots@data <- data.frame(plots@data, ET=extract(ET, plots))
plots$ET[plots$ET > 50000] <- NA

# extract countries
COUNTR.shp@data$ET <- env.extract(COUNTR.shp, ET, fun=mean, fact=25, na.threshold=50000)

rm(ET)

############################
# MEAN ANNUAL TEMPERATURE - BIO 1
############################

# extract plots
ANN_T <- raster("/media/pk33loci/Elements/GIS_data/Bioclim_30s/bio_1.bil")
names(ANN_T) <- "ANN_T"

# extract plots
plots@data <- data.frame(plots@data, ANN_T=extract(ANN_T, plots))
plots@data$ANN_T <- plots@data$ANN_T / 10

# extract countries
ANN_T10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio1.bil")
COUNTR.shp@data$ANN_T <- env.extract(COUNTR.shp, ANN_T10, fun=mean, fact=1)/10

rm(ANN_T)

########################################
# MEAN TEMPERATURE OF WARMEST QUQRTER - BIO 10
########################################

# extract plots
WARM_T <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio10-19_30s_bil/bio_10.bil")
names(WARM_T) <- "WARM_T"

# extract plots
plots@data <- data.frame(plots@data, WARM_T=extract(WARM_T, plots))
plots@data$WARM_T <- plots@data$WARM_T / 10

# extract countries
WARM_T10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio10.bil")/10
COUNTR.shp@data$WARM_T <- env.extract(COUNTR.shp, WARM_T10, fun=mean, fact=1)

rm(WARM_T)

########################################
# ISOTHERMALITY - BIO 3
########################################

ISO_T <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio1-9_30s_bil/bio_3.bil")
names(ISO_T) <- "ISO_T"
plots@data <- data.frame(plots@data, ISO_T=extract(ISO_T, plots))

# extract countries
ISO_T10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio3.bil")
COUNTR.shp@data$ISO_T <- env.extract(COUNTR.shp, ISO_T10, fun=mean, fact=1)

rm(ISO_T)

########################################
# PRECIPITATION IN THE DRIEST QUARTER - BIO 17
########################################

MIN_P <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio10-19_30s_bil/bio_17.bil")
names(MIN_P) <- "MIN_P"

# extract plots
plots@data <- data.frame(plots@data, MIN_P=extract(MIN_P, plots))

# extract countries
MIN_P10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio17.bil")
COUNTR.shp@data$MIN_P <- env.extract(COUNTR.shp, MIN_P10, fun=mean, fact=1)

rm(MIN_P, MIN_P10)

########################################
# PRECIPITATION SEASONALITY - BIO 15
########################################

P_SEAS <- raster("/media/pk33loci/Elements/GIS_data/worldclim_30s_archives/bio10-19_30s_bil/bio_15.bil")
names(P_SEAS) <- "P_SEAS"

# extract plots
plots@data <- data.frame(plots@data, P_SEAS=extract(P_SEAS, plots))

# extract countries
P_SEAS10 <- raster("/media/pk33loci/Elements/GIS_data/worldclim10/bio15.bil")
COUNTR.shp@data$P_SEAS <- env.extract(COUNTR.shp, P_SEAS10, fun=mean, fact=1)

rm(P_SEAS, P_SEAS10)

###############################
# TOPOGRAPHIC HETEROGENEITY  ##
###############################


# extract plots
MIN_ALT <- raster("/media/pk33loci/Elements/GIS_data/Elevation_GMTED2010_USGS/mi30_grd/mi30_grd")
MAX_ALT <- raster("/media/pk33loci/Elements/GIS_data/Elevation_GMTED2010_USGS/mx30_grd/mx30_grd")
# this takes about 4 minutes
ALT_DIF <- MAX_ALT - MIN_ALT
plots@data <- data.frame(plots@data, ALT_DIF=raster::extract(ALT_DIF, plots))


# extract countries
MIN_ALT_CN <- env.extract(COUNTR.shp, MIN_ALT, fun=min, fact=25)
MAX_ALT_CN <- env.extract(COUNTR.shp, MAX_ALT, fun=max, fact=25)
COUNTR.shp@data$ALT_DIF <- MAX_ALT_CN - MIN_ALT_CN

rm(MIN_ALT, MAX_ALT, ALT_DIF)

######################
# INSULARITY       ###
######################

## OLD APPROACH WHERE TRUE ISLANDS AND SHELF ILANDS WERE ALL LUMPED TOGETHER
 
 ## extract plots
 MAINL <- readOGR(dsn = "/media/pk33loci/Elements/GIS_data/Boundaries/GLOBAL_SHORELINE", 
                 layer = "main_landmasses")
 MAINL <- spTransform(MAINL,  CRSobj = WGS84)
 CONTS <- over(x=plots, y=MAINL)
 is.island <- is.na(CONTS$continent == "<NA>")*1
 plots@data$ISL_LS <- ifelse(is.island == 1, "island", "mainland")
 
 
 ## extract countries
 CONTS <- over(SpatialPoints(coordinates(COUNTR.shp), proj4string=CRS(WGS84)), y=MAINL)
 is.island <- is.na(CONTS$continent == "<NA>")*1
 COUNTR.shp@data$ISL_LS <- ifelse(is.island == 1, "island", "mainland")

 
# ----------------------

# UPDATED APPROACH, IN WHICH SHELF ISLANDS ARE TREATED
# AS EFFECTIVELY MAINLANDS

ISLAND <- raster("/media/pk33loci/Elements/GIS_data/ISLANDNESS/rasters/ISLAND_clean.tif")
ALL.LAND <- raster("/media/pk33loci/Elements/GIS_data/ISLANDNESS/rasters/LAND_clean.tif")
LAND_AND_SHELF <- ALL.LAND - ISLAND

# plots
is.mainland.plots <- raster::extract(x = LAND_AND_SHELF , y = plots)
is.island.plots <- ifelse(is.mainland.plots == 1, "mainland", "island")
plots@data$ISL_ST <- is.island.plots

# manually extracted values for COUNTRIES
is.isl.countr <- read.csv("../Data/COUNTRIES/A_Insularity_COUNTRIES.csv")[,c("NAME","ISL_ST")]
COUNTR.shp@data <- dplyr::left_join(COUNTR.shp@data, is.isl.countr, by="NAME")


# YET ANOTHER APPROACH, IN WHICH SHELF ISLANDS ARE TREATED
# AS EFFECTIVELY MAINLANDS, AND DISJUNCT COUNTRIES ARE ALSO FLAGGED

# plots
is.mainland.plots <- raster::extract(x = LAND_AND_SHELF , y = plots)
is.island.plots <- ifelse(is.mainland.plots == 1, "mainland", "island")
plots@data$ISL_DIS <- is.island.plots

# manually extracted values for COUNTRIES
is.isl.countr <- read.csv("../Data/COUNTRIES/A_Insularity_COUNTRIES.csv")[,c("NAME","ISL_DIS")]
COUNTR.shp@data <- dplyr::left_join(COUNTR.shp@data, is.isl.countr, by="NAME")



########################
# ELONGATION         ###
########################

source("0_roundness_or_elongation_metrics_functions.r")

elong.cntr <- vector()

for(i in 1:nrow(COUNTR.shp))
{
  cat("*")
  pol <- COUNTR.shp[i,]  
  elong.cntr[i] <- round(elongation.sample(pol), 3)
}
COUNTR.shp@data <- data.frame(COUNTR.shp@data, ELONG = elong.cntr)


########################
# CLASSICAL REGIONS  ###
########################


# extract plots
REGS <- readOGR(dsn="/media/pk33loci/Elements/GIS_data/Ecoregions_Olsen/terr-ecoregions-TNC", layer="tnc_terr_ecoregions")
proj4string(REGS) <- WGS84
REGS.OVER <- over(x=plots, y=REGS)
plots@data <- data.frame(plots@data, HABITAT=REGS.OVER$WWF_MHTNAM,
                         REALM = REGS.OVER$WWF_REALM2)

# extract countries
REGS.OVER <- over(SpatialPoints(coordinates(COUNTR.shp), proj4string=CRS(WGS84)), y=REGS)

COUNTR.shp@data <- data.frame(COUNTR.shp@data, 
                              HABITAT=REGS.OVER$WWF_MHTNAM,
                              REALM = REGS.OVER$WWF_REALM2)



########################
# PETR's REGIONS     ###
########################

# extract plots
PK_REGS <- read.csv("../Data/COUNTRIES/A_Realm_classification_COUNTRIES.csv")
COUNTR.shp <- sp::merge(x=COUNTR.shp, y=PK_REGS, by="NAME", all.x=TRUE)

# extract countries
PK_REGS.OVER <- over(x=plots, y=COUNTR.shp)
plots@data <- data.frame(plots@data, REALM_PK = PK_REGS.OVER$REALM_3)



################################################################################
# EXPORT THE DATA

write.csv(plots@data, file="../Data/PLOTS/PLOTS_with_environment.csv", row.names = FALSE)

writeOGR(obj = COUNTR.shp, dsn = "../Data/COUNTRIES", 
         layer = "COUNTRIES_with_environment", 
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

















