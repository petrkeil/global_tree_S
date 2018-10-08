################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: April 26 2018
################################################################################

# Description: This is a set of scripts which take the raw data from the 
# ../Data/PLOTS and ../Data/COUNTRIES folders, and put them together in the 
# ../Data folder, from which they will be used for the main analysis.
# Also, some additional cleaning and massaging of the data is done here.

################################################################################

source("0_libraries_functions_settings.r")

PLOTS <- read.csv("../Data/PLOTS/PLOTS_with_environment.csv")
CNTRS <- readOGR(dsn = "../Data/COUNTRIES", layer = "COUNTRIES_with_environment")


################################################################################
# PROCESSING PLOTS
################################################################################

# 1. REMOVE PLOTS WITH DUPLICATED COORDINATES

duplicate.cleaner <- function(PLOTS)
  # this functions finds and deletes duplicates, and gives preference to large
  # plots = if two plots share coordinates, but one of them is larger, the larger
  # is kept preferentially
{
  PLOTS <- PLOTS[order(PLOTS$Area_ha, decreasing=TRUE),]
  latlon <- data.frame(Lat = PLOTS$Lat, Lon = PLOTS$Lon)
  PLOTS <- PLOTS[duplicated(latlon) == FALSE,]
  PLOTS <- PLOTS[order(PLOTS$Numerical_ID),]
  return(PLOTS)
}

PLOTS <- duplicate.cleaner(PLOTS)

# ------------------------------------------------------------------------------

# 2. REMOVE PLOTS THAT DON'T MEET THE STRICTER CRITERIA ON TREE DEFINITION OR DBH

good.tree.id <- PLOTS$Tree_def %in% c("trees")
good.dbh.id <- PLOTS$min_DBH_cm >= 10
good.plots <- (good.tree.id + good.dbh.id) > 0
good.plots[is.na(good.plots)] <- FALSE

PLOTS.subset <- PLOTS[good.plots,]
nrow(PLOTS.subset)


# ------------------------------------------------------------------------------

# 3. SELECT ONLY THE VARIABLES RELEVANT FOR THE ANALYSES

PLOTS <- dplyr::select(PLOTS, Dataset, Loc_ID = Plot_ID, Lat, Lon, 
                       Area_km = Area_ha,
                       S, N, min_DBH_cm, GPP, ET, ANN_T, WARM_T, ISO_T, MIN_P, 
                       P_SEAS, 
                       ALT_DIF, ISLAND, INSULARITY, ELONGATION, 
                       HABITAT, REALM, REALM_PK = REALM_PK)
PLOTS$Area_km <- PLOTS$Area_km/100
PLOTS <- data.frame(PLOTS, DAT_TYPE="Plot")
PLOTS <- na.omit(PLOTS)

# ---------------

PLOTS.subset <- dplyr::select(PLOTS.subset, Dataset, Loc_ID = Plot_ID, Lat, Lon, 
                       Area_km = Area_ha,
                       S, N, min_DBH_cm, GPP, ET, ANN_T, WARM_T, ISO_T, MIN_P, 
                       P_SEAS, 
                       ALT_DIF, ISLAND, INSULARITY, ELONGATION, 
                       HABITAT, REALM, REALM_PK = REALM_PK)
PLOTS.subset$Area_km <- PLOTS.subset$Area_km/100
PLOTS.subset <- data.frame( PLOTS.subset, DAT_TYPE="Plot")
PLOTS.subset <- na.omit(PLOTS.subset)


################################################################################
# PROCESSING COUNTRIES
################################################################################


CNTRS.dat <- dplyr::select(CNTRS@data, Aggregated, Dataset, Loc_ID = NAME, Lat, Lon, Area_km,
                            S, N = TREE_DENS, min_DBH_cm, GPP, ET, ANN_T, WARM_T, 
                            ISO_T, MIN_P, P_SEAS, 
                            ALT_DIF, ISLAND, INSULARITY, ELONGATION, 
                           HABITAT, REALM, REALM_PK = REALM_3)
CNTRS.dat <- data.frame(CNTRS.dat, DAT_TYPE="Country")


# remove Antarctica and Oceania
CNTRS.dat <- CNTRS.dat[CNTRS.dat$REALM!="Antarctic",]
CNTRS.dat <- CNTRS.dat[CNTRS.dat$REALM!="Oceania",]

# remove NA's in the Dataset field
CNTRS.dat <- na.omit(CNTRS.dat)


# fix REALM_PK for China, US and Brazil
CNTRS.dat$REALM_PK <- as.character(CNTRS.dat$REALM_PK )
CNTRS.dat[CNTRS.dat$Loc_ID == "China", 'REALM_PK'] <- "Indo-Malay"
CNTRS.dat[CNTRS.dat$Loc_ID == "United States",'REALM_PK'] <- "Nearctic"
CNTRS.dat[CNTRS.dat$Loc_ID == "Brazil",'REALM_PK'] <- "Neotropic"

# create an indicator vector for the level of aggregation
Aggregated <- as.character(dplyr::select(CNTRS.dat, Aggregated)[,1])

# ------------------------------------------------------------------------------
# DATA WITH FULL DETAIL, DISAGGREGATED LARGE COUNTRIES, AND ALL PLOTS

# remove the large aggregated spatial units (China, US, Brazil)
CNTRS.fine <- CNTRS.dat[CNTRS.dat$Loc_ID %in% c("China", "United States", "Brazil") == FALSE,]
CNTRS.fine <- na.omit(CNTRS.fine)

# remove the Aggregated field
CNTRS.fine <- CNTRS.fine[,-1]

ALL <- rbind(PLOTS, CNTRS.fine)

# ------------------------------------------------------------------------------
# SMALLER DATA WITH CHINA, US AND BRAZIL INTACT, AND WITH ONLY THE BEST PLOTS

# remove the fine dis-aggregated parts of China, US and Brazil
CNTRS.coarse <- CNTRS.dat[CNTRS.dat$Aggregated == "Original country",]

# remove the Aggregated field
CNTRS.coarse <- CNTRS.coarse[,-1]


ALL.subset <- rbind(PLOTS.subset, CNTRS.coarse)



################################################################################
# EXPORT THE DATA
################################################################################

write.csv(ALL, 
          file="../Data/Main_dataset_full_detail.csv", 
          row.names = FALSE)

write.csv(ALL.subset, 
          file="../Data/Main_dataset_subset.csv", 
          row.names = FALSE)






