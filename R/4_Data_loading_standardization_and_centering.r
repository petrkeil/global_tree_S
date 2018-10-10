# load the data
PLT <- read.csv("../Data/Main_dataset_full_detail.csv")


# calculate tree density (note the x+1 step!!)
PLT$Tree_dens <- (PLT$N + 1) / PLT$Area_km 


# select only the variables of interest from the larger data.frame
DAT <- dplyr::select(PLT, S, Area_km, Tree_dens, min_DBH=min_DBH_cm, 
                     GPP, ET, ANN_T, WARM_T, ISO_T, MIN_P, P_SEAS, ALT_DIF, ELONG,
                     ISLAND = ISL_LS,
                     ISL_LS, ISL_ST, ISL_DIS, REALM=REALM_PK, Lat, Lon, DAT_TYPE, Loc_ID) 


# order the data frame by regions and area
DAT <- DAT[order(DAT$REALM),]


# transform area and tree density
DAT$Area_km <- log(DAT$Area_km)
DAT$Tree_dens <- log(DAT$Tree_dens)


# remove NAs
DAT <- DAT[rowSums(is.na(DAT)) == 0,]


# mean and sd of Area (will be used later to bring these to their origina scale)
A.mean <- mean(DAT$Area_km)
A.sd <- sd(DAT$Area_km)


# see column names and their numbers
data.frame(names(DAT))


# means and sd of the rest of the variables
centr <- attributes(scale(DAT[,2:13]))$'scaled:center'
scale <- attributes(scale(DAT[,2:13]))$'scaled:scale'
scale.tab <- data.frame(var=names(centr), centr, scale)
write.csv(scale.tab, file="scale_tab.csv", row.names = FALSE)


# do the actual scaling
DAT[,2:13] <- scale(DAT[,2:13])
