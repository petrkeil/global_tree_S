
source("0_libraries_functions_settings.r")

ISL <- raster::raster("/media/pk33loci/Elements/GIS_data/ISLANDNESS/rasters/ALL_CLASSES_clean.tif")
CTR <- readOGR(dsn = "../Data/COUNTRIES", layer = "COUNTRIES")



pdf("../Figures/Islandness_polygons.pdf", width = 30, height = 100)
    par(mfrow=c(36,10))
    
    for(i in 1:nrow(CTR))
    {
      print(i)
      pol <- CTR[i,]
      rst <- crop(ISL, y = extent(pol))
      plot(pol, main = pol@data$NAME[1])
      plot(rst, axes = FALSE, box=FALSE, legend=FALSE, add=TRUE)
      plot(pol, add=TRUE)
    }
dev.off()