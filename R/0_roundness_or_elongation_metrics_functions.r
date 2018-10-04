################################################################################
# Functions for calculation of the ELONGATION METRICS

# libraries for working with ellipses 
library(cluster) # for the ellipsoid hull funciton
library(car) # another one for ellipses


# ------------------------------------------------------------------------------
# Function that takes a polygon (must be in WGS84 Lat Long system), 
# and calculates its elongation based on the ratio of the first two
# eigenvalues of a Principal Coordinate Analysis of its ***borders***.

# Arguments: pol - the polygon
#            draw - should the polygon be plotted?

elongation.sample <- function(pol, draw = FALSE)
{
  LAT = coordinates(pol)[1,2]
  LON = coordinates(pol)[1,1]
  AED <- paste("+proj=aeqd +lat_0=", 
               LAT, 
               " +lon_0=", 
               LON, 
               " +x_0=100000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
               sep = "")
  pol <- spTransform(pol, CRSobj = CRS(AED))
  
  # convert polygon to lines
  lin <- as(pol, "SpatialLines")
  # sample points around the line at regular intervals
  pts <- spsample(lin, n = 100, type = "regular")
  # distance matrix for the points
  dst <- dist(coordinates(pts))
  # Principal Coordinates Analysis of the distance matrix
  eigs <- cmdscale(dst, k = 2, eig = TRUE)$eig[1:2]
  # Ratio of the first and the second eigenvalue (squared)
  elong <- sqrt(eigs[1]) / sqrt(eigs[2])
  
  if(draw)
  {
    plot(pol, main  = paste(pol@data$NAME[1], ", ", round(elong, 2 ), sep = ""), col="grey")
    points(pts, col="red")
  }
  
  return(elong)
}


# ------------------------------------------------------------------------------
# Function that takes a polygon (must be in WGS84 Lat Long system), 
# and calculates its elongation based on the ratio of the first two
# eigenvalues of a Principal Coordinate Analysis of an ***ellipsoid hull***
# fitted around the borders.

# Arguments: pol - the polygon
#            draw - should the polygon be plotted?

elongation.ellipse <- function(pol, draw = FALSE)
{
  LAT = coordinates(pol)[1,2]
  LON = coordinates(pol)[1,1]
  AED <- paste("+proj=aeqd +lat_0=", 
               LAT, 
               " +lon_0=", 
               LON, 
               " +x_0=100000 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
               sep = "")
  pol <- spTransform(pol, CRSobj = CRS(AED))
  
  # convert polygon to lines
  lin <- as(pol, "SpatialLines")
  # convert lines to points
  pts <- as(lin, "SpatialPoints")
  # fit the ellipsoid hull to the points
  elps <- ellipsoidhull(pts@coords)
  # convert elps to an actual set of points
  pts <- ellipse(center = elps$loc, shape = elps$cov, radius = sqrt(elps$d2), 
                 segments = 100, draw = FALSE, add=FALSE)
  
  # distance matrix for the points
  dst <- dist(pts)
  # Principal Coordinates Analysis of the distance matrix
  eigs <- cmdscale(dst, k = 2, eig = TRUE)$eig[1:2]
  # Ratio of the first and the second eigenvalue (squared)
  elong <- sqrt(eigs[1]) / sqrt(eigs[2])
  
  if(draw)
  {
    plot(pol, main  = paste(pol@data$NAME[1], ", ", round(elong, 2 ), sep = ""), col="grey")
    lines(pts, col="red")
  }
  
  return(elong)
}

#-------------------------------------------------------------------------------
# Should the plots showing the elongation metric for each polygon be plotted?

test.on.all.polygons = FALSE

if(test.on.all.polygons)
{
  CTR <- readOGR(dsn = "../Data/COUNTRIES", layer = "COUNTRIES")
  
  # plot a huge figure with all the polygons and their elongation value 
  pdf("../Figures/Elongations.pdf", width = 30, height = 100)
  par(mfrow=c(36,10))
  for(i in 1:nrow(CTR))
  {
    pol <- CTR[i,]
    elong <- elongation.sample(pol, draw = TRUE)
  }
  dev.off()
  
  # compare the elongation metrics
  res <- list()
  for(i in 1:nrow(CTR))
  {
    print(i)
    pol <- CTR[i,]
    
    res[[i]] <- c(elong.samp = elongation.sample(pol, draw = FALSE),
                  elong.elps = elongation.ellipse(pol, draw = FALSE))
  }
  res <- plyr::ldply(res)
  
  par(mfrow = c(1,2))
  plot(res); abline(a=0, b=1)
  boxplot(res)
  
}
