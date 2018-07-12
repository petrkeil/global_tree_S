################################################################################
# Author: Petr Keil
# Email: pkeil@seznam.cz
# Date: Dec 7 2017
################################################################################

# Description: This is a script that extracts the data on species composition 
# in each country in the World form the Global Tree Search database.
# The main tool used in the script is RSelenium.

################################################################################

# LOAD THE NECESSARY LIBRARIES
library(RSelenium)
library(XML)
library(magrittr)
library(plyr)

# THIS NEEDS TO BE RUN IN THE LINUX TERMINAL -
# sudo docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.0
# sudo docker ps

# START THE DRIVER
remDr <- remoteDriver(port = 4445L)



# ------------------------------------------------------------------------------
# SPECIES LIST

# load the species list
tree.list <- read.csv("../Data/COUNTRIES/global_tree_search_trees.csv")

# remove diploon cuspidatum, since it contains weird characters
tree.list <- tree.list[-17781,]

# split the genera and species to two separate columns
taxa <- as.character(tree.list$Taxon.name) %>% strsplit(split=" ") %>% ldply()  
names(taxa) = c("genus", "species")




# ------------------------------------------------------------------------------
# FUNCTION THAT DOWNLOADS THE DATA FOR A SINGLE TREE

# Arguments:
# genus - character string giving the genus name
# species - character string givint a species name

get.tree <- function(genus, species)
{
  genus <- as.character(genus)
  species <- as.character(species)
  
  # go to the webpage
  remDr$navigate("http://www.bgci.org/global_tree_search.php?sec=globaltreesearch")
  remDr$refresh() # refresh the page
  
  # create R objects from the website elements
  genusElem <- remDr$findElement(using = 'id', value = "genus-field")
  specElem <- remDr$findElement(using = 'id', value = "species-field")
  buttElem <- remDr$findElement(using = 'class', value = "btn_ohoDO")
  
  # fill in the forms with the genus and species names
  genusElem$sendKeysToElement(list(genus))
  specElem$sendKeysToElement(list(species))
  
  # click the search button
  buttElem$clickElement()
  
  # get the output
  out <- remDr$findElement(using = "css", value="td.cell_1O3UaG:nth-child(4)")$getElementText()
  
  return(out)  
}



# ------------------------------------------------------------------------------
# THE MAIN WEB-SCRAPING LOOP

dir.create(path = "../Data/COUNTRIES/scraped_GTS_data")

remDr$open()
for(i in 1:nrow(taxa))
{
  genus <- taxa[i,'genus']
  species <- taxa[i,'species']
  full.name <- paste(genus, species)
  
  print(full.name)
  print(i)
  
  for(j in 1:10)  # try the function 10x before giving up
  {
    out <- tryCatch(get.tree(genus, species), error=function(e) "failed")
    if(out != "failed") break # if success, break and skip to next species
    if(out == "failed") remDr$close(); remDr$open() # if failure, reset connection and try again
  }
  
  write.table(out, file=paste("../Data/COUNTRIES/scraped_GTS_data/", full.name, ".txt", sep=""),
              col.names=FALSE, row.names=FALSE)
}
remDr$close()



################################################################################
# PROCESSING THE SCRAPPED DATA TO A DATA FRAME
################################################################################

fl.path <- "../Data/COUNTRIES/scraped_GTS_data/"

# read the file list
fls <- list.files(fl.path)

# look at one species
spec.no <- 182
filename <- fls[spec.no]


# initial pre-processing steps
all.data <- list()
N <- length(fls)
pb <- txtProgressBar(min = 0, max = N, style = 3)

# loop that goes through the file list and appends the data to a list
for(i in seq_along(fls))
{
  filename <- fls[i]
  species  <-  gsub(pattern=".txt", replacement="", x=filename) 
  
  cntrs <- tryCatch(as.character(read.table(paste(fl.path, filename, sep=""))[1,]),
                    error=function(e) "empty")
  cntrs <- strsplit(cntrs, split="; ")[[1]]
  
  all.data[[species]] <- data.frame(Species = rep(species, times=length(cntrs)),
                                    Country = cntrs)
  setTxtProgressBar(pb, i)
}

# squish the data to a data frame
all.data.table <- ldply(all.data, rbind)
all.data.table <- all.data.table[,2:3]

# export the data
write.csv(all.data.table, "../Data/COUNTRIES/GTS_database.csv", row.names = FALSE)


### This .csv was then manually aligned with names in TM_WORLD_BORDERS-0.3 dataset
# to produce coutnry_list_translation_table.csv

