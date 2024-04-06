# Create tables to summarize the XCI results

### Libraries ###
require(dplyr)

### Functions ###

#' Function to construct the fit path
constructFitPath <- function(donorName) {
  return(file.path(getwd(), donorName, 'inactiveXX',
                   sprintf('%s_fit.rds', donorName)))
}

#' Create a table summarizing the inactiveXX results
#' Include the XCI ratio, whether inactiveXX threw a warning for the donor, 
#' the number of determined and undetermined cells, the XCI ratio, and whether
#' the XCI ratio is skewed or not (False if 0.4 < XCI ratio < 0.6, and True
#' otherwise)
#'
#' @param inactiveXXDataPath where to store the inactiveXX summary table
getInactiveXXTable <- function(inactiveXXDataPath) {
  # List all .rds files
  rdsFiles <- list.files(pattern = '_summary\\.rds$', recursive = TRUE)
  
  # Read and combine all .rds files
  inactiveXXData <- lapply(rdsFiles, function(file) {
    readRDS(file)
  }) %>% bind_rows()
  
  # Add a column to store if the XCI ratio is skewed
  inactiveXXData <- inactiveXXData %>%
    #filter(! inconclusive) %>% # should inconclusive donors be removed here?
    mutate(skewed = donorTau < 0.40 | donorTau > 0.60)
  
  # Save the combined data into one .rds file
  saveRDS(inactiveXXData, file = inactiveXXDataPath)
}

#' Combine the inactiveXX results and the metadata. Add a path to the 
#' inactiveXX results
#' 
#' @param mDatPath where the metadata is stored
#' @param inactiveXXDataPath where the inactiveXX summary results table is 
#'                           stored (this is the table from getInactiveXXTable)
#' @param mDatCombinedPath where to store the metadata and inactiveXX results
#'                         combined table
combineInactiveXXAndMetadata <- function(mDatPath, inactiveXXDataPath, 
                                         mDatCombinedPath) {
  mDat <- readRDS(file = mDatPath)
  inactiveXXData <- readRDS(file = inactiveXXDataPath)
  mDat <- mDat %>%
    merge(inactiveXXData, by.x = 'humanDonorName', by.y = 'donorName') %>%
    distinct(humanDonorName, groupName, .keep_all = TRUE) %>%
    filter(! inconclusive)
  
  # Create a vector to store the results paths
  fitsPath <- vector('list', length = nrow(mDat))
  
  # Loop through each row in the dataframe and add the path to the inactiveX 
  # results 
  for (i in 1:nrow(mDat)) {
    donorName <- mDat$humanDonorName[i]
    fitPath <- constructFitPath(donorName)
    if (file.exists(fitPath)) {
      fitsPath[[i]] <- fitPath
    } else {
      fitsPath[[i]] <- NA
    }
  }
  
  # Add the fit and results path to mDat
  mDat$fitPath <- fitsPath
  
  saveRDS(mDat, file = mDatCombinedPath)
}
