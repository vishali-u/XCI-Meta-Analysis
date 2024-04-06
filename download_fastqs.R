# Second step: Download fastq files from the HCA

#' Download an HCA fastq file
#'
#' @param fastqDest the complete file path for this fastq file
#' @param fastqURL a URL link to the fastq file
downloadFASTQ <- function(fastqDest, fastqURL) {
  # the needed fastq file has already been downloaded
  if (file.exists(fastqDest)) {
    return()
  }
  
  message(sprintf("Downloading %s", fastqDest))
  
  # download the fastq file using curl
  cmd <- sprintf('curl -L "%s" -o "%s" -s', fastqURL, fastqDest)
  system(capture.output(cat(cmd)))
  
  # if the fastq file has not been downloaded, stop running
  if (! file.exists(fastqDest)) {
    stop(sprintf('%s was not downloaded.', fastqDest))
  }
  
  # check if the file is a JSON file which would indicate an error occurred
  mimeType <- system(sprintf("file --mime-type -b '%s'", fastqDest), 
                     intern = TRUE)
  if (grepl('application/json', mimeType)) {
    stop(sprintf('Downloaded file %s is a JSON file', fastqDest))
  }
}

#' Download all the fastq files in a dataset
#' 
#' @param dataDir the parent directory for the directory in which the fastq
#'                files should be stored
#' @param mDat the metadata for this dataset
downloadAll <- function(dataDir, mDat) {
  for (i in 1:nrow(mDat)) {
    
    donorName <- mDat$humanDonorName[i]
    group <- mDat$groupName[i]
    
    # create a directory for the donor
    if (! dir.exists(file.path(dataDir, donorName))) {
      dir.create(file.path(dataDir, donorName))
    }
    
    # create a directory for the group
    if (! dir.exists(file.path(dataDir, donorName, group))) {
      dir.create(file.path(dataDir, donorName, group))
    }
    
    fastqDirectory <- file.path(dataDir, donorName, group, 'fastq')
    if (! dir.exists(fastqDirectory)) {
      dir.create(fastqDirectory)
    }
    
    # download a FASTQ file
    fastqName <- mDat$fileName[i]
    fastqURL <- mDat$file_url[i]
    fastqDest <- file.path(fastqDirectory, fastqName)
    downloadFASTQ(fastqDest = fastqDest, fastqURL = fastqURL)
  }
}
