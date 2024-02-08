# Functions to prepare data for the inactiveXX pipeline by downloading and 
# aligning the fastq files

####### 1. Libraries #######
require(data.table)
require(dplyr)

####### 2. Constants and File Paths#######
# Cellranger
nThreads <- 32
refGenome <- '/home/umaiyal1/scratch/refdata-gex-GRCh38-2020-A'

# STARsolo script paths:
tenX <- '/home/umaiyal1/XCI-Meta-Analysis/star_scripts/10x.sh'
dropSeq <- '/home/umaiyal1/XCI-Meta-Analysis/star_scripts/drop_seq.sh'
inDrop <- '/home/umaiyal1/XCI-Meta-Analysis/star_scripts/in_drop.sh'

####### 3. Functions ######

#' Create a mapped directory for one donor if those directories
#' do not already exist
#' 
#' @param donorName the name of the donor the directories should be made for
createMappedDirectory <- function(dataDir, donorName, group, useCellRanger=FALSE) {
  
  if (useCellRanger) {
    mappedDirectory <- file.path(dataDir, donorName, group)
  } else {
    mappedDirectory <- file.path(dataDir, donorName, group, 'mapped')
    if (! dir.exists(mappedDirectory)) {
      dir.create(mappedDirectory)
    } 
  }
  
  return(mappedDirectory)
}

#' Return TRUE if this directory contains an outs directory, and FALSE otherwise
doesOutsExist <- function(directory) {
  outsPath <- file.path(directory, 'outs')
  return(dir.exists(outsPath))
}

#' Run STARsolo on all the downloaded fastq files
#' Assumes that the sequencing method is either 10x, in-drop, drop-seq,
#' smart-seq2, or fluidigm C1
#' 
#' @param fastqDir the directory holding the fastq files
#' @param mappedDir the destination for the mapped files
#' @param seqMethod the protocol used to prepare the data
#' @return the path to the bam file and the path to the matrices 
runSTARsolo <- function(fastqDir, mappedDir, seqMethod) {
  
  bamName <- paste(mappedDir, 'Aligned.sortedByCoord.out.bam', sep = '/')
    
  if (file.exists(bamName)) {
    message(sprintf('%s already exists.', bamName))
  }
    
  if (grepl('10X', seqMethod, ignore.case = TRUE)) {
    cmd <- sprintf('cd %s && sh %s %s %s', mappedDir, tenX, fastqDir, group)
  } 
    
  else if (grepl('drop-seq', seqMethod, ignore.case = TRUE)) {
    cmd <- sprintf('cd %s && sh %s %s %s', mappedDir, dropSeq, fastqDir, group)
  } 
    
  else if (grepl('inDrop', seqMethod, ignore.case = TRUE)) {
    cmd <- sprintf('cd %s && sh %s %s %s', mappedDir, inDrop, fastqDir, group)
  }
    
  result <- system2(cmd, stdout = TRUE, stderr = TRUE)
  if (result$status != 0) {
    stop(paste("The shell command failed with message:", result$stderr))
  }
  
  message(sprintf('Finished mapping: %s', bamName))
    
  if (! file.exists(bamName)) {
    stop(sprintf('%s was not created.', bamName))
  }
}

#' Run cellranger on the fastq files for a group
#'
runCellRanger <- function(groupName, fastqPath, mappedDest) {
  # If the directory contains the 'outs' directory, assume cellranger has
  # finished running.
  if (doesOutsExist(file.path(mappedDest, 'mapped'))) {
    return()
  }
  
  message(sprintf("Started mapping files at: %s", mappedDest))
  
  cmd <- sprintf("cd %s && cellranger count --id=mapped --transcriptome=%s --fastqs=%s --sample=%s --expect-cells=300 --localcores=%s --localmem=128",
                 mappedDest, refGenome, fastqPath, groupName, nThreads)
  
  system(cmd)
  
  # result <- system2(cmd, stdout = TRUE, stderr = TRUE)
  # if (result$status != 0) {
  #   stop(paste("The shell command failed with message:", result$stderr))
  # }
  
  message(sprintf("Finished mapping. Files at: %s", 
                  file.path(mappedDest, 'mapped')))
}


#' Index one or more BAM file(s)
#'
#' @param bamPaths a vector of bam file paths
indexBAM <- function(bamPaths) {
  
  for (bamPath in bamPaths) {
    indexName <- paste0(bamPath, '.bai')
    
    # Skip if the BAM file has already been indexed
    if (file.exists(indexName)) {
      next
    }
    
    # Index the bam file with samtools index with 4 threads
    cmd <- sprintf('samtools index -@4 %s', bamPath)
    system(cmd)
    
    # Check if the indexed BAM file exists (.bai file extension)
    if (!file.exists(indexName)) {
      stop(sprintf('%s was not created.', indexName))
    } else {
      message(sprintf('Indexed: %s', bamPath))
    }
  }
}

#' Align the fastq files, and index the BAM files if necessary
#' 
#' Assumes there are fewer than 1 billion donors
#'
alignData <- function(mDat, dataDir, useCellRanger = FALSE, 
                      startRow = 0, endRow = 1e9) {
  
  # Ensure startRow and endRow are within the bounds of mDat if running in parallel jobs
  startRow <- max(1, startRow)
  endRow <- min(length(unique(mDat$humanDonorName)), endRow)

  # Filter mDat for the specified range
  for (donor in unique(mDat$humanDonorName)[startRow:endRow]) {
    
    mDatDonor <- mDat[humanDonorName == donor]
    
    for(grp in unique(mDatDonor$groupName)) {
      tgts <- mDatDonor[groupName == grp]
      
      # All files in a group for a donor has the same preparation protocol?
      seqMethod <- tgts$library_preparation_protocol.library_construction_approach[1]
      
      mappedDirectory <- createMappedDirectory(dataDir = dataDir, 
                                               donorName = donor, 
                                               group = grp,
                                               useCellRanger = useCellRanger)
      fastqDirectory <- file.path(dataDir, donor, grp, 'fastq')
      
      if (useCellRanger) {
        runCellRanger(groupName = grp,
                      fastqPath = fastqDirectory,
                      mappedDest = mappedDirectory)

      } # TODO: add STAR
    }
  }
}
