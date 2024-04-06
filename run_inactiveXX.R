# Run inactiveXX on HCA data

### Libraries ###
require(inactiveXX)

####### Constants and Directories #######
useVartrix <- FALSE
nParallel <- 8

# Path to reference genome used to produce BAM files, BAM files, and metadata
# Reference genome must also be indexed
refGenomePath <- '~/scratch/refdata-gex-GRCh38-2020-A/fasta/genome.fa'

####### Functions #######
#' Run inactiveXX on the dataset stored in mDat. Store the results in directory
#' with the resultsPath as a parent directory. Also store a summary table that
#' includes: the donor ID, whether inactiveXX threw a warning for this donor,
#' the XCI ratio, the total number of cells predicted by inactiveXX, the total
#' number of cells that could be determined (as paternal or maternal) by 
#' inactiveXX, and the total number of cells that were undetermined by 
#' inactiveXX.
#' 
#' Assumes all the BAM files for the specified donor exist and are indexed.
#' 
#' @param mDat metadata about the dataset
#' @param donorIndex the index for a single donor
#' @param resultsPath a parent directory to store the inactiveXX results
runInactiveXX <- function(mDat, donorIndex, resultsPath) {
  
  donor <- unique(mDat$humanDonorName)[donorIndex]
  
  # Initialize variables 
  inconclusive <- FALSE # Store whether inactiveXX worked for this donor
  donorTau <- NULL # Store the XCI skew
  
  # Create a directory for this donor in the results folder
  # If a directory for the donor already exists, assume that donor is done
  if (! dir.exists(file.path(resultsPath, donor))) {
    dir.create(file.path(resultsPath, donor))
  }
  
  # Create a directory for the inactiveXX results in this donor's folder
  if (! dir.exists(file.path(resultsPath, donor, 'inactiveXX'))) {
    dir.create(file.path(resultsPath, donor, 'inactiveXX'))
  }
  
  # Get the BAM file paths for this donor
  donorFiles <- mDat[mDat$humanDonorName == donor,]
  bams <- donorFiles$bamFile
  
  donorResultsPath <- file.path(resultsPath, donor, 'inactiveXX')
  
  # Name the bams using the sample name
  names(bams) <- donorFiles$groupName
  
  # Call heterozygous SNPs and save the counts 
  hSNPs <- hetSNPsFromRNA(bams,
                          refGenomePath,
                          nParallel = nParallel,
                          outputs = file.path(donorResultsPath,
                                              sprintf('%s_XCnts.tsv', 
                                                      names(bams))),
                          useVartrix = useVartrix)
  
  # Filter counts
  xCnts <- filterCountsX(hSNPs)
  
  # Run inferInactivX with default parameters except for nParallel and tauDiffWarnOnly
  fit <- withCallingHandlers({
    inferInactiveX(xCnts, nParallel = nParallel, tauDiffWarnOnly = TRUE)
  }, warning = function(w) {
    
    # Since a warning was thrown for this donor, set inconclusive to true
    inconclusive <<- TRUE
    
    # Allow the execution to continue
    invokeRestart("muffleWarning")
  })
  
  donorTau <- fit$tau
  saveRDS(fit, file = file.path(donorResultsPath,
                                sprintf('%s_fit.rds', donor)))
  
  # Summarize the number of determined and undetermined cells
  xiStates <- fit$cellSummary$stateXi
  xiCounts <- table(factor(xiStates, 
                           levels = c('Maternal', 'Paternal', 'Undetermined')))
  xiCounts['Determined'] <- sum(xiCounts['Maternal'], xiCounts['Paternal'])
  
  # Make a summary table for future use
  donorResult <- data.frame(
    donorName = donor,
    inconclusive = inconclusive,
    donorTau = donorTau,
    totalNumberOfCells = length(fit$states),
    totalNumberOfDeterminedCells = xiCounts[['Determined']],
    totalNumberOfUndeterminedCells = xiCounts[['Undetermined']],
    stringsAsFactors = FALSE
  )
  
  saveRDS(donorResult, file = file.path(donorResultsPath,
                                        sprintf('%s_summary.rds', donor)))
  
  # Plot the solution and save the resulting plot
  jpeg(filename = file.path(donorResultsPath, 
                            sprintf('%s_plot.jpeg', donor)), 
       height = 575, width = 668)
  
  plotSolutions(fit)
  dev.off()
}