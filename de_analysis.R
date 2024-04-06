# # Find genes that are differentially expressed between cells with differernt
# # copies of an inactived X chromosome

require(Seurat) 
require(ggplot2)
require(dplyr)

#' Run differential expression for a donor in a dataset.
#' Save a table of differentially expressed genes and a volcano plot for this
#' donor. 
#' 
#' @param mDat the metadata for a dataset
#' @param sampleIndex the index for the donor to run DE on
#' @param resultsPath a location to store DE results
runDifferentialExpression <- function(mDat, sampleIndex, resultsPath) {
  
  # Store the data for the donor
  donorName <- mDat$humanDonorName[sampleIndex]
  sampleName <- mDat$groupName[sampleIndex]
  donorSratPath <- mDat$sratPath[[sampleIndex]]
  srat <- readRDS(file = donorSratPath)
  
  # Remove the cells that have XiState as NA or Undetermined
  # If a cell refers to NA, that means that this cell was filtered out when
  # running inactiveXX??
  srat <- subset(srat, subset = XiState %in% c("Maternal", "Paternal"))
  srat <- SetIdent(srat, value = "XiState")
  
  # TODO: should some we filter out donors that have fewer than 100 cells?
  deGenes <- FindMarkers(srat, ident.1 = 'Maternal', ident.2 = 'Paternal',
                         test.use = 'wilcox')
  
  # identCounts <- table(Idents(srat))
  # if (identCounts[[1]] < 100 || identCounts[[2]] < 100) {
  #   stop("There are fewer than 100 cells for one of the genotypes (maternal
  #         or paternal).")
  # } else {
  #   deGenes <- FindMarkers(srat, ident.1 = 'Maternal', ident.2 = 'Paternal',
  #                          test.use = 'wilcox')
  # }
  
  deGenes$neg_logPval <- - log10(deGenes$p_val)
  deGenes <- na.omit(deGenes)
  
  # Make a separate folder for the DE analysis results
  deAnalysisResults <- file.path(resultsPath, donorName, 'de_analysis')
  if (! dir.exists(deAnalysisResults)) {
    dir.create(deAnalysisResults)
  }
  
  saveRDS(deGenes, 
          file = file.path(deAnalysisResults,
                           sprintf('%s_%s_deGenes.rds', donorName, sampleName)))
  
  # Make a volcano plot for the donor
  # Choose 1.5 as the threshold fold change and 0.05 as the p-value threshold
  log2fcThreshold <- 1.5
  adjPValueThreshold <- 0.05
  
  deGenes$color <- ifelse((deGenes$p_val_adj < adjPValueThreshold),'red', 'grey')
  
  dePlot <- ggplot(data = deGenes, aes(x = avg_log2FC, y = neg_logPval)) +
    geom_point(aes(color = color), alpha = 0.4) +
    scale_color_identity() +
    theme_minimal() +
    theme(legend.position = 'none') +
    labs(title = sprintf('Differential Expression for Donor %s', donorName),
         x = 'Log2(Fold Change)',
         y = '-Log10(P-value)') +
    geom_vline(xintercept = log2fcThreshold, linetype = 'dashed') +
    geom_vline(xintercept = -log2fcThreshold, linetype = 'dashed')
  
  ggsave(file.path(deAnalysisResults, sprintf('%s_dePlot.jpeg', donorName)),
         plot = dePlot, width = 10, height = 8, dpi = 300)
}
