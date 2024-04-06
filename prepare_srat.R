# Fourth step: Prepare a seurat object each donor in a dataset and store it

# TODO: add a check for the files before calling Read10X 
# The file names must be barcodes.tsv genes.tsv and matrix.mtx
# Must rename features.tsv to genes.tsv

### Libraries ###
require(Seurat) 

prepareSeurat <- function(countsPath, xiStates, metadataPath=NULL,
                          sampleID) {
  
  # Load the data and store the non-normalized counts in a Seurat object
  srat.data <- Read10X(data.dir = countsPath)
  srat <- CreateSeuratObject(counts = srat.data)
  
  if (!is.null(metadataPath)) {
    # Load the metadata
    metadata <- read.csv(metadataPath)
    metadata <- metadata[metadata$suspension_uuid == sampleID,]
    
    # Extract the barcodes, count up the occurrences, and make the row names
    metadata$barcode <- sub('^(.*?)-.*', '\\1', metadata$index)
    metadata$barcodeCount <- ave(rep(1, nrow(metadata)),
                                 metadata$barcode, FUN = cumsum)
    metadata$rowName <- paste(metadata$barcode, metadata$barcodeCount, sep = '-')
    rownames(metadata) <- metadata$rowName
    
    # Identify common barcodes between metadata and Seurat object
    commonBarcodes <- intersect(rownames(metadata), names(srat@active.ident))
    metadata <- metadata[commonBarcodes, ]
    
    # Add the metadata
    srat <- AddMetaData(srat, metadata = metadata) 
  }
  
  # Calculate percentage of all counts belonging to a subset of the features 
  # for each cell (as a quality control). Filter results based on this.
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
  srat <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                   percent.mt < 5)
  
  # Normalize expression measurements for each cell (using default
  # normalization methods from the Seurat package)
  srat <- NormalizeData(srat)
  
  # Scale the data and then perform dimensional reduction (pca)
  srat <- FindVariableFeatures(srat, selection.method = 'vst',
                               nfeatures = 2000)
  all.genes <- rownames(srat)
  srat <- ScaleData(srat, features = all.genes)
  features <- VariableFeatures(object = srat)
  srat <- RunPCA(srat, features = features)

  # Cluster the cells
  srat <- FindNeighbors(srat, dims = 1:10)
  srat <- FindClusters(srat, resolution = 0.5)

  # Perform non-linear dimensional reduction
  srat <- RunUMAP(srat, dims = 1:10)
  
  srat <- AddMetaData(srat, metadata = xiStates, col.name = 'XiState')
  
  return(srat)
}
