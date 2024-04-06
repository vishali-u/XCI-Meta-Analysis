# Install all packages needed to run inactiveXX pipeline
# module bcftools must be loaded first (required to install gsl)

# First run: module load apptainer/1.1.8 StdEnv/2020 gcc/9.3.0 bcftools/1.11 r/4.2.1

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

if (! requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo")
}

if (! requireNamespace("checkmate", quietly = TRUE)) {
  install.packages("checkmate")
}

if (! requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

bettermc_url <- 'https://cran.r-project.org/src/contrib/Archive/bettermc/bettermc_1.2.1.tar.gz'
if (! requireNamespace("bettermc", quietly = TRUE)) {
  install.packages(bettermc_url, repos = NULL, type = "source")
}

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (! requireNamespace("zlibbioc", quietly = TRUE)) {
  BiocManager::install("zlibbioc", update = TRUE)
}

if (! requireNamespace("XVector", quietly = TRUE)) {
  BiocManager::install("XVector", update = TRUE)
}

if (! requireNamespace("xopen", quietly = TRUE)) {
  install.packages("xopen")
}

if (! requireNamespace("XML", quietly = TRUE)) {
  install.packages("XML")
}

if (! requireNamespace("waldo", quietly = TRUE)) {
  install.packages("waldo")
}

if (! requireNamespace("viridisLite", quietly = TRUE)) {
  install.packages("viridisLite")
}

if (! requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}

if (! requireNamespace("vipor", quietly = TRUE)) {
  install.packages("vipor")
}

if (! requireNamespace("VGAM", quietly = TRUE)) {
  install.packages("VGAM")
}

if (! requireNamespace("BH", quietly = TRUE)) {
  install.packages("BH")
}

if (! requireNamespace("dqrng", quietly = TRUE)) {
  install.packages("dqrng")
}

if (! requireNamespace("uwot", quietly = TRUE)) {
  install.packages("uwot")
}

if (! requireNamespace("usethis", quietly = TRUE)) {
  install.packages("usethis")
}

if (! requireNamespace("tweenr", quietly = TRUE)) {
  install.packages("tweenr")
}

if (! requireNamespace("tidyselect", quietly = TRUE)) {
  install.packages("tidyselect")
}

if (! requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

if (! requireNamespace("abind", quietly = TRUE)) {
  install.packages("abind")
}

if (! requireNamespace("ade4", quietly = TRUE)) {
  install.packages("ade4")
}

if (! requireNamespace("Cairo", quietly = TRUE)) {
  install.packages("Cairo")
}

if (! requireNamespace("crosstalk", quietly = TRUE)) {
  install.packages("crosstalk")
}

if (! requireNamespace("deldir", quietly = TRUE)) {
  install.packages("deldir")
}

if (! requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi", update = TRUE)
}

if (! requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges", update = TRUE)
}

if (! requireNamespace("IRanges", quietly = TRUE)) {
  BiocManager::install("IRanges", update = TRUE)
}

if (! requireNamespace("askpass", quietly = TRUE)) {
  install.packages("askpass")
}

if (! requireNamespace("backports", quietly = TRUE)) {
  install.packages("backports")
}

if (! requireNamespace("beachmat", quietly = TRUE)) {
  BiocManager::install("beachmat")
}

if (! requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocManager::install("BiocParallel", update = TRUE)
}

if (! requireNamespace("BiocSingular", quietly = TRUE)) {
  BiocManager::install("BiocSingular", update = TRUE)
}

if (! requireNamespace("BiocNeighbors", quietly = TRUE)) {
  BiocManager::install("BiocNeighbors", update = TRUE)
}

if (! requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment", update = TRUE)
}

if (! requireNamespace("BiocGenerics", quietly = TRUE)) {
  BiocManager::install("BiocGenerics", update = TRUE)
}

if (! requireNamespace("littler", quietly = TRUE)) {
  BiocManager::install("littler", update = TRUE, ask = FALSE, force = TRUE)
}

if (! requireNamespace("rstan", quietly = TRUE)) {
  BiocManager::install("rstan", update = TRUE, ask = FALSE, force = TRUE)
}

if (! requireNamespace("StanHeaders", quietly = TRUE)) {
  BiocManager::install("StanHeaders", update = TRUE, ask = FALSE, force = TRUE)
}

if (! requireNamespace("DelayedArray", quietly = TRUE)) {
  BiocManager::install("DelayedArray", update = TRUE)
}

if (! requireNamespace("beachmat", quietly = TRUE)) {
  BiocManager::install("beachmat", update = TRUE)
}

if (! requireNamespace("ggbeeswarm", quietly = TRUE)) {
  install.packages("ggbeeswarm")
}

if (! requireNamespace("dendextend", quietly = TRUE)) {
  install.packages("dendextend")
}

if (! requireNamespace("dendsort", quietly = TRUE)) {
  install.packages("dendsort")
}

if (! requireNamespace("distr", quietly = TRUE)) {
  install.packages("distr")
}

if (! requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}

if (! requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

if (! requireNamespace("listenv", quietly = TRUE)) {
  install.packages("listenv")
}

if (! requireNamespace("lmtest", quietly = TRUE)) {
  install.packages("lmtest")
}

if (! requireNamespace("locfit", quietly = TRUE)) {
  install.packages("locfit")
}

if (! requireNamespace("markdown", quietly = TRUE)) {
  install.packages("markdown")
}

if (! requireNamespace("fastcluster", quietly = TRUE)) {
  install.packages("fastcluster")
}

if (! requireNamespace("future", quietly = TRUE)) {
  install.packages("future")
}

if (! requireNamespace("future.apply", quietly = TRUE)) {
  install.packages("future.apply")
}

if (! requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}

if (! requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

if (! requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}

if (! requireNamespace("goftest", quietly = TRUE)) {
  install.packages("goftest")
}

if (! requireNamespace("gplots", quietly = TRUE)) {
  install.packages("gplots")
}

if (! requireNamespace("graphlayouts", quietly = TRUE)) {
  install.packages("graphlayouts")
}

if (! requireNamespace("gridGraphics", quietly = TRUE)) {
  install.packages("gridGraphics")
}

if (! requireNamespace("grImport", quietly = TRUE)) {
  install.packages("grImport")
}

if (! requireNamespace("grImport2", quietly = TRUE)) {
  install.packages("grImport2")
}

if (! requireNamespace("gtools", quietly = TRUE)) {
  install.packages("gtools")
}

if (! requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}

if (! requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

if (! requireNamespace("ica", quietly = TRUE)) {
  install.packages("ica")
}

if (! requireNamespace("jpeg", quietly = TRUE)) {
  install.packages("jpeg")
}

if (! requireNamespace("markdown", quietly = TRUE)) {
  install.packages("markdown")
}

if (! requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

if (! requireNamespace("ROCR", quietly = TRUE)) {
  install.packages("ROCR")
}

if (! requireNamespace("Rtsne", quietly = TRUE)) {
  install.packages("Rtsne")
}

if (! requireNamespace("scattermore", quietly = TRUE)) {
  install.packages("scattermore")
}

if (! requireNamespace("sctransform", quietly = TRUE)) {
  install.packages("sctransform")
}

if (! requireNamespace("segmented", quietly = TRUE)) {
  install.packages("segmented")
}

if (! requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}

if (! requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}

if (! requireNamespace("SeuratObject", quietly = TRUE)) {
  install.packages("SeuratObject")
}

if (! requireNamespace("sfsmisc", quietly = TRUE)) {
  install.packages("sfsmisc")
}

if (! requireNamespace("spatstat.data", quietly = TRUE)) {
  install.packages("spatstat.data")
}

if (! requireNamespace("spatstat.explore", quietly = TRUE)) {
  install.packages("spatstat.explore")
}

if (! requireNamespace("spatstat.geom", quietly = TRUE)) {
  install.packages("spatstat.geom")
}

if (! requireNamespace("spatstat.random", quietly = TRUE)) {
  install.packages("spatstat.random")
}

if (! requireNamespace("spatstat.sparse", quietly = TRUE)) {
  install.packages("spatstat.sparse")
}

if (! requireNamespace("spatstat.utils", quietly = TRUE)) {
  install.packages("spatstat.utils")
}

if (! requireNamespace("tensor", quietly = TRUE)) {
  install.packages("tensor")
}

if (! requireNamespace("tidygraph", quietly = TRUE)) {
  install.packages("tidygraph")
}

if (! requireNamespace("tiff", quietly = TRUE)) {
  install.packages("tiff")
}

if (! requireNamespace("pbapply", quietly = TRUE)) {
  install.packages("pbapply")
}

if (! requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (! requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

if (! requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}

if (! requireNamespace("polyclip", quietly = TRUE)) {
  install.packages("polycliip")
}

if (! requireNamespace("progressr", quietly = TRUE)) {
  install.packages("progressr")
}

if (! requireNamespace("RANN", quietly = TRUE)) {
  install.packages("RANN")
}

if (! requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

if (! requireNamespace("miloR", quietly = TRUE)) {
  BiocManager::install("miloR")
}

if (! requireNamespace("grr", quietly = TRUE)) {
  install.packages("grr")
}

matrix_utils_url <- 'https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz'
if (! requireNamespace("Matrix.utils", quietly = TRUE)) {
  install.packages(matrix_utils_url, repos = NULL, type = "source")
}

if (! requireNamespace("Rsamtools", quietly = TRUE)) {
  BiocManager::install("Rsamtools")
}

if (! requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

if (! requireNamespace("VariantAnnotation", quietly = TRUE)) {
  BiocManager::install("VariantAnnotation")
}

if (! requireNamespace("GenomicFeatures", quietly = TRUE)) {
  BiocManager::install("GenomicFeatures")
}

if (! requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

if (! requireNamespace("SNPRelate", quietly = TRUE)) {
  BiocManager::install("SNPRelate", update = TRUE)
}

if (! requireNamespace("RcppTOML", quietly = TRUE)) {
  install.packages("RcppTOML")
}

# gsl package requires bcftools to be installed
if (! requireNamespace("gsl", quietly = TRUE)) {
  install.packages("gsl")
}

# Optional package for alleleIntegrator (not needed for inactiveXX but needed
# for one function in alleleIntegrator)
#devtools::install_github('VanLoo-lab/ascat/ASCAT')

devtools::install_github('constantAmateur/alleleIntegrator')

devtools::install_github('constantAmateur/inactiveXX')
