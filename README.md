# XCI-Meta-Analysis

## Main scripts for analysis:
1. `format_metadata.R` prepare the metadata
2. `download_fastqs.R` download the HCA fastq files
3. `align_data.R` align the HCA fastq files
4. `run_inactiveXX.R` run inactiveXX on the BAM files
5. `summarize_inactiveX_results.R` create a filtered dataframe only including females that inactiveXX worked for
6. `prepare_srat.R` filter and normalize raw counts, return a Seurat object with the XCI states added
7. `de_analysis.R` run differential expression analysis

## Other scripts
1. `install_inactiveXX.R` includes most (if not all) packages need to use inactiveXX (on Compute Canada)

## Results
1. supplementary tables:
    a. table 1 - inactiveXX results
    b. table 2 - differential expression analysis results
