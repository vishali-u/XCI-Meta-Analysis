# XCI-Meta-Analysis

## Main scripts for analysis:
1. `format_metadata.R` prepares the metadata for HCA datasets \
    a. Uses `formatMetadataHCA` function which was taken from [inactiveXX paper code](https://github.com/constantAmateur/XiPaperCode/blob/main/prepData.R)
3. `download_fastqs.R` downloads the HCA fastq files
4. `align_data.R` aligns the HCA fastq files
5. `run_inactiveXX.R` runs inactiveXX on the BAM files
6. `summarize_inactiveX_results.R` creates a filtered data frame including females that inactiveXX worked for
7. `prepare_srat.R` filter and normalize raw counts, return a Seurat object with the XCI states added
8. `de_analysis.R` run differential expression analysis

## Other scripts
1. `install_inactiveXX.R` includes most (if not all) packages needed to install inactiveXX on Compute Canada

## Results
1. supplementary tables: \
    a. table 1 - inactiveXX results \
    b. table 2 - differential expression analysis results
