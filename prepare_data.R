# Prepare data for the inactiveXX pipeline by downloading and aligning the 
# fastq files

####### 1. Libraries #######
library(data.table)
library(dplyr)

####### 2. Constants #######






####### 3. Functions? #######
# The formatMetadataHCA function came from XiPaperCode (https://github.com/constantAmateur/XiPaperCode/tree/main)

#' Format HCA metadata
#'
#' HCA metadata is provided in a somewhat standard format.  Check the assumptions we're making and try and auto-generate cellranger compatible file names.  The default assumption here is that the "sequencing input", given by sequencing_input.biomaterial_core.biomaterial_id, is the lowest level of grouping that signifies a collection of files that should be processed together.
#'
#' Note that for almost all metadata sequencing_input.biomaterial_core.biomaterial_id and cell_suspension.biomaterial_core.biomaterial_id are the same thing.  But there are a subset where they are not and the sequencing_input seems the more dependable (and logically is the lower level) of the two.
#'
#' @param mDat data.frame containing HCA metadata.
#' @param dataName Unique name for this data.
#' @param useDonorName Use the \code{donorName} string to name donors.  If false, name numerically.
#' @param channelName What to use to name the channels.
#' @param donorName What to use to name the channels.
#' @return Updated version of mDat
formatMetadataHCA = function(mDat,dataName,useDonorName=FALSE,channelName = mDat$sequencing_input.biomaterial_core.biomaterial_id,donorName = mDat$donor_organism.biomaterial_core.biomaterial_id){
  #Decide on namings 
  if(useDonorName){
    humanDonorName = sprintf('%s_%s',dataName,donorName)
  }else{
    humanDonorName = sprintf('%s%03d',dataName,match(donorName,unique(donorName)))
  }
  mDat$humanDonorName = humanDonorName
  mDat$groupName = channelName
  #We only want the fastq files
  mDat = mDat[mDat$file_type == 'sequence_file' &
                mDat$file_format %in% c('fastq','fastq.gz','fq.gz') &
                mDat$read_index %in% c('read1','read2','index1'),]
  message('Found library prep methods')
  print(table(mDat$library_preparation_protocol.library_construction_approach))
  message('Derived from')
  print(table(mDat$library_preparation_protocol.nucleic_acid_source))
  message("Sequencing input types")
  print(table(mDat$sequencing_input_type))
  #The basic expectation I'm making is that per "sequencing input", there is one donor, specimen, sample, library prep protocol, cell suspension
  x = split(mDat,mDat$groupName)
  #First check whatever we're using for channel/donor is unique
  if(!all(lengths(lapply(split(mDat$groupName,mDat$sequencing_input.biomaterial_core.biomaterial_id),unique))==1))
    stop("Specified channelName is not unique by channel!")
  if(!all(lengths(lapply(split(mDat$humanDonorName,mDat$sequencing_input.biomaterial_core.biomaterial_id),unique))==1))
    stop("Specified donorName is not unique by channel!")
  if(!all(sapply(x,function(e) length(unique(e$donor_organism.biomaterial_core.biomaterial_id)))==1))
    stop("Some sequencing inputs have more than one donor")
  if(!all(sapply(x,function(e) length(unique(e$specimen_from_organism.provenance.document_id)))==1))
    stop("Some sequencing inputs have more than one specimen")
  if(!all(sapply(x,function(e) length(unique(e$sample.biomaterial_core.biomaterial_id)))==1))
    stop("Some sequencing inputs have more than one sample")
  if(!all(sapply(x,function(e) length(unique(e$library_preparation_protocol.library_construction_approach)))==1))
    stop("Some sequencing inputs have more than one library prep protocol")
  if(!all(sapply(x,function(e) length(unique(e$cell_suspension.biomaterial_core.biomaterial_id)))==1))
    stop("Some sequencing inputs have more than one cell suspension")
  #Now check if we can uniquely resolve file name for all
  #First assumption is that within each, sequencing_process.provenance.document_id, there's one set of reads.
  if(!all(sapply(split(mDat$read_index,mDat$sequencing_process.provenance.document_id),function(e) all(table(e)==1)))){
    warning("More than one group of reads per sequencing process.  You'll have to manually build file names.")
    return(mDat)
  }else{
    message("Auto constructing cellranger compatible filenames")
    #Construct name if we can
    for(nom in names(x)){
      tgts = x[[nom]]
      tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                              nom,
                              match(tgts$sequencing_process.provenance.document_id,unique(tgts$sequencing_process.provenance.document_id)),
                              c(read1='R1',read2='R2',index1='I1')[tgts$read_index]
      )
      x[[nom]] = tgts
    }
    mDat = do.call(rbind,x)
  }
  return(mDat)
}






####### 4. Baron et. al, 2016 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/f86f1ab4-1fbb-4510-ae35-3ffd752d4dfc
# Paper: https://pubmed.ncbi.nlm.nih.gov/27667365/
# 1 female control
# in-drop

# Data named after last name of the first author of the project
# Links to file manifest should be updated regularly (otherwise - HTPS error)

prepBaronData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22f86f1ab4-1fbb-4510-ae35-3ffd752d4dfc%22%5D%7D%7D&objectKey=manifests%2Fcc66d913-45ac-5934-9c08-5ab495e34f66.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/baron_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 5. Shrestha et. al, 2021 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/daa371e8-1ec3-43ef-924f-896d901eab6f
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8492318/
# 3 female controls
# 10X 3' v2

prepareShresthaData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22daa371e8-1ec3-43ef-924f-896d901eab6f%22%5D%7D%7D&objectKey=manifests%2Fb074f643-a801-5ec0-bb65-77b5a15049d1.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/shrestha_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 6. Fang et. al, 2019 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/1c6a960d-52ac-44ea-b728-a59c7ab9dc8e
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6573026/
# 1 female control
# drop-seq

prepareFangData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%221c6a960d-52ac-44ea-b728-a59c7ab9dc8e%22%5D%7D%7D&objectKey=manifests%2F8b780eb3-73f4-55eb-bb59-9e295a33a431.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/fang_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 7. Lawlor et. al, 2017 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/c6ad8f9b-d26a-4811-b2ba-93d487978446
# Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5287227/
# 2 female controls
# fluidigm C1-based library preparation

prepareLawlorData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22c6ad8f9b-d26a-4811-b2ba-93d487978446%22%5D%7D%7D&objectKey=manifests%2Fab9ee984-8dfa-513b-b0de-010475b6ace4.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/lawlor_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$library_preparation_protocol.nucleic_acid_source == 'single cell']
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 8. Wang et. al, 2016 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/99101928-d9b1-4aaf-b759-e97958ac7403
# Paper: https://pubmed.ncbi.nlm.nih.gov/27364731/
# 2 female controls
# Smart-seq2

prepareWangData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2299101928-d9b1-4aaf-b759-e97958ac7403%22%5D%7D%7D&objectKey=manifests%2F5a6d2219-ba5f-5f41-9287-683428d6343c.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/wang_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$library_preparation_protocol.nucleic_acid_source == 'single cell']
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 9. Enge et. al, 2017 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/cddab57b-6868-4be4-806f-395ed9dd635a
# Paper: https://pubmed.ncbi.nlm.nih.gov/28965763/
# 2 female controls
# Smart-seq2

prepareEngeData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22cddab57b-6868-4be4-806f-395ed9dd635a%22%5D%7D%7D&objectKey=manifests%2Fb59df8c9-61c2-5e07-b5b6-1175f4ac6892.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/enge_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$library_preparation_protocol.nucleic_acid_source == 'single cell']
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 10. Segerstolpe et. al, 2016 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/ae71be1d-ddd8-4feb-9bed-24c3ddb6e1ad
# Paper: https://europepmc.org/article/MED/27667667
# 1 female control
# Smart-seq2

prepareSegerstolpeData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22ae71be1d-ddd8-4feb-9bed-24c3ddb6e1ad%22%5D%7D%7D&objectKey=manifests%2F648d58bb-9b6a-5af3-861d-f177fd525677.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/segerstolpe_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$library_preparation_protocol.nucleic_acid_source == 'single cell']
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 11. Lee et. al, 2021 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/c5ca43aa-3b2b-4216-8eb3-f57adcbc99a1
# Paper: https://www.biorxiv.org/content/10.1101/2021.04.05.438347v2.full.pdf
# 6 female control
# CITE-seq

prepareLeeData <- function() {
  fileManifest <- 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp30&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22c5ca43aa-3b2b-4216-8eb3-f57adcbc99a1%22%5D%7D%7D&objectKey=manifests%2Ff702f2fe-f247-59f4-9080-4504134b2dcd.ca5182c0-4376-57e5-a328-f00b9a3a3252.tsv'
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/lee_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$library_preparation_protocol.nucleic_acid_source == 'single cell']
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}

####### 12. Quake, 2021 (pancreas) #######
# Project: https://data.humancellatlas.org/explore/projects/10201832-7c73-4033-9b65-3ef13d81656a
# Paper: https://www.biorxiv.org/content/10.1101/2021.07.19.452956v2.full
#  female control
# 

prepareQuakeData <- function() {
  fileManifest <- ''
  fileDest <- '/home/umaiyal1/XCI-Meta-Analysis/pancreas/quake_manifest.tsv'
  download.file(url = fileManifest, destfile = fileDest)
  
  mDat <- fread(fileDest)
  mDat <- mDat[mDat$file_format == 'fastq.gz',]
  mDat <- mDat[mDat$library_preparation_protocol.nucleic_acid_source == 'single cell']
  mDat <- mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
  mDat <- mDat[mDat$donor_organism.sex == 'female',]
  mDat <- mDat[mDat$donor_organism.diseases == 'normal',]
  
  
}








