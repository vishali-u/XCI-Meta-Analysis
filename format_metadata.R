# First step: prepare metadata. Create cellranger filenames. Create donor and
# sample IDs. Create filepaths for mapped files (BAM files and matrices).
# Create a file path for a seurat object for each donor. Save the formatted 
# metadata. 

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
formatMetadataHCA = function(mDat,dataName,useDonorName=FALSE,
                             channelName = mDat$sequencing_input.biomaterial_core.biomaterial_id,
                             donorName = mDat$donor_organism.biomaterial_core.biomaterial_id){
  #Decide on namings 
  if(useDonorName){
    humanDonorName = sprintf('%s_%s',dataName,donorName)
  }else{
    humanDonorName = sprintf('%s%03d',dataName,match(donorName,unique(donorName)))
  }
  mDat$humanDonorName = humanDonorName
  # V: added gsub to substitute spaces with underscores
  mDat$groupName = gsub(' ', '_', channelName)
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

#' Add the location of the bam files and count files to the mDat object
#' and return the object
#' 
#' Assumes dest is a valid location.
#' 
#' @param mDat all the metadata for a dataset
#' @param useCellRanger a boolen set to TRUE if the data was aligned using
#'        cellranger (default is set to FALSE)
#' @param dataDir the parent directory where the fastq and mapped files will
#'                be stored
#' @return metadata including file paths for mapped files and a seurat object
addMappedFiles <- function(mDat, useCellRanger = FALSE, dataDir) {
  
  if (useCellRanger) {
    bamPath <- 'outs/possorted_genome_bam.bam'
    countsPath <- 'outs/filtered_gene_bc_matrices/GRCh38'
  } else {
    # Update once working with STAR
    bamPath <- 'Aligned.sortedByCoord.out.bam'
    countsPath <- 'Solo.out/Gene/filtered'
  }
  
  mDat$bamFile <- file.path(dataDir, mDat$humanDonorName, mDat$groupName, 
                            'mapped', bamPath)
  mDat$countsPath <- file.path(dataDir, mDat$humanDonorName, mDat$groupName, 
                               'mapped', countsPath)
  
  mDat$sratPath <- file.path(dataDir, mDat$humanDonorName, mDat$groupName,
                             sprintf('%s_%s_srat.rds', 
                                     mDat$humanDonorName, 
                                     mDat$groupName))
  
  return(mDat)
}