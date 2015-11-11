#' Load data as OTU table and convert to a phyloseq object.
#'
#' Load data as OTU table and convert to a phyloseq object.
#'
#' @usage amp_load(otutable, metadata)
#'
#' @param otutable (required) A OTU table generated from workflow scripts v.4+.
#' @param metadata (required) A metadata file with sample names in first column.
#' @param refseq Reference sequences for all OTUs.
#' @param rarefy Rarefy all samples to the same sequencing depth.
#' 
#' @return A phyloseq object.
#' 
#' @export
#' @import phyloseq
#' @import dplyr
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_load <- function(otutable, metadata, refseq = NULL, rarefy = NULL){
  otu_counts <- otutable[,1:(ncol(otutable)-7)]
  
  # Remove whitespace from the otutable as this will break the structure of the taxonomy
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  otutable$Kingdom<-trim(as.character(otutable$Kingdom))
  otutable$Phylum<-trim(as.character(otutable$Phylum))
  otutable$Class<-trim(as.character(otutable$Class))
  otutable$Order<-trim(as.character(otutable$Order))
  otutable$Family<-trim(as.character(otutable$Family))
  otutable$Genus<-trim(as.character(otutable$Genus))
  otutable$Species<-trim(as.character(otutable$Species))
  
  taxonomy <- otutable[,(ncol(otutable)-6):ncol(otutable)] %>% as.matrix()
  rownames(metadata) <- metadata[,1]
  
  if (!is.null(refseq)){
    data <- phyloseq(otu_table(otu_counts, taxa_are_rows=T),
                     tax_table(taxonomy),
                     sample_data(metadata),
                     refseq)
  } else {
    data <- phyloseq(otu_table(otu_counts, taxa_are_rows=T),
                     tax_table(taxonomy),
                     sample_data(metadata))
  }
    
  data <- transform_sample_counts(data, function(x) x / 1)
  if(!is.null(rarefy)){data <- rarefy_even_depth(data, sample.size = rarefy, rngseed = 712)}
  return(data)
}
