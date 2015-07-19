#' Export otutable from a phyloseq object
#'
#' Export otutable from a phyloseq object.
#'
#' @usage amp_export_table(data)
#'
#' @param data (required) A phyloseq object.
#' @param file Name of the file containing the exported otutable.
#' @param id Name the samples using any variable in the metadata.
#' @param sort.samples Vector to sort the samples by.
#' 
#' @export
#' @import phyloseq
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export_table <- function(data, file = "exported_otutable.txt", id = NULL, sort.samples = NULL){
  
  samples <- data.frame(otu_table(data)@.Data)
  
  if(!is.null(id)){
    colnames(samples) <- as.character(unlist(sample_data(data)[,id]))
  }
  
  if(!is.null(sort.samples)){
    samples <- samples[,sort.samples]
  }
  
  e_bak <- cbind.data.frame(samples, 
                            data.frame(tax_table(data)@.Data, 
                                       OTU = rownames(tax_table(data))))
  
  e_bak2 <- mutate(e_bak, 
                   sum = rowSums(e_bak[,1:nrow(sample_data(data))])) %>%
    arrange(desc(sum)) %>%
    select(-sum)
  
  write.table(e_bak2, file = file, quote = F, row.names = F, sep = "\t")
}
