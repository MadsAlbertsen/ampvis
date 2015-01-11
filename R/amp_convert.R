#' Convert a phyloseq object to a list of dataframes.
#'
#' Convert a phyloseq object to a list of dataframes.
#'
#' @usage amp_convert(data)
#'
#' @param data (required) A phyloseq object.
#' 
#' @return A list of dataframes.
#' 
#' @export
#' @import phyloseq
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_convert <- function(data){
  data <- list(abund = as.data.frame(otu_table(data)),
               tax = data.frame(tax_table(data), OTU = rownames(tax_table(data))),
               sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  return(data)
}
