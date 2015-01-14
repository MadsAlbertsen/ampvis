#' Export sequences from a phyloseq object
#'
#' Export sequences from a phyloseq object
#'
#' @usage amp_export(data)
#'
#' @param data (required) A phyloseq object.
#' 
#' @export
#' @import phyloseq
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_export <- function(data, file = "exported_sequences.fa"){
  writeXStringSet(refseq(data), file = file)
  
}
