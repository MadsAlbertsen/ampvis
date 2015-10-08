#' Calculate basic statistics for each sample
#'
#' Calculate basic statistics for each sample and combine it with metadata.
#'
#' @usage amp_stats(data)
#'
#' @param data (required) A phyloseq object.
#' @param measure Alpha-diversity measures to be included (default:observed).
#' 
#' @export
#' @import phyloseq
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_stats <- function(data, measures = "Observed"){
  Reads <- sample_sums(data)
  dmeta <- sample_data(data)
  rich <- estimate_richness(data, measures = measures)
  
  combined <- cbind.data.frame(dmeta, Reads, rich) %>% 
    arrange(Reads)
  return(combined)
}
