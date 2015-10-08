#' Calculate rarefaction curve for each sample.
#'
#' Calculate rarefaction curve for each sample using the vegan rarecurve function directly from a phyloseq object.
#'
#' @usage amp_rarecurve(data)
#'
#' @param data (required) A phyloseq object.
#' @param step Step size for sample sizes in rarefaction curves (default: 100).
#' @param ylim vector of y-axis limits.
#' @param xlim vector of x-axis limits.
#' @param label Label rarefaction curves (default: F).
#' 
#' @export
#' @import phyloseq
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rarecurve <- function(data, step = 100, ylim = NULL, xlim = NULL, label = F){
  
  abund = otu_table(data)@.Data %>% as.data.frame()
  
  if (is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, label = label)
  }
  if (!is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, xlim = xlim, label = label)
  }
  if (!is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, label = label)
  }
  if (is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, xlim = xlim, label = label)
  }
}
