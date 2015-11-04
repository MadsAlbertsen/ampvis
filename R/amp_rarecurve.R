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
#' @param color Color lines by metadata.
#' 
#' @export
#' @import phyloseq
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rarecurve <- function(data, step = 100, ylim = NULL, xlim = NULL, label = F, color = NULL){
  
  abund = otu_table(data)@.Data %>% as.data.frame()
  
if (!is.null(color)) {
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  group_vector<-sample_data(data)[,color]@.Data %>% as.data.frame()
  names(group_vector)<-"color_variable"
  group_vector<-as.character(group_vector$color_variable)
  groups<-unique(group_vector)
  n = length(groups)
  cols = gg_color_hue(n)
  
  col_vector<-rep("black",length(group_vector))
  for (i in 1:length(group_vector)){
    col_vector[i]<-cols[match(group_vector[i],groups)]
  }
} else {
  cols="black"
}
  
  if (is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, label = label, col = cols)
  }
  if (!is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, xlim = xlim, label = label, col = cols)
  }
  if (!is.null(ylim) & is.null(xlim)){
    rarecurve(t(abund), step = step, ylim = ylim, label = label, col = cols)
  }
  if (is.null(ylim) & !is.null(xlim)){
    rarecurve(t(abund), step = step, xlim = xlim, label = label, col = cols)
  }
}
