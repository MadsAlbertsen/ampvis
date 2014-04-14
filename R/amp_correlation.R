
#' Generates a ggplot2 style ordinate plot of amplicon data.
#'
#' A nice long description
#'
#' @usage amp_rabund(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param trans Transform the raw counts, currently only supports "none" or "sqrt" (default: "sqrt").
#' @param ordinate.type Either PCA or NMDS (default: "PCA").
#' @param ncomp The number of principal components to extract (default: 5)
#' @param plot.x Variable to plot on the x-axis (default: "PC1")
#' @param plot.y Variable to plot on the y-axis (default: "PC2")
#' @param plot.color Color the points by a sample variable.
#' @param plot.point.size Size of the plotted points (default: 3)
#' @param plot.species Plot loadings as points (default: False)
#' @param plot.nspecies Plot the n most extreme species with their genus classification (default: 0).
#' @param envfit.factor A vector of factor variables from the sample data used for envfit to the model.
#' @param envfit.numeric A vector of numerical variables from the sample data used for envfit to the model.
#' @param envfit.significant The significance treshold for displaying envfit parameters (default: 0.01).
#' @param tax.clean Add best assignment to Genus level classification if none exists.
#' @param output Either plot or complete (default: "plot").
#' 
#' @return A ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @import plyr
#' @import reshape2
#' @import phyloseq
#' @import grid
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

data <- V13f
cor.treshold <- 0.8

amp_ordinate <- function(data, trans = "sqrt"){
  
  ## Load the data
  
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
  
  ## Transform the data
  
  abund1 <- abund
  if (trans == "sqrt"){
    abund1 <- sqrt(abund)
  }
  
  ## Calculate correlation matrix
  
  cor <- cor(t(abund1))
  cor[lower.tri(cor)] <- NA
  cor2 <- melt(cor)
  cor3 <- cor2[!is.na(cor2$value), ]
  cor4 <- cor3[cor3$Var1 != cor3$Var2, ]
  cor5 <- subset(cor4, abs(value) > cor.treshold)

  
  g <- graph.data.frame(cor5, directed = F)
  
  V(g)$label.cex <- 0.7
  V(g)$label.color <- "black"
  V(g)$label.font <- 2
  V(g)$frame.color <- "black"
  V(g)$color <- "white"
  V(g)$size <- 10
  E(g)$width <- 2
  E(g)$color <- (E(g)$value + 1) * 100
  
  rwg <- colorRampPalette(c("red", "white", "green"))
  palette(rwg(100))
  
  plot(g)
  
  
}
