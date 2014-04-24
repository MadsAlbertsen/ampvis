#' Generates a ggplot2 style rank abundance plot of the most abundant OTUs 
#'
#' A nice long description
#'
#' @usage amp_rabund(data)
#'
#' @param data (required) A phyloseq object.
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 20000)
#' @param plot.type Either point or boxplot (default:point).
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
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rabund <- function(data, tax.show = 50, scale.seq = 20000, tax.clean = T, plot.type = "point", plot.log = F, output = "plot"){
  
  ## Extract all data from the phyloseq object
  
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
  
  ## Change Proteobacteria to Class level
  
  if(tax.clean == T){
    tax$Phylum <- as.character(tax$Phylum)
    tax$Class <- as.character(tax$Class)
    for (i in 1:nrow(tax)){
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] == "p__Proteobacteria"){
        tax$Phylum[i] <- tax$Class[i]   
      }
    }
    tax$Phylum <- gsub("p__", "", tax$Phylum)
    tax$Phylum <- gsub("c__", "", tax$Phylum)
    tax$Genus <- gsub("g__", "", tax$Genus)
    
    tax$Genus[is.na(tax$Genus)] <- ""
    tax$Phylum[is.na(tax$Phylum)] <- ""
  }
  
  ## Merge the taxonomic and abundance information
  
  abund2 <- cbind.data.frame(rownames(abund), tax, abund)
  colnames(abund2)[1] <- "OTU"
  
  ## Melt data to long format
  
  abund3 <- melt(abund2, id.var = c("OTU",colnames(tax)), measure.vars=rownames(sample))
  colnames(abund3)[ncol(abund3)] <- "Abundance" 
  colnames(abund3)[ncol(abund3)-1] <- "Sample" 
  
  ## Subset to X most abundant OTUs
  TotalCounts <- ddply(abund3, ~OTU+Genus+Phylum, summarise, Abundance = median(Abundance))
  TotalCounts <- TotalCounts[with(TotalCounts, order(-Abundance)),]
  abund4 <- subset(abund3, abund3$OTU %in% TotalCounts$OTU[1:tax.show])
  abund4$OTU <- factor(abund4$OTU, levels = rev(TotalCounts$OTU[1:tax.show])) 
  
  abund4$Abundance <- abund4$Abundance/scale.seq*100
  
  ## plot the data
  
p <-ggplot(abund4, aes(x = OTU, y = Abundance)) +
  coord_flip() +
  scale_x_discrete(labels = rev(paste(TotalCounts$Genus[1:tax.show], TotalCounts$Phylum[1:tax.show], sep = "; "))) +
  ylab("Abundance (%)")
  if (plot.type == "point"){
    p <- p + geom_point()
  }
  if (plot.type == "boxplot"){
    p <- p + geom_boxplot()
  }
  if (plot.log ==T){
   p <- p + scale_y_log10() 
  }
  
  outlist <- append(outlist, list(dataframe = abund4, plot = p))

  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }

}
