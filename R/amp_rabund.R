#' Generates a ggplot2 style rank abundance plot of the most abundant OTUs 
#'
#' A nice long description
#'
#' @usage amp_rabund(data)
#'
#' @param data (required) A phyloseq object.
#' @param group Group the data based on a sample variable.
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: OTU)
#' @param tax.empty Either "remove" OTUs without taxonomic information or "rename" with OTU ID (default: rename).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 20000)
#' @param plot.type Either point or boxplot (default:point).
#' @param names Two taxonomic levels used for naming, supplied as a vector (default: c("Genus", "Phylum")).
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

amp_rabund <- function(data, group = NULL, tax.show = 50, scale.seq = 20000, tax.clean = T, plot.type = "boxplot", plot.log = F, output = "plot", names = c("Genus", "Phylum"), tax.aggregate = "OTU", tax.empty = "rename"){
  
  ## Extract all data from the phyloseq object
  
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(data.frame(sample_data(data)))
  
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
    tax$Class <- gsub("c__", "", tax$Class)
    tax$Order <- gsub("o__", "", tax$Order)
    tax$Family <- gsub("f__", "", tax$Family)
    tax$Genus <- gsub("g__", "", tax$Genus)
    tax[is.na(tax)] <- ""
    if (!is.null(tax$Species)){tax$Species <- gsub("s__", "", tax$Species)} 
    
    if(tax.empty == "rename"){
      for (i in 1:nrow(tax)){
        if (tax[i,tax.aggregate] == ""){
          tax[i,tax.aggregate] <- rownames(tax)[i]
        }
      }    
    }
    
    if(tax.empty == "remove"){
      tax <- subset(tax, tax[,tax.aggregate] != "")
      abund <- subset(abund, rownames(abund) %in% rownames(tax))
    }
  }
  
  ## Merge the taxonomic and abundance information
  
  abund2 <- cbind.data.frame(rownames(abund), tax, abund)
  colnames(abund2)[1] <- "OTU"
  
  ## Melt data to long format
  abund3 <- melt(abund2, id.var = c("OTU",colnames(tax)), measure.vars=rownames(sample))
  colnames(abund3)[ncol(abund3)] <- "Abundance" 
  colnames(abund3)[ncol(abund3)-1] <- "Sample" 
  colnames(sample)[1] <- "Sample"
  
  ## Add sample data
  temp1 <- merge(x = abund3, y = data.frame(sample), by.x = "Sample", by.y = "Sample")
  
  ## Tax aggregate
  
  temp2 <- ddply(temp1, c(tax.aggregate, group, names[1], names[2], "Sample"), summarise, Abundance = sum(Abundance))
  
  ## Subset to X most abundant "OTUs"
  TotalCounts <- ddply(temp2, c(tax.aggregate,names[1],names[2]), summarise, Abundance = median(Abundance))
  TotalCounts <- TotalCounts[with(TotalCounts, order(-Abundance)),]
  
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){
      tax.show <- nrow(TotalCounts)
    }
  }
  if (!is.numeric(tax.show)){
    tax.show <- nrow(TotalCounts)
  }
  
  abund4 <- subset(temp2, temp2[, tax.aggregate] %in% TotalCounts[1:tax.show, tax.aggregate])
  abund4[,tax.aggregate] <- factor(abund4[, tax.aggregate], levels = rev(TotalCounts[1:tax.show, tax.aggregate])) 
  
  
  abund4$Abundance <- abund4$Abundance/scale.seq*100
  
  ## plot the data
  
  if (is.null(group)){
    p <-ggplot(abund4, aes_string(x = tax.aggregate, y = "Abundance"))   
  }
  if (!is.null(group)){
    p <-ggplot(abund4, aes_string(x = tax.aggregate, y = "Abundance", color = group))   
  }

  p <- p +
  coord_flip() +
  scale_x_discrete(labels = rev(paste(TotalCounts[1:tax.show, names[1]], TotalCounts[1:tax.show, names[2]], sep = "; "))) +
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
