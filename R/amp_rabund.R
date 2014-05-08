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
#' @param name2 A taxonomic level used for naming (default: "Phylum").
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
#' @import data.table
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rabund <- function(data, group = "Sample", tax.show = 50, scale.seq = 20000, tax.clean = T, plot.type = "boxplot", plot.log = F, output = "plot", name2 = "Phylum", tax.aggregate = "OTU", tax.empty = "rename"){
  
  ## Extract all data from the phyloseq object
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(data.frame(sample_data(data)))
  tax <- data.frame(tax, OTU = rownames(tax))
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
  
  
  
  ## Clean the taxonomy
  if(tax.clean == T){
    for ( i in 1:ncol(tax)){
      tax[,i] <- as.character(tax[,i])  
    }
  #### Change Proteobacteria to Class level  
    for (i in 1:nrow(tax)){
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] == "p__Proteobacteria"){
        tax$Phylum[i] <- tax$Class[i]   
      }
    }
  
  #### Remove the greengenes level prefix
    tax$Phylum <- gsub("p__", "", tax$Phylum)
    tax$Phylum <- gsub("c__", "", tax$Phylum)
    tax$Class <- gsub("c__", "", tax$Class)
    tax$Order <- gsub("o__", "", tax$Order)
    tax$Family <- gsub("f__", "", tax$Family)
    tax$Genus <- gsub("g__", "", tax$Genus)
    tax[is.na(tax)] <- ""
    if (!is.null(tax$Species)){tax$Species <- gsub("s__", "", tax$Species)} 
  
  #### Handle empty taxonomic strings
    if(tax.empty == "rename"){  
      t2 <- tax
      a1 <- data.frame(OTU = as.character(t2[,"OTU"]), temp = as.character(t2[,tax.aggregate]), OTU1 = as.character(t2[,"OTU"]))
      a2 <- subset(a1, temp != "")[,1:2]
      colnames(a2)[2] <- tax.aggregate 
      a3 <- subset(a1, temp == "")[,c(1,3)]
      colnames(a3)[2] <- tax.aggregate
      a4 <- rbind(a2, a3)
      a5 <- t2[ , -which(names(t2) %in% tax.aggregate)]      
      a7 <- join(a5, a4, by = "OTU")
      tax <- a7
    }
    if(tax.empty == "remove"){
      tax <- subset(tax, tax[,tax.aggregate] != "")
      abund <- subset(abund, rownames(abund) %in% rownames(tax))
    }
  }
  
  ## Merge the taxonomic and abundance information
  
  abund2 <- cbind.data.frame(tax, abund)
  
  ## Melt data to long format
  abund3 <- melt(abund2, id.var = c(colnames(tax)), measure.vars=rownames(sample))
  colnames(abund3)[ncol(abund3)] <- "Abundance" 
  colnames(abund3)[ncol(abund3)-1] <- "Sample" 
  colnames(sample)[1] <- "Sample"
  
  ## Merge sample data
  temp1 <- join(x = abund3, y = data.frame(sample), by = "Sample")
  
  
  ## Summarise to specific taxonomic levels and groups using data.table for blazing speed
  colnames(temp1)[colnames(temp1) == tax.aggregate] <- "var1"
  colnames(temp1)[colnames(temp1) == name2] <- "var2"
  
  if(group != "Sample"){
    colnames(temp1)[colnames(temp1) == group] <- "var3"
    DT <- data.table(temp1)
    DT2 <- DT[, lapply(.SD, sum, na.rm=TRUE), by=list(var1, var2, var3, Sample), .SDcols=c("Abundance") ]
    temp2 <- data.frame(DT2)
    colnames(temp2)[colnames(temp2) == "var3"] <- group
  }
  if(group == "Sample"){
    DT <- data.table(temp1)
    DT2 <- DT[, lapply(.SD, sum, na.rm=TRUE), by=list(var1, var2, Sample), .SDcols=c("Abundance") ]
    temp2 <- data.frame(DT2)
  }
  
  colnames(temp2)[colnames(temp2) == "var1"] <- tax.aggregate
  colnames(temp2)[colnames(temp2) == "var2"] <- name2
  
  ## Subset to X most abundant "OTUs"
  TotalCounts <- ddply(temp2, c(tax.aggregate,name2), summarise, Abundance = median(Abundance))
  TotalCounts <- TotalCounts[with(TotalCounts, order(-Abundance)),]
  
  ## Make sure we only show a possible number of taxa
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){
      tax.show <- nrow(TotalCounts)
    }
  }
  if (!is.numeric(tax.show)){
    tax.show <- nrow(TotalCounts)
  }
    
  abund4 <- subset(temp2, temp2[, tax.aggregate] %in% TotalCounts[1:tax.show, tax.aggregate])
  abund4 <- droplevels(abund4)
  abund4[,tax.aggregate] <- factor(abund4[, tax.aggregate], levels = rev(TotalCounts[1:tax.show, tax.aggregate])) 
  
  ## Scale to a specific abundance
  abund4$Abundance <- abund4$Abundance/scale.seq*100
  
  ## plot the data
  
  if (group == "Sample"){
    p <-ggplot(abund4, aes_string(x = tax.aggregate, y = "Abundance"))   
  }
  if (group != "Sample"){
    p <-ggplot(abund4, aes_string(x = tax.aggregate, y = "Abundance", color = group))   
  }

  p <- p +
  coord_flip() +
  scale_x_discrete(labels = rev(paste(TotalCounts[1:tax.show, tax.aggregate], TotalCounts[1:tax.show, name2], sep = "; "))) +
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
