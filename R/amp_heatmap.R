#' Generates a ggplot2 style heatmap from amplicon data in phyloseq format 
#'
#' A nice long description
#'
#' @usage amp_headtmap(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param group A variable from the associated sample data to group samples by.
#' @param scale A variable from the associated sample data to scale the abundance by.
#' @param normalise A sample or group to normalise the counts to.
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: Phylum)
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.empty Either "remove" OTUs without taxonomic information or "rename" with OTU ID (default: rename).
#' @param order.x A taxonomy group or vector to order the x-axis by.
#' @param order.y A sample or vector to order the y-axis by.
#' @param plot.numbers Plot the values on the heatmap (default: T)
#' @param plot.breaks A vector of breaks for the abundance legend.
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 20000)
#' @param output To output a plot or the complete data inclusive dataframes (default: plot)
#' 
#' @return A ggplot2 object or a list with the plot and associated dataframes.
#' 
#' @export
#' @import ggplot2
#' @import plyr
#' @import reshape2
#' @import phyloseq
#' @import grid
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_heatmap <- function(data, group = NULL, normalise = NULL, scale = NULL, tax.aggregate = "Phylum", tax.show = 10, tax.clean = T, tax.empty = "rename", order.x = NULL, order.y = NULL, plot.numbers = T, plot.breaks = NULL, scale.seq = 20000, output = "plot"){
  
  ## Extract all data from the phyloseq object
  
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  
  ## Extract group information
  if (!is.null(group)){
    grp <- cbind.data.frame(rownames(sample),sample[,group])
    colnames(grp) <- c("Sample", group)  
  }
  if (is.null(group)){
    grp <- cbind.data.frame(rownames(sample),rownames(sample))
    colnames(grp) <- c("Sample", "Samples")  
  }
  
  ## Scale the data by a select variable

  if (!is.null(scale)){
    variable <- as.numeric(sample[,scale])
    abund <- t(t(abund)*variable)
  }
  
  ## Change Proteobacteria to Class level
  
  if(tax.clean == T){
    for ( i in 1:ncol(tax)){
      tax[,i] <- as.character(tax[,i])  
    }
    
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
  
  ## Aggregate to a specific taxonomic level
  
  abund3 <- melt(abund2, id.var = tax.aggregate, measure.vars=rownames(sample))
  colnames(abund3)[2] <- "Sample"
  colnames(abund3)[3] <- "Abundance"
  
  colnames(abund3)[colnames(abund3) == tax.aggregate] <- "var1"
  DT <- data.table(abund3)
  DT2 <- DT[, lapply(.SD, sum, na.rm=TRUE), by=list(var1, Sample), .SDcols=c("Abundance") ]   
  abund4 <- data.frame(DT2)
  colnames(abund4)[colnames(abund4) == "var1"] <- tax.aggregate
  
  abund4 <- subset(abund4, abund4[,tax.aggregate] != "<NA>")
  
  ## Add group information
  
  abund5 <- join(abund4, grp, by="Sample")
  
  ## Take the average to group level
  
  colnames(abund5)[colnames(abund5) == tax.aggregate] <- "var1"
  colnames(abund5)[colnames(abund5) == colnames(grp[2])] <- "var2"
  DT3 <- data.table(abund5)
  DT4 <- DT3[, lapply(.SD, mean, na.rm=TRUE), by=list(var1, var2), .SDcols=c("Abundance") ]   
  abund6 <- data.frame(DT4)
  colnames(abund6)[colnames(abund6) == "var1"] <- tax.aggregate
  colnames(abund6)[colnames(abund6) == "var2"] <- colnames(grp[2])
  
  ## Find the X most abundant levels
  
  TotalCounts <- ddply(abund6, tax.aggregate, summarise, Abundance = sum(Abundance))
  TotalCounts <- TotalCounts[with(TotalCounts, order(-Abundance)),]
  
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){
      tax.show <- nrow(TotalCounts)
    }
    abund7 <- subset(abund6, abund6[,1] %in% TotalCounts[1:tax.show,tax.aggregate])  
  }
  if (!is.numeric(tax.show)){
    abund7 <- subset(abund6, abund6[,1] %in% tax.show)  
  }

  ## Normalise to a specific group (The Abundance of the group is set as 1)

  if(!is.null(normalise)){
    colnames(abund7) <- c("var1", "var2", "Abundance")
    temp <- dcast(abund7, var1~var2, value.var = "Abundance")
    colnames(temp)[1] <- tax.aggregate
    temp2 <- temp[,-1]  
    temp3 <- temp2/temp2[,normalise]
    temp4 <- cbind.data.frame(temp[,1], temp3)    
    temp5 <- melt(temp4)
    colnames(temp5) <- c(tax.aggregate, colnames(grp[2]), "Abundance")
    abund7 <- temp5
  }
  abund7$Abundance <- round(abund7$Abundance, 1)
  
  
  ## Order.y
  
  if (is.null(order.y)){
    TotalCounts <- ddply(abund7, tax.aggregate, summarise, Abundance = sum(Abundance))
    TotalCounts <- TotalCounts[with(TotalCounts, order(Abundance)),]
    abund7[,1] <- factor(abund7[,1], levels = TotalCounts[,1])
    }
  
  if (!is.null(order.y)){
    if (length(order.y) == 1){
      temp1 <- subset(abund7, abund7[,2] %in% order.y)
      temp1 <- temp1[with(temp1, order(temp1[,3])),]
      abund7[,1] <- factor(abund7[,1], levels = temp1[,1])
    }
    if (length(order.y) > 1){
      abund7[,1] <- factor(abund7[,1], levels = order.y)
    }
  }
  
  ## Order.x
  
  if (!is.null(order.x)){
    if (length(order.x) == 1){
      temp1 <- subset(abund7, abund7[,1] %in% order.x)
      temp1 <- temp1[with(temp1, order(temp1[,3])),]
      abund7[,2] <- factor(abund7[,2], levels = temp1[,2])
    }    
    if (length(order.x) > 1){
      abund7[,2] <- factor(abund7[,2], levels = order.x)
    }
  }

  ## Scale to percentages if not normalised and scaled
  
  if (is.null(scale) & is.null(normalise)){
    abund7[,3] <- round(abund7[,3]/scale.seq*100,1)
  }
  
  ## Make a heatmap style plot
  
  p <- ggplot(abund7, aes_string(x = colnames(grp[2]), y = tax.aggregate, label = formatC("Abundance", format = "f", digits = 1))) + 
    geom_tile(aes(fill = Abundance), colour = "white", size = 0.5) + 
    labs(x = "", y = "", fill = "Abundance") + 
    theme(axis.text.x = element_text(size = 10, hjust = 1, angle = 90)) + 
    theme(axis.text.y = element_text(size = 12)) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    if (plot.numbers == T){
      p <- p + geom_text(size = 4, colour = "grey30")  
    }
  if (is.null(plot.breaks)){
    p <- p +scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), trans = "log10")
  }
  if (!is.null(plot.breaks)){
      p <- p +scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), trans = "log10", breaks=plot.breaks)
    }
  
  ## Define the output 
  
  if (output == "complete"){
    outlist <- list(heatmap = p, sampledata = sample, taxonomy = tax, abundance = as.data.frame(otu_table(data)), data = abund7)
    return(outlist)  
  }
  if (output == "plot"){
    return(p)
  } 
}