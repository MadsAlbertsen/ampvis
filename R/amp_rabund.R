#' Generates a ggplot2 style rank abundance plot of the most abundant OTUs 
#'
#' A nice long description
#'
#' @usage amp_rabund(data)
#'
#' @param data (required) A phyloseq object.
#' @param group Group the data based on a sample variable.
#' @param order.group A vector defining the order of groups.
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: OTU)
#' @param tax.add Additional taxonomic levels to display for each entry (default: Phylum) 
#' @param tax.empty Either "remove" OTUs without taxonomic information or "rename" with OTU ID (default: rename).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 10000)
#' @param plot.type Either point, boxplot or curve (default: boxplot).
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

amp_rabund <- function(data, group = "Sample", order.group = NULL, tax.show = 50, scale.seq = 10000, tax.clean = T, plot.type = "boxplot", plot.log = F, output = "plot", tax.add = NULL, tax.aggregate = "Genus", tax.empty = "rename"){
  
  ## Extract all data from the phyloseq object
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(data.frame(sample_data(data)))
  tax <- data.frame(tax, OTU = rownames(tax))
  if (is.null(tax$Species)){tax$Species <- ""}      
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
  
  ## Temporary fix to handle showing just 1 taxonomic level
  tax.add2 <- "something"
  if (is.null(tax.add)){
    tax.add2 <- NULL
    tax.add <- "Kingdom"
  }
  
  ## Clean up the taxonomy
  for ( i in 1:ncol(tax) ){
    tax[,i] <- as.character(tax[,i])  
  }
  
  ## Change Proteobacteria to Class level  
  if(tax.clean == T){
    for (i in 1:nrow(tax)){
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] == "p__Proteobacteria"){
        tax$Phylum[i] <- tax$Class[i]   
      }
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
  
  ## How to handle empty taxonomic assignments
  if(tax.empty == "rename"){
    tax[tax$Phylum == "","Phylum"] <- "Unclassified"
    for (i in 1:nrow(tax)) {   
        if (tax[i,"Species"] == "") {
          if (tax[i,"Genus"] != "") { rn <- paste("g__", tax[i,"Genus"], "_", tax[i,"OTU"], sep = "") } else{
            if (tax[i,"Family"] != "") { rn <- paste("f__", tax[i,"Family"], "_", tax[i,"OTU"], sep = "") } else{
              if (tax[i,"Order"] != "") { rn <- paste("o__", tax[i,"Order"], "_", tax[i,"OTU"], sep = "") } else{
                if (tax[i,"Class"] != "") { rn <- paste("c__", tax[i,"Class"], "_", tax[i,"OTU"], sep = "") } else{
                  if (tax[i,"Phylum"] != "") { rn <- paste("p__", tax[i,"Phylum"], "_", tax[i,"OTU"], sep = "") }
                }
              }
           }
         }
       }
      tax[i,tax[i,] == ""] <- rn
    }
  }
  
  if(tax.empty == "remove"){
    tax <- subset(tax, tax[,tax.aggregate] != "")
    abund <- subset(abund, rownames(abund) %in% rownames(tax))
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
  colnames(temp1)[colnames(temp1) == tax.add] <- "var2"
  
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
  colnames(temp2)[colnames(temp2) == "var2"] <- tax.add
  
  if (plot.type != "curve"){
  
  ## Subset to X most abundant "OTUs"
  TotalCounts <- ddply(temp2, c(tax.aggregate,tax.add), summarise, Abundance = median(Abundance))
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
    if(!is.null(order.group)){
      abund4[,group] <- factor(abund4[,group], levels = rev(order.group))
    }
    p <-ggplot(abund4, aes_string(x = tax.aggregate, y = "Abundance", color = group))   
  }

  p <- p +
  coord_flip() +  
  ylab("Abundance (%)") + 
  guides(col = guide_legend(reverse = TRUE))
  
  if (!is.null(tax.add2)){
    p <- p + scale_x_discrete(labels = rev(paste(TotalCounts[1:tax.show, tax.aggregate], TotalCounts[1:tax.show, tax.add], sep = "; ")))
  } else {
    p <- p + scale_x_discrete(labels = rev(TotalCounts[1:tax.show, tax.aggregate]))
  }
    
  
  if (plot.type == "point"){
    p <- p + geom_point() 
  }
  if (plot.type == "boxplot"){
    p <- p + geom_boxplot()
  }
  if (plot.log ==T){
   p <- p + scale_y_log10() 
  }
  }
  
  ## If type = curve then generate a second dataframe
  
  if (plot.type == "curve"){
    temp2$Abundance <- temp2$Abundance/scale.seq*100
    
    if (group != "Sample"){      
      temp3 <- ddply(temp2, c(tax.aggregate, group), summarise, Mean = mean(Abundance))
      temp3 <- temp3[with(temp3, order(-Mean)),]
      test <- ddply(temp3, group, transform , dummy = 1)
      TotalCounts <- ddply(test, group, transform, Cumsum = cumsum(Mean), Rank = cumsum(dummy))
      if(!is.null(order.group)){
        TotalCounts[,group] <- factor(TotalCounts[,group], levels = rev(order.group))
      }
      
      p <- ggplot(data = TotalCounts, aes_string(x = "Rank", y = "Cumsum", color = group)) +
        geom_line(size = 2) +
        ylim(0,100) +
        xlab("Rank abundance") +
        ylab("Cummulative read abundance (%)")  
      if (plot.log ==T){
        p <- p + scale_x_log10() 
      }  
    }
    
    if (group == "Sample"){
      TotalCounts <- ddply(temp2, tax.aggregate, summarise, Mean = mean(Abundance))  
      TotalCounts <- TotalCounts[with(TotalCounts, order(-Mean)),]      
      TotalCounts$Cumsum <- cumsum(TotalCounts$Mean)
      TotalCounts$x <- 1:nrow(TotalCounts)
      
      p <- ggplot(data = TotalCounts, aes(x = x, y = Cumsum)) +
        geom_line(size = 2) +
        ylim(0,100) +
        xlab("Rank abundance") +
        ylab("Cummulative read abundance (%)")  
        if (plot.log ==T){
          p <- p + scale_x_log10() 
        }
    } 
    
    abund4 <- TotalCounts
    
  }
  
  outlist <- append(outlist, list(dataframe = abund4, plot = p))

  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }
}
