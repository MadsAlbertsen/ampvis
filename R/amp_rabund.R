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
#' @param tax.empty Either "remove" OTUs without taxonomic information, "rename" with best classification or add the "OTU" name (default: rename).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 10000)
#' @param plot.type Either point, boxplot or curve (default: boxplot).
#' @param plot.flip Flip the axis of the plot (default: F).
#' @param point.size Size of points (default: 2).
#' @param output Either plot or complete (default: "plot").
#' @param sort.by Sort the boxplot by either Median, Mean or Total (default = "Median")
#' 
#' @return A ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import phyloseq
#' @import grid
#' @import data.table
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rabund <- function(data, group = "Sample", order.group = NULL, tax.show = 50, scale.seq = 10000, tax.clean = T, plot.type = "boxplot", plot.log = F, output = "plot", tax.add = NULL, tax.aggregate = "Genus", tax.empty = "best", tax.class = NULL, point.size = 2, plot.flip = F, sort.by = "Median"){
  
  ## Check the input data type and convert to list if it's a phyloseq object
  data <- list(abund = as.data.frame(otu_table(data)@.Data),
               tax = data.frame(tax_table(data)@.Data, OTU = rownames(tax_table(data))),
               sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Extract the data into seperate objects for readability
  abund <- data[["abund"]]  
  tax <- data[["tax"]]
  sample <- data[["sample"]]
  
  ## Make a name variable that can be used instead of tax.aggregate to display multiple levels 
  suppressWarnings(
    if (!is.null(tax.add)){
      if (tax.add != tax.aggregate) {
        tax <- data.frame(tax, Display = apply(tax[,c(tax.add,tax.aggregate)], 1, paste, collapse="; "))
      }
    } else {
      tax <- data.frame(tax, Display = tax[,tax.aggregate])
    }
  )  
  
  # Aggregate to a specific taxonomic level
  abund3 <- cbind.data.frame(Display = tax[,"Display"], abund) %>%
    melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
  
  abund3 <- data.table(abund3)[, Abundance:=sum(Abundance), by=list(Display, Sample)] %>%
    setkey(Display, Sample) %>%
    unique() %>% 
    as.data.frame()
  
  ## Add group information
  suppressWarnings(
    if (group != "Sample"){
      if (length(group) > 1){
        grp <- data.frame(Sample = rownames(sample), Group = apply(sample[,group], 1, paste, collapse = " ")) 
      } else{
        grp <- data.frame(Sample = rownames(sample), Group = sample[,group]) 
      }
      abund3$Group <- grp$Group[match(abund3$Sample, grp$Sample)]
      abund5 <- abund3
    } else{ abund5 <- data.frame(abund3, Group = abund3$Sample)}
  ) 
  
  if (plot.type != "curve"){
    ## Find the X most abundant levels and sort
    TotalCounts <- group_by(abund5, Display) %>%
      summarise(Median = median(Abundance), Total = sum(Abundance), Mean = mean(Abundance))
    if(sort.by == "Median"){TotalCounts %<>% arrange(desc(Median)) %>% as.data.frame()}
    if(sort.by == "Mean"){TotalCounts %<>% arrange(desc(Mean)) %>% as.data.frame()}
    if(sort.by == "Total"){TotalCounts %<>% arrange(desc(Total)) %>% as.data.frame()}
    
    abund5$Display <- factor(abund5$Display, levels = rev(TotalCounts$Display))
    
    ## Make sure we only show a possible number of taxa
    if (!is.numeric(tax.show)){
      tax.show <- nrow(TotalCounts)
    }
    
    ## Subset to X most abundant levels
    if (is.numeric(tax.show)){
      if (tax.show > nrow(TotalCounts)){  
        tax.show <- nrow(TotalCounts)
      }
      abund7 <- subset(abund5, abund5$Display %in% TotalCounts[1:tax.show,"Display"])  
    }
    
    ## Subset to a list of level names
    if (!is.numeric(tax.show)){
      if (tax.show != "all"){
        abund7 <- subset(abund5, abund5$Display %in% tax.show)    
      }
      ### Or just show all  
      if (tax.show == "all"){
        tax.show <- nrow(TotalCounts)  
        abund7 <- subset(abund5, abund5$Display %in% TotalCounts[1:tax.show,"Display"])  
      }
    }
    
    ## Scale to a specific abundance
    abund7$Abundance <- abund7$Abundance/scale.seq*100
    
    ## plot the data
    
    if (group == "Sample"){
      p <-ggplot(abund7, aes(x = Display, y = Abundance))   
    }
    if (group != "Sample"){
      if(!is.null(order.group)){
        abund7$Group <- factor(abund7$Group, levels = rev(order.group))
      }
      p <-ggplot(abund7, aes(x = Display, y = Abundance, color = Group))   
    }
    
    p <- p +  ylab("Read Abundance (%)") + guides(col = guide_legend(reverse = TRUE))
    
    if (plot.flip == F){ p <- p + coord_flip() } else{
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    
    if (plot.type == "point"){ p <- p + geom_point(size = point.size) }
    if (plot.type == "boxplot"){p <- p + geom_boxplot(outlier.size = point.size)}
    if (plot.log ==T){ p <- p + scale_y_log10()}
    
    outlist <- list(plot = p, data = abund7)
  }
  
  ## If type = curve then generate a second dataframe
  
  if (plot.type == "curve"){
    temp3 <- group_by(abund5, Display, Group) %>%
      summarise(Mean = mean(Abundance))
    
    TotalCounts <- temp3[with(temp3, order(-Mean)),] %>%
      group_by(Group) %>%
      mutate(dummy = 1) %>%
      mutate(Cumsum = cumsum(Mean), Rank = cumsum(dummy)) %>%
      as.data.frame()
    
    if(!is.null(order.group)){
      TotalCounts$Group <- factor(TotalCounts$Group, levels = rev(order.group))
    }
    TotalCounts$Cumsum <- TotalCounts$Cumsum/scale.seq * 100
    
    p <- ggplot(data = TotalCounts, aes(x = Rank, y = Cumsum, color = Group)) +
      geom_line(size = 2) +
      ylim(0,100) +
      xlab("Rank abundance") +
      ylab("Cummulative read abundance (%)")  
    if (plot.log ==T){
      p <- p + scale_x_log10() 
    } 
    
    outlist <- list(plot = p, data = TotalCounts)
  }
  
  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }
}
