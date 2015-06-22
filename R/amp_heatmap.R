#' Generate a heatmap from amplicon data
#'
#' Generate a heatmap in ggplot2 format from amplicon data in phyloseq format. Use sample metadata to aggregate sampes and taxonomy to aggregate OTUs.
#'
#' @usage amp_headtmap(data)
#'
#' @param data (required) A phyloseq object including sample data (or a list).
#' @param group A variable from the associated sample data to group samples by.
#' @param scale A variable from the associated sample data to scale the abundance by.
#' @param normalise A specific sample or group to normalise the counts to, or "relative".
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: Phylum)
#' @param tax.add Additional taxonomic levels to display for each entry e.g. "Phylum" (default: none) 
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: best).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param calc Calculate and display mean or max across the groups (default: "mean").
#' @param order.x A taxonomy group or vector to order the x-axis by.
#' @param order.y A sample or vector to order the y-axis by.
#' @param plot.numbers Plot the values on the heatmap (default: T)
#' @param plot.breaks A vector of breaks for the abundance legend.
#' @param plot.colorscale Either sqrt or log10 (default: "sqrt")
#' @param plot.na Wether to color missing values with the lowest color in the scale (default: F).
#' @param plot.text.size The size of the plotted text (default: 4). 
#' @param plot.theme Chose different standard layouts choose from "normal" or "clean" (default: "normal").
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 10000)
#' @param output To output a plot or the complete data inclusive dataframes (default: plot)
#' 
#' @return A ggplot2 object or a list with the ggplot2 object and associated dataframes.
#' 
#' @export
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import phyloseq
#' @import data.table
#' @import grid
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_heatmap <- function(data, group = "Sample", normalise = NULL, scale = NULL, tax.aggregate = "Phylum", tax.add = NULL, tax.show = 10, tax.class = NULL, tax.empty = "best", order.x = NULL, order.y = NULL, plot.numbers = T, plot.breaks = NULL, plot.colorscale = "sqrt", plot.na = F, scale.seq = 10000, output = "plot",plot.text.size = 4, plot.theme = "normal", calc = "mean"){
  
  data <- list(abund = as.data.frame(otu_table(data)@.Data),
               tax = data.frame(tax_table(data)@.Data, OTU = rownames(tax_table(data))),
               sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data)))))
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Extract the data into seperate objects for readability
  abund <- data[["abund"]]  
  tax <- data[["tax"]]
  sample <- data[["sample"]]
  
  ## Scale the data by a selected metadata variable
  if (!is.null(scale)){
    variable <- as.numeric(sample[,scale])
    abund <- t(t(abund)*variable)
  }
  
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
  
  abund3 <- data.table(abund3)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
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
  
  ## Take the average to group level
  
  if (calc == "mean"){
    abund6 <- data.table(abund5)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
              setkey(Display, Group) %>%
              unique() %>% 
              as.data.frame()
  }
  
  if (calc == "max"){
      abund6 <- data.table(abund5)[, Abundance:=max(sum), by=list(Display, Group)] %>%
        setkey(Display, Group) %>%
        unique() %>% 
        as.data.frame()
  }  
  
  
  ## Find the X most abundant levels
  if (calc == "mean"){
      TotalCounts <- group_by(abund6, Display) %>%
        summarise(Abundance = sum(Abundance)) %>%
        arrange(desc(Abundance))
  }
  
  if (calc == "max"){
    TotalCounts <- group_by(abund6, Display) %>%
      summarise(Abundance = max(Abundance)) %>%
      arrange(desc(Abundance))
  }
  
  
  ## Subset to X most abundant levels
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){  
      tax.show <- nrow(TotalCounts)
    }
    abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show])
  }
  
  ## Subset to a list of level names
  if (!is.numeric(tax.show)){
    if (tax.show != "all"){
      abund7 <- filter(abund6, Display %in% tax.show)    
    }
    ### Or just show all  
    if (tax.show == "all"){
      tax.show <- nrow(TotalCounts)  
      abund7 <- filter(abund6, Display %in% TotalCounts$Display[1:tax.show]) 
    }
  }
  abund7 <- as.data.frame(abund7)
  
  ## Normalise to a specific group (The Abundance of the group is set as 1)  
  
  if(!is.null(normalise)){
    if (normalise != "relative"){
      temp <- dcast(abund7, Display~Group, value.var = "Abundance")
      temp1 <- cbind.data.frame(Display = temp$Display, temp[,-1]/temp[,normalise])   
      abund7 <- melt(temp1, id.var = "Display", value.name="Abundance", variable.name="Group")
    }
  } 
  if(!is.null(normalise)){
    if (normalise == "relative"){
      temp <- dcast(abund7, Display~Group, value.var = "Abundance")
      temp1 <- cbind.data.frame(Display = temp[,1], temp[,-1]/apply(as.matrix(temp[,-1]), 1, mean))    
      abund7 <- melt(temp1, id.var = "Display" , value.name="Abundance", variable.name="Group")
    }
  }
  
  ## Order.y
  if (is.null(order.y)){
    abund7$Display <- factor(abund7$Display, levels = rev(TotalCounts$Display))
  }
  if (!is.null(order.y)){
    if (length(order.y) == 1){      
      temp1 <- filter(abund7, Group == order.y) %>%
        group_by(Display) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      
      abund7$Display <- factor(abund7$Display, levels = rev(temp1$Display))
    }
    if (length(order.y) > 1){
      abund7$Display <- factor(abund7$Display, levels = order.y)
    }
  }
  
  ## Order.x
  if (!is.null(order.x)){
    if (length(order.x) == 1){
      temp1 <- filter(abund7, Display == order.x) %>%
        group_by(Group) %>%
        summarise(Mean = mean(Abundance)) %>%
        arrange(desc(Mean))
      abund7$Group <- factor(abund7$Group, levels = as.character(temp1$Group))
        
    }    
    if (length(order.x) > 1){
      abund7$Group <- factor(abund7$Group, levels = order.x)
    }
  }
  
  ## Handle NA values
  if(plot.na == F){ plot.na <- "grey50" }else{ plot.na <-"#EF8A62" }  
  
  ## Scale to percentages if not normalised and scaled
  
  if (is.null(scale) & is.null(normalise)){
    abund7[,3] <- abund7[,3]/scale.seq*100
  }
  
  ## Make a heatmap style plot
  p <- ggplot(abund7, aes_string(x = "Group", y = "Display", label = formatC("Abundance", format = "f", digits = 1))) +     
    geom_tile(aes(fill = Abundance), colour = "white", size = 0.5) +
    theme(axis.text.x = element_text(size = 10, hjust = 1, angle = 90)) + 
    theme(axis.text.y = element_text(size = 12)) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  if (plot.numbers == T){
    abund8 <- abund7
    abund8$Abundance <- round(abund8$Abundance, 1)
    p <- p + geom_text(data = abund8, size = plot.text.size, colour = "grey10")  
  }
  if (is.null(plot.breaks)){
    p <- p +scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), trans = plot.colorscale, na.value=plot.na)
  }
  if (!is.null(plot.breaks)){
    p <- p +scale_fill_gradientn(colours = brewer.pal(3, "RdBu"), trans = plot.colorscale, breaks=plot.breaks, na.value=plot.na)
  }
  if (is.null(normalise)){
    p <- p + labs(x = "", y = "", fill = "% Read\nAbundance")  
  }
  if (!is.null(normalise)){
    p <- p + labs(x = "", y = "", fill = "Relative")  
  }
  
  if(plot.theme == "clean"){
    p <- p + theme(legend.position = "none",
                   axis.text.y = element_text(size = 8, color = "black"),
                   axis.text.x = element_text(size = 8, color = "black"),
                   axis.title = element_blank(),
                   text = element_text(size = 8, color = "black"),
                   axis.ticks.length = unit(1, "mm"),
                   plot.margin = unit(c(0,0,0,0), "mm"),
                   title = element_text(size = 8)
    )
  }
  
  ## Define the output 
  if (output == "complete"){
    outlist <- list(heatmap = p, data = abund7)
    return(outlist)  
  }
  if (output == "plot"){
    return(p)
  }
}
