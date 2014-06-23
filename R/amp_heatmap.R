#' Generates a ggplot2 style heatmap from amplicon data in phyloseq format 
#'
#' A nice long description
#'
#' @usage amp_headtmap(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param group A variable from the associated sample data to group samples by.
#' @param scale A variable from the associated sample data to scale the abundance by.
#' @param normalise A specific sample or group to normalise the counts to, or "relative".
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: Phylum)
#' @param tax.add Additional taxonomic levels to display for each entry e.g. "Phylum" (default: none) 
#' @param tax.show The number of taxa to show or a vector of taxa names (default: 10).
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.empty Either "remove" OTUs without taxonomic information, "rename" with best classification or add the "OTU" name (default: rename).
#' @param order.x A taxonomy group or vector to order the x-axis by.
#' @param order.y A sample or vector to order the y-axis by.
#' @param plot.numbers Plot the values on the heatmap (default: F)
#' @param plot.breaks A vector of breaks for the abundance legend.
#' @param plot.colorscale Either sqrt or log (default: "sqrt")
#' @param plot.na Wether to color missing values with the lowest color in the scale (default: F).
#' @param plot.text.size The size of the plotted text (default: 4). 
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 10000)
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

amp_heatmap <- function(data, group = NULL, normalise = NULL, scale = NULL, tax.aggregate = "Phylum", tax.add = NULL, tax.show = 10, tax.clean = T, tax.empty = "rename", order.x = NULL, order.y = NULL, plot.numbers = T, plot.breaks = NULL, plot.colorscale = "sqrt", plot.na = F, scale.seq = 10000, output = "plot", tax.clean.proteobacteria = T,plot.text.size = 4){
  
  ## Extract all data from the phyloseq object
  abund<-as.data.frame(otu_table(data))
  tax <- as.data.frame(tax_table(data))
  tax <- data.frame(tax, OTU = rownames(tax))
  sample <- suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  if (is.null(tax$Species)){tax$Species <- ""}
  
  if(plot.na == F){ plot.na <- "grey50" }else{ plot.na <-"#EF8A62" }
  
  ## Extract group information
  if (!is.null(group)){
    grp <- cbind.data.frame(Sample = rownames(sample),sample[,group])
    colnames(grp)[2] <- group
  }
  if (is.null(group)){
    grp <- cbind.data.frame(Sample = rownames(sample), Samples = rownames(sample))
  }
  
  ## Scale the data by a select variable
  if (!is.null(scale)){
    variable <- as.numeric(sample[,scale])
    abund <- t(t(abund)*variable)
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
  
  tax$Kingdom <- gsub("k__", "", tax$Kingdom)
  tax$Phylum <- gsub("p__", "", tax$Phylum)
  tax$Phylum <- gsub("c__", "", tax$Phylum)
  tax$Class <- gsub("c__", "", tax$Class)
  tax$Order <- gsub("o__", "", tax$Order)
  tax$Family <- gsub("f__", "", tax$Family)
  tax$Genus <- gsub("g__", "", tax$Genus)
  tax[is.na(tax)] <- ""
  if (!is.null(tax$Species)){tax$Species <- gsub("s__", "", tax$Species)} 
  
  t <- tax
  
  ## How to handle empty taxonomic assignments
  if (tax.empty == "OTU"){
    for (i in 1:nrow(tax)) {
      if (tax[i,"Species"] == "") {tax[i,"Species"] <- tax[i,"OTU"]}
      if (tax[i,"Genus"] == "") {tax[i,"Genus"] <- tax[i,"OTU"]}
      if (tax[i,"Family"] == "") {tax[i,"Family"] <- tax[i,"OTU"]}
      if (tax[i,"Order"] == "") {tax[i,"Order"] <- tax[i,"OTU"]}
      if (tax[i,"Class"] == "") {tax[i,"Class"] <- tax[i,"OTU"]}
      if (tax[i,"Phylum"] == "") {tax[i,"Phylum"] <- tax[i,"OTU"]}
    }
  }
  
  if(tax.empty == "rename"){
    tax[tax$Kingdom == "","Kingdom"] <- "Unclassified"
    for (i in 1:nrow(tax)) {   
      if (tax[i,"Species"] == "") {
        if (tax[i,"Genus"] != "") { rn <- paste("g__", tax[i,"Genus"], "_", tax[i,"OTU"], sep = "") } else{
          if (tax[i,"Family"] != "") { rn <- paste("f__", tax[i,"Family"], "_", tax[i,"OTU"], sep = "") } else{
            if (tax[i,"Order"] != "") { rn <- paste("o__", tax[i,"Order"], "_", tax[i,"OTU"], sep = "") } else{
              if (tax[i,"Class"] != "") { rn <- paste("c__", tax[i,"Class"], "_", tax[i,"OTU"], sep = "") } else{
                if (tax[i,"Phylum"] != "") { rn <- paste("p__", tax[i,"Phylum"], "_", tax[i,"OTU"], sep = "") } else{
                  if (tax[i,"Kingdom"] != "") { rn <- paste("k__", tax[i,"Kingdom"], "_", tax[i,"OTU"], sep = "") } 
                }
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
  
  #Make a name variable that can be used instead of tax.aggregate to display multiple levels 
  dname <- tax.aggregate
  a <- data.frame(tax[,tax.aggregate])
  if (!is.null(tax.add)){
    dname <- c(tax.aggregate, tax.add)
    a <- data.frame(apply(tax[,dname], 1, paste, collapse="; "))
  }
  tax <- cbind(tax, a)   
  colnames(tax)[ncol(tax)] <- "Display"
  tax$Display <- as.character(tax$Display)
  
  ## Merge the taxonomic and abundance information
  abund2 <- cbind.data.frame(tax, abund)
  
  ## Aggregate to a specific taxonomic level
  abund3 <- melt(abund2, id.var = "Display", measure.vars=rownames(sample))
  colnames(abund3)[2:3] <- c("Sample", "Abundance")
  
  DT <- data.table(abund3)
  DT2 <- DT[, lapply(.SD, sum, na.rm=TRUE), by=list(Display, Sample), .SDcols=c("Abundance") ]   
  abund4 <- data.frame(DT2)
  
  ## Add group information
  abund5 <- join(abund4, grp, by="Sample")
  
  ## Take the average to group level
  colnames(abund5)[colnames(abund5) == colnames(grp[2])] <- "var2"
  DT3 <- data.table(abund5)
  DT4 <- DT3[, lapply(.SD, mean, na.rm=TRUE), by=list(Display, var2), .SDcols=c("Abundance") ]   
  abund6 <- data.frame(DT4)
  colnames(abund6)[colnames(abund6) == "var2"] <- colnames(grp[2])
  
  ## Find the X most abundant levels
  TotalCounts <- ddply(abund6, "Display", summarise, Abundance = sum(Abundance))
  TotalCounts <- TotalCounts[with(TotalCounts, order(-Abundance)),]
  
  ## Subset to X most abundant levels
  if (is.numeric(tax.show)){
    if (tax.show > nrow(TotalCounts)){
      tax.show <- nrow(TotalCounts)
    }
    abund7 <- subset(abund6, abund6[,1] %in% TotalCounts[1:tax.show,"Display"])  
  }
  
  ## Subset to a list of level names
  if (!is.numeric(tax.show)){
    if (tax.show != "all"){
      abund7 <- subset(abund6, abund6[,1] %in% tax.show)    
    }
    ### Or just show all  
    if (tax.show == "all"){
      tax.show <- nrow(TotalCounts)  
      abund7 <- subset(abund6, abund6[,1] %in% TotalCounts[1:tax.show,"Display"])  
    }
  }
  
  ## Normalise to a specific group (The Abundance of the group is set as 1)
  
  if(!is.null(normalise)){
    if (normalise != "relative"){
      colnames(abund7) <- c("var1", "var2", "Abundance")
      temp <- dcast(abund7, var1~var2, value.var = "Abundance")
      colnames(temp)[1] <- "Display"
      temp2 <- temp[,-1]  
      temp3 <- temp2/temp2[,normalise]
      temp4 <- cbind.data.frame(temp[,1], temp3)       
      colnames(temp4)[1] <- "Display" 
      temp5 <- melt(temp4, id.var = "Display")
      colnames(temp5) <- c("Display", colnames(grp[2]), "Abundance")
      abund7 <- temp5
    }
  }
  
  if(!is.null(normalise)){
    if (normalise == "relative"){
      colnames(abund7) <- c("var1", "var2", "Abundance")
      temp <- dcast(abund7, var1~var2, value.var = "Abundance")
      colnames(temp)[1] <- "Display"
      temp2 <- temp[,-1]  
      rel <- apply(as.matrix(temp2), 1, mean)
      temp3 <- temp2/rel
      temp4 <- cbind.data.frame(temp[,1], temp3)    
      colnames(temp4)[1] <- "Display" 
      temp5 <- melt(temp4, id.var = "Display")
      colnames(temp5) <- c("Display", colnames(grp[2]), "Abundance")
      abund7 <- temp5    
    }
  }
  
  ## Order.y
  
  if (is.null(order.y)){
    abund7[,1] <- factor(abund7[,1], levels = rev(TotalCounts[,1]))
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
    abund7[,3] <- abund7[,3]/scale.seq*100
  }
  
  ## Make a heatmap style plot
  
  p <- ggplot(abund7, aes_string(x = colnames(grp[2]), y = "Display", label = formatC("Abundance", format = "f", digits = 1))) +     
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
    p <- p + labs(x = "", y = "", fill = "Abundance")  
  }
  if (!is.null(normalise)){
    p <- p + labs(x = "", y = "", fill = "Relative")  
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