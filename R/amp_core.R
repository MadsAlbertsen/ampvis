#' Generates a ggplot2 style core community plots
#'
#' A nice long description
#'
#' @usage amp_core(data)
#'
#' @param data (required) A phyloseq object.
#' @param group Group the data based on a sample variable (default: "Sample").
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: best).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: OTU).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 10000)
#' @param plot.type Either core or frequency (default: frequency).
#' @param weight Weight the frequency by abundance (default: T).
#' @param abund.treshold Treshold for considering something abundant in percent (default: 0.1).
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

amp_core <- function(data, group = "Sample", scale.seq = 10000, tax.class = NULL, tax.empty = "best", plot.type = "frequency", output = "plot",  tax.aggregate = "OTU", weight = T, abund.treshold = 0.1){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Aggregate to a specific taxonomic level
  if (tax.aggregate != "OTU"){ data <- tax_glom(data, taxrank=tax.aggregate) }
    
  ## Extract all data from the phyloseq object
  abund<-as.data.frame(otu_table(data))
  tax <- data.frame(tax_table(data), OTU = rownames(tax_table(data)))
  sample <- suppressWarnings(data.frame(sample_data(data)))
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
    
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
  colnames(temp1)[colnames(temp1) == group] <- "var2"
  DT <- data.table(temp1)
  DT2 <- DT[, lapply(.SD, mean, na.rm=TRUE), by=list(var1, var2), .SDcols=c("Abundance") ]
  temp2 <- data.frame(DT2)
  colnames(temp2)[colnames(temp2) == "var2"] <- group
  colnames(temp2)[colnames(temp2) == "var1"] <- tax.aggregate
  
  temp2$freq <- 1
  temp2 <- subset(temp2, Abundance > 0)
  
  ## Make a nice frequency plot
  
  if(plot.type == "frequency"){
    temp3 <- ddply(temp2, c(tax.aggregate), summarise, Frequency = sum(freq), Mean = mean(Abundance))
      
    if(weight == T){
      p <- ggplot(data = temp3, aes(x = Frequency, weight = Mean / sum(Mean)*100)) +
        ylab("Read abundance (%)") +
        xlab(paste("Frequency (Observed in N ", group, "s)", sep = ""))
        if(max(temp3$Frequency) > 30){ p <- p + geom_bar()}
        if(max(temp3$Frequency) <= 30){ p <- p + geom_bar(binwidth = 1)}
    }
    
    if(weight == F){
      p <- ggplot(data = temp3, aes(x = Frequency)) +
        ylab("Count") +
        xlab(paste("Frequency (Observed in N ", group, "s)", sep=""))
        if(max(temp3$Frequency) > 30){ p <- p + geom_bar()}
        if(max(temp3$Frequency) <= 30){ p <- p + geom_bar(binwidth = 1)}
    } 
  }
  
  if (plot.type == "core") {
    temp2$Abundance <- temp2$Abundance/scale.seq * 100
    temp2$HA <- ifelse(temp2$Abundance > abund.treshold, 1, 0)
    temp3 <- ddply(temp2, tax.aggregate, summarise, Frequency = sum(freq), freq_HA = sum(HA),  Mean = mean(Abundance))
    
    p <- ggplot(data = temp3, aes(x = Frequency, y = freq_HA)) +
      geom_jitter(size = 3, alpha = 0.5) +
      ylab(paste("Highly abundant in N ", group, "s (>", abund.treshold, "%)" , sep="")) +
      xlab(paste("Observed in N ", group, "s", sep=""))
  }
  
  outlist <- append(outlist, list(data = temp3, plot = p))
  
  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }
}
