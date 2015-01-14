#' Generates a ggplot2 style core community plots
#'
#' A nice long description
#'
#' @usage amp_core(data)
#'
#' @param data (required) A phyloseq object (or a list).
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
#' @import dplyr
#' @import reshape2
#' @import phyloseq
#' @import grid
#' @import data.table
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_core <- function(data, group = "Sample", scale.seq = 10000, tax.class = NULL, tax.empty = "best", plot.type = "frequency", output = "plot",  tax.aggregate = "OTU", weight = T, abund.treshold = 0.1){
  
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
  
  # Aggregate to a specific taxonomic level
  abund1 <- cbind.data.frame(Display = tax[,tax.aggregate], abund) %>%
    melt(id.var = "Display", value.name= "Abundance", variable.name = "Sample")
  abund1 <- data.table(abund1)[, sum:=sum(Abundance), by=list(Display, Sample)] %>%
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
      abund1$Group <- grp$Group[match(abund1$Sample, grp$Sample)]
      abund2 <- abund1
    } else{ abund2 <- data.frame(abund1, Group = abund1$Sample)}
  )
  
  ## Take the average to group level
  abund3 <- data.table(abund2)[, Abundance:=mean(sum), by=list(Display, Group)] %>%
    setkey(Display, Group) %>%
    unique() %>%
    filter(Abundance > 0) %>%  
    mutate(freq = 1)
  
  ## Make a nice frequency plot
  
  if(plot.type == "frequency"){
    temp3 <- group_by(abund3, Display) %>%
      summarise(Frequency = sum(freq), Total = sum(sum)) %>%
      mutate(Percent = Total / (length(unique(abund3$Group)) * scale.seq) * 100) %>%
      as.data.frame()
    
    if(weight == T){
      p <- ggplot(data = temp3, aes(x = Frequency, weight = Percent)) +
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
    abund3$Abundance <- abund3$Abundance/scale.seq * 100
    abund3$HA <- ifelse(abund3$Abundance > abund.treshold, 1, 0)
    temp3 <- group_by(abund3, Display) %>%
      summarise(Frequency = sum(freq), freq_A= sum(HA), Abundance = round(mean(Abundance),2)) %>%
      as.data.frame()
    
    p <- ggplot(data = temp3, aes(x = Frequency, y = freq_A)) +
      geom_jitter(size = 3, alpha = 0.5) +
      ylab(paste("Abundant in N ", group, "s (>", abund.treshold, "%)" , sep="")) +
      xlab(paste("Observed in N ", group, "s", sep=""))
  }
  
  if(tax.aggregate == "OTU"){
    colnames(temp3)[1] <- "OTU"
    core <- merge(x = temp3, y = tax, by = "OTU") 
    temp3 <- core
  }
  
  if(output == "complete"){ return(list(data = temp3, plot = p, abund = abund, tax = tax, sample = sample)) }
  if(output == "plot"){ return(p) }
}
