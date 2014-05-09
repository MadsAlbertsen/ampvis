#' Generates a ggplot2 style core community plots
#'
#' A nice long description
#'
#' @usage amp_core(data)
#'
#' @param data (required) A phyloseq object.
#' @param group Group the data based on a sample variable (default: "Sample").
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (defualt: OTU)
#' @param tax.empty Either "remove" OTUs without taxonomic information or "rename" with OTU ID (default: rename).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 20000)
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

amp_core <- function(data, group = "Sample", scale.seq = 20000, tax.clean = T, plot.type = "Frequency", plot.log = F, output = "plot",  tax.aggregate = "OTU", tax.empty = "rename", weight = T, abund.treshold = 0.1){
  
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
    if(tax.empty == "rename" & tax.aggregate != "OTU"){  
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
  colnames(temp1)[colnames(temp1) == group] <- "var2"
  DT <- data.table(temp1)
  DT2 <- DT[, lapply(.SD, sum, na.rm=TRUE), by=list(var1, var2), .SDcols=c("Abundance") ]
  temp2 <- data.frame(DT2)
  colnames(temp2)[colnames(temp2) == "var2"] <- group
  colnames(temp2)[colnames(temp2) == "var1"] <- tax.aggregate
  
  temp2$freq <- 1
  temp2 <- subset(temp2, Abundance > 0)
  
  ## Make a nice frequency plot
  
  if(plot.type == "frequency"){
    temp3 <- ddply(temp2, c(tax.aggregate), summarise, Frequency = sum(freq), Mean = mean(Abundance)*sum(freq))
    if(weight == T){
      p <- ggplot(data = temp3, aes(x = Frequency, weight = Mean / sum(Mean)*100)) +
        geom_bar(binwidth = 1) +        
        ylab("Read abundance (%)") +
        xlab(paste("Frequency (Observed in N", group, "s)"))
    }
    
    if(weight == F){
      p <- ggplot(data = temp3, aes(x = Frequency)) +
        geom_bar(binwidth = 1) +
        scale_x_discrete(breaks = 1:max(temp3$Frequency)) +
        ylab("Count") +
        xlab(paste("Frequency (Observed in N ", group, "s)", sep=""))
    } 
  }
  
  if (plot.type == "core") {
    temp2$Abundance <- temp2$Abundance/scale.seq * 100
    temp2$HA <- ifelse(temp2$Abundance > abund.treshold, 1, 0)
    temp3 <- ddply(temp2, c(tax.aggregate), summarise, Frequency = sum(freq), freq_HA = sum(HA),  Mean = mean(Abundance)*sum(freq))
    
    p <- ggplot(data = temp3, aes(x = Frequency, y = freq_HA)) +
      geom_jitter(size = 3, alpha = 0.5) +
      ylab(paste("Highly abundant in N ", group, "s (>", abund.treshold, "%)" , sep="")) +
      xlab(paste("Observed in N ", group, "s", sep=""))
  }
  
  outlist <- append(outlist, list(data = temp3, plot = p))
  
  if(output == "complete"){ return(outlist) }
  if(output == "plot"){ return(p) }
}
