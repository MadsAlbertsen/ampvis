#' Tests if there is a significant difference in abundance between selected conditions
#'
#' A nice long description
#'
#' @usage amp_test(data, group)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param group (required) A variable from the associated sample data to group samples by.
#' @param tax.aggregate Group data at specific taxonomic level (defaul: "Genus").
#' @param tax.clean Replace the phylum Proteobacteria with the respective Classes instead (default: T).
#' @param tax.empty Either "remove" OTUs without taxonomic information or "rename" with OTU ID (default: rename).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 20000)
#' @param sig.pval Significance cutoff for plotting the data (default: 0.01)
#' 
#' @return A p-value for each comparison.
#' 
#' @export
#' @import ggplot2
#' @import plyr
#' @import reshape2
#' @import phyloseq
#' @import grid
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_test <- function(data, group, tax.aggregate = "Genus", tax.clean = T, tax.empty = "rename", scale.seq = 20000, sig.pval = 0.01){
  
  ## Extract all data from the phyloseq object 
  abund <- as.data.frame(otu_table(data))
  tax <- as.data.frame(tax_table(data))
  tax <- cbind.data.frame(tax, rownames(tax))
  colnames(tax)[ncol(tax)] <- "OTU"
  sample <- suppressWarnings(data.frame(sample_data(data)))
  
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
  
  ## Combine it to a single data frame
  q <- cbind.data.frame(tax, abund)
  q2 <- melt(q, id.vars=colnames(tax))
  colnames(q2)[(ncol(q2)-1):ncol(q2)] <- c("SampleID","Abundance")
  colnames(sample)[1] <- "SampleID"
  
  ## Aggregate to a specific level
  q3 <- ddply(q2, c(tax.aggregate,"SampleID"), summarise, Abundance = sum(Abundance))
  
  ## Merge with sample data
  q4 <- merge(x = sample, y = q3)
  
  ## Test significance
  t4 <- droplevels(q4)
  colnames(t4)[colnames(t4) == group] <- "group"
  
  ### If the group is factor data
  if (is.factor(sample[,group])){
    
    ### If the group contains 2 levels then use a simple t-test
    if (length(levels(sample[,group]))==2){
      res.ttest <- ddply(t4, tax.aggregate, summarise, pval = t.test(Abundance~group)$p.value)    
      sig <- subset(res.ttest, pval <= sig.pval)
      q5 <- subset(q4, q4[,tax.aggregate] %in% sig[, 1])
      q5$Abundance <- q5$Abundance/scale.seq*100
      p <- ggplot(q5, aes_string(x = tax.aggregate, y = "Abundance", color = group)) +
        geom_boxplot() +
        coord_flip() +
        ylab("Abundance (%)")
    }
    
    ### If the group contains more that 2 levels then use something like anova
    if (length(levels(sample[,group]))>2){
      
      #res.tukey <- 
      #### Need to extract the pvalue of each comparison and then convert it into a table somehow
      
      #fit = lm(formula = t4$Abundance ~ qt[,group])
      #fit.anova <- anova(fit)
      #fit.aov <- aov(fit)
      #fit.bartlett <- bartlett.test(qt$Abundance, qt[,group])
      #fit.tukeyHSD <- TukeyHSD(fit.aov)
    }
    
  }  
  
  outlist <- list(result.ttest = res.ttest, plot = p)
  return(outlist)
}
