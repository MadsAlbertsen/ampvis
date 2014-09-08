#' Used for cleaning and renaming taxonomy in phyloseq objects
#'
#' Used internally in other ampvis functions.
#'
#' @usage amp_rename(data)
#'
#' @param data (required) A phyloseq object.
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.empty Either "remove" OTUs without taxonomic information at X level, with "best" classification or add the "OTU" name (default: best).
#' @param tax.level The taxonomic level to remove OTUs with empty taxonomy, only used when tax.empty = "remove" (default: Genus).
#' 
#' @return A phyloseq object with cleaned and renamed taxonomy.
#' 
#' @export
#' @import phyloseq
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rename <- function(data, tax.class = NULL, tax.empty = "best", tax.level = "Genus"){

  abund <- as.data.frame(otu_table(data))
  tax <- as.data.frame(tax_table(data))
  
  ## First make sure that all entires are strings
  for ( i in 1:ncol(tax) ){
    tax[,i] <- as.character(tax[,i])  
  }
  
  ## Change a specific phylum to class level
  if(!is.null(tax.class)){
    for (i in 1:nrow(tax)){
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] == tax.class){
        tax$Phylum[i] <- tax$Class[i]   
      }
    }
  }
  
  ## Remove the underscore classifier from the data  
  tax$Kingdom <- gsub("k__", "", tax$Kingdom)
  tax$Phylum <- gsub("p__", "", tax$Phylum)
  tax$Phylum <- gsub("c__", "", tax$Phylum)
  tax$Class <- gsub("c__", "", tax$Class)
  tax$Order <- gsub("o__", "", tax$Order)
  tax$Family <- gsub("f__", "", tax$Family)
  tax$Genus <- gsub("g__", "", tax$Genus)
  
  ## Check if there is a species level otherwise add it for consistency
  if (!is.null(tax$Species)){
    tax$Species <- gsub("s__", "", tax$Species)
  } else {
    tax$Species <- ""
  }
  
  tax[is.na(tax)] <- ""
    
  ## How to handle empty taxonomic assignments
  if (tax.empty == "OTU"){
    for (i in 1:nrow(tax)) {
      if (tax[i,"Species"] == "") {tax[i,"Species"] <- rownames(tax)[i]}
      if (tax[i,"Genus"] == "") {tax[i,"Genus"] <- rownames(tax)[i]}
      if (tax[i,"Family"] == "") {tax[i,"Family"] <- rownames(tax)[i]}
      if (tax[i,"Order"] == "") {tax[i,"Order"] <- rownames(tax)[i]}
      if (tax[i,"Class"] == "") {tax[i,"Class"] <- rownames(tax)[i]}
      if (tax[i,"Phylum"] == "") {tax[i,"Phylum"] <- rownames(tax)[i]}
    }
  }
  
  if(tax.empty == "best"){
    tax[tax$Kingdom == "","Kingdom"] <- "Unclassified"
    for (i in 1:nrow(tax)) {   
      if (tax[i,"Species"] == "") {
        if (tax[i,"Genus"] != "") { rn <- paste("g__", tax[i,"Genus"], "_", rownames(tax)[i], sep = "") } else{
          if (tax[i,"Family"] != "") { rn <- paste("f__", tax[i,"Family"], "_", rownames(tax)[i], sep = "") } else{
            if (tax[i,"Order"] != "") { rn <- paste("o__", tax[i,"Order"], "_", rownames(tax)[i], sep = "") } else{
              if (tax[i,"Class"] != "") { rn <- paste("c__", tax[i,"Class"], "_", rownames(tax)[i], sep = "") } else{
                if (tax[i,"Phylum"] != "") { rn <- paste("p__", tax[i,"Phylum"], "_", rownames(tax)[i], sep = "") } else{
                  if (tax[i,"Kingdom"] != "") { rn <- paste("k__", tax[i,"Kingdom"], "_", rownames(tax)[i], sep = "") } 
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
    tax <- subset(tax, tax[,tax.level] != "")
    abund <- subset(abund, rownames(abund) %in% rownames(tax))
  }
  
  data@tax_table <- tax_table(as.matrix(tax))
  data@otu_table <- otu_table(as.matrix(abund), taxa_are_rows = T)
  return(data)
}
