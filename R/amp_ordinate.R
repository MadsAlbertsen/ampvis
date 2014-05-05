#' Generates a ggplot2 style ordinate plot of amplicon data.
#'
#' A nice long description
#'
#' @usage amp_ordinate(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param trans Transform the raw counts, currently supports "none" or "sqrt" (default: "sqrt").
#' @param ordinate.type Either PCA or NMDS (default: "PCA").
#' @param ncomp The number of principal components to extract if using PCA (default: 5)
#' @param plot.x Variable to plot on the x-axis if using PCA (default: "PC1")
#' @param plot.y Variable to plot on the y-axis if using PCA (default: "PC2")
#' @param plot.color Color the points by a sample variable.
#' @param plot.point.size Size of the plotted sample points (default: 3)
#' @param plot.species Plot loadings as points (default: F)
#' @param plot.nspecies Plot the n most extreme species with their genus classification (default: 0).
#' @param plot.label Label points using a sample variable.
#' @param plot.group Uses plot.color and groups samples by either "chull" or "centroid".
#' @param plot.group.label Add label to the groups (default: F).
#' @param envfit.factor A vector of factor variables from the sample data used for envfit to the model.
#' @param envfit.numeric A vector of numerical variables from the sample data used for envfit to the model.
#' @param envfit.significant The significance treshold for displaying envfit parameters (default: 0.01).
#' @param envfit.resize Scale the size of the numeric arrows (default: 1).
#' @param tax.clean Add best assignment to Genus level classification if none exists.
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
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_ordinate <- function(data, trans = "sqrt", ordinate.type = "PCA", ncomp = 5, plot.x = "PC1", plot.y = "PC2", plot.color = NULL, plot.point.size = 3, plot.species = F, plot.nspecies = NULL, plot.label = NULL, plot.group = NULL, plot.group.label = F, envfit.factor = NULL, envfit.numeric = NULL, envfit.significant = 0.001, envfit.resize = 1, tax.clean =T, output = "plot"){
  
  ## Load the data and extract relevant dataframes from the phyloseq object
  
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(as.data.frame(sample_data(data)))
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
  
  ## Transform the data
  
  abund1 <- abund
  if (trans == "sqrt"){
    abund1 <- sqrt(abund)
  }
  
  ## Clean taxonomic names by assigning the best classification to genus
  
  if (tax.clean == T){
    t <- tax
    t$Kingdom <- as.character(t$Kingdom)
    t$Phylum <- as.character(t$Phylum)
    t$Class <- as.character(t$Class)
    t$Order <- as.character(t$Order)
    t$Family <- as.character(t$Family)
    t$Genus <- as.character(t$Genus)
    t[is.na(t$Phylum),"Phylum"] <- "p__Unclassified"
    a <- which( colnames(t)=="Class")
    for (i in 1:nrow(t)){
      if (t[i,a-1] == "p__"){t[i,a-1] <- "p__Unclassified"}
      if (is.na(t[i,a])){t[i,a] <- t[i, a-1]}
      if (t[i,a] == "c__"){t[i,a] <- t[i, a-1]}
      if (is.na(t[i,a+1])){t[i,a+1] <- t[i, a]}
      if (t[i,a+1] == "o__"){t[i,a+1] <- t[i, a]}
      if (is.na(t[i,a+2])){t[i,a+2] <- t[i, a+1]}
      if (t[i,a+2] == "f__"){t[i,a+2] <- t[i, a+1]}
      if (is.na(t[i,a+3])){t[i,a+3] <- t[i, a+2]}  
      if (t[i,a+3] == "g__"){t[i,a+3] <- t[i, a+2]}
    }  
    tax <-t  
  }
  
  ## Calculate NMDS
  
  if (ordinate.type == "NMDS"){
      plot.x = "MDS1"
      plot.y = "MDS2"
      
      model <- metaMDS(t(abund1))
      combined <- cbind.data.frame(model$points, sample) 
    
      loadings <- cbind.data.frame(model$species, tax)
      loadings$MDS1 <- loadings$MDS1*attributes(model$species)$shrinkage[1]
      loadings$MDS2 <- loadings$MDS2*attributes(model$species)$shrinkage[2]
      OTU <- gsub("denovo", "", rownames(loadings))
      loadings <- cbind.data.frame(loadings, OTU)
      
      species <- cbind(loadings, loadings[,plot.x]^2 + loadings[,plot.y]^2)
      colnames(species)[ncol(species)] <- "dist"
      species <- species[with(species, order(-dist)), ]
      
      outlist <- append(outlist, list(nmds.model = model, nmds.scores = combined, nmds.loadings = species))
  }
  
  ## Calculate PCA
  
  if (ordinate.type == "PCA"){
    
    model <- rda(t(abund1))
    exp <- round(model$CA$eig/model$CA$tot.chi*100,1)
    
    scores <- scores(model, choices = 1:ncomp)$sites
    combined <- cbind.data.frame(sample, scores)
      
    loadings <- cbind.data.frame(scores(model, choices = 1:ncomp)$species, tax)
    OTU <- gsub("denovo", "", rownames(loadings))
    loadings <- cbind.data.frame(loadings, OTU)
    
    species <- cbind(loadings, loadings[,plot.x]^2 + loadings[,plot.y]^2)
    colnames(species)[ncol(species)] <- "dist"
    species <- species[with(species, order(-dist)), ]
    
    outlist <- append(outlist, list(pca.model = model, pca.scores = combined, loadings = species))
  }
    
    ## Fit environmental factors using vegans envfit function
  
    if(!is.null(envfit.factor)){   
      s3  <- suppressWarnings(as.data.frame(as.matrix(sample[, envfit.factor])))
      if (ordinate.type == "PCA"){ef.f <- envfit(model, s3, permutations = 999, choices=c(plot.x, plot.y))}
      if (ordinate.type == "NMDS"){ef.f <- envfit(model, s3, permutations = 999)}
      temp <- cbind.data.frame(rownames(ef.f$factors$centroids),ef.f$factors$var.id, ef.f$factors$centroids)
      colnames(temp)[1:2] <- c("Name","Variable")
      temp1 <- cbind.data.frame(names(ef.f$factors$pvals), ef.f$factors$pvals)
      colnames(temp1) <- c("Variable", "pval")
      temp2 <- merge(temp, temp1)
      f.sig <- subset(temp2, pval <= envfit.significant)   
      if (ordinate.type == "NMDS"){colnames(f.sig)[3:4] <- c("MDS1", "MDS2")}
      outlist <- append(outlist, list(eff.model = ef.f))
    }
    
     ## Fit environmental numeric data using vegans envfit function
  
    if (!is.null(envfit.numeric)){
      if (ordinate.type == "PCA") {ef.n <- envfit(model, sample[, envfit.numeric], permutations = 999, choices=c(plot.x, plot.y))}
      if (ordinate.type == "NMDS") {ef.n <- envfit(model, sample[, envfit.numeric], permutations = 999)}
      temp <- cbind.data.frame(rownames(ef.n$vectors$arrows), ef.n$vectors$arrows*sqrt(ef.n$vectors$r), ef.n$vectors$pvals)
      colnames(temp)[c(1,4)] <- c("Name","pval")      
      n.sig <- subset(temp, pval <= envfit.significant)          
      if (ordinate.type == "NMDS"){colnames(n.sig)[2:3] <- c("MDS1", "MDS2")}
      outlist <- append(outlist, list(efn.model = ef.n))
    }
  
    ## Plot  
  
    ### Plot: Basic plot either with or without colors
    if(is.null(plot.color)){
      p <- ggplot(combined, aes_string(x = plot.x, y = plot.y))  
    }
    if(!is.null(plot.color)){
      p <- ggplot(combined, aes_string(x = plot.x, y = plot.y, color = plot.color))  
    }
  
    ### Plot: Add loadings as points
    if(plot.species == T){
      p <- p + geom_point(data = species, color = "grey")        
    }
  
    ### Plot: Add samples as points
    p <- p + geom_point(size = plot.point.size)
    
  
    ### Plot: If PCA add explained variance to the axis
    if(ordinate.type == "PCA"){
      p <- p + xlab(paste(plot.x, " [",exp[plot.x],"%]", sep = "")) +
        ylab(paste(plot.y, " [",exp[plot.y],"%]", sep = ""))  
    }    
    
    ### Plot: Plot the names of the n most extreme species
    if(!is.null(plot.nspecies)){
      p <- p + geom_text(data = species[1:plot.nspecies,], aes_string(x = plot.x, y = plot.y, label = "Genus"), colour = "black", size = 4, hjust = -0.05, vjust = 1)
    }
    
    ### Plot: Environmental factors
    if((!is.null(envfit.factor))){
      if (nrow(f.sig) != 0){
      p <- p + geom_text(data = f.sig, aes_string(x = plot.x, y = plot.y, label = "Name"), colour = "darkred", size = 4, hjust = -0.05, vjust = 1)
      }
    }
    
    ### Plot: Environmental numeric data
    if(!is.null(envfit.numeric)){
      if (nrow(n.sig) != 0){
      n.sig[, plot.x] <- n.sig[, plot.x]*envfit.resize
      n.sig[, plot.y] <- n.sig[, plot.y]*envfit.resize
      n.sig2 <- n.sig
      n.sig2[, plot.x] <- n.sig[, plot.x]*1.1
      n.sig2[, plot.y] <- n.sig[, plot.y]*1.1

      p <- p + geom_segment(data=n.sig, aes_string(x = 0, xend = plot.x, y = 0, yend = plot.y),arrow = arrow(length = unit(0.3, "cm")), colour = "darkred", size = 1) +
        geom_text(data=n.sig2, aes_string(x = plot.x, y = plot.y, label = "Name"), colour = "darkred")
      }
    }  
  
    ###Plot: Group the data either by centroid or chull
    if (!is.null(plot.group) & !is.null(plot.color)){
      ts <- data.frame(group = combined[,plot.color], x = combined[,plot.x], y = combined[,plot.y])
      os <- ddply(ts, ~group, summarize, cx = mean(x), cy = mean(y))
      os2 <- merge(combined, os, by.x=plot.color, by.y = "group")
      if (plot.group == "centroid"){
        p <- p + geom_segment(data=os2, aes_string(x = plot.x, xend = "cx", y = plot.y, yend = "cy", color = plot.color), size = 1)
      }
      
      if (plot.group == "chull"){
        find_hull <- function(df) df[chull(df[,plot.x], df[,plot.y]), ]
        hulls <- ddply(combined, plot.color, find_hull)
        p <- p + geom_polygon(data=hulls, aes_string(fill = plot.color), alpha = 0.2)
      }
      
      if (plot.group.label == T){
        p <- p + 
          geom_text(data=os2, aes_string(x = "cx", y = "cy", label = plot.color), size = 4, color = "black") 
      
      }
    }
  
  ### Plot: Label all samples using sample data
  
  if (!is.null(plot.label)){
    p <- p + geom_text(aes_string(label=plot.label), size = 3, vjust = 1.5, color = "grey40")
  }  
  
  outlist <- append(outlist, list(plot = p))
  
  ## Export the data
  
  if(output == "plot"){
    return(outlist$plot)
  }
  if(output == "complete"){
    return(outlist)
  }
  
}