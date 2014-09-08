#' Generates a ggplot2 style ordinate plot of amplicon data.
#'
#' A nice long description
#'
#' @usage amp_ordinate(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param scale A variable from the associated sample data to scale the abundance by.
#' @param trans Transform the raw counts, currently supports "none" or "sqrt" (default: "sqrt").
#' @param ordinate.type Either PCA or NMDS (default: "PCA").
#' @param constrain Constrain the PCA by a sample variable.
#' @param ncomp The number of principal components to extract if using PCA (default: 5)
#' @param plot.x Variable to plot on the x-axis if using PCA (default: "PC1")
#' @param plot.y Variable to plot on the y-axis if using PCA (default: "PC2")
#' @param plot.color Color the points by a sample variable.
#' @param plot.color.order Order the groups used for coloring by a vector.
#' @param plot.point.size Size of the plotted sample points (default: 3)
#' @param plot.species Plot loadings as points (default: F)
#' @param plot.nspecies Plot the n most extreme species with their taxonomic classification (default: 0).
#' @param plot.nspecies.tax Taxonomic level used in plot.nspecies (default: "Genus").
#' @param plot.label Label points using a sample variable.
#' @param plot.group Uses plot.color and groups samples by either "chull" or "centroid".
#' @param plot.group.label Add label based on the centroid of the specified group.
#' @param plot.group.label.size Text size of the labels.
#' @param trajectory Connects points based on a sample variable e.g. date.
#' @param trajectory.group Split the trajectory by a group.
#' @param envfit.factor A vector of factor variables from the sample data used for envfit to the model.
#' @param envfit.numeric A vector of numerical variables from the sample data used for envfit to the model.
#' @param envfit.significant The significance treshold for displaying envfit parameters (default: 0.01).
#' @param envfit.resize Scale the size of the numeric arrows (default: 1).
#' @param tax.empty Option to add "best" classification or just the "OTU" name to each "OTU" (default: best).
#' @param scale.species Rescale the plotted loadings to maximise visability (default: F).
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

amp_ordinate <- function(data, scale = NULL, trans = "sqrt", ordinate.type = "PCA", ncomp = 5, plot.x = "PC1", plot.y = "PC2", plot.color = NULL, plot.color.order = NULL, plot.point.size = 3, plot.species = F, plot.nspecies = NULL, plot.nspecies.tax = "Genus", plot.label = NULL, plot.group = NULL, plot.group.label = NULL, envfit.factor = NULL, envfit.numeric = NULL, envfit.significant = 0.001, envfit.resize = 1, tax.empty ="best", output = "plot", constrain = NULL, scale.species = F, trajectory = NULL, trajectory.group = trajectory, plot.group.label.size = 4){
  
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.empty = tax.empty)  
  
  ## Load the data and extract relevant dataframes from the phyloseq object
  abund<-as.data.frame(otu_table(data))
  tax<-as.data.frame(tax_table(data))
  sample <- suppressWarnings(as.data.frame(sample_data(data)))
  
  outlist <- list(abundance = abund, taxonomy = tax, sampledata = sample)
  
  ## Scale the data
  if (!is.null(scale)){
    variable <- unlist(sample[,scale])
    abund <- t(t(abund)*variable)
  }
  
  ## Transform the data
  abund1 <- abund
  if (trans == "sqrt"){
    abund1 <- sqrt(abund)
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
    
    if(scale.species == T){
      maxx <- max(abs(scores[,plot.x]))/max(abs(loadings[,plot.x]))
      loadings[, plot.x] <- loadings[, plot.x] * maxx * 0.8
      
      maxy <- max(abs(scores[,plot.y]))/max(abs(loadings[,plot.y]))
      loadings[, plot.y] <- loadings[, plot.y] * maxy * 0.8 
    }
    
    species <- cbind(loadings, loadings[,plot.x]^2 + loadings[,plot.y]^2)
    colnames(species)[ncol(species)] <- "dist"
    species <- species[with(species, order(-dist)), ]
    
    outlist <- append(outlist, list(nmds.model = model, nmds.scores = combined, nmds.loadings = species))
  }
  
  ## Calculate PCA
  
  if (ordinate.type == "PCA"){
    
    if(is.null(constrain)){
      model <- rda(t(abund1))  
      exp <- round(model$CA$eig/model$CA$tot.chi*100,1)
    }
    
    if(!is.null(constrain)){
      constrain1 <- sample[, constrain]
      colnames(constrain1) <- "constrain"
      model <- rda(t(abund1) ~ constrain1$constrain)  
      plot.x <- "RDA1"
      if (model$CCA$rank > 1){
        plot.y <- "RDA2"
      }
      exp <- round(model$CA$eig/model$tot.chi*100,1)
      expCCA <- round(model$CCA$eig/model$CA$tot.chi*100,1)
      exp <- c(exp, expCCA)
    }
    
    scores <- scores(model, choices = 1:ncomp)$sites
    combined <- cbind.data.frame(sample, scores)
    
    loadings <- cbind.data.frame(scores(model, choices = 1:ncomp)$species, tax)
    OTU <- gsub("denovo", "", rownames(loadings))
    loadings <- cbind.data.frame(loadings, OTU)
    
    if(scale.species == T){
      maxx <- max(abs(scores[,plot.x]))/max(abs(loadings[,plot.x]))
      loadings[, plot.x] <- loadings[, plot.x] * maxx * 0.8
      
      maxy <- max(abs(scores[,plot.y]))/max(abs(loadings[,plot.y]))
      loadings[, plot.y] <- loadings[, plot.y] * maxy * 0.8
      
    }
    
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
  
  ## Order the colors
  
  if (!is.null(plot.color.order)) {
    combined[,plot.color] <- factor(combined[,plot.color], levels = plot.color.order)
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
  
  ###Plot: Add trajectory based on e.g. data
  if (!is.null(trajectory)){
    traj <- combined[order(combined[,trajectory]),]
    p <- p + geom_path(data = traj, aes_string(group = trajectory.group))
  }
  
  ### Plot: Add samples as points
  p <- p + geom_point(size = plot.point.size)
  
  ### Plot: If PCA add explained variance to the axis
  if(ordinate.type == "PCA"){
    p <- p + xlab(paste(plot.x, " [",exp[plot.x],"%]", sep = "")) +
      ylab(paste(plot.y, " [",exp[plot.y],"%]", sep = ""))  
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
  }
  
  if (!is.null(plot.group.label)){
    ts <- data.frame(group = combined[,plot.group.label], x = combined[,plot.x], y = combined[,plot.y])
    os <- ddply(ts, ~group, summarize, cx = mean(x), cy = mean(y))
    os2 <- merge(combined, os, by.x=plot.group.label, by.y = "group")
    os3<- os2[!duplicated(os2[,plot.group.label]),]
    p <- p + geom_text(data=os3, aes_string(x = "cx", y = "cy", label = plot.group.label), size = plot.group.label.size , color = "black", fontface = 2) 
  }
  
  ### Plot: Plot the names of the n most extreme species
  if(!is.null(plot.nspecies)){
      p <- p + geom_text(data = species[1:plot.nspecies,], aes_string(x = plot.x, y = plot.y, label = plot.nspecies.tax), colour = "black", size = 4)  
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