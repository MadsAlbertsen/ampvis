#' Generates correlation networks from amplicon data.
#'
#' A nice long description
#'
#' @usage amp_correlation(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param trans Transform the raw counts by "sqrt".
#' @param treshold.cor Absolute correlation treshold (default: 0.8).
#' @param treshold.count Minimum count in order to use a pair of observations for test of correlation (default: 10)
#' @param treshold.pval P-value treshold for correlations (default: 0.01).
#' @param treshold.abundance Average abundance treshold (default: 0) 
#' @param treshold.obs Minimum number of observations in order to calculate correlations (default: 10).
#' @param ignore.zero Remove observations where 1 obeservation has 0 counts from the correlation test (default: F)
#' @param ignore.global.zero Remove observations where 1 obeservation has 0 counts from the initial correlation test (default: T)
#' @param label Plot taxonomic classifications (e.g. "Genus") instead of points.
#' @param scale.abundance Scale the size of nodes by abundance (default: F).
#' @param scale.seq The number of sequences in the pre-filtered samples (default: 20000)
#' @param tax.empty Either "remove" OTUs without taxonomic information, add "best" classification or add the "OTU" name (default: best).
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.aggregate The taxonomic level that the data should be aggregated to (default: OTU).
#' @param scale.size Scale the size of the plotted objects (default: 1)
#' @param output Either plot or complete (default: "plot").
#' @param alpha Transparency of the plotted points (default: 1).
#' 
#' @return A ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @import dplyr
#' @import igraph
#' @import reshape2
#' @import phyloseq
#' @import grid
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}


amp_correlation <- function(data, trans = "none", treshold.cor = 0.8, treshold.pval = 0.01, tax.clean = T, ignore.zero = F, ignore.global.zero = T, scale.seq = 20000, treshold.abundance = 0, treshold.obs = 10, label = NULL, scale.abundance = F, output = "plot", treshold.count = 10, scale.size = 1, tax.class = NULL, tax.empty = "best", tax.aggregate = "OTU", alpha = 1){
  
  print("Loading, cleaning and transforming the data")
  ## Clean up the taxonomy
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Aggregate to a specific taxonomic level
  if (tax.aggregate != "OTU"){ data <- tax_glom(data, taxrank=tax.aggregate) }
  
  abund<-as.data.frame(otu_table(data))
  tax <- data.frame(tax_table(data), OTU = rownames(tax_table(data)))
  sample <- suppressWarnings(data.frame(sample_data(data)))
  
  ## Transform the data
  
  abund1 <- abund
  if (trans == "sqrt"){
    abund1 <- sqrt(abund)
    treshold.count <- sqrt(treshold.count)
    treshold.abundance <- sqrt(treshold.abundance)
  }
  
  ## Calculate correlation matrix
  print("Calculating the correlation matrix")
  if (ignore.global.zero == T){
    abund1[abund1 == 0] <- NA
  }
  
  cor <- cor(t(abund1), use = "pairwise.complete.obs")
  cor[lower.tri(cor)] <- NA
  cor2 <- melt(cor)
  cor3 <- cor2[!is.na(cor2$value), ]
  cor4 <- cor3[cor3$Var1 != cor3$Var2, ]
  
  ## Subset based on correlation
  
  cor5 <- subset(cor4, abs(value) >= treshold.cor)
  cor6 <- cbind.data.frame(cor5, pval = rep(1, nrow(cor5)), nobs = rep(0, nrow(cor5)))  
  
  ## Make formal correlation test of the subset of points
  
  print(paste("Testing", nrow(cor6), "of", nrow(cor4) ,"correlations for significance", sep = " "))
  tcount <- 0
  pmt <- proc.time()
  
  for (i in 1:nrow(cor6)){
    x <- as.numeric(abund1[rownames(abund1)==cor6$Var1[i],])
    y <- as.numeric(abund1[rownames(abund1)==cor6$Var2[i],])
    c <- cbind(x,y)
    if (ignore.zero == T){
      c[c==0] <- NA
    } else {      
      c[is.na(c)] <- 0
    }
    c <- subset(c, !is.na(x) | !is.na(y)) 
    c2 <- subset(c, c[,1]+c[,2] >= treshold.count)
    if (nrow(c2) >= treshold.obs){
      cor6$pval[i] <- cor.test(c2[,1],c2[,2])$p.value*nrow(abund)  
      cor6$value[i] <- cor.test(c2[,1],c2[,2])$estimate 
    }  
    if (!is.null(nrow(c2))) { cor6$nobs[i] <- nrow(c2)}
    tcount <- tcount +1
    if(tcount == 1000){
      t <- as.numeric(proc.time()[3] - pmt[3])
      print(paste("Tested", i, "correlations in", t, "s", sep = " "))
      tcount <- 0
    }
  }
  
  cor7 <- subset(cor6, pval <= treshold.pval & abs(value) >= treshold.cor)
  
  print(paste(nrow(cor7), "of", nrow(cor6) ,"correlations were significant", sep = " "))
  
  ## Format the data to graph format using igraph
  
  print("Generating a correlation plot")
    
  g <- graph.data.frame(cor7, directed = F)
  
  ## Calculate a graph   
  
  t <- layout.fruchterman.reingold(g)
  
  ## Extract the individual points in the graph and add abundance 
    
  gpoints <- data.frame( "OTU" = V(g)$name, "x" = t[,1], "y" = t[,2])
  gpoints1 <-merge(gpoints, tax)
  
  TotalSum <- apply(abund, 1, function(x) mean(x)/scale.seq*100)
  TotalSum1 <- data.frame("OTU" = names(TotalSum), "Abundance" = TotalSum)
  
  gpoints2 <- merge(TotalSum1, gpoints1)
  
  ## Subset based on average abundance
  
  gpoints3 <- subset(gpoints2, Abundance >= treshold.abundance)
  
  ## Extract links between scaffolds
  
  seq <- merge(cor7, gpoints3[,c(1,3,4)], by.x = "Var1", by.y = "OTU")
  seg1 <- merge(seq, gpoints3[,c(1,3,4)], by.x = "Var2", by.y = "OTU")
  colnames(seg1) <- c("OTU1","OTU2","Correlation","pval", "nobs", "x", "y", "xend", "yend")
  
  ## Make a ggplot2 object
  
  p <- ggplot(data = gpoints3, aes(x = x, y = y)) +
    geom_segment(data=seg1, aes(x=x, y=y, xend=xend, yend=yend, color = Correlation), size = 1) +
    geom_point(size = 4*scale.size, alpha = alpha) +
    scale_color_continuous(low="red", high = "green", limits=c(-1,1), name = "Correlation (r)")
  
  
  if (!is.null(label)){
    p <- ggplot(data = gpoints3, aes_string(x = "x", y = "y", label = label)) +
      geom_segment(data=seg1, aes(x=x, y=y, xend=xend, yend=yend, label = NA, color = Correlation), size = 1) +
      geom_text(size = 4*scale.size, color = "black") +
      scale_color_continuous(low="red", high = "green", limits=c(-1,1), name = "Correlation (r)")
  }
  
  if (scale.abundance == T & is.null(label)){
    p <- ggplot(data = gpoints3, aes(x = x, y = y, size = Abundance)) +
      scale_size_area(name = "Abundance (%)", max_size=20*scale.size, breaks = c(0.01, 0.05, 0.1, 0.5, 1)) +
      geom_segment(data=seg1, aes(x=x, y=y, xend=xend, yend=yend, color = Correlation), size = 1) +
      geom_point(alpha = alpha) +
      scale_color_continuous(low="red", high = "green", limits=c(-1,1), name = "Correlation (r)")
  }

  if (scale.abundance == T & !is.null(label)){
    p <- ggplot(data = gpoints3, aes_string(x = "x", y = "y", size = "Abundance", label = label)) +
      scale_size_area(name = "Abundance (%)", max_size=10*scale.size, breaks = c(0.01, 0.1, 1)) +
      geom_segment(data=seg1, aes(x=x, y=y, xend=xend, yend=yend, color = Correlation, label = NA), size = 1) +
      geom_text(color = "black") +
      scale_color_continuous(low="red", high = "green", limits=c(-1,1), name = "Correlation (r)")
  }  
  
  ## Remove background lines
  
  p <- p + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank())
  
  colnames(cor7) <- c("x", "y", "R", "pval", "nobs")
  cor7 <- arrange(cor7, pval)
  
  if (output == "complete") { return(list(plot = p, cor = cor7, abund = abund, tax = tax, metadata = sample))} else {return(p)}
}