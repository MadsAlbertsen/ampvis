#' Tests if there is a significant difference in abundance between selected conditions
#'
#' A nice long description
#'
#' @usage amp_test_species(data, design)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param group (required) The group to test against.
#' @param sig Significance treshold (default: 0.01).
#' @param fold Log2fold filter default for displaying significant results (default: 0)
#' @param tax.aggregate Group data at specific taxonomic level (default: "OTU").
#' @param tax.class Converts a specific phyla to class level instead (e.g. "p__Proteobacteria").
#' @param tax.empty Either "remove" OTUs without taxonomic information at X level, with "best" classification or add the "OTU" name (default: best).
#' @param tax.display Display additional taxonomic levels in the plot output e.g. "Genus".
#' @param label Label the significant entries with tax.display (default:F).
#' @param plot.type Either "boxplot" or "point" (default: point)
#' @param plot.show Display the X most significant results.
#' @param plot.point.size The size of the plotted points.
#' 
#' @return A p-value for each comparison.
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

amp_test_species <- function(data, group, tax.aggregate = "OTU", tax.display = NULL, test = "Wald", fitType = "parametric", sig = 0.01, fold = 0, tax.class = NULL, tax.empty = "best", scale.seq = 10000, label = F, plot.type = "point", plot.show = NULL, plot.point.size = 2){
  
  
  ## Clean the taxonomy  
  data <- amp_rename(data = data, tax.class = tax.class, tax.empty = tax.empty, tax.level = tax.aggregate)
  
  ## Aggregate to a specific taxonomic level
  if (tax.aggregate != "OTU"){ data <- tax_glom(data, taxrank=tax.aggregate) }
  
  ## Convert to DESeq2 object
  groupF <- as.formula(paste("~", group, sep=""))
  data_deseq = phyloseq_to_deseq2(physeq=data, design=groupF)
  
  ## Test for significant differential abundance
  data_deseq_test = DESeq(data_deseq, test=test, fitType=fitType)
  
  ## Extract the results
  res = results(data_deseq_test, cooksCutoff = FALSE)  
  res_tax = cbind(as.data.frame(res), as.matrix(tax_table(data)[rownames(res), ]), OTU = rownames(res))
  
  res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))
  res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]
  
  if (nrow(res_tax_sig) > 1){  
  
  ## Plot the data
  ### MA plot
  res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
  res_tax$Significant[is.na(res_tax$Significant)] <- "No"
  
  p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) + 
    geom_point(size = plot.point.size) +
    scale_x_log10() +
    scale_color_manual(values=c("black", "red")) +
    labs(x = "Mean abundance", y = "Log2 fold change")
  
  if(label == T){
    if (!is.null(tax.display)){
      rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
    }  else {
      rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
    }
    p1 <- p1 + geom_text(data = subset(res_tax_label, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
      
  }
  
    
  ### Points plot of significant differential abundant entries
  res_tax_sig_abund = cbind(as.matrix(tax_table(data)[rownames(res_tax_sig)]), as.data.frame(otu_table(data)[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"])    
  
  res_tax_sig_abund_long <- melt(res_tax_sig_abund, id.vars=c(colnames(tax_table(data)),"OTU","padj"),value.name="Count", variable.name="Sample")
  
  colnames(sample_data(data))[1] <- "Sample"
  metadata <- suppressWarnings(as.matrix(sample_data(data))[,c("Sample",group)])
    
  point_df <- merge(x = res_tax_sig_abund_long, y = metadata, by = "Sample")
  point_df$Abundance <- point_df$Count / scale.seq * 100  
  
  if (!is.null(tax.display)){
    point_df <- data.frame(point_df, Display = apply(point_df[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  }  else {
    point_df <- data.frame(point_df, Display = point_df[,tax.aggregate])
  }
    
  point_df$Display <- factor(point_df$Display, levels = rev(as.character(unique(point_df$Display))))
  
  if(!is.null(plot.show)){
    point_df <- subset(point_df, Display %in% as.character(unique(point_df$Display)[1:plot.show]))
  }
  
  
  p2 <-ggplot(data = point_df, aes_string(x = "Display", y = "Abundance", color = group)) +
    labs(x = "", y = "Read Abundance (%)") +
    coord_flip()  
  
  if (plot.type == "point"){
    p2 <- p2 + geom_jitter(position = position_jitter(width = .05), size = plot.point.size)
  } else{
    p2 <- p2 + geom_boxplot(outlier.size=1)
  }
  
  out <- list(results = res, plot_MA = p1, sig_res = res_tax_sig, plot_sig = p2 , sig_res_plot_data = point_df)
  
  return(out)
  } else{
    print("No significant")
    return(res)
  }
}
