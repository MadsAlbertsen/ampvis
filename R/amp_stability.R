#' Calculate similarity over time
#'
#' A nice long description.
#'
#' @usage amp_stability(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param date (required) The name of the variable that stores the date information.
#' @param group Split the dataset into selected groups using a metadata variable.
#' @param method Dissimilarity index from vegdist (default: "bray").
#' @param color Color the plots based on 1 of the grouping variables (default: unique colors).
#' @param order Vector to order the groups by,
#' @param plot.type Either "time" or "delta" (default: "time").
#' @param plot.theme Chose different standard layouts choose from "normal" or "clean" (default: "clean").
#' @param output Either "plot" or "complete" (default: "plot").
#' 
#' @export
#' @import phyloseq
#' @import dplyr
#' @import vegan
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_stability <- function(data, date, group = NULL, plot.type = "time", plot.theme = "clean", output = "plot", method = "bray", color = "Plant1", order = NULL){
  
  abund = t(as.data.frame(otu_table(data)@.Data))
  sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  
  if (is.null(group)){
    sample$tgroup <- "A"
    group <- "tgroup"
  }
  
  if (length(group) > 1){
    sample$newGroup <- do.call(paste, c(as.list(sample[,group]), sep=" "))
    oldGroup <- unique(cbind.data.frame(sample[,group], group = sample$newGroup))
    group <- "newGroup"
  }
  
  betad <- vegdist(x = abund, method = method, binary = F)
  bt <- melt(as.matrix(betad)) 
  
  id <- colnames(sample)[1]
  
  bt2 <- merge(x = bt, y = sample[,c(id, date, group)], by.x = "Var1", by.y = id)
  
  colnames(bt2)[4:5] <- c("D1", "Plant1")
  bt2 <- merge(x = bt2, y = sample[,c(id,date, group)], by.x = "Var2", by.y = id)
  colnames(bt2)[6:7] <- c("D2", "Plant2")
  bt2$D1 <- as.Date(bt2$D1)
  bt2$D2 <- as.Date(bt2$D2)
  
  if (plot.type == "delta"){
    bt3 <- mutate(bt2, DeltaDays = as.numeric(abs(D1-D2))) %>%
      filter(Plant1 == Plant2) %>%
      filter(Var1 != Var2) %>%
      mutate(group = Plant1) %>%
      mutate(Similarity = 1-value)
    
    if (group == "newGroup"){ bt3 <- merge(bt3, oldGroup)}
    
    if(!is.null(order)){
      bt3[,"Plant1"] <- factor(bt3[,"Plant1"], levels = order)  
    }
    
    
    p <- ggplot(bt3, aes_string(x = "DeltaDays", y = "Similarity", color = color, group = "Plant1")) +
      geom_point() +
      xlab("Time between samples (Days)") +
      ylab(paste("Similarity (", method, ")", sep = "")) +
      scale_color_discrete(name = "") +
      ylim(0, 1) +
      xlim(0,max(bt3$DeltaDays)) +
      geom_smooth(size = 1, se = F)     
  }
  
  if (plot.type == "time"){
    bt3 <- mutate(bt2, DeltaDays = as.numeric(D2-D1)) %>%
      filter(DeltaDays > 0) %>%
      filter(Plant1 == Plant2) %>%
      group_by(Plant1, D1) %>%
      mutate(MinDays = min(DeltaDays)) %>%
      filter(DeltaDays == MinDays) %>%
      group_by(Plant1, D1, Plant2, D2) %>%
      summarise(value = mean(value)) %>%
      arrange(Plant1, D1) %>% 
      mutate(group = Plant1) %>%
      mutate(Similarity = 1-value) %>%
      as.data.frame()
    
    if (group == "newGroup"){ bt3 <- merge(bt3, oldGroup) }
    
    if(!is.null(order)){
      bt3[,"Plant1"] <- factor(bt3[,"Plant1"], levels = order)  
    }
    
      p <- ggplot(bt3, aes_string(x = "D2", y = "Similarity", color = color, group = "Plant1")) +
      geom_point() +
      geom_line() +
      ylim(0,1) +
      xlab("Date") +
      ylab(paste("Similarity (", method, ")", sep = "")) +
      scale_color_discrete(name = "")
  }
  
  if (plot.theme == "clean"){
    p <- p +
      theme(axis.text = element_text(color = "black"),
            text = element_text(color = "black"),
            plot.margin = unit(c(0,0,0,0), "mm"),
            plot.margin = unit(c(0,0,0,0), "mm"),
            panel.grid = element_blank(),
            legend.key = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            legend.position = c(0.9,0.9)
      )
  }
  
  if(output == "plot"){
    return(p)  
  } else{
    list <- list(plot = p, data = bt3)
    return(list)
  }
  
}
