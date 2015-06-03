#' Calculate similarity over time
#'
#' A nice long description.
#'
#' @usage amp_stability(data)
#'
#' @param data (required) A phyloseq object including sample data.
#' @param date (required) The name of the variable that stores the date information.
#' @param group Split the dataset into selected groups using a metadata variable.
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

amp_stability <- function(data, date, group = NULL, plot.type = "time", plot.theme = "clean", output = "plot"){
  
  abund = t(as.data.frame(otu_table(data)@.Data))
  sample = suppressWarnings(as.data.frame(as.matrix(sample_data(data))))
  
  if (is.null(group)){
    sample$tgroup <- "A"
    group <- "tgroup"
  }
  
  
  betad <- vegdist(x = abund, method = "bray", binary = F)
  bt <- melt(as.matrix(betad)) 
  
  id <- colnames(sample)[1]
  
  bt2 <- merge(x = bt, y = sample[,c(id,date, group)], by.x = "Var1", by.y = id)
  
  colnames(bt2)[4:5] <- c("D1", "Plant1")
  bt2 <- merge(x = bt2, y = sample[,c(id,date, group)], by.x = "Var2", by.y = id)
  colnames(bt2)[6:7] <- c("D2", "Plant2")
  bt2$D1 <- as.Date(bt2$D1)
  bt2$D2 <- as.Date(bt2$D2)
  
  if (plot.type == "delta"){
    bt3 <- mutate(bt2, DeltaDays = as.numeric(abs(D1-D2))) %>%
      mutate(Keep = ifelse(Plant1 == Plant2 & DeltaDays != 0, "Yes", "No")) %>%
      subset(Keep == "Yes")
    
    p <- ggplot(bt3, aes(x = DeltaDays, y = 1-value, color = Plant1)) +
      geom_point() +
      xlab("Time between samples (Days)") +
      ylab("Similarity") +
      scale_color_discrete(name = "") +
      ylim(0, 1) +
      xlim(0,max(bt3$DeltaDays)) +
      geom_smooth(size = 1, se = F)     
  }
  
  if (plot.type == "time"){
    bt3 <- mutate(bt2, DeltaDays = as.numeric(D2-D1)) %>%
      filter(DeltaDays > 0) %>%
      group_by(Plant1, D1) %>%
      mutate(MinDays = min(DeltaDays)) %>%
      mutate(Keep = ifelse(Plant1 == Plant2 & DeltaDays == MinDays, "Yes", "No")) %>%
      subset(Keep == "Yes") %>%
      arrange(Plant1, D1) %>% 
      as.data.frame()
    
    p <- ggplot(bt3, aes(x = D2, y = 1-value, color = Plant1)) +
      geom_point() +
      geom_line() +
      ylim(0,1) +
      xlab("Date") +
      ylab("Similarity") +
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
