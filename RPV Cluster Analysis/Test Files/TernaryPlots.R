# This script will plot a ternary plot showing ratio of Ni:Mn:Si for each cluster
# Written by Ben Jenkins Dec 2018

TernaryPlotFunction <- function(EdgeRemoval){
  require(ggtern)
  if(EdgeRemoval == "N"){
  
  TernaryPlot <<- ggtern(ClusterInfoForBarChart %>%
                           spread(Element, `Concentration (at.%)`) %>%
                           mutate(Edge = if_else(parse_number(as.character(X)) %in% EdgeClusters,"Edge",NA_character_)),
                        aes(Ni, Si, Mn)) + 
    geom_point(color = "red") + theme_bw() +
    theme(text = element_text(size=20, colour = "black"), 
          panel.grid.major = element_line(colour = "black"),
          axis.text = element_text(size=20, colour = "black", face="bold"),
          plot.title = element_text(hjust = 0.5))+ 
    theme(legend.position="right") + 
    guides(colour = guide_legend(ncol = 1),
           shape = guide_legend(ncol = 2)) +
    scale_shape_manual(values=seq(0,15))
  print(TernaryPlot)
  }else{
    TernaryPlot <<- ggtern(ClusterInfoForBarChart %>%
                             spread(Element, `Concentration (at.%)`) %>%
                             mutate(Edge = if_else(parse_number(as.character(X)) %in% EdgeClusters,"Edge",NA_character_)) %>%
                             filter(is.na(Edge)),
                           aes(Ni, Si, Mn)) + 
      geom_point(color = "red") + theme_bw() +
      theme(text = element_text(size=20, colour = "black"), 
            panel.grid.major = element_line(colour = "black"),
            axis.text = element_text(size=20, colour = "black", face="bold"),
            plot.title = element_text(hjust = 0.5))+ 
      theme(legend.position="right") + 
      guides(colour = guide_legend(ncol = 1),
             shape = guide_legend(ncol = 2)) +
      scale_shape_manual(values=seq(0,15))
    print(TernaryPlot)
  }
  detach("package:ggtern", unload = TRUE)
  detach("package:tidyverse", unload = TRUE)
  detach("package:ggplot2", unload = TRUE)
  library(tidyverse)
}


