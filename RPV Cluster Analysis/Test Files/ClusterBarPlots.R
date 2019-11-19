# This script will plot a bar chart showing the composition of each cluster
# The Fe content can be plotted alongside the cluster compositions to determine how Fe content varies with size
# Written by Ben Jenkins Dec 2018

ClusterBarChartPlot <- function(FeLinePlot,
                                XAxisLabels,
                                EdgeRemoval){

#### Allows cluster size to avoid overlaps ####  
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

if(EdgeRemoval == "N"){
  
}else{
  ClusterInfoForBarChart <- ClusterInfoForBarChart %>%
    spread(Element, `Concentration (at.%)`) %>%
    mutate(Edge = if_else(parse_number(as.character(X)) %in% EdgeClusters,"Edge",NA_character_))
  ClusterInfoForBarChart <- ClusterInfoForBarChart %>%
    filter(is.na(Edge)) %>%
    gather(Element, `Concentration (at.%)`, Fe, Si, Cu, Ni, Mn, P, Other)
}

#### Cluster Plot ####
ClusterBarChart <- ggplot(ClusterInfoForBarChart %>%
                            filter(Element != "Fe"),
                          aes(x=factor(TotalAtoms),
                              `Concentration (at.%)`,
                              fill=Element)) +
  geom_bar(stat = "identity")

if(FeLinePlot == "Y"){
  ClusterBarChart <- ClusterBarChart +
    geom_line(data=ClusterInfoForBarChart %>% filter(Element == "Fe"),
              aes(x=factor(TotalAtoms),
                  `Concentration (at.%)`,
                  group=1,
                  colour=Element),
              size=1.0)
} else{
}

ClusterBarChart <- ClusterBarChart +
  theme_bw() + #Graph Theme
  scale_y_continuous(limits = c(0, 100.00001),
                     breaks=seq(0,100,10),
                     labels = fmt_dcimals(0)) + #y-axis limits and label distribution
  scale_x_discrete(labels = function(x) round(as.numeric(x), digits=0),
                   breaks = unique((ClusterInfoForBarChart %>% 
                                      filter(Element == "Cu") %>%
                                      arrange(TotalAtoms))$TotalAtoms)[seq(1,1000,XAxisLabels)]) +
  scale_fill_manual(values = c("Si" = "#CCCCCC", "Fe" = "#FF00FF", "Cu" = "#FF6600",
                               "Ni" = "#00CC00", "Mn" = "#CCCC00", "P" = "#FF00FF", "Other" = "#000000"),
                    breaks=c("Si", "", "Cu", "Ni", "Mn", "P", "Other"),
                    labels = c("Si", "", "Cu", "Ni", "Mn", "P", "Other")) + #colour bar regions
  scale_colour_manual(values = c("Si" = "#CCCCCC", "Fe" = "#FF00FF", "Cu" = "#FF6600", "Ni" = "#00CC00",
                                 "Mn" = "#CCCC00", "P" = "#FF00FF", "Other" = "#000000")) + #colour lines
  xlab('Uncorrected Cluster Size (number of ions)') + #X-axis title
  ylab('Concentration (at.%)') + # y-axis title
  theme(text = element_text(face="bold", size=15, colour="black")) +
  theme(plot.title = element_text(face="bold", size = 15, hjust=0.5)) + #Centre title
  theme(axis.text.x = element_text(face="bold", size=12, colour="black", angle=90, vjust =0.5)) + #x-axis label size and font
  theme(axis.text.y = element_text(face="bold", size=12, colour="black", angle=90, vjust =0.5)) + #x-axis label size and font
  scale_y_continuous(limits = c(0, 100.00001), breaks=seq(0,100,10)) + #y-axis limits and label distribution
  theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1))

print(ClusterBarChart)
}
