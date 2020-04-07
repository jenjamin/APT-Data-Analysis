# This function takes the IVAS cluster search csv output
# It separates ions into their constituent elements
# For RPV steel users, all Fm ions within the cluster are assumed to be Ni
# It then converts the cluster information into a new dataframe called "ClusterAtomCount"
# Written by Ben Jenkins Dec 2018

if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Rdisop")
}

CompositionCalculator <- function(DataframeName){
  ClusterInformation <- DataframeName  %>%
    filter(!grepl("Matrix",X)) %>%
    select("X",contains("Ions"), contains("Count")) %>%
    rename_(.dots=setNames(names(.), (gsub(".Count", "", names(.)))))
  
  ClusterInfoAllElements <<- ClusterInformation
  
  #### Converts ion counts to individual atom counts ####
  require(PeriodicTable)
  require(InterpretMSSpectrum)
  data("periodicTable")
  ElementsList <- c(unique(periodicTable$symb))
  ClusterInformationConverted <- ClusterInfoAllElements
  for(Ions in colnames(ClusterInfoAllElements %>% select(-"X",-contains("Ions")))){
    if(Ions %in% unique(periodicTable$symb) == FALSE){
      SplitIon <- CountChemicalElements(Ions)
      i <- 1
      for(DecomposedElements in names(SplitIon[])){
        ClusterInformationConverted <- ClusterInformationConverted %>%
          mutate(!!paste((DecomposedElements), sep = "_") := (!!rlang::sym(paste((DecomposedElements), sep = "_"))) + 
                   ((!!rlang::sym(paste((Ions), sep = "_"))) * as.numeric(SplitIon[i]))) 
        i < i + 1 
      }
      ClusterInformationConverted <<- ClusterInformationConverted %>% select(-Ions)
    }
  }
  #rm(DecomposedElements, i, Ions, SplitIon, ElementsList) # Remove vars from global envir
  
  
  ClusterInformationConvertedCount <<- ClusterInformationConverted %>%
    mutate(TotalAtoms = rowSums(select(.,-X, -Solute.Ions, -Ranged.Ions, -Total.Ions))) %>%
    select(-Total.Ions)
  
  #### Create dataframe ready to plot bar charts ####  
  ClusterInfoForBarChart <<- ClusterInformationConvertedCount %>%
    mutate(Other = TotalAtoms - (C + Zr + Nb + Fe + O + H)) %>%
    select(X, Solute.Ions, Ranged.Ions, C, Zr, Nb, Fe, O, H, Other, TotalAtoms) %>%
    gather(Element, Count, -X, -Solute.Ions, -Ranged.Ions, -TotalAtoms) %>%
    mutate(`Concentration (at.%)` = 100*Count/TotalAtoms,
           TotalAtoms = TotalAtoms + (parse_number(as.character(X))/10000))
}

ClusterBarChartPlot <- function(XAxisLabels){
  
  #### Allows cluster size to avoid overlaps ####  
  fmt_dcimals <- function(decimals=0){
    function(x) format(x,nsmall = decimals,scientific = FALSE)
  }
  
  #### Cluster Plot ####
  ClusterBarChart <- ggplot(ClusterInfoForBarChart,
                            aes(x=factor(TotalAtoms),
                                `Concentration (at.%)`,
                                fill=Element)) +
    geom_bar(stat = "identity")
  
  
  ClusterBarChart <- ClusterBarChart +
    theme_classic() + #Graph Theme
    scale_y_continuous(limits = c(0, 100.00001),
                       breaks=seq(0,100,10),
                       labels = fmt_dcimals(0)) + #y-axis limits and label distribution
    scale_x_discrete(labels = function(x) round(as.numeric(x), digits=0),
                     breaks = unique((ClusterInfoForBarChart %>% 
                                        filter(Element == "Zr") %>%
                                        arrange(TotalAtoms))$TotalAtoms)[seq(1,1000,XAxisLabels)]) +
    scale_fill_manual(values = c("C" = "#660033", "Fe" = "#FF00FF", "H" = "#00ff00",
                                 "Nb" = "#cc6600", "O" = "#00ffff", "Zr" = "#660066", "Other" = "#000000")) + #colour bar regions
    scale_colour_manual(values = c("C" = "#660033", "Fe" = "#FF00FF", "H" = "#00ff00",
                                   "Nb" = "#cc6600", "O" = "#00ffff", "Zr" = "#660066", "Other" = "#000000")) + #colour lines
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
