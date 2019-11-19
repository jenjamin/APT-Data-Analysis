# This script will calulcate cluster size, Volume Fraction, and Number Density 
# Size calculations assume clusters are spherical in shape
# Written by Ben Jenkins Dec 2018

ClusterCalculations <- function(AtomicDensity, DetectionEfficiency){
  
  #Dataframe with estimated cluster size (including and excluding Fe)
  EdgeClustersIdentity <- ClusterInformationConvertedCount %>%
                                              filter(grepl("Cluster", X))%>%
                                  mutate(Edge = if_else(parse_number(as.character(X)) %in% EdgeClusters,"Edge",NA_character_))
  EdgeClustersIdentity <- EdgeClustersIdentity %>%
    mutate(`Estimated Size (with Fe)` = (TotalAtoms/DetectionEfficiency)/AtomicDensity,
           `Estimated Size (No Fe)`= (NonFeAtoms/DetectionEfficiency)/AtomicDensity,
           `Estimated Diameter (with Fe)` = 2*((3*`Estimated Size (with Fe)`)/(4*pi))^(1/3),
           `Estimated Diameter (No Fe)` = 2*((3*`Estimated Size (No Fe)`)/(4*pi))^(1/3))
  
  TotalNumberAtomsInDataset <- sum(ClusterInformationConvertedCount$TotalAtoms)
  
  VolumeOfTip <- (10^-27)*(TotalNumberAtomsInDataset/DetectionEfficiency)/AtomicDensity
  
  `Volume Fraction (%)` <- 100 *
    sum(EdgeClustersIdentity$NonFeAtoms)/
    TotalNumberAtomsInDataset
  
  EdgeClustersRemoved <- EdgeClustersIdentity %>% filter(is.na(Edge))
  
  NumberDensity <- (nrow(EdgeClustersIdentity) - 0.5 * nrow(EdgeClustersIdentity %>% filter(is.na(Edge))))/VolumeOfTip
  
  `Average Volume No Fe (nm3)` <- mean(EdgeClustersRemoved$`Estimated Size (No Fe)`)
  `SE Volume No Fe (nm3)` <- sd(EdgeClustersRemoved$`Estimated Size (No Fe)`)/
    sqrt(length(EdgeClustersRemoved$`Estimated Size (No Fe)`))
  `Average Diameter No Fe (nm)`<- mean(EdgeClustersRemoved$`Estimated Diameter (No Fe)`)  
  `SE Diameter No Fe (nm)` <- sd(EdgeClustersRemoved$`Estimated Diameter (No Fe)`)/
    sqrt(length(EdgeClustersRemoved$`Estimated Diameter (No Fe)`))
  
  `Average Volume with Fe (nm3)` <- mean(EdgeClustersRemoved$`Estimated Size (with Fe)`)
  `SE Volume with Fe (nm3)` <- sd(EdgeClustersRemoved$`Estimated Size (with Fe)`)/
    sqrt(length(EdgeClustersRemoved$`Estimated Size (with Fe)`))
  `Average Diameter with Fe (nm)`<- mean(EdgeClustersRemoved$`Estimated Diameter (with Fe)`)  
  `SE Diameter with Fe (nm)` <- sd(EdgeClustersRemoved$`Estimated Diameter (with Fe)`)/
    sqrt(length(EdgeClustersRemoved$`Estimated Diameter (with Fe)`))
  
  ClusterCalculationsSummary <<- cbind(NumberDensity,
                                      `Volume Fraction (%)`,
                                      `Average Volume No Fe (nm3)`,
                                      `SE Volume No Fe (nm3)`,
                                      `Average Diameter No Fe (nm)`,
                                      `SE Diameter No Fe (nm)`,
                                      `Average Volume with Fe (nm3)`,
                                      `SE Volume with Fe (nm3)`,
                                      `Average Diameter with Fe (nm)`,
                                      `SE Diameter with Fe (nm)`)
  print(ClusterCalculationsSummary)

}

