setwd(dirname(parent.frame(2)$ofile))

source("Scripts/Edge Cluster Detection.R")
source("Scripts/ClusterCompositions.R")

OriginalPosFile= "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\R5083_08811-v01.pos"

#### NbFe Clusters ####
clusterCSVFilePathNbFe = "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\R5083_08811-v01 - Top-Level ROI - Cluster Analysis (Nb, Fe).csv"
clusterIndexedPosNbFe = "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\R5083_08811-v01.cluster.indexed.pos"
clusterPosNbFe = "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\R5083_08811-v01.cluster.pos"
#### Fe Clusters ####
clusterCSVFilePathFe = "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\R5083_08811-v01.matrix - Top-Level ROI - Cluster Analysis (Fe).csv"
clusterIndexedPosFe = "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\Fe Clusters.cluster.indexed.pos"
clusterPosFe = "D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\Fe Clusters.cluster.pos"

#### NbFe Clusters ####
findEdgeClustersConvex(
  posFileName = OriginalPosFile,
  clusterStatsFile = clusterCSVFilePathNbFe,
  clusterIndexPosFile = clusterIndexedPosNbFe,
  clusterPosFile = clusterPosNbFe,
  SamplingFraction = 0.005,
  NNDMultiplier = 2
)
# Calculate Number Density and Vf
DatasetVolume = (sum(read.csv(clusterCSVFilePathNbFe, skip = 10) %>% select(Ranged.Ions))/(42.9*0.52))*(10^-27)
NumberDensityNbFe = (nrow(read.csv(clusterCSVFilePathNbFe, skip = 10)) - (0.5 * length(EdgeClusters)))/DatasetVolume
VfNbFe = 100 * sum(read.csv(clusterCSVFilePathNbFe, skip = 10) %>%
                   filter(!grepl("Matrix",X)) %>%
                   select(Ranged.Ions))/
  sum(read.csv(clusterCSVFilePathNbFe, skip = 10) %>% select(Ranged.Ions))

CompositionCalculator(read.csv(clusterCSVFilePathNbFe, skip = 10) %>%
                        filter(!grepl("Matrix",X)) %>% 
                        filter(!parse_number(as.character(X)) %in% EdgeClusters))
rm(ClusterInfoAllElements, ClusterInformationConverted, ClusterInformationConvertedCount, periodicTable)

MeanClusterCompositionNbFe <- ClusterInfoForBarChart %>% 
  group_by(Element) %>%
  summarise(Mean = mean(`Concentration (at.%)`),
            StdErr = sd(`Concentration (at.%)`)/sqrt(n()))
ClusterBarChartPlot(1)

#### Fe Clusters ####
findEdgeClustersConvex(
  posFileName = OriginalPosFile,
  clusterStatsFile = clusterCSVFilePathFe,
  clusterIndexPosFile = clusterIndexedPosFe,
  clusterPosFile = clusterPosFe,
  SamplingFraction = 0.005,
  NNDMultiplier = 2
)
# Calculate Number Density and Vf
NumberDensityFe = (nrow(read.csv(clusterCSVFilePathFe, skip = 10)) - (0.5 * length(EdgeClusters)))/DatasetVolume
VfFe = 100 * sum(read.csv(clusterCSVFilePathFe, skip = 10) %>%
                   filter(!grepl("Matrix",X)) %>%
                   select(Ranged.Ions))/
  sum(read.csv(clusterCSVFilePathNbFe, skip = 10) %>% select(Ranged.Ions))

CompositionCalculator(read.csv(clusterCSVFilePathFe, skip = 10) %>%
                               filter(!grepl("Matrix",X)) %>% 
                               filter(!parse_number(as.character(X)) %in% EdgeClusters))
rm(ClusterInfoAllElements, ClusterInformationConverted, ClusterInformationConvertedCount, periodicTable)

MeanClusterCompositionFe <- ClusterInfoForBarChart %>% 
  group_by(Element) %>%
  summarise(Mean = mean(`Concentration (at.%)`),
            StdErr = sd(`Concentration (at.%)`)/sqrt(n()))
ClusterBarChartPlot(1)

write.table(MeanClusterCompositionNbFe, "clipboard", sep="\t", row.names=FALSE)
