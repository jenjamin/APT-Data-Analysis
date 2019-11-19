# This script allows you to import a cluster csv file
# Cluster composition bar charts can then be plotted
# Ternary plots for ratios of Mn Ni and Si can be plotted
# Cluster Sizes, Volume Fractions, Number Densities can then be calcultaed

setwd(dirname(parent.frame(2)$ofile))

require(tidyverse)
source("Test Files//CSVConverter.R")
source("Test Files//Edge Cluster Detection.R")
source("Test Files//ClusterBarPlots.R")
source("Test Files//TernaryPlots.R")
source("Test Files//ClusterCalculator.R")

# Type the paths to your cluster analysis csv, pos file, and indexed range file
# Ensure that "/" in file path are re-written as "//" or "/"
`Cluster Analysis CSV` <- c("Test Files//006_full - Top-Level ROI - Cluster Analysis (Ni, Cu).csv")
`Original POS File Location` <- c("Test Files//006_full.pos")
`Indexed Cluster Pos File` <- c("Test Files//006_full.cluster.indexed.pos")
`Cluster Pos File` <- c("Test Files//006_full.cluster.pos")
`Atomic Density of Material` <- 84.59 # Atomic density of BCC Fe is ~ 85.49 atoms/nm3
`Detection Efficiency of LEAP` <- 0.37 # Detection efficiency of reflectron LEAP 3000 ~ 0.37, reflectron LEAP 5000 ~ 0.52
SamplingFraction <- 0.005 # 
AlphaNND <- 2

# Imports csv
ClusterCSVImport <- read.csv(paste0(`Cluster Analysis CSV`),
                             skip = 10)

# Run this function to generate dataframes for further analyses
CSVInputConverter(ClusterCSVImport)

# This function will determine which clusters are on the edge of your dataset
# You need to decide whether to include these in cluster composition calculations
# They will not be included in cluster size calculations
# This function may take a minute or two to run. Lower sampling fractions should increase speed.
findEdgeClustersConvex(`Original POS File Location`,
                       `Cluster Analysis CSV`,
                       `Indexed Cluster Pos File`,
                       `Cluster Pos File`,
                       SamplingFraction,
                       AlphaNND)

# View Alpha shape that is generated
plot(AlphaHullShape, indexAlpha = 1)
rm(AlphaHullShape)

# Run this function to plot bar chart of cluster compositions - object will be named ClusterBarChart
# Set first object to "Y" if you want Fe concentration plotting - else set to "N"
# Second option allows you to choose how often labels appear on x-axis
ClusterBarChartPlot(FeLinePlot = "Y",
                    XAxisLabels = 5,
                    EdgeRemoval = "Y")

# Run this function to plot a ternary plot showing ratio of Ni:Mn:Si for each cluster
# Option for plotting composition of edge clusters: "N" will plot edge clusters, "Y" will remove them from analyses
TernaryPlotFunction(EdgeRemoval = "Y")

# Perform calculations for cluster size (edge clusters not included), volume fraction, number density (edge clusters = 0.5)
# Option 1 are the atomic density of the material you are studying (85.49 atoms/nm3 for BCC Fe)
# Option 2 detection efficiency of intstrument used (0.37 for LEAP 3000, 0.52 for LEAP 5000)
ClusterCalculations(AtomicDensity = `Atomic Density of Material`,
                    DetectionEfficiency = `Detection Efficiency of LEAP`)

rm(periodicTable, ClusterCSVImport, ClusterInfoForBarChart, ClusterInfoAllElements, ClusterInformationConvertedCount)

# Run this if you want to copy cluster calculations summary table to your clipboard (can then paste into excel etc.)
write.table(ClusterCalculationsSummary, 
            "clipboard",
            sep="\t", row.names=FALSE, col.names=TRUE)
