# Script for analysing a grain boundary containing pos file
# Will determine start and end of interface, GB composition and width
# Will return a .pos file of just the GB region

#### Read intial .pos file ####
PosLocation <- c("Simulated GB.pos")
source("readposR.R")
PosFileInput <- readposR(PosLocation)

##### Combine .pos file with .rrng file to get chemical information ####
RangeFileLocation <- c("Range For Simulation.rrng")

source("Combine Pos and Range File.R")
PosFileRanger(PosFileInput,
              RangeFileLocation)

rm(Ion, PosFileInput, PosFileRanger, PosLocation, RangeFileLocation, readposR)

#### Generate 1D concentration profile ####
source("PosOneDimensionalPlot.R")
OneDCountFunc(RangedPos, 40, "z")

# Calculate extent of GB region

# Export .pos of GB region

# Contour plots etc. of GB region 

source("readposR.R")
RangeFilePath <- c("Range For Simulation.rrng")
PosFile <- sample_n(readposR("Simulated GB.pos"),1500)
