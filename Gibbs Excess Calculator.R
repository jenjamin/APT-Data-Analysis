# This script will take a one dimensional concentration profile from IVAS
# It will then allow you to plot 1D concentration profiles of individual elements
# Ensure that the csv file is of:
# 1 - fixed sample count profile
# 2 - atomic count and decomposed complex ions
# Written by Ben Jenkins - Feb 2019

setwd(dirname(parent.frame(2)$ofile))

#Area of ROI in nm
Area <- 225

# Detection Efiiciency of LEAP
DetectionEfficiency <- 0.37 

# CSV file path for 1D conc profile
CSVFilePath <- c("D:/DPhilY1/Thesis/ThesisAPTImages/NeutronIrrad/Run7727/GB ROI 15x15nm.csv")

source("Determine Gibbs Excess.R")
