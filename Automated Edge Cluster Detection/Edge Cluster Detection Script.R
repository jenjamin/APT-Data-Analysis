# Script to identifying edge clusters in APT datasets
# Required input of POS File, Cluster Range File, Cluster Analysis CSV
# Written by Ben Jenkins - Feb 2019

setwd(dirname(parent.frame(2)$ofile))

#Import: POS File, Cluster Range File, Cluster Analysis CSV in that order
FilePaths <- c("Test Files/006_full.pos",
               "Test Files/006_full.cluster.indexed.rrng",
               "Test Files/006_full - Top-Level ROI - Cluster Analysis (Cu, Ni).csv")

#### Parameters For Calculating Alpha Value####
# How much one would like to sample the POS file by
SamplingFraction <- 0.005 

source("AlphaHull.R")

# View Alpha shape that is generated
plot(AlphaHullShape, indexAlpha = 1)