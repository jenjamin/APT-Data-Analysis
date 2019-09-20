StartTime <- Sys.time()
#This function allows one to import a pos file into R with x, y, z, and mass-to-charge-ratio information
readposR <- function(posFileName) {
  # size of floats = 4 bytes
  sizeOfFloat = 4
  # get file size
  fileInfo <- file.info(posFileName)
  fileSize <- fileInfo$size / sizeOfFloat
  to.read = file(posFileName, "rb")
  posFile <- readBin(to.read,
                     double(),
                     size = 4,
                     n = fileSize,
                     endian = "big")
  close(to.read)
  posFile <-
    t(matrix(posFile, nrow = 4, dimnames = list(c("x", "y", "z", "m"))))
  posFile <- as.data.frame(posFile)
  return(posFile)
}

#This function creates a new range file (in the same location as the original range file)
#The new range file will have edge clusters coloured red and non-edge clusters coloured blue
EdgeClusterRangeFileWriter <- function(RangeFilePath) {
  ClusterRangeFileImp <-
    read.delim(RangeFilePath)
  
  ClusterRangeFile <- ClusterRangeFileImp %>%
    filter(grepl("Color", `X.Ions.`))  %>%
    separate(`X.Ions.`, c("X.Ions.", "Colour"), sep = "Color:")
  
  ClusterEdgeColour <-
    cbind(ClusterRangeFile, EdgeForExcelClusterImport)
  
  
  ClusterEdgeColour <- ClusterEdgeColour %>% 
    rowwise() %>%
    mutate(Colour2 = if_else(
      is.na(Edge),
      paste0(sample(c(
        "00", "05", "10", "15", "20", "30", "35", "40", "45", "50", "55", "60"
      ), 1),
      sample(c(
        "00", "05", "10", "15", "20", "30", "35", "40", "45", "50", "55", "60"
      ), 1),
      "FF"),
      paste0("FF",
             sample(c(
               "00", "05", "10", "15", "20", "30", "35", "40", "45", "50", "55", "60"
             ), 1),
             sample(c(
               "00", "05", "10", "15", "20", "30", "35", "40", "45", "50", "55", "60"
             ), 1))
    )) %>%
    mutate(
      Colour2 = paste0(" Color:",Colour2)
    ) %>%
    select(X.Ions., Colour2) %>%
    unite(X.Ions., `X.Ions.`, Colour2, sep ="")
  
  RangeFileExp <- rbind(ClusterRangeFileImp %>%
                          filter(!grepl("Color", `X.Ions.`)),
                        ClusterEdgeColour)
  
  colnames(RangeFileExp) <- c("[Ions]")
  
  write.table(
    RangeFileExp,
    file = gsub(".rrng", ".EdgeClusterRangeFile.rrng", RangeFilePath),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
}

library(tidyverse)
library("geometry", lib.loc = "~/R/win-library/3.4")
library(spatstat)

#### Get file paths of pos, range file and cluster information ####
posFilePath <- FilePaths[1]
IndexedRangeFilePath <- FilePaths[2]
ClusterImportPath <- FilePaths[3]

#pos file path
posFile <-
  readposR(posFilePath)

#cluster file path
ClusterImport <-
  read.csv(ClusterImportPath,
           skip = 10
  )


#### Cluster Import ####
ClusterImport <-
  ClusterImport %>% select(
    X,
    Solute.Ions,
    Ranged.Ions,
    Total.Ions,
    Center_x..nm..Ranged,
    Center_y..nm..Ranged,
    Center_z..nm..Ranged,
    Extent_x..nm..Ranged,
    Extent_y..nm..Ranged,
    Extent_z..nm..Ranged
  ) %>%
  mutate(
    x_max = Center_x..nm..Ranged + Extent_x..nm..Ranged,
    x_min = Center_x..nm..Ranged - Extent_x..nm..Ranged,
    y_max = Center_y..nm..Ranged + Extent_y..nm..Ranged,
    y_min = Center_y..nm..Ranged - Extent_y..nm..Ranged,
    z_max = Center_z..nm..Ranged + Extent_z..nm..Ranged,
    z_min = Center_z..nm..Ranged - Extent_z..nm..Ranged
  ) %>%
  filter(grepl("Cluster", X))

#### Find extreme co-ords for each cluster (modelling as cuboids) ####
ClusterImport <- ClusterImport %>%
  mutate(
    ClusterID = as.numeric(gsub("Cluster ", "", X)),
    Point1 = paste0(x_max, ",", y_max, ",", z_max),
    Point2 = paste0(x_max, ",", y_max, ",", z_min),
    Point3 = paste0(x_max, ",", y_min, ",", z_max),
    Point4 = paste0(x_max, ",", y_min, ",", z_min),
    Point5 = paste0(x_min, ",", y_max, ",", z_max),
    Point6 = paste0(x_min, ",", y_max, ",", z_min),
    Point7 = paste0(x_min, ",", y_min, ",", z_max),
    Point8 = paste0(x_min, ",", y_min, ",", z_min)
  ) %>%
  select(ClusterID,
         Point1,
         Point2,
         Point3,
         Point4,
         Point5,
         Point6,
         Point7,
         Point8) %>%
  gather(Point, value, -ClusterID) %>%
  arrange(ClusterID) %>%
  separate(value, c("x", "y", "z"), sep = ",") %>%
  mutate(x = as.numeric(x),
         y = as.numeric(y),
         z = as.numeric(z))

ClusterLocationMatrix <-
  as.matrix(ClusterImport %>% select(x, y, z))

#### Filter pos file to improve speed ####
PosFraction <- floor(nrow(posFile) * SamplingFraction)
FilterPosFile <- posFile %>% sample_n(PosFraction)
AlphaValue <- 2*round(ceiling(100*max(nndist(FilterPosFile %>% select(x,y,z), k=1))),2)/100
print(paste0(
  "Alpha Value: ",
  AlphaValue,
  " Sampling fraction: ",
  SamplingFraction
))

#### Using alpha-shape 3d ####
library("alphashape3d", lib.loc = "~/R/win-library/3.4")
AlphaHullShape <-
  ashape3d(as.matrix(FilterPosFile[, 1:3]), alpha =  c(AlphaValue, 2*AlphaValue))
#plot(AlphaHullShape, indexAlpha = 1)

#### Ensuring there is only one connected volume ####
comp <- components_ashape3d(AlphaHullShape, indexAlpha = "all")
NumberVolumes <- as.data.frame(table(comp[[1]]))


#### Determing which clusters are edge clusters and returning values ####
if (nrow(NumberVolumes) != 1) {
  print(paste0("Error. Multiple volumes created"))
} else{
  print(paste0("Alpha value has created one volume"))
  
  # determining which clusters are edge clusters
  ClustersInHull <-
    inashape3d(AlphaHullShape, indexAlpha = 1, ClusterLocationMatrix)
  
  # getting names of clusters that are edge
  ClusterLocation2 <- cbind(ClusterImport, ClustersInHull) %>%
    filter(ClustersInHull == "FALSE")
  
  TotalEdgeClusters <-
    as.data.frame(unique(ClusterLocation2$ClusterID))
  TotalEdgeClusters <-
    TotalEdgeClusters %>% transmute(ID = gsub("Cluster ", "", unique(ClusterLocation2$ClusterID)))
  TotalEdgeClusters <- sort(TotalEdgeClusters$ID)
  
  EdgeForExcelClusterImport <-
    ClusterImport %>% select(ClusterID) %>% distinct(ClusterID)
  EdgeForExcelClusterImport <-
    EdgeForExcelClusterImport %>% mutate(Edge = ifelse(
      match(EdgeForExcelClusterImport$ClusterID, TotalEdgeClusters),
      "Edge",
      ""
    ))
  write.table(EdgeForExcelClusterImport,
              "clipboard",
              sep = "\t",
              row.names = FALSE)
  # write.table(
  #   EdgeForExcelClusterImport %>% filter(Edge == "Edge") %>% select(ClusterID),
  #   "clipboard",
  #   sep = "\t",
  #   row.names = FALSE
  # )
  nrow(EdgeForExcelClusterImport %>% filter(Edge == "Edge"))
}

print(paste0(
  "Total number of clusters detected: ",
  nrow(EdgeForExcelClusterImport %>% filter(Edge == "Edge"))
))

EdgeClusterRangeFileWriter(
  IndexedRangeFilePath)


EndTime <- Sys.time()
TotalTime <- EndTime - StartTime
print(TotalTime)
