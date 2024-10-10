#Function for reading .apt file to dataframe
readaptR <- function(posFileName) {

  # Sourcing necessary scripts
  source("Scripts For Grain Boundary Analysis/ScriptsForReadAPT/APTBranchesDict.R")
  source("Scripts For Grain Boundary Analysis/ScriptsForReadAPT/APTFileHeader2.R")
  source("Scripts For Grain Boundary Analysis/ScriptsForReadAPT/APTSectionHeaderAuto.R")
  source("Scripts For Grain Boundary Analysis/ScriptsForReadAPT/PARAPROBE_Transcoder2.R")
  
  # Save information into list called apt
  apt <- PARAPROBE_Transcoder2(posFileName)
  
  posFile <- data.frame(
    x = c(apt$Position[1,]),
    y = c(apt$Position[2,]),
    z = c(apt$Position[3,]),
    m = c(apt$Mass[1,])
  )
  return(posFile)
}