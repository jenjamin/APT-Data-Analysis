#### Create class for APTBranchesDict####
APTBranchesDict <- setRefClass(
  "APTBranchesDict",
  fields = list(dict = "list"),
  
  methods = list(
    initialize = function() {
      # Initialize the 'dict' list with all the properties
      dict <<- list()
      
      # Set the keywords
      dict$keyword <<- c('tof', 'pulse', 'freq', 'tElapsed', 'erate', 'tstage', 'TargetErate', 
                         'TargetFlux', 'pulseDelta', 'Pres', 'VAnodeMon', 'Temp', 'AmbTemp', 'FractureGuard', 
                         'Vref', 'Noise', 'Uniformity', 'xstage', 'ystage', 'zstage', 'z', 'tofc', 'Mass', 'tofb', 
                         'xs', 'ys', 'zs', 'rTip', 'zApex', 'zSphereCorr', 'XDet_mm', 'YDet_mm', 'Multiplicity', 
                         'Vap', 'Detector Coordinates', 'Position')
      
      # Index of sections
      dict$kwnsect <<- list(1:length(dict$keyword))
      
      # Initialize the other properties with corresponding values
      dict$iHeaderSize <<- rep(148, length(dict$keyword))
      dict$iHeaderVersion <<- rep(2, length(dict$keyword))
      dict$iSectionVersion <<- rep(1, length(dict$keyword))
      dict$eRelationshipType <<- rep(1, length(dict$keyword))
      dict$eRecordType <<- rep(1, length(dict$keyword))
      dict$eRecordDataType <<- c(3, 3, 3, 3, 3, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 3, 1, 1, 1, 
                                 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3)
      dict$iDatatypeSize <<- c(32, 32, 32, 32, 32, 16, 32, 32, 16, 32, 
                               32, 32, 32, 16, 32, 32, 32, 32, 32, 32, 
                               64, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
                               32, 32, 32, 32, 32, 32)
      dict$iRecordSize <<- c(4, 4, 4, 4, 4, 2, 4, 4, 2, 4, 
                             4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 
                             8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                             4, 4, 4, 4, 8, 12)
      
      # Calculate iElements based on iRecordSize and iDatatypeSize
      dict$iElements <<- sapply(1:length(dict$keyword), function(i) dict$iRecordSize[i] / (dict$iDatatypeSize[i] / 8))
    }
  )
)
