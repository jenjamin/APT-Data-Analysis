#### Create class for APTSectionHeaderAuto####
# APTSectionHeaderAuto class definition in R
APTSectionHeaderAuto <- setRefClass("APTSectionHeaderAuto",
                                    fields = list(
                                      healthy = "logical",
                                      cSignature = "character",
                                      iHeaderSize = "integer",
                                      iHeaderVersion = "integer",
                                      wcSectionType = "character",
                                      iSectionVersion = "integer",
                                      eRelationshipType = "integer",
                                      eRecordType = "integer",
                                      eRecordDataType = "integer",
                                      iDataTypeSize = "integer",
                                      iRecordSize = "integer",
                                      wcDataUnit = "character",
                                      llRecordCount = "numeric",
                                      llByteCount = "numeric",
                                      iElements = "numeric"
                                    ),
                                    
                                    methods = list(
                                      initialize = function(fid) {
                                        # Constructor
                                        healthy <<- TRUE
                                        cSignature <<- ""
                                        iHeaderSize <<- integer(0)
                                        iHeaderVersion <<- integer(0)
                                        wcSectionType <<- "Failure"
                                        iSectionVersion <<- integer(0)
                                        eRelationshipType <<- integer(0)
                                        eRecordType <<- integer(0)
                                        eRecordDataType <<- integer(0)
                                        iDataTypeSize <<- integer(0)
                                        iRecordSize <<- integer(0)
                                        wcDataUnit <<- character(0)
                                        llRecordCount <<- 0
                                        llByteCount <<- 0
                                        
                                        # Read from file
                                        cSignature <<- rawToChar(readBin(fid, "raw", n = 4))
                                        if (cSignature != "SEC") {
                                          healthy <<- FALSE
                                          cat("Section is not a valid SEC!\n")
                                          wcSectionType <<- "Failure"
                                          return()
                                        }
                                        
                                        iHeaderSize <<- readBin(fid, "integer", n = 1)
                                        iHeaderVersion <<- readBin(fid, "integer", n = 1)
                                        if (iHeaderVersion != 2) {
                                          healthy <<- FALSE
                                          cat("Section header has a different version!\n")
                                          wcSectionType <<- "Failure"
                                          return()
                                        }
                                        
                                        # Read section type
                                        tmp <- paste(rawToChar(readBin(fid, what = "raw", n = 64), multiple = TRUE),collapse = "")
                                        if (sum(utf8ToInt(tmp) > 255) > 0) {
                                          healthy <<- FALSE
                                          cat("UTF16 conversion of wcSectionType did not work!\n")
                                          wcSectionType <<- "Failure"
                                          return()
                                        }
                                        
                                        wcSectionType <<- tmp
                                        iSectionVersion <<- readBin(fid, "integer", n = 1)
                                        eRelationshipType <<- readBin(fid, "integer", n = 1)
                                        if (eRelationshipType != 1) {
                                          healthy <<- FALSE
                                          cat("eRelationshipType != ONE_TO_ONE but that is the only currently supported!\n")
                                          wcSectionType <<- "Failure"
                                          return()
                                        }
                                        
                                        eRecordType <<- readBin(fid, "integer", n = 1)
                                        if (eRecordType != 1) {
                                          healthy <<- FALSE
                                          cat("eRecordType != FIXED_SIZE but that is the only currently supported!\n")
                                          wcSectionType <<- "Failure"
                                          return()
                                        }
                                        
                                        eRecordDataType <<- readBin(fid, "integer", n = 1)
                                        iDataTypeSize <<- readBin(fid, "integer", n = 1)
                                        iRecordSize <<- readBin(fid, "integer", n = 1)
                                        wcDataUnit <<- paste(rawToChar(readBin(fid, what = "raw", n = 32), multiple = TRUE),collapse = "")
                                        llRecordCount <<- readBin(fid, "numeric", n = 1)
                                        llByteCount <<- readBin(fid, "numeric", n = 1)
                                        iElements <<- iRecordSize / (iDataTypeSize / 8)
                                        
                                        cat(sprintf("Reading *.APT section __%s__ successful\n", wcSectionType))
                                      },
                                      
                                      print_info = function() {
                                        cat(sprintf("cSignature: %s\n", cSignature))
                                        cat(sprintf("iHeaderSize: %d\n", iHeaderSize))
                                        cat(sprintf("iHeaderVersion: %d\n", iHeaderVersion))
                                        cat(sprintf("wcSectionType: %s\n", wcSectionType))
                                        cat(sprintf("iSectionVersion: %d\n", iSectionVersion))
                                        cat(sprintf("eRelationshipType: %d\n", eRelationshipType))
                                        cat(sprintf("eRecordType: %d\n", eRecordType))
                                        cat(sprintf("eRecordDataType: %d\n", eRecordDataType))
                                        cat(sprintf("iDataTypeSize: %d\n", iDataTypeSize))
                                        cat(sprintf("iRecordSize: %d\n", iRecordSize))
                                        cat(sprintf("wcDataUnit: %s\n", wcDataUnit))
                                        cat(sprintf("llRecordCount: %f\n", llRecordCount))
                                        cat(sprintf("llByteCount: %f\n", llByteCount))
                                        cat(sprintf("iElements: %f\n", iElements))
                                      }
                                    )
)
