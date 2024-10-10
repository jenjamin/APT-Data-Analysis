#PARAPROBE_Transcoder2
setRefClass("Data",
            fields = c(
              "healthy",
              "aptfn",
              "aptbranches",
              "header",
              "idtfyd_sections",
              "tipbox",
              "tof",
              "pulse",
              "freq",
              "tElapsed",
              "erate",
              "tstage",
              "TargetErate",
              "TargetFlux",
              "pulseDelta",
              "Pres",
              "VAnodeMon",
              "Temp",
              "AmbTemp",
              "FractureGuard",
              "Vref",
              "Noise",
              "Uniformity",
              "xstage",
              "ystage",
              "zstage",
              "z",
              "tofc",
              "Mass",
              "tofb",
              "xs",
              "ys",
              "zs",
              "rTip",
              "zApex",
              "zSphereCorr",
              "XDet_mm",
              "YDet_mm",
              "Multiplicity",
              "Vap",
              "DetectorCoordinates",
              "Position",
              "known_sections"
            ),
            methods = 
              PARAPROBE_Transcoder2 <- function(fn) {
                # Constructor for the class object/instance
                obj <- list()
                
                if (!missing(fn)) {
                  obj$healthy <- TRUE
                  obj$aptfn <- fn
                  
                  # Define all *.APT branches possible (assuming you have a similar function in R)
                  obj$aptbranches <- APTBranchesDict()
                  
                  # Open the file for read
                  fid <- file(obj$aptfn, "rb")
                  
                  # Read the header of the file
                  obj$header <- APTFileHeader2(fid)
                  print(obj$header)
                  
                  if (obj$header$healthy == FALSE) {
                    close(fid)
                    return(obj)
                  }
                  
                  # Auto-detect all sections in the file
                  obj$idtfyd_sections <- list()
                  for (i in seq_along(obj$aptbranches$dict$keyword)) {
                    obj$idtfyd_sections[[i]] <- APTSectionHeaderAuto(fid)
                    
                    if (obj$idtfyd_sections[[i]]$healthy == FALSE) {
                      cat(paste0("Section ", i, " failed!\n"))
                      close(fid)
                      return(obj)
                    }
                    
                    # Print the section header
                    print(obj$idtfyd_sections[[i]])
                    
                    # Load specific data based on what the section encodes
                    ni <- obj$header$llIonCount
                    nj <- obj$idtfyd_sections[[i]]$iElements
                    
                    cat(paste0("ni: ", ni, " nj: ", nj, "\n"))
                    
                    # Switch cases to read specific data
                    obj <- switch(
                      obj$idtfyd_sections[[i]]$wcSectionType,
                      'Failure' = {
                        cat(paste('Autodetect the', i, 'APT header failed!\n'))
                        return()
                      },
                      'tof' = {
                        cat("tof\n")
                        obj$tof <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$tof), c(nj, ni))) stop("Reading tof failed!")
                        obj
                      },
                      'pulse' = {
                        cat("pulse\n")
                        obj$pulse <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$pulse), c(nj, ni))) stop("Reading pulse failed!")
                        obj
                      },
                      'freq' = {
                        cat("freq\n")
                        obj$freq <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$freq), c(nj, ni))) stop("Reading freq failed!")
                        obj
                      },
                      'tElapsed' = {
                        cat("tElapsed\n")
                        obj$tElapsed <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$tElapsed), c(nj, ni))) stop("Reading tElapsed failed!")
                        obj
                      },
                      'erate' = {
                        cat("erate\n")
                        obj$erate <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$erate), c(nj, ni))) stop("Reading erate failed!")
                        obj
                      },
                      'tstage' = {
                        cat("tstage\n")
                        obj$tstage <- matrix(readBin(fid, "integer", size = 2, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$tstage), c(nj, ni))) stop("Reading tstage failed!")
                        obj
                      },
                      'pulseDelta' = {
                        cat("pulseDelta\n")
                        obj$pulseDelta <- matrix(readBin(fid, "integer", size = 2, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$pulseDelta), c(nj, ni))) stop("Reading pulseDelta failed!")
                        obj
                      },
                      'Pres' = {
                        cat("Pres\n")
                        obj$Pres <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$Pres), c(nj, ni))) stop("Reading Pres failed!")
                        obj
                      },
                      'VAnodeMon' = {
                        cat("VAnodeMon\n")
                        obj$VAnodeMon <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$VAnodeMon), c(nj, ni))) stop("Reading VAnodeMon failed!")
                        obj
                      },
                      'Temp' = {
                        cat("Temp\n")
                        obj$Temp <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$Temp), c(nj, ni))) stop("Reading Temp failed!")
                        obj
                      },
                      'AmbTemp' = {
                        cat("AmbTemp\n")
                        obj$AmbTemp <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$AmbTemp), c(nj, ni))) stop("Reading AmbTemp failed!")
                        obj
                      },
                      'FractureGuard' = {
                        cat("FractureGuard\n")
                        obj$FractureGuard <- matrix(readBin(fid, "integer", size = 2, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$FractureGuard), c(nj, ni))) stop("Reading FractureGuard failed!")
                        obj
                      },
                      'Vref' = {
                        cat("Vref\n")
                        obj$Vref <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$Vref), c(nj, ni))) stop("Reading Vref failed!")
                        obj
                      },
                      'Noise' = {
                        cat("Noise\n")
                        obj$Noise <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$Noise), c(nj, ni))) stop("Reading Noise failed!")
                        obj
                      },
                      'Uniformity' = {
                        cat("Uniformity\n")
                        obj$Uniformity <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$Uniformity), c(nj, ni))) stop("Reading Uniformity failed!")
                        obj
                      },
                      'xstage' = {
                        cat("xstage\n")
                        obj$xstage <- matrix(readBin(fid, "integer", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$xstage), c(nj, ni))) stop("Reading xstage failed!")
                        obj
                      },
                      'ystage' = {
                        cat("ystage\n")
                        obj$ystage <- matrix(readBin(fid, "integer", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$ystage), c(nj, ni))) stop("Reading ystage failed!")
                        obj
                      },
                      'zstage' = {
                        cat("zstage\n")
                        obj$zstage <- matrix(readBin(fid, "integer", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$zstage), c(nj, ni))) stop("Reading zstage failed!")
                        obj
                      },
                      'z' = {
                        cat("z\n")
                        obj$z <- matrix(readBin(fid, "integer64", size = 8, n = ni * nj, endian = "little"), ncol = ni)
                        if (!identical(dim(obj$z), c(nj, ni))) stop("Reading z failed!")
                        obj
                      },
                      'Mass' = {
                        cat("Mass\n")
                        obj$Mass <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$Mass), c(nj, ni))) stop("Reading Mass failed!")
                        obj
                      },
                      'tofb' = {
                        cat("tofb\n")
                        obj$tofb <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$tofb), c(nj, ni))) stop("Reading tofb failed!")
                        obj
                      },
                      'xs' = {
                        cat("xs\n")
                        obj$xs <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$xs), c(nj, ni))) stop("Reading xs failed!")
                        obj
                      },
                      'ys' = {
                        cat("ys\n")
                        obj$ys <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$ys), c(nj, ni))) stop("Reading ys failed!")
                        obj
                      },
                      'zs' = {
                        cat("zs\n")
                        obj$zs <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$zs), c(nj, ni))) stop("Reading zs failed!")
                        obj
                      },
                      'rTip' = {
                        cat("rTip\n")
                        obj$rTip <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$rTip), c(nj, ni))) stop("Reading rTip failed!")
                        obj
                      },
                      'zApex' = {
                        cat("zApex\n")
                        obj$zApex <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$zApex), c(nj, ni))) stop("Reading zApex failed!")
                        obj
                      },
                      'zSphereCorr' = {
                        cat("zSphereCorr\n")
                        obj$zSphereCorr <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), nrow = nj, ncol = ni, byrow = TRUE)
                        if (!all.equal(dim(obj$zSphereCorr), c(nj, ni))) stop("Reading zSphereCorr failed!")
                        obj
                      },
                      'XDet_mm' = {
                        cat("XDet_mm\n")
                        obj$XDet_mm <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$XDet_mm), c(nj, ni))) stop("Reading XDet_mm failed!")
                        obj
                      },
                      'YDet_mm' = {
                        cat("YDet_mm\n")
                        obj$YDet_mm <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$YDet_mm), c(nj, ni))) stop("Reading YDet_mm failed!")
                        obj
                      },
                      'Multiplicity' = {
                        cat("Multiplicity\n")
                        obj$Multiplicity <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), nrow = nj, ncol = ni, byrow = TRUE)
                        if (!all.equal(dim(obj$Multiplicity), c(nj, ni))) stop("Reading Multiplicity failed!")
                        obj
                      },
                      'Vap' = {
                        cat("Vap\n")
                        obj$Vap <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), ncol = ni)
                        if (!all.equal(dim(obj$Vap), c(nj, ni))) stop("Reading Vap failed!")
                        obj
                      },
                      'DetectorCoordinates' = {
                        cat("DetectorCoordinates\n")
                        obj$DetectorCoordinates <- matrix(readBin(fid, "numeric", size = 4, n = ni * nj, endian = "little"), nrow = nj, ncol = ni, byrow = TRUE)
                        if (!all.equal(dim(obj$DetectorCoordinates), c(nj, ni))) stop("Reading DetectorCoordinates failed!")
                        obj
                      },
                      'Position' = {
                        cat('Position\n')
                        # first read the preceding header with bounds
                        tipbox <- matrix(readBin(fid, what = "numeric", n = 6, size = 4, endian = "little"), nrow = 2, ncol = 3)
                        if (!all(dim(tipbox) == c(2, 3))) {
                          cat('Reading Position box bounds failed!\n')
                          break
                        }
                        Position <- matrix(readBin(fid, what = "numeric", n = ni * nj, size = 4, endian = "little"), nrow = nj, ncol = ni, byrow = TRUE)
                        if (!all(dim(Position) == c(nj, ni))) {
                          cat('Reading Position failed!\n')
                          break
                        }
                        obj$tipbox <- tipbox
                        obj$Position <- Position
                        obj
                      },
                      {
                        cat("Stop reading sections\n")
                        break
                      }
                    )
                  }
                  cat("Done, closing the *.APT file\n")
                  close(fid)
                }
                
                return(obj)
              }
            
            # Dummy functions for APTBranchesDict, APTFileHeader2, APTSectionHeaderAuto
            # These need to be implemented or replaced with equivalent logic
            # APTBranchesDict <- function() {
            #   list(dict = list(keyword = c('tof', 'pulse', 'freq', 'Mass')))
            # }
            # 
            # APTFileHeader2 <- function(fid) {
            #   header <- list(
            #     healthy = TRUE, 
            #     llIonCount = 100  # dummy value
            #   )
            #   return(header)
            # }
            # 
            # APTSectionHeaderAuto <- function(fid) {
            #   section <- list(
            #     healthy = TRUE,
            #     iElements = 2,  # dummy value
            #     wcSectionType = "tof"  # example section type
            #   )
            #   return(section)
            # }
)