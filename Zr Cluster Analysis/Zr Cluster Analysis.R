read.pos <- function(posFileName){
  # size of floats = 4 bytes
  sizeOfFloat = 4
  # get file size
  fileInfo<-file.info(posFileName)
  fileSize<-fileInfo$size/sizeOfFloat
  to.read = file(posFileName, "rb")
  posFile<-readBin(to.read,double(),size=4,n=fileSize,endian="big")
  close(to.read)
  posFile<-t(matrix(posFile,nrow=4,dimnames=list(c("x","y","z","m"))))
  posFile<-as.data.frame(posFile)
  return(posFile)
}

a = read.pos("D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\Fe Clusters.cluster.indexed.pos")
b = read.pos("D:\\Post Doc\\MIDAS 3 Material\\X2 RXA Cycle 4\\Run 8811\\X2 RXA Cycle 4 - Run 8811\\recons\\recon-v01\\default\\Fe Clusters.cluster.pos")

