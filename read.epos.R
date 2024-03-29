library(tidyverse)
read.eposint <- function(posFileName){
  # size of floats = 4 bytes
  sizeOfFloat = 4
  # get file size
  fileInfo<-file.info(posFileName)
  fileSize<-fileInfo$size/sizeOfFloat
  to.read = file(posFileName, "rb")
  eposFile<-readBin(to.read,integer(),size=4,n=fileSize,endian="big")
  close(to.read)
  eposFile<-t(matrix(eposFile,nrow=11,dimnames=list(c("x","y","z","m","tof","vdc","vp","xd","yd","dP","nm"))))
  eposFile<-as.data.frame(eposFile)
  return(eposFile)
}
read.eposdbl <- function(posFileName){
  # size of floats = 4 bytes
  sizeOfFloat = 4
  # get file size
  fileInfo<-file.info(posFileName)
  fileSize<-fileInfo$size/sizeOfFloat
  to.read = file(posFileName, "rb")
  eposFile<-readBin(to.read,double(),size=4,n=fileSize,endian="big")
  close(to.read)
  eposFile<-t(matrix(eposFile,nrow=11,dimnames=list(c("x","y","z","m","tof","vdc","vp","xd","yd","dP","nm"))))
  eposFile<-as.data.frame(eposFile)
  return(eposFile)
}

epos <- cbind(read.eposdbl("D:/WorkForOthers/Zhao T91/T91 Oxide - Run 8884/T91 Oxide - Run 8884/recons/recon-v01/default/R5083_08884-v01.epos") %>% 
                select(-dP, -nm),
              read.eposint("D:/WorkForOthers/Zhao T91/T91 Oxide - Run 8884/T91 Oxide - Run 8884/recons/recon-v01/default/R5083_08884-v01.epos") %>% 
                select(dP, nm))

CrRich <- epos %>%
  filter(-0.48-5 < x & x < -0.48+5 &
           7.75-5 < y & y < 7.75+5 &
           82.32-5 < z & z < 82.32+5)

FeRich <- epos %>%
  filter(2.55-5 < x & x < 2.55+5 &
           7.75-5 < y & y < 7.75+5 &
           317.82-5 < z & z < 317.82+5)

rm(epos)
gc()

# Number of multiples = 2 * number of nm = 0
TotalMultiples <- as.data.frame(table(epos$nm))
2*(TotalMultiples %>%
  filter(Var1 == 0) %>%
    select(Freq))

CrMultiples <- as.data.frame(table(CrRich$nm))
2*(CrMultiples %>%
     filter(Var1 == 0) %>%
     select(Freq))


FeMultiples <- as.data.frame(table(FeRich$nm))
2*(FeMultiples %>%
     filter(Var1 == 0) %>%
     select(Freq))


#### Plot of multiples vs depth (z) ####

MultipleDF <- as.data.frame(matrix(ncol = 3, nrow = 0))
i = 0
while (i < max(epos$z)) {
  DF <- epos %>% filter(z > i & z < i + 10)
  TotalMultiples <- as.data.frame(table(DF$nm))
  TotalNumberEvents = 2*(TotalMultiples %>%
                         filter(Var1 == 0) %>%
                         select(Freq)) +
  (TotalMultiples %>%
       filter(Var1 == 1) %>%
       select(Freq))
  TotalNumberMultiples <- 2*(TotalMultiples %>%
                             filter(Var1 == 0) %>%
                             select(Freq))
  Depth <- mean(DF$z)
  DepthVsMutiples <- tibble(Depth, as.numeric(TotalNumberEvents), as.numeric(TotalNumberMultiples))
  MultipleDF <<- rbind(MultipleDF, DepthVsMutiples)
  print(i)
  i <- i + 10
}
colnames(MultipleDF) <- c("Depth", "Events", "Multiple Events")
gc()
MultipleDF <- MultipleDF %>% mutate('Multiple Rate' = 100 * `Multiple Events`/Events)

ggplot(MultipleDF,
       aes(Depth, `Multiple Rate`)) +
  geom_point() +
  theme_classic(16) +
  labs(x = "Depth (nm)",
       y = "Multiple Rate (%)")
