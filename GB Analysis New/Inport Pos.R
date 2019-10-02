#### Inputs .pos file of GB region that has been determined by the analysis ####

setwd(dirname(parent.frame(2)$ofile))

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

PosFile <- readposR("GB Containing ROI.pos")

#### Input RangeFile and Tidy ####
RangeFile <- read.delim("Test Range.rrng", colClasses = "character", header = FALSE)
RowsToSkip <- as.numeric(gsub("Number=", "", RangeFile[2,])) + 5
Ranges <- RangeFile %>% slice(RowsToSkip:n())
rm(RowsToSkip, RangeFile)

RangesDF <- data.frame()

i = 0
for(i in unique(str_count(Ranges$V1, ":") - 2)){
  
  Elements <- c()
  for(j in seq(1,i,1)){
    Elements <- append(Elements,paste("Element",j))
  }
  
  ColumnNames <- c("Start", "End", "Volume",
                   Elements, "Color")
  
  Ranges %>%
    mutate(NumberIons = str_count(V1, ":") - 2) %>%
    filter(NumberIons == i) %>%
    separate(V1,
             ColumnNames,
             sep = " ")
  
  RangesDF <- bind_rows(RangesDF, 
                    Ranges %>%
                      mutate(NumberIons = str_count(V1, ":") - 2) %>%
                      filter(NumberIons == i) %>%
                      separate(V1,
                               ColumnNames,
                               sep = " ")
                    )

  }
   
rm(ColumnNames, Elements,i, j) 

#### Creat R-friendly range file ####

RangesDF2 <- cbind(
  RangesDF %>%
    mutate(Start = as.numeric(str_extract(Start,"[^=]+$")),
           End = as.numeric(str_extract(End, "[^=]+$")),
           Volume = gsub("Vol:", "", Volume),
           Color = gsub("Color:", "", Color)) %>%
  select(-NumberIons, -`Element 1`, -`Element 2`),
  RangesDF %>%
    select_if((grepl("Element",names(.)))) %>%
    unite("Ion") %>%
    mutate(Ion = paste(gsub("1|Name|:|NA|_| ","",Ion)))
)

rm(Ranges, RangesDF)

#### Plot Mass Spec ####
ggplot(PosFile %>%
         filter(m < 30 & m > 25)) +
  geom_histogram(aes(m), binwidth = 0.01) + 
  scale_y_continuous(trans = log10_trans()) 

PosFile %>%
  mutate(Ion = if_else(, ,))

