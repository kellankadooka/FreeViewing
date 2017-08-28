#this script produces a csv containing the MEAN or Median of the locations for every frame
rm(list=ls())

setwd("~/Desktop/Data Sandbox/Variability")  #set wd 
fileNames <- Sys.glob("*.csv")

for (fileName in fileNames) {  #for loop
  
sheet <- read.csv(fileName, header=FALSE) #read each file
  
sheet$V3[sheet$V3 <=0] <-  NA  #remove non screen loactions
sheet$V3[sheet$V3 >= 1280] = NA
sheet$V4[sheet$V4 <= 0] = NA
sheet$V4[sheet$V4 >= 1024] = NA

notnasum <- function(y) {    #function to figure out how many data points are available for each frame
    sum(!is.na(y))
}  
        
refsheet <- sheet %>%    #create reference sheet indacting number of available data points for each from 
  group_by(V5) %>% 
  summarise_each(funs(notnasum))

sheet <-left_join(sheet,refsheet, by = 'V5')  #join ref sheet and orginal

sheet$V3.x[sheet$V3.y <= 5] <- NA  #remove if less than 5 points
sheet$V4.x[sheet$V4.y <= 5] <- NA

meanwithoutNA <- function(z) {    #function for mean of locations
   median(z, na.rm=TRUE)
}

sheetnew <- aggregate(sheet[, 3:4], list(sheet$V5), FUN = "meanwithoutNA")  #make new sheet based on values by frame

emptydf  <- data.frame(matrix(ncol = 1, nrow = 3600))  #create empty df that is length of video
colnames(emptydf)  <- 'Group.1'
emptydf$Group.1<-seq.int(nrow(emptydf))

sheetnew  <- left_join(emptydf,sheetnew, by = "Group.1", copy = TRUE)  #join sheet placing NA for frames without info 

sheetnew$Group.1 <- NULL   #remove index
sheetnew[is.na(sheetnew)] <- -999   #net NA to -999
write.table(sheetnew,    #write files
            fileName, 
            sep = ",",
            col.names = FALSE)
}
