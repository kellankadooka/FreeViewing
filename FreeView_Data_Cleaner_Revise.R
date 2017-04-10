#for each video, if consecutive length == 20 of -9999 == true (print 'video' False)

#testpage
setwd("~/Desktop/CLEANED_DATA/datacleaner") #set wd

rm(list=ls()) #Clear Enivronement 
library(plyr) #set of packages
library(dplyr)
library(tidyr)
#library(reshape2)
library(car)
library(ggplot2)
library(nlme)
library(psych)
library(stringr)

folders <- list.files("~/Desktop/CLEANED_DATA/datacleaner") #set folder of files

for (folder in folders) {  
  df <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/",folder,".csv", sep = "")), header=TRUE, na.string= ".") #read csv

  df <- df[- grep("PLAY_SOUND_VIDEO", df$SAMPLE_MESSAGE),] #remove PLAY_SOUND_VIDEO
  
  df  <- df %>% fill(SAMPLE_MESSAGE) #fill sample message column with value from above
  
  df_clean  <- df[grepl("Video Frame number", df$SAMPLE_MESSAGE),] #remove rows unless it contains "Video Frame number"
  
  df_clean$SAMPLE_MESSAGE <- gsub("Video Frame number = ", "", df_clean$SAMPLE_MESSAGE) #remove video frame number string
  
  #adult should contain adult fixation to remove, if statment removes if present
  if ("adult_fixationX.xvd" %in% df_clean$VIDEO_NAME == TRUE) {
    df_clean <- df_clean[- grep("adult_fixationX.xvd", df_clean$VIDEO_NAME),]
  }
  
  #df_clean <- df_clean[- grep("adult_fixationX.xvd", df_clean$VIDEO_NAME),] #remove rows that are for distracter video
  ####CHILD VERSION
  #########
  if ("rattle.avi" %in% df_clean$VIDEO_NAME  == TRUE){
    df_clean <- df_clean[- grep("rattle.avi", df_clean$VIDEO_NAME),] #removes distacter video rows
  }
  
  if ("ball.xvd" %in% df_clean$VIDEO_NAME == TRUE) {
    df_clean <- df_clean[- grep("ballX.xvd", df_clean$VIDEO_NAME),]
  }
  
  if ("kitten.avi" %in% df_clean$VIDEO_NAME == TRUE){
    df_clean <- df_clean[- grep("kitten.avi", df_clean$VIDEO_NAME),]
  }
  
  #df_clean <- df_clean[- grep("rattle.avi", df_clean$VIDEO_NAME),] #removes distacter video rows
  #df_clean <- df_clean[- grep("ballX.xvd", df_clean$VIDEO_NAME),]
  #df_clean <- df_clean[- grep("kitten.avi", df_clean$VIDEO_NAME),]
  #########
  df_clean <- df_clean[order(df_clean$identifier),] # order df by video identifier
  
  colnames(df_clean)[colnames(df_clean) == 'RECORDING_SESSION_LABEL'] <- 'ID' #rename row
  
  colnames(df_clean)[colnames(df_clean) == 'RIGHT_GAZE_X'] <- 'xpos' #rename row
  
  colnames(df_clean)[colnames(df_clean) == 'RIGHT_GAZE_Y'] <- 'ypos' #rename row
  
  colnames(df_clean)[colnames(df_clean) == 'SAMPLE_MESSAGE'] <- 'vf' #rename row
  
  df_clean$X <- NULL #remove index column
  df_clean$VIDEO_NAME <- NULL #remove video name column
  df_clean$ID <- gsub("fv0", "", df_clean$ID) #replace participant number with integer by replacing fv0
  df_clean$ID. <- NULL # remove ID column (will be replaced with )
  is.na(df_clean$xpos) <- df_clean$xpos < 0
  is.na(df_clean$xpos) <- df_clean$xpos > 1280
  is.na(df_clean$ypos) <- df_clean$ypos < 0
  is.na(df_clean$ypos) <- df_clean$ypos > 1024
  df_clean[is.na(df_clean)] <- -9999 #replace NA coordinate with -9999
  
  listofdf <- split(df_clean, df_clean$identifier) #split by video identifier
  
  
  for (i in seq_along(listofdf)) {
    filename = paste(i,".csv")
    write.table(listofdf[[i]], file = (paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/",filename, sep = "")), sep = ",", col.names = FALSE)
  } #write tables for each video of each participant


  
  #####
  #reintroduce files
  # each section reads the file, subtraccts the timestamp from the first timestamp value, removes the index then writes csv file
  parnum <- str_sub(folder, start = -3)
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/1 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file = (paste("1_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, ) #Filename in the form "Video number_ participant ID.csv" ex 1_001.csv
  
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/2 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file =  (paste("2_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, )
  
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/3 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file =  (paste("3_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, )
  
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/4 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file =  (paste("4_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, )
  
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/5 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file =  (paste("5_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, )
  
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/6 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file =  (paste("6_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, )
  
  df_eachvideo <- read.csv((paste("~/Desktop/CLEANED_DATA/datacleaner/",folder,"/7 .csv", sep = "")), header=FALSE)
  value <- df_eachvideo$V6[1]
  df_eachvideo$V6 <- df_eachvideo$V6 - value
  df_eachvideo$V1 <- NULL
  write.table( df_eachvideo, file =  (paste("7_",parnum,".csv", sep ='')), sep = ",", col.names = FALSE, )
  #
  
}



