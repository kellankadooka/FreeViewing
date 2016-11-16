setwd("~/Desktop/CLEANED_DATA/Child/fv036") #set wd

rm(list=ls()) #Clear Enivronement

#####Files read in should be csv file with added index column, for a single participant

df <- read.csv("fv036.csv", header=TRUE, na.string= ".") #

df <- df[- grep("PLAY_SOUND_VIDEO", df$SAMPLE_MESSAGE),] #remove PLAY_SOUND_VIDEO

df  <- df %>% fill(SAMPLE_MESSAGE) #fill sample message column with value from above

df_clean  <- df[grepl("Video Frame number", df$SAMPLE_MESSAGE),] #remove rows unless it contains "Video Frame number"

df_clean$SAMPLE_MESSAGE <- gsub("Video Frame number = ", "", df_clean$SAMPLE_MESSAGE) #remove video frame number string

#df_clean <- df_clean[- grep("adult_fixationX.xvd", df_clean$VIDEO_NAME),] #remove rows that are for distracter video

####CHILD VERSION
#########
df_clean <- df_clean[- grep("rattle.avi", df_clean$VIDEO_NAME),] #removes distacter video rows
df_clean <- df_clean[- grep("ballX.xvd", df_clean$VIDEO_NAME),]
df_clean <- df_clean[- grep("kitten.avi", df_clean$VIDEO_NAME),]
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
df_clean[is.na(df_clean)] <- -9999 #replace NA coordinate with -9999

listofdf <- split(df_clean, df_clean$identifier) #split by video identifier


for (i in seq_along(listofdf)) {
  filename = paste(i,".csv")
  write.table(listofdf[[i]], filename, sep = ",", col.names = FALSE)
}   ###write file for each video of the participant



#####
#reintroduce files
# each section reads the file, subtraccts the timestamp from the first timestamp value, removes the index then writes csv file

df_eachvideo <- read.csv("1 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "1_036.csv", sep = ",", col.names = FALSE, ) #Filename in the form "Video number_ participant ID.csv" ex 1_001.csv

df_eachvideo <- read.csv("2 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "2_036.csv", sep = ",", col.names = FALSE, )

df_eachvideo <- read.csv("3 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "3_036.csv", sep = ",", col.names = FALSE, )

df_eachvideo <- read.csv("4 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "4_036.csv", sep = ",", col.names = FALSE, )

df_eachvideo <- read.csv("5 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "5_036.csv", sep = ",", col.names = FALSE, )

df_eachvideo <- read.csv("6 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "6_036.csv", sep = ",", col.names = FALSE, )

df_eachvideo <- read.csv("7 .csv", header=FALSE)
value <- df_eachvideo$V6[1]
df_eachvideo$V6 <- df_eachvideo$V6 - value
df_eachvideo$V1 <- NULL
write.table( df_eachvideo, file = "7_036.csv", sep = ",", col.names = FALSE, )
#





