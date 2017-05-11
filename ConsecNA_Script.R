#Takes all files in the directory (all cleaned files) and determines whether there is X amount of consecutive seconds of missing data
setwd("~/Desktop/CLEANED_DATA/cleandataholder")
rm(list=ls()) #Clear Enivronement

files <- list.files("~/Desktop/CLEANED_DATA/cleandataholder")

for (file in files) {  
  df <- read.csv((paste("~/Desktop/CLEANED_DATA/cleandataholder/",file, sep = "")), header=FALSE, na.string= ".") #


navalue <- c("-9999")

for (i in navalue) df[[i]] <- grepl(i, df$V3)

if (with(rle(df$"-9999"), sum(lengths[values] >= 4800)) >=1) { # change value of 9600 to find appropriate length missing
  print(paste(file,"true"))
}
}


#else {
#print(paste(file,"false"))



