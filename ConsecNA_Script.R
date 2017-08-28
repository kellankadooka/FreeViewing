#Takes all files in the directory (all cleaned files) and determines whether there is X amount of consecutive seconds of missing data
setwd("~/Desktop/FreeViewing/data")
rm(list=ls()) #Clear Enivronement

data <- list.files("~/Desktop/FreeViewing/data")

for (file in files) {  
  df <- read.csv((paste("~/Desktop/CLEANED_DATA/cleandataholder/",file, sep = "")), header=FALSE, na.string= ".") #


navalue <- c("-9999")

# to find consecutive missing data
for (i in navalue) df[[i]] <- grepl(i, df$V3)

if (with(rle(df$"-9999"), sum(lengths[values] >= 4800)) >=1) { # change value of 9600 to find appropriate length missing
  print(paste(file,"true"))
}
}


for (file in data) {  
  df <- read.csv((paste("~/Desktop/FreeViewing/data/",file, sep = "")), header=FALSE, na.string= "-9999") #
  vector  <- colMeans(is.na(df))
  if (vector[3] > .50) {
    print(paste(file,"over .5 missing"))
    } else {
      print(paste(file,"looks good"))
  }
  }
  
  navalue <- c("-9999")

colMeans(is.na(df))

#else {
#print(paste(file,"false"))

df <- read.csv("1_006.csv", header= FALSE, na.string = "-9999")
if (colMeans(is.na(df)) > .01) {
  print(paste("06 ah shit"))
}
