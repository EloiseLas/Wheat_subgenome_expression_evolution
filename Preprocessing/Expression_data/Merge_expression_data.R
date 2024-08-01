library(dplyr)
library(readr)
library(data.table)
library(stringr)

## Merging all the condition in one file:
files <- list.files(path="Data/4565", full.names = TRUE)
all.df <- lapply(files, function(i){read.table(i, sep = "\t", header=FALSE, row.names = 1, skip=1, col.names=c("Geneid",str_sub(str_sub(i,11), end = -5)))})
data<-bind_cols(all.df)

write.table(data,"Data/expression_data_merged.csv",sep=",")


