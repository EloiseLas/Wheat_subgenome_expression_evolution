library(plyr)
library(dplyr)
library(tidyverse)

## Load the data:
expression_data<-read.table("Data/expression_data_TPM_no_UQ_batch_normalized_new.csv", sep = ",", header=TRUE, row.names = 1)
expression_data<-result

## Compute stdev of the genes:
expression_data$std<-apply(expression_data, 1, sd, na.rm=TRUE)

## Filter genes with 0 variance:
std_0<-subset(expression_data, expression_data$std==0)

expression_data_filtered<-subset(expression_data, expression_data$std!=0)
expression_data_filtered<-expression_data_filtered[ ,  !names(expression_data_filtered) %in% 
            c("std")]

write.csv(expression_data_filtered,
          file = "Data/expression_data_TPM_no_UQ_batch_filtered_new.csv")
