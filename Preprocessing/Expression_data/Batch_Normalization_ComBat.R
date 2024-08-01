library(dplyr)
library(plyr)
library(tidyr)
library(sva)
library(gt)
library(gtExtras)

## Load the data:
expression_data_UQ_normalized<-read.table("Data/expression_data_TPM_UQ_normalized.csv", sep = ",", header=TRUE, row.names = 1)
# expression_data_UQ_normalized<-dat.log10
expression_data_UQ_normalized<-t(expression_data_UQ_normalized)

metadata<-read.table("Data/run-study-attrib-biosample.txt", sep = "\t", header=FALSE, skip = 1, fill = TRUE)
metadata<-select(metadata, V1, V2)

## Merge the data:
all <- merge(expression_data_UQ_normalized, metadata, by.x = "row.names", by.y = "V1", all.x=TRUE, all.y=FALSE)  # merge by row names (by=0 or by="row.names")
rownames(all) <- all$Row.names
all <- subset(all, select = -c(Row.names))

## Save a table with the conditions and their batches:
conditions_batches_all<-bind_cols(list(conditions=rownames(all),batches=all$V2))
conditions_batches_all[conditions_batches_all=="SRP035348"] <- NA
write.csv(conditions_batches_all,
          file = "Data/conditions_batches_all.csv")

## Filter the data:
count(all, "V2")

selected_data<-all[!(all$V2 %in% "SRP035348"),]
selected_data<-selected_data %>% drop_na(V2)

unselected_data<-all[all$V2 %in% "SRP035348",]
unselected_data<-t(subset(unselected_data, select = -c(V2)))
na_data<-all[which(is.na(all$V2)),]
na_data<-t(subset(na_data, select = -c(V2)))

count(selected_data, "V2")

batch <- selected_data[['V2']]

expression_data_selected<-subset(selected_data, select = -c(V2))

## Save a table with the selected conditions and their batches:
conditions_batches<-bind_cols(list(conditions=rownames(selected_data),batches=selected_data$V2))
write.csv(conditions_batches,
          file = "Data/conditions_batches.csv")

expression_data_selected<-t(expression_data_selected)

## Apply ComBat:
corrected.data = ComBat(dat=expression_data_selected, batch=batch)

result <- round(corrected.data, digits=6)
write.csv(result,
          file = "Data/expression_data_TPM_UQ_batch_normalized.csv")

## Plots:

result<-read.table("Data/expression_data_TPM_UQ_batch_normalized_new.csv", sep = ",", header=TRUE, row.names = 1)
result<-read.table("Data/expression_data_TPM_UQ_batch_normalized.csv", sep = ",", header=TRUE, row.names = 1)

x <- melt(result[,1:20])

plt <- ggplot(data = x, aes(x = variable, y = value))
plt + geom_boxplot() + theme_minimal() + labs(x = "Title", y = "x")

## About the studies and tissue:
rownames(metadata_all) <- metadata_all$V1
studies_all<-metadata_all[rownames(expression_data_UQ_normalized),]
studies<-studies_all[!duplicated(studies_all[,2]), ]
write.csv(studies,
          file = "Data/studies_info.csv")

studies_clean<-studies[,2:3] %>% drop_na()
studies_clean$tissue<-c("NA", "Lodicule","Spikelet","Leaf","Leaf","Caryopsis","Caryopsis","NA","Grain","Leaf","Seed","Upper Leaf", "Leaf","NA","Leaf","Spike and rachis","Leaf","Meiocytes","NA","Leaf","Leaf","Leaf","Grain (immature)","Leaf","Leaf","Endosperm","Spike","Spike","Leaf","Crown","Stem","Spikelet","Leaf","Coleoptile","Lemma","Grain","Head","Leaf","Grain")

studies_clean %>%
  gt() %>%
  cols_label(
    V2 = 'Study accession',
    V3 = 'Study Title',
    tissue = 'Tissue'
  ) %>%
  gt_theme_538() %>% 
  tab_style(
    style = list(
      cell_fill(color = "lightblue"),
      "font-variant: small-caps;"
    ),
    locations = cells_body(columns = V2)
  )

tb %>% gtsave(filename = "Output/Table_studies.png")

condition_batch<-studies_all[,1:2]
batch_tissue<-data.frame(cbind(studies_clean$V2,studies_clean$tissue))
condition_batch_tissue<-merge(condition_batch,batch_tissue,by.x = "V2",by.y = "X1", all.x = TRUE)
colnames(condition_batch_tissue)<-c("Batch","Condition","Tissue")
write.csv(condition_batch_tissue,
          file = "Data/tissues_info.csv")                    
