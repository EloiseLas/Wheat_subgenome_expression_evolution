library(ade4)
library(adegraphics)
library(ggplot2)
library(umap)
library(plyr)
library(dplyr)

## Load the data before and after ComBat:
expression_data<-t(read.table("Data/expression_data_TPM_UQ_normalized.csv", sep = ",", header=TRUE, row.names = 1))
metadata_batch<-read.table("Data/conditions_batches.csv", sep = ",", header=TRUE, row.names = 1)
metadata<-read.table("Data/conditions_batches_all.csv", sep = ",", header=TRUE, row.names = 1)
tissue_data<-read.table("Data/tissues_info.csv", sep = ",", header=TRUE, row.names = 1)
expression_data_batch<-t(read.table("Data/expression_data_TPM_UQ_batch_filtered.csv", sep = ",", header=TRUE, row.names = 1))

merged<-merge(expression_data_batch,tissue_data, by.x=0, by.y="Condition", all.x=TRUE)
tissue_ordered_batch<-data.frame(cbind(merged$Row.names,merged$Tissue))
colnames(tissue_ordered_batch)<-c("condition","tissue")
tissue_ordered_batch <- tissue_ordered_batch %>% mutate_all(na_if,"")

merged<-merge(expression_data,tissue_data, by.x=0, by.y="Condition", all.x=TRUE)
tissue_ordered<-data.frame(cbind(merged$Row.names,merged$Tissue))
colnames(tissue_ordered)<-c("condition","tissue")
tissue_ordered <- tissue_ordered %>% mutate_all(na_if,"")

## PCA on data before ComBat:

# ade4:
acp <- dudi.pca(expression_data, center=TRUE, scale=TRUE, scann = TRUE)
summary(acp)

plot(acp$l1)

a1<-ggplot(acp$li[,1:2], aes(Axis1, Axis2, colour = metadata$batches)) + 
  geom_point(size=3) + 
  ggtitle("PCA") +
  xlab("PC1 (21.975%)") + ylab("PC2 (8.187%)") +
  labs(colour= "Batch") + theme(text = element_text(size=23))

expression_data_batch <- data.frame(scale(expression_data_batch,acp$cent,acp$norm))

acp_batch_l1<-as.matrix(expression_data_batch)%*%as.matrix(acp$c1)

ggplot(acp_batch_l1[,1:2], aes(CS1, CS2, colour = metadata_batch$batches)) + 
  geom_point() + 
  ggtitle("PCA") +
  xlab("PC1 (21.975%)") + ylab("PC2 (8.187%)") +
  labs(colour= "Batch")


# prcomp:
pca <- prcomp(expression_data)

summary(pca)$importance[,1:2]

ggplot(pca$x[,1:2], aes(PC1, PC2, colour = metadata$batches)) + 
  geom_point() + 
  ggtitle("PCA") +
  xlab("PC1 (24.689%)") + ylab("PC2 (12.648%)") +
  labs(colour= "Condition") +
  theme(text = element_text(size=18))

# Project the data after ComBat in the same space:
pca_batch <- scale(expression_data_batch, pca$center, pca$scale) %*% pca$rotation 

ggplot(pca_batch[,1:2], aes(PC1, PC2, colour = metadata_batch$batches)) + 
  geom_point() + 
  ggtitle("PCA") +
  xlab("PC1 (24.689%)") + ylab("PC2 (12.648%)") +
  labs(colour= "Condition")

## PCA of the batch normalized data:
acp <- dudi.pca(expression_data_batch, center=TRUE, scale=TRUE, scann = TRUE)
summary(acp)

a2<-ggplot(acp$l1[,1:2], aes(RS1, RS2, colour = metadata_batch$batches)) + 
  geom_point(size=3) + 
  ggtitle("PCA") +
  xlab("PC1 (11.867%)") + ylab("PC2 (8.259%)") +
  labs(colour= "Batch")+
  theme(text = element_text(size=23)) + guides(colour="none")

aplot::plot_list(a1, a2, ncol=2, labels = c("Before ComBat","After Combat"),tag_size=20)

means<-colMeans(expression_data_batch)
hist(means)

# With tissue information:

ggplot(acp$l1[,1:2], aes(RS1, RS2, colour = tissue_ordered$tissue)) + 
  geom_point() + 
  ggtitle("PCA") +
  xlab("PC1 (11.867%)") + ylab("PC2 (8.259%)") +
  labs(colour= "Tissues")

## Umap on data before ComBat:

cond.umap <- umap(expression_data)
layout <- cond.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, tissue_ordered$tissue) 

ggplot(final[,1:2], aes(X1, X2, colour = tissue_ordered$tissue)) + 
  geom_point() + 
  ggtitle("UMAP") +
  xlab("First Component") + ylab("Second Component") +
  labs(colour= "Tissues")

ggplot(final[,1:2], aes(X1, X2, colour = metadata$batches)) + 
  geom_point() + 
  ggtitle("UMAP") +
  xlab("First Component") + ylab("Second Component") +
  labs(colour= "Batches")

## Umap on data after ComBat:

cond.umap_batch <- umap(expression_data_batch)
layout <- cond.umap_batch[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, tissue_ordered_batch$tissue) 

ggplot(final[,1:2], aes(X1, X2, colour = tissue_ordered_batch$tissue)) + 
  geom_point(size=3) +
  xlab("First Component") + ylab("Second Component") +
  labs(colour= "Tissues") +
  theme(text = element_text(size=23))

ggplot(final[,1:2], aes(X1, X2, colour = metadata_batch$batches)) + 
  geom_point() + 
  ggtitle("UMAP") +
  xlab("First Component") + ylab("Second Component") +
  labs(colour= "Batches")


