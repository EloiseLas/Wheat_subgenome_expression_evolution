library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(forcats)

## Load the data:
iso_table<-read.table("Data/Genes_iso.csv", sep = ",", header=TRUE)[,2:3]

homeologs<-read.table("Data/homoeolog_Traes_EGI.csv", sep = ",", header=TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("X")]
homeologs<-merge(homeologs,iso_table, by.x="EGI_LOC_A", by.y="ID", all.x=TRUE)
names(homeologs)[names(homeologs) == 'Isoform_number'] <- 'Isoform_number_A'
homeologs<-merge(homeologs,iso_table, by.x="EGI_LOC_B", by.y="ID", all.x=TRUE)
names(homeologs)[names(homeologs) == 'Isoform_number'] <- 'Isoform_number_B'
homeologs<-merge(homeologs,iso_table, by.x="EGI_LOC_D", by.y="ID", all.x=TRUE)
names(homeologs)[names(homeologs) == 'Isoform_number'] <- 'Isoform_number_D'

## Compute isoform number std between homeologs:
homeologs <- transform(homeologs, stdev=apply(homeologs[,8:10], 1, sd, na.rm=TRUE))

homeologs$length_diff <- homeologs$stdev>200

plot(homeologs$stdev)

homeologs_clean<-homeologs %>% drop_na()

write.csv(homeologs_clean,
          file = "Output/homoeolog_isoforms_detailed.csv")

## Statistical tests:
mean(homeologs$Isoform_number_A, na.rm=TRUE)
mean(homeologs$Isoform_number_B, na.rm=TRUE)
mean(homeologs$Isoform_number_D, na.rm=TRUE)

iso<-c(homeologs_clean$Isoform_number_A,homeologs_clean$Isoform_number_B,homeologs_clean$Isoform_number_D)
homeo=as.factor(c(rep('A',length(homeologs_clean$Isoform_number_A)),rep('B',length(homeologs_clean$Isoform_number_B)),rep('D',length(homeologs_clean$Isoform_number_D))))
lm=lm(iso~homeo)
anova(lm)

kruskal.test(iso, homeo) 
pairwise.wilcox.test(iso, homeo)

## Plots:
homeologs_clean<-homeologs_clean %>% arrange(desc(stdev))
homeologs_clean$ID<-as.factor(homeologs_clean$ID)
plot(homeologs_clean$stdev, pch=21, col="lightgreen",  bg="lightgreen",cex=2,xlab="Triads ordered by their Isoform Number Standard Deviation", ylab="Isoform Number Standard Deviation",cex.lab=1.5,
     cex.axis=1.5)
ggplot(homeologs_clean) +
  aes(x = ID, y = stdev) + 
  geom_point() + 
  theme(text = element_text(size=15)) +
  xlab('Triads ID') + ylab('Isoform Number Standard Deviation in the triads') +
  scale_x_discrete(expand = c(0,0))

p1<-ggplot(homeologs_clean, aes(x=Isoform_number_A))+
  geom_histogram(color="black",fill="#FF5733", bins = 25) + xlab("Isoform number in A") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 1, 
             color = "black", size=1)

p2<-ggplot(homeologs_clean, aes(x=Isoform_number_B))+
  geom_histogram(color="black",fill="#0EA5DE",bins = 25) +
  xlab("Isoform number in B") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 1, 
             color = "black", size=1)

p3<-ggplot(homeologs_clean, aes(x=Isoform_number_D))+
  geom_histogram(color="black",fill="#6C16B8",bins = 25) +
  xlab("Isoform number in D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 1, 
             color = "black", size=1)

p1<-ggplot(homeologs_clean, aes(x=diffAB))+
  geom_histogram(color="black",fill="#FF5733", bins = 50) + xlab("Isoform number difference between A and B") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

p2<-ggplot(homeologs_clean, aes(x=diffBD))+
  geom_histogram(color="black",fill="#0EA5DE",bins = 50) +
  xlab("Isoform number difference between B and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

p3<-ggplot(homeologs_clean, aes(x=diffAD))+
  geom_histogram(color="black",fill="#6C16B8",bins = 50) +
  xlab("Isoform number difference between A and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

aplot::plot_list(p1, p2, p3, labels=c('AB','BD','AD'), nrow=3)

## Link between isoform number and expression level:

expression_level<-read.table("Output/homeologs_pvalue.csv", sep = ",", header=TRUE, row.names = 1)
merged<-merge(expression_level,homeologs,by="ID")
all_iso<-c(merged$Isoform_number_A,merged$Isoform_number_B,merged$Isoform_number_D)
all_means<-c(merged$means_A,merged$means_B,merged$means_D)
data<-data.frame(unlist(cbind(all_iso,all_means)))
data$all_iso <- as.numeric(data$all_iso)
data$all_means <- as.numeric(data$all_means)
ggplot(data) + geom_point(aes(x = all_iso, y = all_means),size=3) +
  xlab("Gene isoform number") + ylab("Gene mean expression level") +
  theme(text = element_text(size=15))
