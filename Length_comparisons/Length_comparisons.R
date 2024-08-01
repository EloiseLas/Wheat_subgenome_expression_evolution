library(dplyr)
library(plyr)
library(tidyr)
library(RColorBrewer)

## Load the data:

length_table<-read.table("Data/Genes_length.csv", sep = ",", header=TRUE)[,2:4]

homeologs<-read.table("Data/homoeolog_Traes_EGI.csv", sep = ",", header=TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("X")]
homeologs<-merge(homeologs,length_table, by.x="EGI_LOC_A", by.y="ID", all.x=TRUE)
names(homeologs)[names(homeologs) == 'Length'] <- 'length_A'
names(homeologs)[names(homeologs) == 'High_std'] <- 'High_std_A'
homeologs<-merge(homeologs,length_table, by.x="EGI_LOC_B", by.y="ID", all.x=TRUE)
names(homeologs)[names(homeologs) == 'Length'] <- 'length_B'
names(homeologs)[names(homeologs) == 'High_std'] <- 'High_std_B'
homeologs<-merge(homeologs,length_table, by.x="EGI_LOC_D", by.y="ID", all.x=TRUE)
names(homeologs)[names(homeologs) == 'Length'] <- 'length_D'
names(homeologs)[names(homeologs) == 'High_std'] <- 'High_std_D'

## Compute length diff and stdev between homeologs:
homeologs <- transform(homeologs, stdev=apply(cbind(homeologs[,8],homeologs[,10],homeologs[,12]), 1, sd, na.rm=TRUE))
homeologs$diffAB <- abs(homeologs$length_A-homeologs$length_B)
homeologs$diffBD <- abs(homeologs$length_D-homeologs$length_B)
homeologs$diffAD <- abs(homeologs$length_A-homeologs$length_D)

homeologs$length_diff <- homeologs$stdev>200

homeologs_clean<-homeologs %>% drop_na()

write.csv(homeologs_clean,
          file = "Output/homoeolog_length_detailed.csv")

homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("length_A", "length_B", "length_D","stdev")]

write.csv(homeologs,
          file = "Output/homoeolog_Traes_EGI_length.csv")

## Statistical tests:

mean(homeologs$length_A, na.rm=TRUE)
mean(homeologs$length_B, na.rm=TRUE)
mean(homeologs$length_D, na.rm=TRUE)

homeologs_clean<-homeologs %>% drop_na()

lgth<-c(homeologs_clean$length_A,homeologs_clean$length_B,homeologs_clean$length_D)
homeo=as.factor(c(rep('A',length(homeologs_clean$length_A)),rep('B',length(homeologs_clean$length_B)),rep('D',length(homeologs_clean$length_D))))
lm=lm(lgth~homeo)
anova(lm)
pvalue[[i]]<-anova(lm)$"Pr(>F)"[1]

kruskal.test(lgth, homeo) 
pairwise.wilcox.test(lgth, homeo)

## Plots:
plot(homeologs$stdev)

homeologs_graph<-homeologs %>% drop_na(stdev)
homeologs_graph<-homeologs_graph %>% arrange(desc(stdev))
homeologs_graph$ID<-as.factor(homeologs_graph$ID)
png(filename = "Output/Plot_length_stdev_big.png",
    width = 5500, height = 3500, units = "px",
    bg = "white", res = 200)
plot(homeologs_graph$stdev, pch=21, col="lightpink",  bg="lightpink",cex=2,xlab="Triads ordered by their Length Standard Deviation", ylab="Length Standard Deviation",cex.lab=2,
     cex.axis=1.5)
dev.off()

ggplot(homeologs_graph) +
  aes(x = ID, y = stdev, color = stdev > 200) + 
  geom_point() + 
  scale_colour_manual(name = 'Standard Deviation > 200', values = setNames(c('#FF5733','black'),c(T, F))) +
  theme(text = element_text(size=13)) +
  xlab('Triads ID') + ylab('Length Standard Deviation in the triads')


p1<-ggplot(homeologs_clean, aes(x=length_A))+
  geom_histogram(color="black",fill="#FF5733", bins = 50) + xlab("Homeologs length in A") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

p2<-ggplot(homeologs_clean, aes(x=length_B))+
  geom_histogram(color="black",fill="#0EA5DE",bins = 50) +
  xlab("Homeologs length in B") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 0.6, 
             color = "black", size=1)

p3<-ggplot(homeologs_clean, aes(x=length_D))+
  geom_histogram(color="black",fill="#6C16B8",bins = 50) +
  xlab("Homeologs length in D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 0.6, 
             color = "black", size=1)

p1<-ggplot(homeologs_clean, aes(x=diffAB))+
  geom_histogram(color="black",fill="#FF5733", bins = 50) + xlab("Length difference between A and B") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

p2<-ggplot(homeologs_clean, aes(x=diffBD))+
  geom_histogram(color="black",fill="#0EA5DE",bins = 50) +
  xlab("Length difference between B and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

p3<-ggplot(homeologs_clean, aes(x=diffAD))+
  geom_histogram(color="black",fill="#6C16B8",bins = 50) +
  xlab("Length differencebetween A and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 200, 
             color = "black", size=1)

aplot::plot_list(p1, p2, p3, labels=c('AB','BD','AD'), nrow=3)

## Link between length isoform number and expression:

iso_table<-read.table("Data/Genes_iso.csv", sep = ",", header=TRUE)[,2:3]
iso_len<-merge(iso_table,length_table, by.x="ID", by.y="ID", all.x=TRUE)
plot(iso_len$Length,iso_len$Isoform_number)
ggplot(iso_len) + geom_point(aes(x = Length, y = Isoform_number)) +
  xlab('Gene Length') + ylab('Gene Isoform Number') + theme(text = element_text(size=17)) 

expression_data_UQ_batch_normalized<-read.table("Data/expression_data_TPM_UQ_batch_filtered_new.csv", sep = ",", header=TRUE, row.names = 1)
means<-rowMeans(expression_data_UQ_batch_normalized)
means<-data.frame(means)
data<-merge(iso_len,means,by.x="ID", by.y=0)
cols <- rev(brewer.pal(11, 'RdYlBu'))

ggplot(data) + geom_point(aes(x = log10(Length), y = Isoform_number, fill=means), size=4, alpha=0.7, pch=21) +
  xlab('log10(Gene Length)') + ylab('Gene Isoform Number') + theme(text = element_text(size=20)) + scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 2.1) +
  guides(fill=guide_legend(title="Mean Expression Level"))
