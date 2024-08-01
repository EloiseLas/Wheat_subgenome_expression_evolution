library(ggplot2)

## Annotation pieplot:
data=read.table("Output/homeologs_annotation_similarity_count.csv", h=T, sep=",")
pie(data$A,data$category)

data=read.table('Output/annotation_category_counts.csv', h=T, sep=",")
pie(data$ID,data$category)

bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) +
  theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie

## Modules pieplot:
data=read.table("Output/triads_modules_comparisons_TPM_v4.csv", h=T, sep=",")
mod_cat_group <- data %>% group_by(precise_category) %>% tally()
bp<- ggplot(mod_cat_group, aes(x="", y=n, fill=precise_category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c( "#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30")) + 
  theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie

## Expression level pieplot:
data=read.table('Output/homeo_pvalue_category_counts_new.csv', h=T, sep=",")
# Barplot
bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + 
  theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie

rownames(data)<-data$category
new_data<-data.frame(Category=character(),Number=integer())
new_data[1,]<-c("A_high",(data["A+",2]+data["F_ABD",2]+data["F_ADB",2]))
new_data[2,]<-c("A_low",(data["A-",2]+data["F_BDA",2]+data["F_DBA",2]))
new_data[3,]<-c("B_high",(data["B+",2]+data["F_BAD",2]+data["F_BDA",2]))
new_data[4,]<-c("B_low",(data["B-",2]+data["F_ADB",2]+data["F_DAB",2]))
new_data[5,]<-c("D_high",(data["D+",2]+data["F_DBA",2]+data["F_DAB",2]))
new_data[6,]<-c("D_low",(data["D-",2]+data["F_BAD",2]+data["F_ABD",2]))
new_data[7,]<-c("S",data["S",2])
new_data$Number<-as.numeric(new_data$Number)

bp<- ggplot(new_data, aes(x=Category, y=Number, fill=Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#FF5733", "#FAA897", "#0EA5DE", "#24D2F9", "#6C16B8", "#CAA0EE")) +
  theme(text = element_text(size=23)) + xlab("Simplified Expression Level Categories") + guides(fill="none")
bp

new_data<-data.frame(Category=character(),Number=integer())
new_data[1,]<-c("A_high",(data["A+",2]+data["F_ABD",2]+data["F_DAB",2]+data["B-",2]+data["F_ADB",2]))
new_data[2,]<-c("B_high",(data["B+",2]+data["F_BAD",2]+data["F_BDA",2]+data["F_DBA",2]+data["A-",2]))
new_data[3,]<-c("S",(data["S",2]+data["D-",2]+data["D+",2]))
new_data$Number<-as.numeric(new_data$Number)

bp<- ggplot(new_data, aes(x=Category, y=Number, fill=Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#FF5733", "#0EA5DE", "#F9EA30")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10000)) +
  theme(text = element_text(size=15))
bp

data=read.table('Output/homeo_pvalue_category_counts_nolengthdiff.csv', h=T, sep=",")

bp<- ggplot(data, aes(x="", y=D, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#FAA897", "#0EA5DE", "#24D2F9", "#6C16B8", "#CAA0EE", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + ggtitle("Expression level categories for same length homeologs")
pie

data=read.table('Output/homeo_pvalue_category_counts_lengthdiff.csv', h=T, sep=",")

bp<- ggplot(data, aes(x="", y=D, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#FAA897", "#0EA5DE", "#24D2F9", "#6C16B8", "#CAA0EE", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + ggtitle("Expression level categories for different length homeologs")
pie


## Length pieplot:

data=read.table('Output/length_category_counts.csv', h=T, sep=",")
pie(data$ID,data$category)
data[,1:2]

bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")

pie

## Isoform number pieplot:

data=read.table('Output/isoform_category_counts.csv', h=T, sep=",")
pie(data$ID,data$category)

bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie

## Expression profiles cor pieplot:

data=read.table('Output/cor_category_counts.csv', h=T, sep=",")
pie(data$ID,data$category)

bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30")) + ggtitle("Expression profiles categories") +
  theme(text = element_text(size=15))
pie

data=read.table('Output/homeo_cor_category.csv', h=T, sep=",")
subset_F<-subset(data,category=='F')
hist(subset_F$cor_AB)
hist(subset_F$cor_BD)
hist(subset_F$cor_AD)

mean(subset_F$cor_AB)
mean(subset_F$cor_BD)
mean(subset_F$cor_AD)

subset_A<-subset(data,category=='A')
hist(subset_A$cor_AB)
hist(subset_A$cor_AD)
hist(subset_A$cor_BD)

mean(subset_A$cor_AB)
mean(subset_A$cor_AD)

subset_D<-subset(data,category=='D')
hist(subset_D$cor_AD)

mean(subset_D$cor_AD)
mean(subset_D$cor_BD)


## Expression profiles mr pieplot:

data=read.table('Output/logmr_category_counts.csv', h=T, sep=",")
pie(data$ID,data$category)

bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30")) + ggtitle("Expression profiles categories") +
  theme(text = element_text(size=15))
pie

## Expression profile both pieplot:

data=read.table('Output/profile_both_category_counts.csv', h=T, sep=",")
pie(data$ID,data$category)

bp<- ggplot(data, aes(x="", y=ID, fill=category)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30")) + theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie


## Chromosomes pieplot :

chr_cat<-read.table("Output/homeo_chr_category.csv", sep = ",", header=TRUE, row.names = 1)
chr_cat_group <- chr_cat %>% group_by(category) %>% tally()
chr_cat_group$category <- as.factor(chr_cat_group$category)
cols <- c("A" = "#FF5733", "B" = "#24D2F9", "D" = "#CAA0EE", "F" = "#6CCC78", "S"="#F9EA30")
bp<- ggplot(chr_cat_group, aes(x="", y=n, fill=category)) +
  geom_bar(width = 1, stat = "identity") 
pie3 <- bp + coord_polar("y", start=0) + scale_fill_manual(values=cols) + theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie3

