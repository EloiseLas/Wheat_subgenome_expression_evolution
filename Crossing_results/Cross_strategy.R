library(ade4)
library(adegraphics)
library(ca)
library(plyr)
library(dplyr)
library(tidyverse)
library(factoextra)
library(RColorBrewer)
library(hrbrthemes)
library(paletteer)
library(ggthemes)
library(explor)
library(AnnotationHub)
library(stringr)
library(clusterProfiler)
library(GO.db)
library(ggplot2)
library(enrichplot)
library(pheatmap)
library(AnnotationHub)
library(GOSemSim)
library(AnnotationForge)
library(fpc)
library(vcd)
library(gridExtra)
library(emmeans)
library(ggmosaic)
library(qgraph)
library(gt)
library(gtExtras)
library(openintro)


homeo_ID<-read.table("Data/homoeolog_Traes_EGI.csv", sep = ",", header=TRUE, row.names = 8)
homeo_ID<-homeo_ID[ ,  !names(homeo_ID) %in% 
            c("X")]

## Load all the categorical variables:

modules<-read.table("Output/triads_modules_comparisons_TPM_v4.csv", sep = ",", header=TRUE, row.names = 1)
modules_clean<-data.frame(cbind(modules$ID,modules$precise_category))
colnames(modules_clean)<-c("ID","modules")

expression_cat_n<-read.table("Output/homeo_pvalue_category_detail_wilcox_new.csv", sep = ",", header=TRUE, row.names = 1)
expression_clean_n<-data.frame(cbind(expression_cat_n$ID,expression_cat_n$category))
colnames(expression_clean_n)<-c("ID","expression_level")

expression_cat<-read.table("Output/homeo_pvalue_category_detail_wilcox.csv", sep = ",", header=TRUE, row.names = 1)
expression_clean<-data.frame(cbind(expression_cat$ID,expression_cat$category))
colnames(expression_clean)<-c("ID","expression_level")

annotation_cat<-read.table("Output/homeo_annotation_GO_category_precise2.csv", sep = ",", header=TRUE, row.names = 1)
annotation_clean<-data.frame(cbind(annotation_cat$ID,annotation_cat$category))
colnames(annotation_clean)<-c("ID","annotation")

length_cat<-read.table("Output/homeo_length_category.csv", sep = ",", header=TRUE, row.names = 1)
length_clean<-data.frame(cbind(length_cat$ID,length_cat$category,length_cat$High_std_flag))
colnames(length_clean)<-c("ID","length","high_std")

iso_cat<-read.table("Output/homeo_iso_category_detail.csv", sep = ",", header=TRUE, row.names = 1)
iso_clean<-data.frame(cbind(iso_cat$ID,iso_cat$category))
colnames(iso_clean)<-c("ID","isoform")

profile_cor_cat<-read.table("Output/homeo_cor_category.csv", sep = ",", header=TRUE, row.names = 1)
profile_cor_clean<-data.frame(cbind(profile_cor_cat$ID,profile_cor_cat$category))
colnames(profile_cor_clean)<-c("ID","profile_cor")

profile_mr_cat<-read.table("Output/homeo_logmr_category.csv", sep = ",", header=TRUE, row.names = 1)
profile_mr_clean<-data.frame(cbind(profile_mr_cat$ID,profile_mr_cat$category))
colnames(profile_mr_clean)<-c("ID","profile_mr")


profile_both_cat<-read.table('Output/homeo_profile_both_category.csv', sep = ",", header=TRUE, row.names = 1)
profile_both_clean<-data.frame(cbind(profile_both_cat$ID,profile_both_cat$category))
colnames(profile_both_clean)<-c("ID","profile")

global_expression_cat_n<-read.table("Output/homeo_global_expression_level_category_new.csv", sep = ",", header=TRUE, row.names = 1)
global_expression_clean_n<-data.frame(cbind(global_expression_cat_n$V1,global_expression_cat_n$category))
colnames(global_expression_clean_n)<-c("ID","global_expression_level")

global_expression_cat<-read.table("Output/homeo_global_expression_level_category.csv", sep = ",", header=TRUE, row.names = 1)
global_expression_clean<-data.frame(cbind(global_expression_cat$V1,global_expression_cat$category))
colnames(global_expression_clean)<-c("ID","global_expression_level")

chr_cat<-read.table("Output/homeo_chr_category.csv", sep = ",", header=TRUE, row.names = 1)
chr_clean<-data.frame(cbind(chr_cat$ID,chr_cat$category))
colnames(chr_clean)<-c("ID","chromosome")

all<-join_all(list(expression_clean_n,length_clean,iso_clean,profile_both_clean,global_expression_clean_n,chr_clean),by="ID")

rownames(all)<-all$ID
all<-all[ ,  !names(all) %in% 
                        c("ID")]

str(all)
all$modules <- as.factor(all$modules)
all$expression_level <- as.factor(all$expression_level)
all$length <- as.factor(all$length)
all$isoform <- as.factor(all$isoform)
all$profile_cor <- as.factor(all$profile_cor)
all$profile_mr <- as.factor(all$profile_mr)
all$profile <- as.factor(all$profile)
all$global_expression_level <- as.factor(all$global_expression_level)
all$chromosome <- as.factor(all$chromosome)
all$high_std <- as.factor(all$high_std)

all_group <- all %>% group_by(expression_level,profile_cor,global_expression_level,length,isoform) %>% tally()
all_group <- all_group[order(all_group$n,decreasing = TRUE),]

test <- all %>% group_by(global_expression_level) %>% tally()

write.csv(all_group,
          file = "Output/Cross_strategy_top_categories.csv")

## Cross the categorical variables:

# Expression x Modules:
merged<-merge(modules,expression_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$precise_category,merged$category.y))

tab<-table(data$X2,data$X3)

plot(ca(tab))

mosaicplot(tab,shade = TRUE)

coa=dudi.coa(unclass(tab),scannf = TRUE)

scatter(coa,method=2,posieig = "bottomright")

ic <- inertia.dudi(coa,row.inertia = T, col.inertia = T)
names(ic)
summary(ic)

plot(ic,contrib="abs")

chisq.test(tab)

# Expression profile x Modules:
merged<-merge(modules,profile_both_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$precise_category,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

all_cor<-c(merged$cor_AB,merged$cor_BD,merged$cor_AD)
all_mr<-c(merged$AB_log_mr,merged$BD_log_mr,merged$AD_log_mr)
all_mod<-c(merged$AB_comp,merged$BD_comp,merged$AD_comp)
data<-data.frame(unlist(cbind(all_cor,all_mr,all_mod)))
data$all_cor <- as.numeric(data$all_cor)
data$all_mr <- as.numeric(data$all_mr)
data$all_mod  <- as.numeric(data$all_mod)
ggplot(data) + geom_point(aes(x = all_cor, y = all_mr, color= all_mod == 1 )) +
  scale_colour_manual(name = 'AB_pv < 0.05', values = setNames(c('#4DBC8B','black'),c(T, F))) +
  xlab("Spearman Correlation Coefficient (SCC)") + ylab("Log Mutual Rank (MR)") +
  guides(color=guide_legend(title="Same module: ")) +
  theme(text = element_text(size=15))

# Annotation x Modules:

merged<-merge(modules,annotation_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$precise_category,merged$category.y))

filter(merged, category.y=="A_diff")
filter(merged, category.y=="B_diff")
filter(merged, category.y=="D_diff")
filter(merged, category.y=="All_diff")

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

# Annotation x Global Expression level:

merged<-merge(global_expression_cat,annotation_cat, by.x="V1", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$global_expression_level,merged$category))

filter(merged, category.y=="A_diff")
filter(merged, category.y=="B_diff")
filter(merged, category.y=="D_diff")
filter(merged, category.y=="All_diff")

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

plot(merged$means,merged$cor_AB)
plot(merged$means,merged$cor_BD)
plot(merged$means,merged$cor_AD)

boxplot(merged$means ~ merged$category.y , 
        ylab="Mean expression level of the triad" , xlab="Annotation Category")

# Annotation x Expression Profile cor:

merged<-merge(profile_cor_cat,annotation_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

filter(merged, category.y=="A_diff")
filter(merged, category.y=="B_diff")
filter(merged, category.y=="D_diff")
filter(merged, category.y=="All_diff")

data_clean<-filter(data, X3!="Unknown")

tab<-table(data_clean$X2,data_clean$X3)

mosaicplot(tab,shade = TRUE)

# Annotation x Expression Profile:

merged<-merge(profile_both_cat,annotation_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))


filter(merged, category.y=="A_diff")
filter(merged, category.y=="B_diff")
filter(merged, category.y=="D_diff")
filter(merged, category.y=="All_diff")
filter(merged, category.y=="Same_GO", category.x=="A" )

data_clean<-filter(data, X3!="Unknown")

tab<-table(data$X2,data$X3)

contTable(tab)

mosaicplot(tab,shade = TRUE)

# Expression Profile cor x Expression Profile mr:

merged<-merge(profile_mr_cat,profile_cor_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

data_clean<-filter(data, X3!="Unknown")

tab<-table(data$X2,data$X3)
plot(tab)

mosaicplot(tab,shade = TRUE)


# Length x Expression:

merged<-merge(expression_cat,length_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

chisq.test(tab)

# Length x Expression:

merged<-merge(expression_cat,expression_cat_n, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

chisq.test(tab)

# Expression profile cor x Expression:

merged<-merge(expression_cat_n,profile_both_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

# Length x Expression profile cor :

merged<-merge(profile_cor_cat,length_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)


# Chr x Expression profile cor :

merged<-merge(profile_both_cat,chr_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

data<-data.frame(cbind(merged$ID,merged$category.x,merged$special))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

# Chr x Expression level :

merged<-merge(expression_cat,chr_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

data<-data.frame(cbind(merged$ID,merged$category.x,merged$special))

sub<-filter(merged,special=="A_4" | special=="A_5")

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

# Isoform number x Expression level :

merged<-merge(iso_cat,expression_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

Diff<-subset(merged,merged$category.x=="A+" & merged$category.x=="A-" & merged$category.y=="A+" & merged$category.y=="A-")

iso<-merged$Isoform_number_A
means<-merged$means_A
data<-data.frame(cbind(iso,means))
data<-data %>% drop_na()
data$iso<-as.numeric(data$iso)
data$means<-as.numeric(data$means)

ggplot(data) +
  aes(x = iso, y = means) + 
  geom_point() + 
  xlab('Global Expression Level') + ylab('-log10(pv)')

# Isoform number x Expression profile cor :

merged<-merge(iso_cat,profile_cor_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

# Isoform number x Length :

merged<-merge(iso_cat,length_cat, by.x="ID", by.y="ID")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

# Global expression level x Expression level :

merged<-merge(expression_cat_n,global_expression_cat_n, by.x="ID", by.y="V1")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

means<-merged$means
length<-merged$category.x
pairwise.wilcox.test(means, length)

ggplot(merged, aes(x=category.x, y=means, fill = category.x)) + 
  geom_boxplot() +
  annotate("text", x = 1, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 2, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 3, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 4, y=2, label = "b", size = unit(5, "pt")) +
  annotate("text", x = 5, y=3.3, label = "d", size = unit(5, "pt")) +
  scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

# Global expression level old vs new:

merged<-merge(global_expression_cat_n,global_expression_cat, by.x="V1", by.y="V1")
plot(merged$means.x,merged$means.y)

# Global expression level x Expression profile :

merged<-merge(profile_cor_cat,global_expression_cat, by.x="ID", by.y="V1")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

plot(merged$means,merged$cor_AB)
plot(merged$means,merged$cor_BD)
plot(merged$means,merged$cor_AD)

data<-data.frame(cbind(merged$cor_AB,merged$means))

sim_matrix<-1-cor(t(data))
hclust_obj <- hclust(dist(sim_matrix), method="complete")
plot(hclust_obj)
data$Clusters<-as.factor(cutree(hclust_obj, k = 5))
Clusters<-lapply(Clusters, function(i){paste("c",i,sep="")})
Clusters<-t(data.frame(Clusters))
Clusters<-data.frame(Clusters)
colnames(Clusters)<-"Clusters"

km.res<-kmeans(data, 3, iter.max = 10, nstart = 25)
data$Clusters<-as.factor(km.res$cluster)

db <- fpc::dbscan(data, 0.03, MinPts = 10, scale = FALSE)
data$Clusters<-as.factor(db$cluster)

ggplot(data) +
  aes(x = X2, y = X1, color=Clusters) + 
  geom_point() + 
  xlab('Global Expression Level') + ylab('-log10(pv)')

ggplot(data) +
  aes(x = X4, y = X2, color=Clusters) + 
  geom_point() + 
  xlab('Global Expression Level') + ylab('-log10(pv)')

chisq.test(tab)

boxplot(merged$means ~ merged$category.x , 
        ylab="Mean expression level of the triad" , xlab="Annotation Category")

means<-merged$means
expression_pro<-merged$category.x
lm=lm(means~expression_pro)
plot(lm)
hist(lm$residuals)
anova(lm)
summary(lm)
effect=emmeans(lm,~expression_l)
summary(pairs(effect,adjust="tukey"))
pairwise.wilcox.test(means, expression_pro, p.adjust.method="hochberg")

ggplot(merged, aes(x=expression_pro, y=means, fill = expression_pro)) + 
  geom_boxplot() +
  annotate("text", x = 1, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 2, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 3, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 4, y=2, label = "b", size = unit(5, "pt")) +
  annotate("text", x = 5, y=3.3, label = "d", size = unit(5, "pt")) +
  scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

# Global expression level x Expression profile :

merged<-merge(profile_both_cat,global_expression_cat_n, by.x="ID", by.y="V1")

sub<-filter(merged, cor_AD<0.3, AD_log_mr<1)

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

plot(merged$means,merged$AB_log_mr)
plot(merged$means,AD)
plot(merged$means,log(merged$BD_mr))

chisq.test(tab)
means<-merged$means
expression_pro<-merged$category.x
pairwise.wilcox.test(means, expression_pro)

ggplot(merged, aes(x=category.x, y=means, fill = category.x)) + 
  geom_violin() +
  annotate("text", x = 1, y=3, label = "a,d", size = unit(9, "pt")) +
  annotate("text", x = 2, y=3, label = "a", size = unit(9, "pt")) +
  annotate("text", x = 3, y=3, label = "d", size = unit(9, "pt")) +
  annotate("text", x = 4, y=2.6, label = "b", size = unit(9, "pt")) +
  annotate("text", x = 5, y=3.3, label = "c", size = unit(9, "pt")) +
  scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=23))

all_cor<-c(merged$cor_AB,merged$cor_BD,merged$cor_AD)
all_mr<-c(merged$AB_log_mr,merged$BD_log_mr,merged$AD_log_mr)
all_global_cat<-c(merged$category.y,merged$category.y,merged$category.y)
all_cat<-c(merged$category.x,merged$category.x,merged$category.x)
all_means<-c(merged$means,merged$means,merged$means)
data<-data.frame(unlist(cbind(all_cor,all_mr,all_global_cat,all_cat,all_means)))
data$all_cor <- as.numeric(data$all_cor)
data$all_mr <- as.numeric(data$all_mr)
data$all_means <- as.numeric(data$all_means)
data$all_global_cat <- factor(data$all_global_cat)
data$all_cat <- factor(data$all_cat)
ggplot(data) + geom_point(aes(x = all_cor, y = all_mr, color= all_global_cat),size=2) +
  xlab("Spearman Correlation Coefficient (SCC)") + ylab("Log Mutual Rank (MR)") +
  guides(color=guide_legend(title="Global expression level of the triad")) +
  theme(text = element_text(size=23))
ggplot(data) + geom_point(aes(x = all_cor, y = all_mr, color= all_cat, alpha = all_means),size=2) +
  xlab("Spearman Correlation Coefficient (SCC)") + ylab("Log Mutual Rank (MR)") +
  guides(color=guide_legend(title="Expression profile category")) +
  scale_alpha("Global expression level of the triad",range = c(0.3, 1)) +
  scale_colour_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE", "#6CCC78", "#F9EA30"))+
  theme(text = element_text(size=23))


# Global expression level x Expression profile :

merged<-merge(profile_mr_cat,global_expression_cat, by.x="ID", by.y="V1")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

plot(merged$means,log(merged$AB_mr))
plot(merged$means,log(merged$AD_mr))
plot(merged$means,log(merged$BD_mr))

chisq.test(tab)

# Global expression level x Length :

merged<-merge(length_cat,global_expression_cat, by.x="ID", by.y="V1")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

chisq.test(tab)

means<-merged$means
length<-merged$category.x
pairwise.wilcox.test(means, length)

ggplot(merged, aes(x=category.x, y=means, fill = category.x)) + 
  geom_boxplot() +
  annotate("text", x = 1, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 2, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 3, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 4, y=2, label = "b", size = unit(5, "pt")) +
  annotate("text", x = 5, y=3.3, label = "d", size = unit(5, "pt")) +
  scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

# Global expression level x Isoform :

merged<-merge(iso_cat,global_expression_cat, by.x="ID", by.y="V1")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

chisq.test(tab)

means<-merged$means
iso<-merged$category.x
pairwise.wilcox.test(means, iso)

ggplot(merged, aes(x=category.x, y=means, fill = category.x)) + 
  geom_boxplot() +
  annotate("text", x = 1, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 2, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 3, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 4, y=2, label = "b", size = unit(5, "pt")) +
  annotate("text", x = 5, y=3.3, label = "d", size = unit(5, "pt")) +
  scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

# Global expression level x chr :

merged<-merge(chr_cat,global_expression_cat, by.x="ID", by.y="V1")

data<-data.frame(cbind(merged$ID,merged$category.x,merged$category.y))

tab<-table(data$X2,data$X3)

mosaicplot(tab,shade = TRUE)

chisq.test(tab)

means<-merged$means
length<-merged$category.x
pairwise.wilcox.test(means, length)

ggplot(merged, aes(x=category.x, y=means, fill = category.x)) + 
  geom_boxplot() +
  annotate("text", x = 1, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 2, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 3, y=2, label = "a", size = unit(5, "pt")) +
  annotate("text", x = 4, y=2, label = "b", size = unit(5, "pt")) +
  annotate("text", x = 5, y=3.3, label = "d", size = unit(5, "pt")) +
  scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

# Chi-sq Heatmap:

all_clean<-all[ ,  !names(all) %in% c("high_std","ID","global_expression_level")]
table<-outer(all_clean, all_clean, Vectorize(\(x, y) -log10(chisq.test(table(x, y))$p.value)))
table[sapply(table, is.infinite)] <- NA
rownames(table)<-c("Expression Level", "Length","Isoform","Expression Profile", "Chromosome")
colnames(table)<-c("Expression Level", "Length","Isoform","Expression Profile", "Chromosome")
pheatmap(table,cluster_rows = FALSE, breaks=c(0,10,20,30,40,50,60,70,80), color = colorRampPalette(c("white", "#01A71D"))(10), cluster_cols = FALSE,display_numbers=TRUE,fontsize=20,legend_breaks =c(0,10,20,30,40,50,60,70,80),fontsize_number=27,legend_labels = c("0","10", "20", "30", "40","50","60","70", "- log10(pvalue)"))

qgraph(table, layout='spring',labels=c("Level", "Length","Isoform","Profile", "Chr"))

assoc(~length + isoform, data=all,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

assoc(~profile + isoform, data=all,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

assoc(~profile + length, data=all,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

assoc(~profile + expression_level, data=all,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

assoc(~chromosome + profile, data=all,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

assoc(~chromosome + length, data=all,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

## Working on interesting subsets:

stats <- function(sub) { # create a function with the name my_function
  glob_expression_level_group <- sub %>% group_by(global_expression_level) %>% tally()
  isoform_group <- sub %>% group_by(isoform) %>% tally()
  length_group <- sub %>% group_by(length) %>% tally()
  flag_group <- sub %>% group_by(high_std) %>% tally()
  chr_group <- sub %>% group_by(chromosome) %>% tally()
  glob_expression_level_group$global_expression_level <- as.factor(glob_expression_level_group$global_expression_level)
  isoform_group$isoform <- as.factor(isoform_group$isoform)
  length_group$length <- as.factor(length_group$length)
  flag_group$high_std <- as.factor(flag_group$high_std)
  chr_group$chromosome <- as.factor(chr_group$chromosome)
  cols <- c("A" = "#FF5733", "B" = "#24D2F9", "D" = "#CAA0EE", "F" = "#6CCC78", "S"="#F9EA30")
  bp<- ggplot(glob_expression_level_group, aes(x="", y=n, fill=global_expression_level)) +
    geom_bar(width = 1, stat = "identity")
  pie1 <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE")) + ggtitle("Global Expression Level Categories")
  
  bp<- ggplot(isoform_group, aes(x="", y=n, fill=isoform)) +
    geom_bar(width = 1, stat = "identity")
  pie2 <- bp + coord_polar("y", start=0) + scale_fill_manual(values=cols) + ggtitle("Isoform categories")
  
  bp<- ggplot(length_group, aes(x="", y=n, fill=length)) +
    geom_bar(width = 1, stat = "identity")
  pie3 <- bp + coord_polar("y", start=0) + scale_fill_manual(values=cols) + ggtitle("Length categories")
  
  bp<- ggplot(flag_group, aes(x="", y=n, fill=high_std)) +
    geom_bar(width = 1, stat = "identity")
  pie4 <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#24D2F9","#FF5733")) + ggtitle("High Standard Deviation")
  
  bp<- ggplot(chr_group, aes(x="", y=n, fill=chromosome)) +
    geom_bar(width = 1, stat = "identity")
  pie5 <- bp + coord_polar("y", start=0) + scale_fill_manual(values=cols) + ggtitle("Chromosome Category")
  
  cowplot::plot_grid(pie1, pie2, pie3, pie4, pie5, ncol=2, nrow = 3,  labels=LETTERS[1:4])
}

# Same expression profile and level:
Same_expression_profile_level<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
stats(Same_expression_profile_level)

# Same expression profile and different levels:
Same_expression_profile_different_expression_level<-subset(all,all$profile_cor=="S" & all$expression_level!="S")
stats(Same_expression_profile_different_expression_level)
expression_level_group <- Same_expression_profile_different_expression_level %>% group_by(expression_level) %>% tally()
bp<- ggplot(expression_level_group, aes(x="", y=n, fill=expression_level)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + ggtitle("Expression Level Categories")
pie

# Same expression level:
Same_expression_level<-subset(all,all$expression_level=="S")
Same_expression_level_clean<-Same_expression_level[ ,  !names(Same_expression_level) %in% 
                                                          c("high_std","expression_level","global_expression_level")]
table<-data.frame(outer(Same_expression_level_clean, Same_expression_level_clean, Vectorize(\(x, y) chisq.test(table(x, y))$statistic)))
table[sapply(table, is.infinite)] <- NA
rownames(table)<-c("Length","Isoform","Expression profile","Chromosome")
colnames(table)<-c("Length","Isoform","Expression profile","Chromosome")
pheatmap(table,cluster_rows = FALSE, breaks=c(0,10,20,30,40,50,60,70,80), color = colorRampPalette(c("white", "#01A71D"))(10), cluster_cols = FALSE,display_numbers=TRUE,fontsize=20,legend_breaks =c(0,10,20,30,40,50,60,70,80),fontsize_number=27,legend_labels = c("0","10", "20", "30", "40","50","60","70", "- log10(pvalue)"))

# Same expression profile:
Same_expression_profile<-subset(all,all$profile=="S")
write.csv(Same_expression_profile,
          file = "Output/Cross_strategy_same_expression_profile.csv")
stats(Same_expression_profile)
expression_level_group <- Same_expression_profile %>% group_by(expression_level) %>% tally()
bp<- ggplot(expression_level_group, aes(x="", y=n, fill=expression_level)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FAA897", "#FF5733", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + theme_minimal(base_size = 25) + guides(fill=guide_legend(title="Categories")) + ylab("") + xlab("")
pie

# Mosaic plots:
par(mfrow=c(2,2))
# Expression level x length:
tab<-table(Same_expression_profile$expression_level,Same_expression_profile$length)
test<-chisq.test(tab)
mosaicplot(~expression_level + length, data=Same_expression_profile,shade = TRUE, main='Expression Level Categories',xlab="",cex.axis=1.2,ylab="Length Categories", sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))
# Expression level x iso:
tab<-table(Same_expression_profile$expression_level,Same_expression_profile$isoform)
test<-chisq.test(tab)
mosaicplot(~expression_level + isoform, data=Same_expression_profile,shade = TRUE, main='Expression Level Categories',xlab="",cex.axis=1.2,ylab="Isoform Number Categories",sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))
# Expression level x chr:
tab<-table(Same_expression_profile$expression_level,Same_expression_profile$chromosome)
test<-chisq.test(tab)
mosaicplot(~expression_level + chromosome, data=Same_expression_profile,shade = TRUE, main='',xlab="",cex.axis=1.2,ylab="Chromosome Categories",sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))
# Expression level x Global expression level:
tab<-table(Same_expression_profile$expression_level,Same_expression_profile$global_expression_level)
test<-chisq.test(tab)
mosaicplot(~expression_level + global_expression_level, data=Same_expression_profile,shade = TRUE, main='',xlab="",cex.axis=1.2,ylab="Global Expression Level Categories",sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))

# Smaller mosaic plot:
par(mfrow=c(1,2))
# Expression level x length:
tab<-table(Same_expression_profile$expression_level,Same_expression_profile$length)
test<-chisq.test(tab)
mosaicplot(~expression_level + length, data=Same_expression_profile,shade = TRUE, main='Expression Level Categories',xlab="",cex.axis=1.2,ylab="Length Categories", sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))
# Expression level x iso:
tab<-table(Same_expression_profile$expression_level,Same_expression_profile$isoform)
test<-chisq.test(tab)
mosaicplot(~expression_level + isoform, data=Same_expression_profile,shade = TRUE, main='Expression Level Categories',xlab="",cex.axis=1.2,ylab="Isoform Number Categories",sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))
# iso x length:
tab<-table(Same_expression_profile$isoform,Same_expression_profile$lenth)
test<-chisq.test(tab)
mosaicplot(~isoform + length, data=Same_expression_profile,shade = TRUE, main='',xlab="",cex.axis=1.2,ylab="Chromosome Categories",sub=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))

assoc(~isoform + length, data=Same_expression_profile,
      direction = c("v", "h"),
      shade = TRUE,
      rot_labels=c(90,90,0,0)
)

assoc(~isoform + expression_level, data=Same_expression_profile,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

assoc(~length + expression_level, data=Same_expression_profile,
      direction = c("v", "h"),
      shade = TRUE,
      compress = TRUE
)

# ggplot(data = Same_expression_profile) + 
#   geom_mosaic(aes(x = product(isoform), fill = isoform),
#               divider = c("vspine")) +   
#   labs(x="Number of Planets in the Star System", 
#        y = "Discovery Method", 
#        title = "Size of Planetary System found by Discovery Method") +
#   facet_grid(~length) +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())

# Link with the average expression level of the triad:
merged<-merge(Same_expression_profile,global_expression_cat_n, by.x=0, by.y="V1")
boxplot(merged$means ~ merged$expression_level , 
        ylab="Mean expression level of the triad" , xlab="Annotation Category")
means<-merged$means
expression_l<-merged$expression_level
lm=lm(means~expression_l)
plot(lm)
hist(lm$residuals)
anova(lm)
summary(lm)
effect=emmeans(lm,~expression_l)
summary(pairs(effect,adjust="tukey"))
pairwise.wilcox.test(means, expression_l, p.adjust.method="holm")

Same_expression_profile_clean<-Same_expression_profile[ ,  !names(Same_expression_profile) %in% 
            c("high_std","profile","global_expression_level")]
table<-data.frame(outer(Same_expression_profile_clean, Same_expression_profile_clean, Vectorize(\(x, y) -log10(chisq.test(table(x, y))$p.value))))
table[sapply(table, is.infinite)] <- NA
rownames(table)<-c("Expression level","Length","Isoform","Chromosome")
colnames(table)<-c("Expression level","Length","Isoform","Chromosome")
p1<-pheatmap(table,cluster_rows = FALSE, breaks=c(0,10,20,30,40,50,60,70,80), color = colorRampPalette(c("white", "#01A71D"))(10), cluster_cols = FALSE,display_numbers=TRUE,fontsize=14,legend_breaks =c(0,10,20,30,40,50,60,70,80),fontsize_number=14,legend_labels = c("0","10", "20", "30", "40","50","60","70", "- log10(pvalue)"))
pheatmap(table,cluster_rows = FALSE, breaks=c(0,10,20,30,40,50,60,70,80), color = colorRampPalette(c("white", "#01A71D"))(10), cluster_cols = FALSE,display_numbers=TRUE,fontsize=20,legend_breaks =c(0,10,20,30,40,50,60,70,80),fontsize_number=27,legend_labels = c("0","10", "20", "30", "40","50","60","70", "- log10(pvalue)"))

# With old data:
ggplot(merged, aes(x=expression_l, y=means, fill = expression_l)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  annotate("text", x = 1, y=3, label = "a, e", size = unit(8, "pt")) +
  annotate("text", x = 2, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 3, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 4, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 5, y=3, label = "a, e, f", size = unit(8, "pt")) +
  annotate("text", x = 6, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 7, y=3.3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 8, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 9, y=3, label = "b, e", size = unit(8, "pt")) +
  annotate("text", x = 10, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 11, y=3, label = "b, f", size = unit(8, "pt")) +
  annotate("text", x = 12, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 13, y=3, label = "d", size = unit(8, "pt")) +
  scale_fill_manual(values=c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Level Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=23))

# With new data:
ggplot(merged, aes(x=expression_l, y=means, fill = expression_l)) + 
  geom_violin() + geom_boxplot(width=0.1) +
  annotate("text", x = 1, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 2, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 3, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 4, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 5, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 6, y=3, label = "a", size = unit(8, "pt")) +
  annotate("text", x = 7, y=3.3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 8, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 9, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 10, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 11, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 12, y=3, label = "b", size = unit(8, "pt")) +
  annotate("text", x = 13, y=3, label = "d", size = unit(8, "pt")) +
  scale_fill_manual(values=c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Level Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=23))

# Different expression profile:
Different_expression_profile<-subset(all,all$profile!="S")
stats(Different_expression_profile)
expression_level_group <- Different_expression_profile %>% group_by(expression_level) %>% tally()
bp<- ggplot(expression_level_group, aes(x="", y=n, fill=expression_level)) +
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) + ggtitle("Expression Level Categories")
pie
# Expression level x length:
tab<-table(Different_expression_profile$expression_level,Different_expression_profile$length)
mosaicplot(tab,shade = TRUE)
chisq.test(tab)
# Expression level x iso:
tab<-table(Different_expression_profile$expression_level,Different_expression_profile$isoform)
mosaicplot(tab,shade = TRUE)
chisq.test(tab)
# Expression level x chr:
tab<-table(Different_expression_profile$expression_level,Different_expression_profile$chromosome)
mosaicplot(tab,shade = TRUE)
chisq.test(tab)
# Expression level x Global expression level:
tab<-table(Different_expression_profile$expression_level,Different_expression_profile$global_expression_level)
mosaicplot(tab,shade = TRUE)
chisq.test(tab)
# Expression level x expression profile:
Different_expression_profile$profile<-factor(Different_expression_profile$profile)
tab<-table(Different_expression_profile$expression_level,Different_expression_profile$profile)
test<-chisq.test(tab)
mosaicplot(~expression_level + profile, data=Different_expression_profile,shade = TRUE,cex.axis=1.2,ylab="Expression Profile Categories",xlab="Expression Level Categories",main=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))


Different_expression_profile_clean<-Different_expression_profile[ ,  !names(Different_expression_profile) %in% 
                                                          c("high_std")]
Different_expression_profile_clean$profile_cor<-factor(Different_expression_profile_clean$profile_cor)
table<-outer(Different_expression_profile_clean, Different_expression_profile_clean, Vectorize(\(x, y) chisq.test(table(x, y))$statistic))
p2<-pheatmap(table,cluster_rows = FALSE, breaks=c(0,100,1000,10000,100000), color = colorRampPalette(c("white", "#01A71D"))(5), cluster_cols = FALSE,display_numbers=TRUE,main="Triads with different expression profiles")

z <- grid.arrange(p1[[4]],p2[[4]],nrow=1)

merged<-merge(Different_expression_profile,global_expression_cat, by.x=0, by.y="V1")

ggplot(merged, aes(x=expression_level, y=means, fill = expression_level)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30")) +
  xlab('Expression Level Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))


## Enrichment analysis:

# On expression categories for triads with similar expression profiles:

expression_data<-read.table("Data/expression_data_TPM_UQ_batch_filtered.csv", sep = ",", header=TRUE, row.names = 1)
geneUniverse_EGI<-rownames(expression_data)
geneUniverse_EGI<-unlist(lapply(geneUniverse_EGI,function(i){str_sub(i,4)}))

annotation_data<-read.table("Data/iwgsc_refseqv2.1_functional_annotation.csv", sep = ",", header=TRUE, row.names = NULL)
annotation_data_all<-subset(annotation_data, f.type=="Gene Ontology")
annotation_data_BP<-subset(annotation_data_all, f.domain=="BP")
annotation_data_CC<-subset(annotation_data_all, f.domain=="CC")
annotation_data_MF<-subset(annotation_data_all, f.domain=="MF")
geneUniverse<-unlist(annotation_data_all$g2.identifier)
annotation_data_all<-cbind(annotation_data_all$f.name,annotation_data_all$g2.identifier)
annotation_data_BP<-cbind(annotation_data_BP$f.name,annotation_data_BP$g2.identifier)
annotation_data_CC<-cbind(annotation_data_CC$f.name,annotation_data_CC$g2.identifier)
annotation_data_MF<-cbind(annotation_data_MF$f.name,annotation_data_MF$g2.identifier)

columns(GO.db)
go <- keys(GO.db, keytype="GOID")
term_description <- select(GO.db, columns=c("GOID","TERM"), keys=go, keytype="GOID")

all_res<-list()
all_genes<-list()
all_genes_EGI<-list()

cats<-c("S", "A+","A-","B+","B-","D+","D-","F_ABD","F_ADB","F_BAD","F_BDA","F_DAB","F_DBA")
#cats<-c("S","D-", "A+","B-","B+","F_BAD","F_ABD","F_ADB","A-","F_BDA","D+","F_DBA","F_DAB")
#cats<-c("S", "A+","A-","B+","B-","D+","D-")

for (i in 1:13){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  size <- length(genes)
  
  # ans.go <- enricher(
  #    gene = genes,
  #    pvalueCutoff = 0.1,
  #    pAdjustMethod = "BH",
  #    universe = geneUniverse,
  #    minGSSize = 10,
  #    maxGSSize = 500,
  #    qvalueCutoff = 0.5,
  #    TERM2GENE=annotation_data_BP,
  #    TERM2NAME = term_description
  # )
  #  
  # tab.go <- as.data.frame(ans.go)
  # # tab.go<- subset(tab.go, Count>5)
  # tab.go_clean<-data.frame(tab.go$ID,-log10(tab.go$p.adjust))
  # colnames(tab.go_clean)<-c("ID",cats[[i]])
  # all_res[[i]]<-tab.go_clean
  
  
  # ans.go2 <- pairwise_termsim(ans.go)
  # p1 <- treeplot(ans.go2) + ggtitle(paste(paste("Category",cats[i]),paste("size =", size)))
  # ggsave(p1, filename = paste(paste("Output/GO_cross_cat/Same_profile_",cats[i],sep = ""),"_CC_treeplot.png",sep = ""), device = "png")
  # 
  # write.csv(tab.go,
  #           file = paste(paste("Output/GO/",categories[i],sep = ""),"_CC_results_TPM.csv",sep = ""))
  
  # ans.kegg <- enrichKEGG(gene = genes_EGI,
  #                        organism = 'taes',
  #                        universe = geneUniverse_EGI,
  #                        pvalueCutoff = 0.05)
  # tab.kegg <- as.data.frame(ans.kegg)
  # tab.kegg<- subset(tab.kegg, Count>5)
  # tab.kegg[1:5, 1:6]
  
  # write.csv(tab.kegg,
  #           file = paste(paste("Output/KEGG/",categories[i],sep = ""),"_results.csv",sep = ""))
  
  # p1 <- barplot(ans.go)
  # ggsave(p1, filename = paste(paste("Output/GO/",categories[i],sep = ""),"_CC_plot_TPM.png",sep = ""), device = "png")
  
  # p2 <- barplot(ans.kegg, showCategory=10) + ggtitle(paste(paste("Category",cats[i]),paste("size =", size)))
  # ggsave(p2, filename = paste(paste("Output/KEGG_cross_cat/Same_profile_",cats[i],sep = ""),"_barplot.png",sep = ""), device = "png")
}

res_merge<-join_all(all_res,by="ID",type = 'full')
res_merge[is.na(res_merge)] = 0
GO_point<-gsub(":", ".", res_merge$ID)
rownames(res_merge)<-GO_point
res_merge<-res_merge[ ,  !names(res_merge) %in% 
            c("ID")]

res_merge_encode<-res_merge %>%
  mutate(across(everything(), 
                ~ case_when((. > 1) & (. <= 2) ~ 1,　(. > 2) & (. <= 3) ~ 2, . > 3 ~ 3, TRUE ~ 0)))

install.packages("./org.Taestivum.eg.db",repos = NULL, type="source")

taGO <- godata("org.Taestivum.eg.db", ont="BP", keytype = "GID")
GO<-gsub("[.]", ":", rownames(res_merge))
sim_matrix<-mgoSim(GO, GO, semData=taGO, measure="Wang", combine=NULL)
hclust_obj <- hclust(dist(sim_matrix), method="complete")
plot(hclust_obj)
Clusters<-as.factor(cutree(hclust_obj, k = 5))
Clusters<-lapply(Clusters, function(i){paste("c",i,sep="")})
Clusters<-t(data.frame(Clusters))
Clusters<-data.frame(Clusters)
colnames(Clusters)<-"Clusters"

ann_colors<-list(
  Clusters = c(c1="#FFC300",c2="#D95F02",c3="#7570B3",c4="#E7298A",c5="#66A61E"))

#  c6="#0EA5DE", c7="#FF5733", c8="#2341DB", c9="#23DB90" 


pheatmap(res_merge_encode, cluster_cols=TRUE, color = colorRampPalette(c("black", "white"))(50), cluster_rows=TRUE,clustering_method="complete",show_rownames=FALSE)

pheatmap(res_merge_encode, cluster_cols=TRUE, cluster_rows= hclust_obj, show_colnames=TRUE, show_rownames=FALSE, cellwidth=20, annotation_row = Clusters,annotation_colors=ann_colors)

names(all_genes)<-cats
names(all_genes_EGI)<-cats

ck <- compareCluster(geneCluster = all_genes, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category",cex_category=10)

ck_2 <- compareCluster(geneCluster = all_genes_EGI, fun = enrichKEGG, organism = 'taes')
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck_2)

cnetplot(ck_2,node_label="category",cex_category=10)

# For same expression profile, same expression levels vs different expression levels:

list<-list()

same<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
list[[1]]<-genes_same

diff<-subset(all,all$profile_cor=="S" & all$expression_level!="S")
homeo<-homeo_ID[row.names(diff),]
genes_diff<-c(homeo[,1],homeo[,2],homeo[,3])
genes_diff<-unlist(lapply(genes_diff,function(i){str_sub(i,end=-3)}))
list[[2]]<-genes_diff

names(list)<-c("Same","Diff")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category")

# Same expression profile vs different expression profile:

list<-list()

same<-subset(all,all$profile=="S")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
list[[1]]<-genes_same

diff<-subset(all,all$profile!="S")
homeo<-homeo_ID[row.names(diff),]
genes_diff<-c(homeo[,1],homeo[,2],homeo[,3])
genes_diff<-unlist(lapply(genes_diff,function(i){str_sub(i,end=-3)}))
list[[2]]<-genes_diff

names(list)<-c("Same","Diff")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category")

# For same expression profile, A high vs rest:

list<-list()

cats<-c("A+","F_ABD","F_ADB")
all_genes<-list()

for (i in 1:3){
  sub<-subset(all,all$profile_cor=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
}

list[[1]]<-unlist(all_genes, recursive = FALSE)

cats<-c("S","D-","B-","B+","F_BAD","A-","F_BDA","D+","F_DBA","F_DAB")
all_genes<-list()

for (i in 1:4){
  sub<-subset(all,all$profile_cor=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
}

list[[2]]<-unlist(all_genes, recursive = FALSE)

names(list)<-c("A_high","Rest")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category")

# For same expression profile, profile B high vs rest:

list<-list()

cats<-c("B+","F_BAD","F_BDA")
all_genes<-list()

for (i in 1:3){
  sub<-subset(all,all$profile_cor=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
}

list[[1]]<-unlist(all_genes, recursive = FALSE)

cats<-c("S","D-", "A+","B-","F_ABD","F_ADB","A-","D+","F_DBA","F_DAB")
all_genes<-list()

for (i in 1:4){
  sub<-subset(all,all$profile_cor=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
}

list[[2]]<-unlist(all_genes, recursive = FALSE)

names(list)<-c("B_high","Rest")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category")

# For same expression profile, D high vs rest:

list<-list()

cats<-c("D+","F_DAB","F_DBA")
all_genes<-list()

for (i in 1:3){
  sub<-subset(all,all$profile_cor=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
}

list[[1]]<-unlist(all_genes, recursive = FALSE)

cats<-c("S","D-", "A+","B-","B+","F_BAD","F_ABD","F_ADB","A-","F_BDA")
all_genes<-list()

for (i in 1:4){
  sub<-subset(all,all$profile_cor=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  all_genes[[i]]<-genes
}

list[[2]]<-unlist(all_genes, recursive = FALSE)

names(list)<-c("D_high","Rest")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category")

# For same expression profile, on expression level simplified categories:

list<-list()
list_EGI<-list()

cats<-c("A+","F_ABD","F_ADB")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list[[1]]<-unlist(all_genes, recursive = FALSE)
list_EGI[[1]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("B+","F_BAD","F_BDA")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list[[2]]<-unlist(all_genes, recursive = FALSE)
list_EGI[[2]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("D+","F_DAB","F_DBA")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list[[3]]<-unlist(all_genes, recursive = FALSE)
list_EGI[[3]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("A-","F_BDA","F_DBA")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list[[4]]<-unlist(all_genes, recursive = FALSE)
list_EGI[[4]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("B-","F_ADB","F_DAB")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list[[5]]<-unlist(all_genes, recursive = FALSE)
list_EGI[[5]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("D-","F_ABD","F_BAD")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile=="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list[[6]]<-unlist(all_genes, recursive = FALSE)
list_EGI[[6]]<-unlist(all_genes_EGI, recursive = FALSE)

same<-subset(all,all$profile=="S" & all$expression_level=="S")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
genes_same_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
genes_same_EGI<-unlist(lapply(genes_same_EGI,function(i){str_sub(i,4)}))

list[[7]]<-genes_same
list_EGI[[7]]<-genes_same_EGI

names(list)<-c("A_high","B_high","D_high","A_low","B_low","D_low","Same")
names(list_EGI)<-c("A_high","B_high","D_high","A_low","B_low","D_low","Same")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID", ont="BP")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)
res<-ck@compareClusterResult
mat<-as.data.frame.matrix(xtabs(-log10(p.adjust)~Description+Cluster,res,addNA=TRUE))
mat_filter<-mat[rowSums(mat > 3) >= 1, ]
png(filename = "Output/Check/GO_comparison_heatmap_BP_new.png",
    width = 3500, height = 5600, units = "px", pointsize = 300,
    bg = "white", res = 200)
pheatmap(mat_filter, cluster_cols=TRUE, cluster_rows=TRUE, show_colnames=TRUE,fontsize = 15)
dev.off()

dotplot(ck,font.size = 15)

cnetplot(ck,node_label="category",showCategory=5,cex_category=10,cex_label_category=1.5,colorEdge=TRUE)

ck <- compareCluster(geneCluster = list_EGI, fun = enrichKEGG, organism = 'taes')
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)
browseKEGG(ck, 'taes01250')

dotplot(ck,font.size = 15,showCategory=3)

cnetplot(ck,node_label="category",showCategory = 4,cex_category=10,cex_label_category=1.7,colorEdge=TRUE)

# For different expression profile, on expression level simplified categories:

list_diff<-list()
list_EGI_diff<-list()

cats<-c("A+","F_ABD","F_ADB")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile!="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list_diff[[1]]<-unlist(all_genes, recursive = FALSE)
list_EGI_diff[[1]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("B+","F_BAD","F_BDA")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile!="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list_diff[[2]]<-unlist(all_genes, recursive = FALSE)
list_EGI_diff[[2]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("D+","F_DAB","F_DBA")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile!="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list_diff[[3]]<-unlist(all_genes, recursive = FALSE)
list_EGI_diff[[3]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("A-","F_BDA","F_DBA")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile!="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list_diff[[4]]<-unlist(all_genes, recursive = FALSE)
list_EGI_diff[[4]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("B-","F_ADB","F_DAB")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile!="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list_diff[[5]]<-unlist(all_genes, recursive = FALSE)
list_EGI_diff[[5]]<-unlist(all_genes_EGI, recursive = FALSE)

cats<-c("D-","F_ABD","F_BAD")
all_genes<-list()
all_genes_EGI<-list()

for (i in 1:3){
  sub<-subset(all,all$profile!="S" & all$expression_level==cats[[i]])
  # sub<-subset(all,all$profile_cor=="S" & all$expression_level=="S")
  homeo<-homeo_ID[row.names(sub),]
  genes<-c(homeo[,1],homeo[,2],homeo[,3])
  genes<-unlist(lapply(genes,function(i){str_sub(i,end=-3)}))
  genes_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
  genes_EGI<-unlist(lapply(genes_EGI,function(i){str_sub(i,4)}))
  all_genes_EGI[[i]]<-genes_EGI
  all_genes[[i]]<-genes
}

list_diff[[6]]<-unlist(all_genes, recursive = FALSE)
list_EGI_diff[[6]]<-unlist(all_genes_EGI, recursive = FALSE)

same<-subset(all,all$profile!="S" & all$expression_level=="S")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
genes_same_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
genes_same_EGI<-unlist(lapply(genes_same_EGI,function(i){str_sub(i,4)}))

list_diff[[7]]<-genes_same
list_EGI_diff[[7]]<-genes_same_EGI

names(list_diff)<-c("A_high_diff","B_high_diff","D_high_diff","A_low_diff","B_low_diff","D_low_diff","Same_diff")
names(list_EGI_diff)<-c("A_high_diff","B_high_diff","D_high_diff","A_low_diff","B_low_diff","D_low_diff","Same_diff")

ck <- compareCluster(geneCluster = list_diff, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID", ont="BP")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)
res<-ck@compareClusterResult
mat<-as.data.frame.matrix(xtabs(-log10(p.adjust)~Description+Cluster,res,addNA=TRUE))
mat_filter<-mat[rowSums(mat > 3) >= 1, ]
png(filename = "Output/GO_comparison_heatmap_BP_big.png",
    width = 3500, height = 5600, units = "px", pointsize = 300,
    bg = "white", res = 200)
pheatmap(mat_filter, cluster_cols=TRUE, cluster_rows=TRUE, show_colnames=TRUE,fontsize = 16)
dev.off()

dotplot(ck,font.size = 15,showCategory=2)

cnetplot(ck,node_label="category",showCategory=2,cex_category=10,cex_label_category=1.5,colorEdge=TRUE)

ck <- compareCluster(geneCluster = list_EGI, fun = enrichKEGG, organism = 'taes')
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)
browseKEGG(ck, 'taes01250')

dotplot(ck,font.size = 15,showCategory=3)

cnetplot(ck,node_label="category",showCategory = 3,cex_category=10,cex_label_category=1.7,colorEdge=TRUE)


both_list<-c(list,list_diff)
ck <- compareCluster(geneCluster = both_list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID", ont="BP")
dotplot(ck,font.size = 15,showCategory=3)
res<-ck@compareClusterResult
mat<-as.data.frame.matrix(xtabs(-log10(p.adjust)~Description+Cluster,res,addNA=TRUE))
mat_filter<-mat[rowSums(mat > 3.9) >= 1, ]
pheatmap(mat_filter, cluster_cols=FALSE, cluster_rows=TRUE, show_colnames=TRUE)

# For different expression profile, on expression profile categories::

list_pro_diff<-list()
list_EGI_pro_diff<-list()

same<-subset(all,all$profile=="A")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
genes_same_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
genes_same_EGI<-unlist(lapply(genes_same_EGI,function(i){str_sub(i,4)}))
list_pro_diff[[1]]<-genes_same
list_EGI_pro_diff[[1]]<-genes_same_EGI

same<-subset(all,all$profile=="B")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
genes_same_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
genes_same_EGI<-unlist(lapply(genes_same_EGI,function(i){str_sub(i,4)}))
list_pro_diff[[2]]<-genes_same
list_EGI_pro_diff[[2]]<-genes_same_EGI

same<-subset(all,all$profile=="D")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
genes_same_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
genes_same_EGI<-unlist(lapply(genes_same_EGI,function(i){str_sub(i,4)}))
list_pro_diff[[3]]<-genes_same
list_EGI_pro_diff[[3]]<-genes_same_EGI

same<-subset(all,all$profile=="F")
homeo<-homeo_ID[row.names(same),]
genes_same<-c(homeo[,1],homeo[,2],homeo[,3])
genes_same<-unlist(lapply(genes_same,function(i){str_sub(i,end=-3)}))
genes_same_EGI<-c(homeo[,4],homeo[,5],homeo[,6])
genes_same_EGI<-unlist(lapply(genes_same_EGI,function(i){str_sub(i,4)}))
list_pro_diff[[4]]<-genes_same
list_EGI_pro_diff[[4]]<-genes_same_EGI

names(list_pro_diff)<-c("A","B","D","F")
names(list_EGI_pro_diff)<-c("A","B","D","F")

ck <- compareCluster(geneCluster = list_pro_diff, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID", ont="BP")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)
res<-ck@compareClusterResult
mat<-as.data.frame.matrix(xtabs(-log10(p.adjust)~Description+Cluster,res,addNA=TRUE))
mat_filter<-mat[rowSums(mat > 2) >= 1, ]
png(filename = "Output/GO_comparison_heatmap_expression_profile_cat_big.png",
    width = 3500, height = 5600, units = "px", pointsize = 400,
    bg = "white", res = 200)
pheatmap(mat_filter, cluster_cols=TRUE, cluster_rows=TRUE, show_colnames=TRUE,fontsize = 16)
dev.off()

cnetplot(ck,node_label="category",showCategory = 3,cex_category=10,cex_label_category=1.7,colorEdge=TRUE)
dotplot(ck,font.size = 15,showCategory=2)

## MCA:

all_mca<-join_all(list(profile_both_clean,expression_clean,length_clean[,1:2],iso_clean,chr_clean),by="ID")
rownames(all_mca)<-all_mca$ID
all_mca<-all_mca[ ,  !names(all_mca) %in% 
            c("ID")]
all_mca$expression_level <- as.factor(all_mca$expression_level)
all_mca$length <- as.factor(all_mca$length)
all_mca$isoform <- as.factor(all_mca$isoform)
all_mca$profile <- as.factor(all_mca$profile)
all_mca$global_expression_level <- as.factor(all_mca$global_expression_level)
all_mca$chromosome <- as.factor(all_mca$chromosome)

mca <- dudi.acm(all_mca, scannf = TRUE)
summary(mca)

plot(mca)
s.label(mca$li[,1:2])

plot(mca$li[,1:2])
s.label(mca$co,add=TRUE)
s.label(mca$co)

plot(mca$li[,2:3])
s.label(mca$co[,2:3],add=TRUE)
s.label(mca$co[,2:3])

s.corcircle(mca$co, 1, 2, clabel = 0.7)
boxplot(mca,3)
boxplot(mca,2)
boxplot(mca,1)
fviz_contrib(mca, choice = "var", axes = 1)
fviz_contrib(mca, choice = "var", axes = 2)
fviz_contrib(mca, choice = "var", axes = 3)

scatter(mca, col = c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30"))

score(mca, type = "boxplot")

# On subsets:

Same_expression_profile_mca<-Same_expression_profile[ ,  !names(Same_expression_profile) %in% 
                    c("profile_cor","global_expression_level","high_std")]

mca <- dudi.acm(Same_expression_profile_mca, scannf = TRUE)
summary(mca)

plot(mca)
s.label(mca$li[,1:2])

plot(mca$li[,1:2])
s.label(mca$co,add=TRUE)
s.label(mca$co)

plot(mca$li[,2:3])
s.label(mca$co[,2:3],add=TRUE)
s.label(mca$co[,2:3])

s.corcircle(mca$co, 1, 2, clabel = 0.7)
boxplot(mca,3)
fviz_contrib(mca, choice = "var", axes = 1)
fviz_contrib(mca, choice = "var", axes = 2)
fviz_contrib(mca, choice = "var", axes = 3)

scatter(mca, col = c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE", "#6C16B8", "#01A71D", "#28EC49", "#C7FDD0", "#6CCC78", "#4DBC8B", "#24F99B", "#F9EA30"))


## Try with quantitative values:

expression_pv<-read.table("Output/homeologs_pvalue.csv", sep = ",", header=TRUE, row.names = 1)
expression_pv_clean<-data.frame(cbind(expression_pv$ID,expression_pv$AB_pv,expression_pv$BD_pv,expression_pv$AD_pv))
colnames(expression_pv_clean)<-c("ID","AB_pv","BD_pv","AD_pv")

expression_logpv<-expression_pv_clean
expression_logpv$AB_logpv<--log10(expression_logpv$AB_pv)
expression_logpv$BD_logpv<--log10(expression_logpv$BD_pv)
expression_logpv$AD_logpv<--log10(expression_logpv$AD_pv)
expression_logpv<-expression_logpv[ ,  !names(expression_logpv) %in% 
                c("AB_pv","BD_pv","AD_pv")]

length_diff<-read.table("Output/homeo_length_diff.csv", sep = ",", header=TRUE, row.names = 1)
length_diff_clean<-data.frame(cbind(length_diff$ID,length_diff$AB_length_diff,length_diff$BD_length_diff,length_diff$AD_length_diff))
colnames(length_diff_clean)<-c("ID","AB_length_diff","BD_length_diff","AD_length_diff")

iso_diff<-read.table("Output/homeo_iso_diff.csv", sep = ",", header=TRUE, row.names = 1)
iso_diff_clean<-data.frame(cbind(iso_diff$ID,iso_diff$AB_iso_diff,iso_diff$BD_iso_diff,iso_diff$AD_iso_diff))
colnames(iso_diff_clean)<-c("ID","AB_iso_diff","BD_iso_diff","AD_iso_diff")

profile_cor<-read.table("Output/homeo_cor_category.csv", sep = ",", header=TRUE, row.names = 1)
profile_cor_coef_clean<-data.frame(cbind(profile_cor$ID,profile_cor$cor_AB,profile_cor$cor_BD,profile_cor$cor_AD))
colnames(profile_cor_coef_clean)<-c("ID","AB_cor","BD_cor","AD_cor")

profile_mr<-read.table("Output/homeo_mr_category.csv", sep = ",", header=TRUE, row.names = 1)
profile_mutualrank_clean<-data.frame(cbind(profile_mr$ID,profile_mr$AB_mr,profile_mr$BD_mr,profile_mr$AD_mr))
colnames(profile_mutualrank_clean)<-c("ID","AB_mr","BD_mr","AD_mr")

all_q<-join_all(list(expression_logpv,length_diff_clean,iso_diff_clean,profile_mutualrank_clean),by="ID")

rownames(all_q)<-all_q$ID
all_q<-all_q[ ,  !names(all_q) %in% 
            c("ID")]

all_q<-all_q %>% drop_na()

acp <- dudi.pca(all_q, center=TRUE, scale=TRUE, scann = TRUE)
summary(acp)
biplot(acp)
s.label(acp$l1)
plot(acp$l1[,1:2])
s.label(acp$co,add=TRUE)

explor(acp)

class_mr<-all$profile_mr[!is.na(all$profile_mr)]
class_exp<-all$expression_level[!is.na(all$expression_level)]

s.arrow(acp$co)

s.class(acp$l1[,1:2],class_mr,col=TRUE)

## Enrichement analysis on translocated homeologs:

list<-list()

A_4<-filter(chr_cat,special=="A_4")
A_4<-A_4 %>% drop_na()
genes_A4<-c(A_4$A,A_4$B,A_4$D)
genes_A4<-unlist(lapply(genes_A4,function(i){str_sub(i,end=-3)}))
list[[1]]<-genes_A4

A_5<-filter(chr_cat,special=="A_5")
A_5<-A_5 %>% drop_na()
genes_A5<-c(A_5$A,A_5$B,A_5$D)
genes_A5<-unlist(lapply(genes_A5,function(i){str_sub(i,end=-3)}))
list[[2]]<-genes_A5

N<-filter(chr_cat,special=="N")
N<-N %>% drop_na()
genes_N<-c(N$A,N$B,N$D)
genes_N<-unlist(lapply(genes_N,function(i){str_sub(i,end=-3)}))
list[[3]]<-genes_N

names(list)<-c("A_4","A_5","N")

ck <- compareCluster(geneCluster = list, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
head(ck)

dotplot(ck)

cnetplot(ck,node_label="category",cex_category=10)

## Summary plot for expression level vs Expression profile vs Global expression level:

all_test<-all

all_test<-all_test %>%
  mutate(test=case_when(
    expression_level == "S" & profile != "S" ~ "Diff_profile_Same_level",
    expression_level != "S" & profile != "S" ~ "Diff_profile_Diff_level",
    expression_level == "S" & profile == "S" ~ "Same_profile_Same_level",
    expression_level != "S" & profile == "S" ~ "Same_profile_Diff_level",
    TRUE ~ "other"))

all_test2<-all

all_test2<-all_test2 %>%
  mutate(test_profile=case_when(
    profile != "S" ~ "Different_expression_profile",
    profile == "S" ~ "Same_expression_profile",
    TRUE ~ "other"))

all_test2<-all_test2 %>%
  mutate(test_level=case_when(
    expression_level != "S" ~ "Different_expression_level",
    expression_level == "S" ~ "Same_expression_level",
    TRUE ~ "other"))

mosaicplot(~global_expression_level + test, data=all_test,shade = TRUE,cex.axis=1.2,ylab="Expression Profile Categories",xlab="Expression Level Categories",main=paste(paste("chi-value = ",round(test$statistic,4)),paste("p-value = ",round(test$p.value,4))))
tab<-table(all_test$global_expression_level,all_test$test)

merged<-merge(all_test,global_expression_cat, by.x=0, by.y="V1")

ggplot(merged, aes(x=test, y=means, fill = test)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#FF5733", "#FAA897", "#24D2F9", "#0EA5DE", "#CAA0EE")) +
  xlab('Expression Level Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

merged2<-merge(all_test2,global_expression_cat_n, by.x=0, by.y="V1")
merged2<-merged2 %>% drop_na()
diff_mean<-mean(filter(merged2,test_profile=="Different_expression_profile")$means)
same_mean<-mean(filter(merged2,test_profile=="Same_expression_profile")$means)

merged2<-merged2 %>%
  mutate(profile_mean=case_when(
    test_profile=="Different_expression_profile" ~ diff_mean,
    test_profile=="Same_expression_profile" ~ same_mean,
    TRUE ~ 0))


grp <- merged2 %>% group_by(test_profile,test_level) %>% tally()

ggplot(merged2) + aes(x = test_level, y = means) +
  geom_boxplot(colour = "black", fill = "steelblue",
               alpha = 0.2, outlier.shape = NA) +
  geom_jitter(colour = "steelblue", alpha = 0.4) +
  facet_grid(. ~ test_profile) + geom_hline(data = merged2, aes(yintercept = profile_mean), color="#FF5733")


dat_text <- data.frame(
  label = c("1 206","4 068", "1 920","9 380"),
  test_profile= c("Different_expression_profile", "Different_expression_profile","Same_expression_profile","Same_expression_profile"),
  test_level= c("Different_expression_level", "Same_expression_level","Different_expression_level", "Same_expression_level"),
  means = c(3,2.6,3.2,2.6)
) 

ggplot(merged2) + aes(x = test_level, y = means) +
    geom_violin(colour = "steelblue", fill="steelblue", alpha = 0.4) +
    facet_grid(. ~ test_profile) + geom_hline(data = merged2, aes(yintercept = profile_mean), color="#FF5733") +
    geom_boxplot(width=0.1) + xlab("Categories") + ylab("Mean Expression Level of the triad") + 
    geom_text(data = dat_text,label=dat_text$label,size =6) +
  theme(text = element_text(size=17))

diff_profile_diff_exp<-filter(merged2,test_profile=="Different_expression_profile",test_level=="Different_expression_level")
diff_profile_same_exp<-filter(merged2,test_profile=="Different_expression_profile",test_level=="Same_expression_level")
same_profile_diff_exp<-filter(merged2,test_profile=="Same_expression_profile",test_level=="Different_expression_level")
same_profile_same_exp<-filter(merged2,test_profile=="Same_expression_profile",test_level=="Same_expression_level")
diff_profile<-filter(merged2,test_profile=="Different_expression_profile")
same_profile<-filter(merged2,test_profile=="Same_expression_profile")
wilcox.test(diff_profile_diff_exp$means, diff_profile_same_exp$means)
wilcox.test(same_profile_diff_exp$means, same_profile_same_exp$means)
wilcox.test(same_profile$means, diff_profile$means)
