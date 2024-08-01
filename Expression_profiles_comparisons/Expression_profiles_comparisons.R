library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
library(emmeans)
library(stringr)
library(ade4)
library(adegraphics)
library(ggfortify)
library(lattice)
library(plotly)
library(bigstatsr)

## Load the data and compute the correlation matrix:
homeologs<-read.table("Output/homoeolog_Traes_EGI_length.csv", sep = ",", header=TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("X")]

expression_data<-read.table("Data/expression_data_TPM_UQ_batch_filtered.csv", sep = ",", header=TRUE, row.names = 1)
expression_data_gene_norm <- sweep(expression_data, 1, rowMeans(expression_data))
sum(is.na(expression_data_gene_norm))
sum(is.na(row.names(expression_data_gene_norm)))

homeo<-c(homeologs[,3],homeologs[,2],homeologs[,1])
homeo<-homeo[!is.na(homeo)]
sum(homeo=="NA")
all_genes<-row.names(expression_data_gene_norm)
sum(is.na(all_genes))
sum(all_genes=="NA")
diff<-setdiff(homeo,row.names(expression_data_gene_norm))

expression_data_gene_norm_filtered<-expression_data_gene_norm[homeo,]
expression_data_gene_norm_filtered_na <- expression_data_gene_norm_filtered[rowSums(is.na(expression_data_gene_norm_filtered)) > 0,]
expression_data_gene_norm_filtered<-expression_data_gene_norm_filtered %>% drop_na()
cor_matrix_filtered<-cor(t(expression_data_gene_norm_filtered),method = "spearman")

cor_mean<-lapply(cor_matrix_filtered,function(i){mean(order(i,decreasing = TRUE)[1:100])})
exp_mean<-rowMeans(expression_data_gene_norm_filtered)
plot(exp_mean,cor_mean)

sum(is.na(expression_data_gene_norm_filtered))
sum(is.na(cor_matrix_filtered[1:10000,1:10000]))

expression_data_gene_norm_fbm<-as_FBM(t(expression_data_gene_norm))
cor_matrix<-big_cor(expression_data_gene_norm_fbm)

cor_matrix_filtered[1,1:20]

## Using Spearman Correlation Coefficient:

N<-nrow(homeologs)
N<-8
cor_AB<-list()
cor_BD<-list()
cor_AD<-list()
cor<-list()
AB_pv<-list()
AD_pv<-list()
BD_pv<-list()


for (i in 1:N){
  print(i)
  hA<-homeologs[i,3]
  hB<-homeologs[i,2]
  hD<-homeologs[i,1]
  h<-c(hA,hB,hD)
  cor_AB[[i]]<-"NA"
  cor_BD[[i]]<-"NA"
  cor_AD[[i]]<-"NA"
  cor[[i]]<-"NA"
  AB_pv[[i]]<-"NA"
  AD_pv[[i]]<-"NA"
  BD_pv[[i]]<-"NA"
  if (sum(is.na(h))<=1){
    hA_expression_level<-t(expression_data_gene_norm[hA,])
    hB_expression_level<-t(expression_data_gene_norm[hB,])
    hD_expression_level<-t(expression_data_gene_norm[hD,])
    hA_expression_level_na<-sum(is.na(hA_expression_level))>0
    hB_expression_level_na<-sum(is.na(hB_expression_level))>0
    hD_expression_level_na<-sum(is.na(hD_expression_level))>0
    h_na<-c(hA_expression_level_na,hB_expression_level_na,hD_expression_level_na)
    if (sum(h_na)<=1){
      if (!hA_expression_level_na & !hB_expression_level_na){
        #qqnorm(hA_expression_level)
        test<-cor.test(hA_expression_level, hB_expression_level,  method = "spearman", na.action="na.omit")
        cor_AB[[i]]<-test$estimate
        AB_pv[[i]]<-test$p.value
      }
      if (!hD_expression_level_na & !hB_expression_level_na){
        # qqnorm(hA_expression_level)
        test<-cor.test(hB_expression_level, hD_expression_level,  method = "spearman", na.action="na.omit")
        cor_BD[[i]]<-test$estimate
        BD_pv[[i]]<-test$p.value
      }
      if (!hA_expression_level_na & !hD_expression_level_na){
        # qqnorm(hA_expression_level)
        test<-cor.test(hA_expression_level, hD_expression_level,  method = "spearman", na.action="na.omit")
        cor_AD[[i]]<-test$estimate
        AD_pv[[i]]<-test$p.value
      }
    }
  }
}

cor_AB<-as.numeric(cor_AB)
cor_BD<-as.numeric(cor_BD)
cor_AD<-as.numeric(cor_AD)
AB_pv<-as.numeric(AB_pv)
AD_pv<-as.numeric(AD_pv)
BD_pv<-as.numeric(BD_pv)

all_cor<-c(cor_AB,cor_BD,cor_AD)
all_pv<-c(AB_pv,BD_pv,AD_pv)

# Plots:
hist(all_cor)
hist(all_pv)

hist(cor_AB)
hist(cor_BD)
hist(cor_AD)

res<-data.frame(cbind(homeologs,cor_AB,cor_BD,cor_AD))

p1<-ggplot(res, aes(x=cor_AB))+
  geom_histogram(color="black",fill="#FF5733", bins = 50) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.2, 1.1)) + xlab("Correlation coefficients between A and B") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 0.6, 
               color = "black", size=1)

p2<-ggplot(res, aes(x=cor_BD))+
  geom_histogram(color="black",fill="#0EA5DE",bins = 50)+
  scale_x_continuous(expand = c(0, 0), limits = c(-0.2, 1.1)) + xlab("Correlation coefficients between B and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 0.6, 
             color = "black", size=1)

p3<-ggplot(res, aes(x=cor_AD))+
  geom_histogram(color="black",fill="#6C16B8",bins = 50) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.2, 1.1)) + xlab("Correlation coefficients between A and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 0.6, 
             color = "black", size=1)

aplot::plot_list(p1, p2, p3, labels=c('AB','BD','AD'), nrow=3)

hist(means_A)
hist(means_B)
hist(means_D)

res<-read.table("Output/homeo_expression_profiles_cor.csv", sep = ",", header=TRUE)

cor_list<-c(cor_AB,cor_BD,cor_AD,cor_BD_rand)
homeo<-c(rep('AB',length(cor_AB)),rep('BD',length(cor_BD)),rep('AD',length(cor_AD)),rep('Random',length(cor_BD_rand)))
data<-data.frame(unlist(cbind(cor_list, homeo)))
data$cor_list <- as.numeric(data$cor_list)
data$homeo <- factor(data$homeo)
data<-data %>% drop_na()
mu <- ddply(data, "homeo", summarise, grp.mean=mean(cor_list))
p<-ggplot(data, aes(x=cor_list, fill=homeo)) +
  geom_density(alpha=0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=homeo),
             linetype="dashed", size = 1.1 ) +
  xlab("Spearman Correlation Coefficient") + ylab("Density") + guides(fill=guide_legend(title="Comparison"),color=guide_legend(title="Comparison")) +
  theme(text = element_text(size=18))
p

# Statistical tests:
homeo<-c(rep('AB',length(cor_AB)),rep('BD',length(cor_BD)),rep('AD',length(cor_AD)),rep('Random',length(all_cor_rand)))
cor<-c(cor_AB,cor_BD,cor_AD,all_cor_rand)
data<-data.frame(cbind(homeo,cor))
data<-data %>% drop_na()
homeo<-as.factor(data$homeo)
cor<-as.numeric(data$cor)
lm<-lm(cor~homeo)
summary(lm)
anova(lm)


mean(cor_AB,na.rm = TRUE)
mean(cor_AD,na.rm = TRUE)
mean(cor_BD,na.rm = TRUE)

kruskal.test(cor, homeo) 
pairwise.wilcox.test(cor, homeo)

sum(is.na(res))

write.csv(res,
          file = "Output/homeo_expression_profiles_cor.csv")


## Using Mutual Rank:

N<-nrow(homeologs)
AB_mr<-list()
AD_mr<-list()
BD_mr<-list()


for (i in 1:N){
  print(i)
  hA<-homeologs[i,3]
  hB<-homeologs[i,2]
  hD<-homeologs[i,1]
  h<-c(hA,hB,hD)
  AB_mr[[i]]<-"NA"
  AD_mr[[i]]<-"NA"
  BD_mr[[i]]<-"NA"
  if (sum(is.na(h))<=1){
    hA_cor_list<-tryCatch({row.names(cor_matrix_filtered)[order(cor_matrix_filtered[hA,],decreasing = TRUE)]},error=function(e){NA})
    hB_cor_list<-tryCatch({row.names(cor_matrix_filtered)[order(cor_matrix_filtered[hB,],decreasing = TRUE)]},error=function(e){NA})
    hD_cor_list<-tryCatch({row.names(cor_matrix_filtered)[order(cor_matrix_filtered[hD,],decreasing = TRUE)]},error=function(e){NA})
    hA_cor_list_na<-sum(!is.na(hA_cor_list))>10
    hB_cor_list_na<-sum(!is.na(hB_cor_list))>10
    hD_cor_list_na<-sum(!is.na(hD_cor_list))>10
    if (hA_cor_list_na & hB_cor_list_na){
        rankB<-match(hB,hA_cor_list)
        rankA<-match(hA,hB_cor_list)
        AB_mr[[i]]<-sqrt(rankA*rankB)
    }
    if (hD_cor_list_na & hB_cor_list_na){
      rankB<-match(hB,hD_cor_list)
      rankD<-match(hD,hB_cor_list)
      BD_mr[[i]]<-sqrt(rankB*rankD)
    }
    if (hA_cor_list_na & hD_cor_list_na){
        rankD<-match(hD,hA_cor_list)
        rankA<-match(hA,hD_cor_list)
        AD_mr[[i]]<-sqrt(rankA*rankD)
    }
  }
}

AB_mr<-as.numeric(AB_mr)
AD_mr<-as.numeric(AD_mr)
BD_mr<-as.numeric(BD_mr)

all_mr<-c(AB_mr,BD_mr,AD_mr)
all_log_mr<-c(log10(AB_mr),log10(BD_mr),log10(AD_mr))
hist(all_mr)

res_mr<-data.frame(cbind(homeologs,AB_mr,BD_mr,AD_mr))

sum(is.na(res_mr))

AB_log_mr<-log10(AB_mr)
BD_log_mr<-log10(BD_mr)
AD_log_mr<-log10(AD_mr)

# Plots:

hist(AB_log_mr)
hist(BD_log_mr)
hist(AD_log_mr)

res_mr_log<-data.frame(cbind(homeologs,AB_log_mr,BD_log_mr,AD_log_mr))

p1<-ggplot(res_mr_log, aes(x=AB_log_mr))+
  geom_histogram(color="black",fill="#FF5733",bins = 50) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) + xlab("Mutual Rank between A and B") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 1, 
             color = "black", size=1)

p2<-ggplot(res_mr_log, aes(x=BD_log_mr))+
  geom_histogram(color="black",fill="#0EA5DE",bins = 50)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) + xlab("Mutual Rank between B and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 1, 
             color = "black", size=1)

p3<-ggplot(res_mr_log, aes(x=AD_log_mr))+
  geom_histogram(color="black",fill="#6C16B8",bins = 50) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) + xlab("Mutual Rank between A and D") + theme(text = element_text(size=15)) +
  geom_vline(xintercept = 1, 
             color = "black", size=1)

aplot::plot_list(p1, p2, p3, labels=c('AB','BD','AD'), nrow=3)


res_mr_log<-read.table("Output/homeo_expression_profiles_log_mutualrank.csv", sep = ",", header=TRUE)

mr_list<-c(res_mr_log$AB_log_mr,res_mr_log$BD_log_mr,res_mr_log$AD_log_mr)
homeo<-c(rep('AB',length(res_mr_log$AB_log_mr)),rep('BD',length(res_mr_log$BD_log_mr)),rep('AD',length(res_mr_log$AD_log_mr)))
data<-data.frame(unlist(cbind(mr_list, homeo)))
data$mr_list <- as.numeric(data$mr_list)
data$homeo <- factor(data$homeo)
data<-data %>% drop_na()
mu <- ddply(data, "homeo", summarise, grp.mean=mean(mr_list))
p<-ggplot(data, aes(x=mr_list, fill=homeo)) +
  geom_density(alpha=0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=homeo),
             linetype="dashed", size = 1.1 ) +
  xlab("Log Mutual Rank") + ylab("Density") + guides(fill=guide_legend(title="Comparison"),color=guide_legend(title="Comparison")) +
  theme(text = element_text(size=18))
p


# Statistical tests:
homeo<-c(rep('AB',length(AB_log_mr)),rep('BD',length(BD_log_mr)),rep('AD',length(AD_log_mr)))
mr<-c(AB_log_mr,BD_log_mr,AD_log_mr)
data<-data.frame(cbind(homeo,mr))
data<-data %>% drop_na()
homeo<-as.factor(data$homeo)
mr<-as.numeric(data$mr)
lm<-lm(mr~homeo)
summary(lm)
anova(lm)

mean(AB_log_mr,na.rm = TRUE)
mean(BD_log_mr,na.rm = TRUE)
mean(AD_log_mr,na.rm = TRUE)

kruskal.test(mr, homeo) 
pairwise.wilcox.test(mr, homeo)

res_all<-data.frame(cbind(homeologs,cor_AB,cor_BD,cor_AD,AB_log_mr,BD_log_mr,AD_log_mr))

write.csv(res_mr,
          file = "Output/homeo_expression_profiles_mutualrank.csv")

write.csv(res_mr_log,
          file = "Output/homeo_expression_profiles_log_mutualrank.csv")

write.csv(res_all,
          file = "Output/homeo_expression_profiles_cor_mutualrank.csv")


## SCC and MR comparison:

plot(all_cor,all_log_mr,xlab = "Spearman Correlation Coefficient", ylab = "Mutual Rank")

all<-data.frame(cbind(homeologs,cor,mr))

all<-all %>% drop_na()

ggplot(all) + geom_point(aes(x = cor, y = mr, colour = cor > 0.6 & mr < 1 )) +
  scale_colour_manual(name = 'SCC > 0.6 and log MR < 0.7', values = setNames(c('red','black'),c(T, F))) +
  xlab("Spearman Correlation Coefficient (SCC)") + ylab("Log Mutual Rank (MR)") +
  theme(text = element_text(size=15)) +
  geom_vline(xintercept = 0.6, 
             color = "red", size=1) +
  geom_hline(yintercept = 1, 
             color = "red", size=1)


## MR on random genes:

x <- 1:50000
N<-10000
N<-nrow(homeologs)
AB_mr_rand<-list()
AD_mr_rand<-list()
BD_mr_rand<-list()

for (i in 1:N){
  print(i)
  h<-sample(x,3)
  hA<-row.names(cor_matrix_filtered)[h[1]]
  hB<-row.names(cor_matrix_filtered)[h[2]]
  hD<-row.names(cor_matrix_filtered)[h[3]]
    hA_cor_list<-row.names(cor_matrix_filtered)[order(cor_matrix_filtered[hA,],decreasing = TRUE)]
    hB_cor_list<-row.names(cor_matrix_filtered)[order(cor_matrix_filtered[hB,],decreasing = TRUE)]
    hD_cor_list<-row.names(cor_matrix_filtered)[order(cor_matrix_filtered[hD,],decreasing = TRUE)]
      rankB<-match(hB,hA_cor_list)
      rankA<-match(hA,hB_cor_list)
      AB_mr_rand[[i]]<-sqrt(rankA*rankB)
      rankB<-match(hB,hD_cor_list)
      rankD<-match(hD,hB_cor_list)
      BD_mr_rand[[i]]<-sqrt(rankB*rankD)
      rankD<-match(hD,hA_cor_list)
      rankA<-match(hA,hD_cor_list)
      AD_mr_rand[[i]]<-sqrt(rankA*rankD)
}

AB_mr_rand<-as.numeric(AB_mr_rand)
AD_mr_rand<-as.numeric(AD_mr_rand)
BD_mr_rand<-as.numeric(BD_mr_rand)

all_mr_rand<-c(AB_mr_rand,BD_mr_rand,AD_mr_rand)
all_log_mr_rand<-c(log10(AB_mr_rand),log10(BD_mr_rand),log10(AD_mr_rand))
hist(all_mr_rand)

## SCC on random genes:

x <- 1:100000

N<-nrow(homeologs)
N<-8
cor_AB_rand<-list()
cor_BD_rand<-list()
cor_AD_rand<-list()
cor_rand<-list()
AB_pv_rand<-list()
AD_pv_rand<-list()
BD_pv_rand<-list()


for (i in 1:N){
  print(i)
  h<-sample(x,3)
  hA<-row.names(expression_data_gene_norm)[h[1]]
  hB<-row.names(expression_data_gene_norm)[h[2]]
  hD<-row.names(expression_data_gene_norm)[h[3]]
  h<-c(hA,hB,hD)
  cor_AB_rand[[i]]<-"NA"
  cor_BD_rand[[i]]<-"NA"
  cor_AD_rand[[i]]<-"NA"
  cor_rand[[i]]<-"NA"
  AB_pv_rand[[i]]<-"NA"
  AD_pv_rand[[i]]<-"NA"
  BD_pv_rand[[i]]<-"NA"
  if (sum(is.na(h))<=1){
    hA_expression_level<-t(expression_data_gene_norm[hA,])
    hB_expression_level<-t(expression_data_gene_norm[hB,])
    hD_expression_level<-t(expression_data_gene_norm[hD,])
    hA_expression_level_na<-sum(is.na(hA_expression_level))>0
    hB_expression_level_na<-sum(is.na(hB_expression_level))>0
    hD_expression_level_na<-sum(is.na(hD_expression_level))>0
    h_na<-c(hA_expression_level_na,hB_expression_level_na,hD_expression_level_na)
    if (sum(h_na)<=1){
      if (!hA_expression_level_na & !hB_expression_level_na){
        #qqnorm(hA_expression_level)
        test<-cor.test(hA_expression_level, hB_expression_level,  method = "spearman", na.action="na.omit")
        cor_AB_rand[[i]]<-test$estimate
        AB_pv_rand[[i]]<-test$p.value
      }
      if (!hD_expression_level_na & !hB_expression_level_na){
        # qqnorm(hA_expression_level)
        test<-cor.test(hB_expression_level, hD_expression_level,  method = "spearman", na.action="na.omit")
        cor_BD_rand[[i]]<-test$estimate
        BD_pv_rand[[i]]<-test$p.value
      }
      if (!hA_expression_level_na & !hD_expression_level_na){
        # qqnorm(hA_expression_level)
        test<-cor.test(hA_expression_level, hD_expression_level,  method = "spearman", na.action="na.omit")
        cor_AD_rand[[i]]<-test$estimate
        AD_pv_rand[[i]]<-test$p.value
      }
    }
  }
}

cor_AB_rand<-as.numeric(cor_AB_rand)
cor_BD_rand<-as.numeric(cor_BD_rand)
cor_AD_rand<-as.numeric(cor_AD_rand)
AB_pv_rand<-as.numeric(AB_pv_rand)
AD_pv_rand<-as.numeric(AD_pv_rand)
BD_pv_rand<-as.numeric(BD_pv_rand)

all_cor_rand<-c(cor_AB_rand,cor_BD_rand,cor_AD_rand)
all_pv_rand<-c(AB_pv,BD_pv,AD_pv)
hist(all_cor_rand)
hist(all_pv_rand)

hist(cor_AB_rand)
hist(cor_BD_rand)
hist(cor_AD_rand)

mean(cor_AB_rand)
mean(cor_BD_rand)
mean(cor_AD_rand)
mean(all_cor_rand)

# QQplot:

plot(all_cor[order(all_cor)],all_cor_h[order(all_cor_h)])
abline(0,1)

plot(cor_AB[order(cor_AB)],cor_AB_h[order(cor_AB_h)])
abline(0,1)

plot(cor_BD[order(cor_BD)],cor_BD_h[order(cor_BD_h)])
abline(0,1)

plot(cor_AD[order(cor_AD)],cor_AD_h[order(cor_AD_h)])
abline(0,1)