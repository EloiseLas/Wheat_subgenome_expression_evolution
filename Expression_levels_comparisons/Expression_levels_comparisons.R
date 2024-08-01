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
library(gridExtra) # for tableGrob
library(patchwork)

expression_data_merged<-read.table("Data/expression_data_merged.csv", sep = ",", header=TRUE)
expression_data_normalized<-read.table("Data/expression_data_normalized.csv", sep = ",", header=TRUE, row.names = 1)
expression_data_batch_normalized<-read.table("Data/expression_data_batch_normalized.csv", sep = ",", header=TRUE, row.names = 1)
expression_data_UQ_batch_normalized<-read.table("Data/expression_data_TPM_UQ_batch_filtered_new.csv", sep = ",", header=TRUE, row.names = 1)

## Expression profile observation:

#LOC100146089
#LOC123046553
#LOC123190289

A<-expression_data_UQ_batch_normalized["LOC100146089",]
B<-expression_data_UQ_batch_normalized["LOC123046553",]
D<-expression_data_UQ_batch_normalized["LOC123190289",]

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

tbl<-data.frame(SCC=c(sprintf("%.3f",0.6395956), sprintf("%.3f",0.40287748), sprintf("%.3f",0.40299060)),logMR = c(sprintf("%.3f",0.3010300), sprintf("%.3f",2.6399475), sprintf("%.3f",3.1842615)))
rownames(tbl)<-c("AB","BD","AD")
tbl1 <- tableGrob(tbl,theme = ttheme_default(base_size = 20))

a1 <- aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

plot_layout <- 
  "AD
   AD
   AD"  

wrap_plots(a1, tbl1, design = plot_layout, widths=c(0.8,0.2))


# 61 D Same_GO, 0.6395956, 0.40287748, 0.40299060, 0.3010300, 2.6399475, 3.1842615
#LOC123179837
#LOC123113145
#LOC123185427

A<-expression_data_UQ_batch_normalized["LOC123185427",]
B<-expression_data_UQ_batch_normalized["LOC123113145",]
D<-expression_data_UQ_batch_normalized["LOC123179837",]

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

tbl<-data.frame(SCC=c(sprintf("%.3f",0.6395956), sprintf("%.3f",0.40287748), sprintf("%.3f",0.40299060)),logMR = c(sprintf("%.3f",0.3010300), sprintf("%.3f",2.6399475), sprintf("%.3f",3.1842615)))
rownames(tbl)<-c("AB","BD","AD")
tbl1 <- tableGrob(tbl,theme = ttheme_default(base_size = 20))

a1 <- aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

plot_layout <- 
  "AD
   AD
   AD"  

wrap_plots(a1, tbl1, design = plot_layout, widths=c(0.8,0.2))


# 16991 F Low, 0.5288706 0.5738464 0.4496742 0.3890756 0.3010300 0.7385606
#LOC123107469
#LOC100037627
#LOC100037628

A<-expression_data_UQ_batch_normalized["LOC123107469",]
B<-expression_data_UQ_batch_normalized["LOC100037627",]
D<-expression_data_UQ_batch_normalized["LOC100037628",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13))+
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

# 3440 S Low

#LOC123184283
#LOC123040041
#LOC123048248

A<-expression_data_UQ_batch_normalized["LOC123184283",]
B<-expression_data_UQ_batch_normalized["LOC123040041",]
D<-expression_data_UQ_batch_normalized["LOC123048248",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 2))

aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

# 14607 S low 0.8014917 0.8664862 0.8659533 0.4771213 0.3010300 0.3890756 goood

#LOC123101968
#LOC123110085
#LOC123119084


A<-expression_data_UQ_batch_normalized["LOC123101968",]
B<-expression_data_UQ_batch_normalized["LOC123110085",]
D<-expression_data_UQ_batch_normalized["LOC123119084",]


plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('conditions') + ylab('expression') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('conditions') + ylab('expression') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('conditions') + ylab('expression') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

tbl<-data.frame(SCC=c(0.8014917, 0.8664862, 0.8659533),logMR = c(0.4771213, 0.3010300, 0.3890756))
rownames(tbl)<-c("AB","BD","AD")
tbl1 <- tableGrob(tbl,theme = ttheme_default(base_size = 20))

a1 <- aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

plot_layout <- 
  "AD
   AD
   AD"  

wrap_plots(a1, tbl1, design = plot_layout, widths=c(0.7,0.3))

# 6687 

#0.8273721, 0.8590279, 0.8375179
#0.4771213, 0.3010300, 0.3890756
#LOC123054467
#LOC123046601
#LOC123190332

A<-expression_data_UQ_batch_normalized["LOC123190332",]
B<-expression_data_UQ_batch_normalized["LOC123046601",]
D<-expression_data_UQ_batch_normalized["LOC123054467",]

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

a3<-aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

tbl<-data.frame(SCC=c(sprintf("%.3f",0.8273721), sprintf("%.3f",0.8590279), sprintf("%.3f",0.8375179)),logMR = c(sprintf("%.3f",0.4771213), sprintf("%.3f",0.3010300), sprintf("%.3f",0.3890756)))
rownames(tbl)<-c("AB","BD","AD")
tbl3 <- tableGrob(tbl,theme = ttheme_default(base_size = 20))


plot_layout <- 
  "AD
   AD
   AD"  

wrap_plots(a3, tbl3, design = plot_layout, widths=c(0.8,0.2))

# 7986

#LOC123073913
#LOC123064754
#LOC123057779

A<-expression_data_UQ_batch_normalized["LOC123057779",]
B<-expression_data_UQ_batch_normalized["LOC123064754",]
D<-expression_data_UQ_batch_normalized["LOC123073913",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

a4<-aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)
tbl<-data.frame(SCC=c(0.59475172, 0.5842469, 0.68100396),logMR = c(0.3890756, 0.4771213, 0.3010300))
rownames(tbl)<-c("AB","BD","AD")
tbl4 <- tableGrob(tbl,theme = ttheme_default(base_size = 20))

plot_layout <- 
  "AD
   AD
   AD"  

wrap_plots(a4, tbl4, design = plot_layout, widths=c(0.7,0.3))

# 52
#LOC123162243
#LOC123158088
#LOC123185119

A<-expression_data_UQ_batch_normalized["LOC123162243",]
B<-expression_data_UQ_batch_normalized["LOC123158088",]
D<-expression_data_UQ_batch_normalized["LOC123185119",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

# 87 0.1775509 0.49279909 0.58974320 4.0170784 2.6122665 0.3010300

#LOC123179881
#LOC123075751
#LOC123185785

A<-expression_data_UQ_batch_normalized["LOC123185785",]
B<-expression_data_UQ_batch_normalized["LOC123075751",]
D<-expression_data_UQ_batch_normalized["LOC123179881",]

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('Conditions') + ylab('Expression') +
  theme(text = element_text(size=20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

tbl<-data.frame(SCC=c(sprintf("%.3f", 0.1775509), sprintf("%.3f", 0.49279909), sprintf("%.3f", 0.58974320)),logMR = c(sprintf("%.3f", 4.0170784), sprintf("%.3f", 2.6122665), sprintf("%.3f", 0.3010300)))
rownames(tbl)<-c("AB","BD","AD")
tbl2 <- tableGrob(tbl,theme = ttheme_default(base_size = 20))

a2 <- aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

plot_layout <- 
  "AD
   AD
   AD"  

wrap_plots(a2, tbl2, design = plot_layout, widths=c(0.8,0.2))

# 146

#LOC123163100
#LOC123085148
#LOC123187040

A<-expression_data_UQ_batch_normalized["LOC123163100",]
B<-expression_data_UQ_batch_normalized["LOC123085148",]
D<-expression_data_UQ_batch_normalized["LOC123187040",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

# 8448

#LOC123077626 LOC123065074 LOC123058063
#0.4268148 0.2971992 0.2952766 0.3010300 0.7385606 0.5395906

A<-expression_data_UQ_batch_normalized["LOC123058063",]
B<-expression_data_UQ_batch_normalized["LOC123065074",]
D<-expression_data_UQ_batch_normalized["LOC123077626",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))


p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

# 9616 F Low
LOC123058494
LOC123065576
LOC101290678

A<-expression_data_UQ_batch_normalized["LOC123058494",]
B<-expression_data_UQ_batch_normalized["LOC123065576",]
D<-expression_data_UQ_batch_normalized["LOC101290678",]

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))

# 7517 A High 0.5521763 0.8597944 0.4552160 1.9715246 0.3010300 3.2850002

#LOC123191075
#LOC123047441
#LOC101290664

A<-expression_data_UQ_batch_normalized["LOC123191075",]
B<-expression_data_UQ_batch_normalized["LOC123047441",]
D<-expression_data_UQ_batch_normalized["LOC101290664",]

data<-data.frame(cbind(c(1:1436),unlist(A),unlist(B),unlist(D)))

p1 <- ggplot(data) + geom_point(aes(x = X1, y = X2), color = "#FF5733") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p2 <- ggplot(data) + geom_point(aes(x = X1, y = X3), color = "#0EA5DE") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

p3 <- ggplot(data) + geom_point(aes(x = X1, y = X4), color ="#6C16B8") +
  xlab('') + ylab('') +
  theme(text = element_text(size=13)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.9, 3))

aplot::plot_list(p1, p2, p3, labels=c('A','B','D'), nrow=3)

# 3463 B Low

#LOC123187043
#LOC123040100
#LOC123048271

A<-data.frame(t(expression_data_UQ_batch_normalized["LOC123187043",]))
B<-data.frame(t(expression_data_UQ_batch_normalized["LOC123040100",]))
D<-data.frame(t(expression_data_UQ_batch_normalized["LOC123048271",]))

plot(unlist(A))
plot(unlist(B))
plot(unlist(D))


p <- ggplot(A, aes(sample = LOC123187043)) + stat_qq(size=3) + stat_qq_line() + xlab("Normal theoretical quantiles") + ylab("Observed quantiles") + theme(text = element_text(size=25))

h<-ggplot(A, aes(x=LOC123187043))+
  geom_histogram(color="black",fill="#FF5733",bins = 50) +
  xlab("Normalized expression") + theme(text = element_text(size=25))

aplot::plot_list(p, h, labels=c('A','B'), nrow=1)

expression=unlist(c(A,B,D))

homeo=as.factor(c(rep('A',dim(A)[1]),rep('B',dim(B)[1]),rep('D',dim(D)[1])))
lm=lm(expression~homeo)

par(mfrow=c(2,2))
plot(lm)

# Other triad:

D<-data.frame(t(expression_data_UQ_batch_normalized["LOC100038334",]))
B<-data.frame(t(expression_data_UQ_batch_normalized["LOC543127",]))
A<-data.frame(t(expression_data_UQ_batch_normalized["LOC123060101",]))

expression=unlist(c(A,B,D))

hist(as.numeric(A))
hist(as.numeric(B))
hist(as.numeric(D))

p <- ggplot(A, aes(sample = LOC123060101)) + stat_qq(size=3) + stat_qq_line() + xlab("Normal theoretical quantiles") + ylab("Observed quantiles") + theme(text = element_text(size=25))


h<-ggplot(A, aes(x=LOC123060101))+
  geom_histogram(color="black",fill="#FF5733",bins = 50) +
  xlab("Normalized expression") + theme(text = element_text(size=25))


aplot::plot_list(p, h, labels=c('A','B'), nrow=1)

homeo=as.factor(c(rep('A',length(A)),rep('B',length(B)),rep('D',length(D))))
lm=lm(expression~homeo)

# Hypothesis check:

par(mfrow=c(1,1))
plot(lm)
plot(expression,lm$residuals)

# Pas de points aberrants: ok !
library(outliers)
res=residuals(lm)
grubbs.test(res, type=10) 

# n mesures sont indépendantes: non...
library(car)
durbinWatsonTest(lm)

# Moyenne a 0: ok !
mean(lm$residuals)

# Homoscédasticité: non...
bartlett.test(expression,homeo)
leveneTest(expression ~ homeo)

# Normalite: ok...
shapiro.test(lm$residuals)
hist(lm$residuals)


anova(lm)
lm=aov(expression~homeo-1)
test<-TukeyHSD(lm)


test<-pairwise.wilcox.test(expression, homeo)
str(test)
test$p.value



## Expression level comparison:

#expression_data_UQ_batch_normalized_sp <- t(apply(expression_data_UQ_batch_normalized, 1, function(i) i/sum(i)*100))
#expression_data_UQ_batch_normalized_sp <- data.frame(expression_data_UQ_batch_normalized_sp)

homeologs<-read.table("Output/homoeolog_Traes_EGI_length.csv", sep = ",", header=TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("X")]

N<-nrow(homeologs)
N<-1000
pvalue<-list()
means_A<-list()
means_B<-list()
means_D<-list()
means<-list()
# log_means_A<-list()
# log_means_B<-list()
# log_means_D<-list()
# log_means<-list()
AB_pv<-list()
AD_pv<-list()
BD_pv<-list()
AB_mean<-list()
AD_mean<-list()
BD_mean<-list()
not_found_in_expression_data<-list()

for (i in 1:N){
  print(i)
  hA<-homeologs[i,3]
  hB<-homeologs[i,2]
  hD<-homeologs[i,1]
  h<-c(hA,hB,hD)
  if ((sum(is.na(h))==3) | (sum(is.na(h))==2)){
    pvalue[[i]]<-"NA"
    means_A[[i]]<-"NA"
    means_B[[i]]<-"NA"
    means_D[[i]]<-"NA"
    means[[i]]<-"NA"
    # log_means_A[[i]]<-"NA"
    # log_means_B[[i]]<-"NA"
    # log_means_D[[i]]<-"NA"
    # log_means[[i]]<-"NA"
    AB_pv[[i]]<-"NA"
    AD_pv[[i]]<-"NA"
    BD_pv[[i]]<-"NA"
    AB_mean[[i]]<-"NA"
    AD_mean[[i]]<-"NA"
    BD_mean[[i]]<-"NA"
  } else {
    hA_expression_level<-t(expression_data_UQ_batch_normalized[hA,])
    hB_expression_level<-t(expression_data_UQ_batch_normalized[hB,])
    hD_expression_level<-t(expression_data_UQ_batch_normalized[hD,])
    hA_expression_level_na<-sum(is.na(hA_expression_level))>0
    hB_expression_level_na<-sum(is.na(hB_expression_level))>0
    hD_expression_level_na<-sum(is.na(hD_expression_level))>0
    h_na<-c(hA_expression_level_na,hB_expression_level_na,hD_expression_level_na)
    means_A[[i]]<-mean(hA_expression_level)
    means_B[[i]]<-mean(hB_expression_level)
    means_D[[i]]<-mean(hD_expression_level)
    # log_means_A[[i]]<-mean(log(hA_expression_level+1))
    # log_means_B[[i]]<-mean(log(hB_expression_level+1))
    # log_means_D[[i]]<-mean(log(hD_expression_level+1))
    expression=c(hA_expression_level,hB_expression_level,hD_expression_level)
    means[[i]]<-mean(expression,na.rm=TRUE)
    AB_mean[[i]]<-mean(c(hA_expression_level,hB_expression_level),na.rm=TRUE)
    AD_mean[[i]]<-mean(c(hA_expression_level,hD_expression_level),na.rm=TRUE)
    BD_mean[[i]]<-mean(c(hD_expression_level,hB_expression_level),na.rm=TRUE)
    # log_means[[i]]<-mean(log(expression+1))
    if ((sum(h_na))>=1){
        pvalue[[i]]<-"NA"
        AB_pv[[i]]<-"NA"
        AD_pv[[i]]<-"NA"
        BD_pv[[i]]<-"NA"
        # not_found_in_expression_data<-c(not_found_in_expression_data, c(hA,hB,hD))
    } else {
        homeo=as.factor(c(rep('A',length(hA_expression_level)),rep('B',length(hB_expression_level)),rep('D',length(hD_expression_level))))
        lm=lm(expression~homeo)
        pvalue[[i]]<-anova(lm)$"Pr(>F)"[1]
        
        # effect=emmeans(lm,~homeo)
        # AB_pv[[i]]<-summary(pairs(effect,adjust="tukey"))$p.value[1]
        # AD_pv[[i]]<-summary(pairs(effect,adjust="tukey"))$p.value[2]
        # BD_pv[[i]]<-summary(pairs(effect,adjust="tukey"))$p.value[3]
        
        # lm=aov(expression~homeo-1)
        # test<-TukeyHSD(lm)
        # AB_pv[[i]]<-test$homeo[,4][1]
        # AD_pv[[i]]<-test$homeo[,4][2]
        # BD_pv[[i]]<-test$homeo[,4][3]
        
        test<-pairwise.wilcox.test(expression, homeo)
        AB_pv[[i]]<-test$p.value[1,1]
        AD_pv[[i]]<-test$p.value[2,1]
        BD_pv[[i]]<-test$p.value[2,2]
    }
  }
}

pvalue<-as.numeric(pvalue)
means<-as.numeric(means)
means_A<-as.numeric(means_A)
means_B<-as.numeric(means_B)
means_D<-as.numeric(means_D)
# log_means<-as.numeric(log_means)
# log_means_A<-as.numeric(log_means_A)
# log_means_B<-as.numeric(log_means_B)
# log_means_D<-as.numeric(log_means_D)
AB_pv<-as.numeric(AB_pv)
AD_pv<-as.numeric(AD_pv)
BD_pv<-as.numeric(BD_pv)
AB_mean<-as.numeric(unlist(AB_mean))
AD_mean<-as.numeric(unlist(AD_mean))
BD_mean<-as.numeric(unlist(BD_mean))

# How many comparisons are significant ?
table(pvalue<=0.5)

hist(pvalue,breaks=50)
abline(v=0.05, col = "red")

## PLots:

res<-data.frame(cbind(means_A, means_B, means_D, pvalue, AB_pv, AD_pv, BD_pv))
res<-res %>% drop_na()
res<-res[!(res$means_A %in% "NA"),]
res<-res[!(res$means_B %in% "NA"),]
res<-res[!(res$means_D %in% "NA"),]

# res4<-data.frame(cbind(log_means_A, log_means_B, log_means_D, pvalue, AB_pv, AD_pv, BD_pv))
# res4<-res4 %>% drop_na()
# res4<-res4[!(res4$log_means_A %in% "NA"),]
# res4<-res4[!(res4$log_means_B %in% "NA"),]
# res4<-res4[!(res4$log_means_D %in% "NA"),]

logpv<--log(pvalue)

res2<-data.frame(cbind(logpv,means))
res2<-res2 %>% drop_na()

res3<-data.frame(cbind(logpv,log_means))
res3<-res3 %>% drop_na()


res9<-data.frame(cbind(c(-log10(AB_pv),-log10(AD_pv),-log10(BD_pv)),c(AB_mean,AD_mean,BD_mean)))
res9<-res9 %>% drop_na()

test<- res9[!duplicated(res9), ]

# log_AB_diff<-log_means_A-log_means_B
# log_AD_diff<-log_means_A-log_means_D
# log_BD_diff<-log_means_B-log_means_D

AB_diff<-means_A-means_B
AD_diff<-means_A-means_D
BD_diff<-means_B-means_D

# res5<-data.frame(cbind(log_AB_diff, log_AD_diff, log_BD_diff, AB_pv, AD_pv, BD_pv, log_means))
# res5<-res5 %>% drop_na()

res6<-data.frame(cbind(AB_diff, AD_diff, BD_diff, AB_pv, AD_pv, BD_pv, means))
res6<-res6 %>% drop_na()

ggplot(res2) +
  aes(x = means, y = logpv) + 
  geom_point() + 
  xlab('Global Expression Level') + ylab('-log10(pv)')

ggplot(res3) +
  aes(x = log_means, y = logpv) + 
  geom_point() + 
  scale_x_log10()

ggplot(res9) +
  aes(x = X2, y = X1) + 
  geom_point() + 
  xlab('Global Expression Level') + ylab('-log10(pv)')

plot(res9$X2,res9$X1)
plot(AD_mean,-log10(AD_pv))
plot(AB_mean,-log10(AB_pv))
plot(BD_mean,-log10(BD_pv))

# With 2 by 2 p-value:

p1<-ggplot(res) + geom_point(aes(x = means_A, y = means_B, colour = AB_pv < 0.05), size=2, alpha=0.7) +
  scale_colour_manual(name = 'AB_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('Average expression level of the homoelog in A') + ylab('Average expression level of the homoelog in B') +
  theme(text = element_text(size=23)) + guides(colour="none")

p2<-ggplot(res) + geom_point(aes(x = means_A, y = means_D, colour = AD_pv < 0.05), size=2, alpha=0.7) +
  scale_colour_manual(name = 'AD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('Average expression level of the homoelog in A') + ylab('Average expression level of the homoelog in D') +
  theme(text = element_text(size=23)) + guides(colour="none")

p3<-ggplot(res) + geom_point(aes(x = means_D, y = means_B, colour = BD_pv < 0.05), size=2, alpha=0.7) +
  scale_colour_manual(name = 'p-value < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('Average expression level of the homoelog in D') + ylab('Average expression level of the homoelog in B') +
  theme(text = element_text(size=23))

aplot::plot_list(p1, p2, p3, labels=c('A','B','C'), nrow=1)

# With 2 by 2 p-value log:

ggplot(res4) + geom_point(aes(x = log_means_A, y = log_means_B, colour = AB_pv < 0.05)) +
  scale_colour_manual(name = 'AB_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means_A') + ylab('means_B')

ggplot(res4) + geom_point(aes(x = log_means_A, y = log_means_B, colour = AB_pv < 0.05)) +
  scale_colour_manual(name = 'AB_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means_A') + ylab('means_B') +
  scale_x_log10() + scale_y_log10()

ggplot(res4) + geom_point(aes(x = log_means_A, y = log_means_D, colour = AD_pv < 0.05)) +
  scale_colour_manual(name = 'AD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('log_means_A') + ylab('log_means_D')

ggplot(res4) + geom_point(aes(x = log_means_D, y = log_means_B, colour = BD_pv < 0.05)) +
  scale_colour_manual(name = 'BD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('log_means_D') + ylab('log_means_B')

# With global p-value:

ggplot(res) + geom_point(aes(x = means_A, y = means_B, colour = pvalue < 0.05)) +
  scale_colour_manual(name = 'pvalue < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means_A') + ylab('means_B')

ggplot(res) + geom_point(aes(x = means_A, y = means_D, colour = pvalue < 0.05)) +
  scale_colour_manual(name = 'pvalue < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means_A') + ylab('means_D')

ggplot(res) + geom_point(aes(x = means_D, y = means_B, colour = pvalue < 0.05)) +
  scale_colour_manual(name = 'pvalue < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means_D') + ylab('means_B')

# Expression log difference:

plot(res5$log_means,res5$log_AB_diff,xlim=c(0,10), ylim=c(-8,8))

ggplot(res5) +
  aes(x = log_means, y = log_AB_diff) + 
  geom_point() + 
  scale_x_log10 ()

ggplot(res5) + geom_point(aes(x = log_means, y = log_AB_diff, colour = AB_pv < 0.05)) +
  scale_colour_manual(name = 'AB_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('log_means') + ylab('log_AB_diff') + 
  scale_x_log10 ()

ggplot(res5) + geom_point(aes(x = log_means, y = log_AD_diff, colour = AD_pv < 0.05)) +
  scale_colour_manual(name = 'AD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('log_means') + ylab('log_AD_diff') + 
  scale_x_log10 ()

ggplot(res5) + geom_point(aes(x = log_means, y = log_BD_diff, colour = BD_pv < 0.05)) +
  scale_colour_manual(name = 'BD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('log_means') + ylab('log_BD_diff') + 
  scale_x_log10 ()

# Expression difference:

ggplot(res6) + geom_point(aes(x = means, y = AB_diff, colour = AB_pv < 0.05)) +
  scale_colour_manual(name = 'AB_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means') + ylab('AB_diff') + 
  scale_x_log10 ()

ggplot(res6) + geom_point(aes(x = means, y = AD_diff, colour = AD_pv < 0.05)) +
  scale_colour_manual(name = 'AD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means') + ylab('AD_diff') + 
  scale_x_log10 ()

ggplot(res6) + geom_point(aes(x = means, y = BD_diff, colour = BD_pv < 0.05)) +
  scale_colour_manual(name = 'BD_pv < 0.05', values = setNames(c('red','black'),c(T, F))) +
  xlab('means') + ylab('BD_diff') + 
  scale_x_log10 ()

# Some expression pattern plot:

homeo_to_test1=c(5,19,43,272)
homeo_to_test3=c(10536,7321,14389,22827)
homeo_to_test2=c(31,384,468,750)
par(mfrow=c(2,2))
for (i in homeo_to_test3){
  triads<-homeologs[i,4:6]
  hA_expression_level<-t(expression_data_UQ_batch_normalized[triads[1,1],])
  hB_expression_level<-t(expression_data_UQ_batch_normalized[triads[1,2],])
  hD_expression_level<-t(expression_data_UQ_batch_normalized[triads[1,3],])
  expression=c(hA_expression_level,hB_expression_level,hD_expression_level)
  conditions=c(colnames(expression_data_UQ_batch_normalized),colnames(expression_data_UQ_batch_normalized),colnames(expression_data_UQ_batch_normalized))
  homeo=as.factor(c(rep('A',length(hA_expression_level)),rep('B',length(hB_expression_level)),rep('D',length(hD_expression_level))))
  interaction.plot(homeo, conditions, expression)
  #boxplot(expression~homeo,col=c("darkblue","darkgreen","darkred"),ylab="expression",xlab="homeologs",sub=i)
}

# PCA ?

data<-data.frame(cbind(means_A, means_B, means_D))
data<-data %>% drop_na()

acp <- dudi.pca(data, center=TRUE, scale=TRUE, scann = TRUE)
summary(acp)

s.label(acp$l1)
s.arrow(acp$co)

biplot(acp)
scatter(acp, posieig="bottomright")
scatter(acp, xax=2, yax=3)

pca <- prcomp(data, scale. = FALSE)
summary(pca)
biplot(pca)
autoplot(pca,loadings = TRUE,loadings.label=1)

plot_ly(x=means_A, y=means_B, z=means_D, type="scatter3d", mode="markers")

# Dataframe to categorize:

homeo_pvalue<-data.frame(cbind(homeologs,AB_pv, AD_pv, BD_pv,means_A, means_B, means_D))
homeo_pvalue<-homeo_pvalue %>% drop_na()
write.csv(homeo_pvalue,
          file = "Output/homeologs_pvalue.csv")

# Selection of triads with significantly different expression profiles:

homeo_pvalue<-data.frame(cbind(homeologs,pvalue))
homeo_pvalue<-homeo_pvalue %>% drop_na()
diff_expression<-homeo_pvalue[(homeo_pvalue$pvalue<0.025),]
write.csv(diff_expression,
          file = "Output/differently_expressed_homeo_UQ_batch_norm.csv")

diff_expression<-read.table("Output/differently_expressed_homeo_UQ_batch_norm.csv", sep = ",", header=TRUE, row.names = 1)

diff_homeologs<-c(diff_expression[,1],diff_expression[,2],diff_expression[,3])
diff_homeologs <- lapply(diff_homeologs, function(i){str_sub(i, end=-3)})
lapply(diff_homeologs, write, "Output/diff_homeologs_ID2.txt", append=TRUE)
write(paste(diff_homeologs, collapse = ','), 'test.txt')

# Selection of triads with really similar expression profiles:

homeo_pvalue<-data.frame(cbind(homeologs,pvalue,means))
homeo_pvalue<-homeo_pvalue %>% drop_na()
sim_expression<-homeo_pvalue[(homeo_pvalue$pvalue>0.8),]
write.csv(sim_expression,
          file = "Output/similarly_expressed_homeo_UQ_batch_norm.csv")

# Comparison of triad and duplet mean expression level:

homeo_mean<-data.frame(cbind(homeologs,means_A, means_B, means_D, means))
homeo_mean$na<-rowSums(is.na(homeo_mean[,12:14]))
duplet<-filter(homeo_mean, na==1)
triad<-subset(homeo_mean, na==0)
uno<-subset(homeo_mean, na==2)
duplet_means<-duplet$means
triad_means<-triad$means
uno_means<-uno$means
t.test(duplet_means, triad_means, var.eq=T)
t.test(uno_means, duplet_means, var.eq=T)

# Check independency between length and mean expression:

homeo_mean<-c(means_A, means_B, means_D)
homeo_length<-c(homeologs$length_A, homeologs$length_A, homeologs$length_D)
plot(homeo_length,homeo_mean, type="p",main="Average expression level of the homeologs in function of their length", xlab="Length", ylab="Average Expression")

# Global expression level category:

homeo_global_mean<-data.frame(cbind(homeologs$ID, means))
homeo_global_mean<-homeo_global_mean %>% drop_na()

hist(means,col = 'skyblue3', breaks = 50)

# homeo_global_mean <- homeo_global_mean %>%
#   mutate(category=case_when(
#     means>=1 ~ "High",
#     means<1 & means>=0 ~ "Medium",
#     means<0 ~ "Low",
#     TRUE ~ "other"
#   ))

homeo_global_mean <- homeo_global_mean %>% 
  mutate(category = cut_number(means,
                               n = 3,
                               right = F,
                               labels = c("Low","Medium","High")))

write.csv(homeo_global_mean,
          file = "Output/homeo_global_expression_level_category.csv")

# Comparison of expression average between the 3 subgenomes:

homeo<-c(rep('A',length(means_A)),rep('B',length(means_B)),rep('D',length(means_D)))
mean_list<-c(means_A, means_B, means_D)
data<-data.frame(cbind(homeo,means))
data<-data %>% drop_na()
homeo<-as.factor(data$homeo)
means<-as.numeric(data$means)
lm<-lm(means~homeo)
summary(lm)
anova(lm)

hist(means_A)
hist(means_B)
hist(means_D)

mean(means_A, na.rm = TRUE)
mean(means_B, na.rm = TRUE)
mean(means_D, na.rm = TRUE)

kruskal.test(means, homeo) 
pairwise.wilcox.test(means, homeo, p.adjust.method = "BH")

# Means stdev plot:

data<-data.frame(cbind(means_A, means_B, means_D, means))
data$stdev<-apply(data[,1:3], 1, sd, na.rm=TRUE)
ggplot(data) + geom_point(aes(x = means, y = stdev)) +
  xlab('means') + ylab('std')

# Plot mean for each homeo depending on global mean:

global_mean_list<-c(means, means, means)
data<-data.frame(unlist(cbind(mean_list, homeo, global_mean_list)))
data$mean_list <- as.numeric(data$mean_list)
data$global_mean_list <- as.numeric(data$global_mean_list)
data$homeo <- factor(data$homeo)
data<-data %>% drop_na()
ggplot(data) + geom_point(aes(x = global_mean_list, group=homeo, y = mean_list, color = homeo)) +
  xlab('Global mean for the triad') + ylab('Mean for each subgenome') + 
  geom_smooth(method='lm', formula= y~x,aes(x = global_mean_list, group=homeo, y = mean_list, color = homeo))

ggplot(data) +
  xlab('Global mean for the triad') + ylab('Mean for each subgenome') + 
  geom_smooth(method='lm', formula= y~x,aes(x = global_mean_list, group=homeo, y = mean_list, color = homeo))



