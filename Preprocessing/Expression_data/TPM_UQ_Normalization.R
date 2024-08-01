library("edgeR")
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)


## Load data:
expression_data_merged<-read.table("Data/expression_data_merged.csv", sep = ",", header=TRUE)
length_table<-read.table("Data/Genes_length.csv", sep = ",", header=TRUE)[,2:3]
expression_data_len<-merge(expression_data_merged,length_table, by.x=0, by.y="ID", all.x=TRUE)
sum(is.na(expression_data_len$Length))
expression_data_len<-expression_data_len %>% drop_na()
len<-expression_data_len$Length
row.names(expression_data_len)<-expression_data_len$Row.names
expression_data<-expression_data_len[,2:1774]

## TPM:
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

expression_data_tpm <- tpm3(expression_data,len)

## Plots TPM:
data<-data.frame(mean=rowMeans(expression_data),length=expression_data_len$Length)
p1<-ggplot(data) +
  aes(x = length, y = mean) + 
  geom_point(size=3) + 
  theme(text = element_text(size=23)) +
  xlab('Genes Length') + ylab('Genes Mean Expression Level') + ggtitle("Before TPM")

data<-data.frame(mean=rowMeans(expression_data_tpm),length=expression_data_len$Length)
p2<-ggplot(data) +
  aes(x = length, y = mean) + 
  geom_point(size=3) + 
  theme(text = element_text(size=23)) +
  xlab('Genes Length') + ylab('') + ggtitle("After TPM")

cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])

cond1<-unname(expression_data_tpm[,"SRR4426270"])
cond2<-unname(expression_data_tpm[,"SRR10990695"])
s1<-sum(log10(cond1+1))
s2<-sum(log10(cond2+1))
expression<-c(cond1,cond2)
condition<-c(rep('SRR4426270',length(cond1)),rep('SRR10990695',length(cond2)))
data<-data.frame(cbind(condition,expression))
data<-data %>% drop_na()
data$condition<-factor(data$condition)
data$expression<-as.numeric(data$expression)
str(data)
boxplot(log10(data$expression+1)~data$condition)

x <- melt(expression_data_tpm[,1:20])

plt <- ggplot(data = x, aes(x = Var2, y = log10(value+1)))
plt + geom_boxplot() + theme_minimal() + labs(x = "Title", y = "x")

## UQ:

# Code by Alexis Vandenbon

# get the UQ of the NON-ZERO values
uqs <- apply(expression_data_tpm, 2, function(x) quantile(x[x!=0],0.75))
mean.uq <- exp(mean(log(uqs))) # the geometric mean
correction.factors <- uqs/mean.uq
correction.factors <- mean.uq/uqs # Modified UQ

# make a DGEList object from the counts
dat <- DGEList(counts=expression_data_tpm)
dat$samples[,3] <- correction.factors

# return the normalized tag counts
dat<-cpm(dat)

# set a small pseudo count
pseudo <- c(dat[dat > 0], 0.01)
#pseudo <- c(expression_data_tpm[expression_data_tpm > 0], 0.01)

# and convert to log10 values
#dat.log10 <- log10(expression_data_tpm+(10*pseudo))
dat.log10 <- log10(dat+(10*pseudo))
dat.log10 <- log10(dat+1)


#dat.log10 <- log10(dat+0.01)

# return result
write.csv(dat.log10,
          file = "Data/expression_data_TPM_no_UQ_normalized.csv")

## Plots UQ:

new_uqs <- apply(dat, 2, function(x) quantile(x[x!=0],0.75))

cond1<-unname(dat[,"SRR4426270"])
cond2<-unname(dat[,"SRR10990695"])
s1<-sum(log10(cond1+1))
s2<-sum(log10(cond2+1))
expression<-c(cond1,cond2)
condition<-c(rep('SRR4426270',length(cond1)),rep('SRR10990695',length(cond2)))
data<-data.frame(cbind(condition,expression))
data<-data %>% drop_na()
data$condition<-factor(data$condition)
data$expression<-as.numeric(data$expression)
str(data)
boxplot(log10(data$expression+1)~data$condition)

x <- melt(dat.log10[,1:20])

plt <- ggplot(data = x, aes(x = Var2, y = value))
plt + geom_boxplot() + theme_minimal() + labs(x = "Title", y = "x")


