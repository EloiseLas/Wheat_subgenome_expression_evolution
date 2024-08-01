library(AnnotationHub)
library(GOSemSim)
library(AnnotationForge)
library(stringr)
library(dplyr)
library(plyr)
library(tidyr)
library(tidyverse)

## With AnnotationHub
hub <- AnnotationHub()
q <- query(hub, "Triticum Aestivum")
id <- q$ah_id[length(q)]
Taes <- hub[[id]]

taGO <- godata(Taes, ont="BP")
# As usual not enough annotation

## Makes an organism package:

annotation_data<-read.table("Data/iwgsc_refseqv2.1_functional_annotation.csv", sep = ",", header=TRUE, row.names = NULL)
annotation_data_all<-subset(annotation_data, f.type=="Gene Ontology")
annotation_data_all<-data.frame(cbind(annotation_data_all$g2.identifier,annotation_data_all$f.name, rep("IEA",length(annotation_data_all$f.name))))
colnames(annotation_data_all)<-c("GID","GO","EVIDENCE")
annotation_data_BP<-subset(annotation_data_all, f.domain=="BP")
annotation_data_BP<-data.frame(cbind(annotation_data_BP$g2.identifier,annotation_data_BP$f.name, rep("IEA",length(annotation_data_BP$f.name))))
colnames(annotation_data_BP)<-c("GID","GO","EVIDENCE")
annotation_data_CC<-subset(annotation_data_all, f.domain=="CC")
annotation_data_CC<-data.frame(cbind(annotation_data_CC$g2.identifier,annotation_data_CC$f.name, rep("IEA",length(annotation_data_CC$f.name))))
colnames(annotation_data_CC)<-c("GID","GO","EVIDENCE")
annotation_data_MF<-subset(annotation_data_all, f.domain=="MF")
annotation_data_MF<-data.frame(cbind(annotation_data_MF$g2.identifier,annotation_data_MF$f.name, rep("IEA",length(annotation_data_MF$f.name))))
colnames(annotation_data_MF)<-c("GID","GO","EVIDENCE")

gene_info<-read.table("Data/Traes_2_EGI.csv", sep = ",", header=TRUE, row.names = NULL)
gene_info<-data.frame(cbind(gene_info$Traes_long,gene_info$EGI))
colnames(gene_info) <- c("GID","SYMBOL")
gene_info<-gene_info %>% drop_na()

makeOrgPackage(gene_info=gene_info, go=annotation_data_all,
               version="0.1",
               maintainer="Some One <so@someplace.org>",
               author="Some One <so@someplace.org>",
               outputDir = ".",
               tax_id="4565",
               genus="Triticum",
               species="aestivum",
               goTable="go")

install.packages("./org.Taestivum.eg.db",repos = NULL, type="source")

taGO <- godata("org.Taestivum.eg.db", ont="BP", keytype = "GID")

# Example:

# finchFile <- system.file("extdata","finch_info.txt",
#                          package="AnnotationForge")
# finch <- read.table(finchFile,sep="\t")
# ## Now prepare some data.frames
# fSym <- finch[,c(2,3,9)]
# fSym <- fSym[fSym[,2]!="-",]
# fSym <- fSym[fSym[,3]!="-",]
# colnames(fSym) <- c("GID","SYMBOL","GENENAME")
# 
# fChr <- finch[,c(2,7)]
# fChr <- fChr[fChr[,2]!="-",]
# colnames(fChr) <- c("GID","CHROMOSOME")
# 
# finchGOFile <- system.file("extdata","GO_finch.txt",
#                            package="AnnotationForge")
# fGO <- read.table(finchGOFile,sep="\t")
# fGO <- fGO[fGO[,2]!="",]
# fGO <- fGO[fGO[,3]!="",]
# colnames(fGO) <- c("GID","GO","EVIDENCE")
# 
# ## Then call the function
# makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
#                version="0.1",
#                maintainer="Some One <so@someplace.org>",
#                author="Some One <so@someplace.org>",
#                outputDir = ".",
#                tax_id="59729",
#                genus="Taeniopygia",
#                species="guttata",
#                goTable="go")
# 
# ## then you can call install.packages based on the return value
# install.packages("./org.Tguttata.eg.db",repos = NULL, type="source")
# 
# tgGO <- godata("org.Tguttata.eg.db", ont="BP", keytype = "GENEID")

## GO similarity measurement between homeologs annotations:

homeologs<-read.table("Output/homoeolog_Traes_EGI_length.csv", sep = ",", header=TRUE)
homeologs<-homeologs[,5:8]

N<-nrow(homeologs)
N<-10
scores_AB<-list()
scores_BD<-list()
scores_AD<-list()

for (i in 1:N){
  print(i)
  hA<-strsplit(homeologs[i,3], '[.]')[[1]][1]
  hB<-strsplit(homeologs[i,2], '[.]')[[1]][1]
  hD<-strsplit(homeologs[i,1], '[.]')[[1]][1]
  res<-GOSemSim::geneSim(hA, hB, semData=taGO, measure="Wang", combine="BMA")
  if (is.atomic(res)){
    scores_AB[[i]]="NA"
  } else {
    scores_AB[[i]]=res$geneSim
  }
  res<-GOSemSim::geneSim(hD, hB, semData=taGO, measure="Wang", combine="BMA")
  if (is.atomic(res)){
    scores_BD[[i]]="NA"
  } else {
    scores_BD[[i]]=res$geneSim
  }
  res<-GOSemSim::geneSim(hA, hD, semData=taGO, measure="Wang", combine="BMA")
  if (is.atomic(res)){
    scores_AD[[i]]="NA"
  } else {
    scores_AD[[i]]=res$geneSim
  }
}

sum(is.na(scores_AB))
sum(is.na(scores_AD))
sum(is.na(scores_BD))

scores_AB<-as.numeric(scores_AB)
scores_AD<-as.numeric(scores_AD)
scores_BD<-as.numeric(scores_BD)

hist(scores_BD)

res<-cbind(homeologs,scores_AB,scores_AD,scores_BD)
res$cat<-rowMeans(res[,5:7])
same_GO<-filter(res,scores_AB!="NA",scores_BD!="NA",scores_AD!="NA")

write.csv(res,
          file = "Output/homeo_annotation_GO_cat2.csv")

write.csv(same_GO,
          file = "Output/homeo_same_GO2.csv")

GOSemSim::geneSim("TraesCS1B03G0213500", "TraesCS1B03G0216400", semData=taGO, measure="Wang", combine="BMA")

# Test on triad 60:
TraesCS1D03G0028500.1
TraesCS5B03G0933900.1
TraesCS1A03G0037900.1
GOSemSim::geneSim("TraesCS1D03G0028500", "TraesCS5B03G0933900", semData=taGO, measure="Wang", combine="BMA")
GOSemSim::geneSim("TraesCS5B03G0933900", "TraesCS1A03G0037900", semData=taGO, measure="Wang", combine="BMA")
