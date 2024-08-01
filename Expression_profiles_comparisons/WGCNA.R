library(WGCNA)
library(tidyverse)     
library(magrittr)      
library(DESeq2)
library(genefilter)
library(dplyr)
library(pheatmap)
library(dplyr)
library(AnnotationHub)
library(stringr)
library(clusterProfiler)
library(GO.db)
library(ggplot2)
library(enrichplot)
library(reshape2)

## Part 1 Version 1, on raw data: 

expression_data_merged<-read.table("Data/expression_data_merged.csv", sep = ",", header=TRUE, row.names = 1)

# Normalize Counts with DESeq:

# The next DESeq2 functions need the values to be converted to integers
expression_data_merged <- round(expression_data_merged) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50)

de_input = as.matrix(expression_data_merged)
de_input[1:5,1:10]

tissue<-read.table("Data/tissue_formated.csv", sep = ",", header=TRUE, col.names=c("condition","tissue"))
meta_df <- data.frame(condition = names(expression_data_merged))
meta_df <- merge(meta_df, tissue, by="condition", all.x=TRUE, all.y=FALSE)  # merge by row names (by=0 or by="row.names")

dds <- DESeqDataSetFromMatrix(round(de_input),meta_df,design = ~ 1)

dds <- DESeq(dds)

# Option 1 with filtering:
vsd <- varianceStabilizingTransformation(dds)

wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)
q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

# Option 2, without filtering:
dds_norm <- vst(dds)
expr_normalized <- assay(dds_norm)
expr_normalized_d <- data.frame(expr_normalized)
row.names(expr_normalized_d)<-row.names(expression_data_merged)
write.csv(expr_normalized,
            file = "Data/expression_data_normalized.csv")

input_mat = t(expr_normalized)

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

## Part 1 Version 2, on UQ + Batch correction + TPM normalized data:

expression_data_UQ_batch_normalized<-read.table("Data/expression_data_TPM_UQ_batch_filtered.csv", sep = ",", header=TRUE, row.names = 1)
expression_data_gene_norm <- sweep(expression_data_UQ_batch_normalized, 1, rowMeans(expression_data_UQ_batch_normalized))
# expression_data_gene_norm <- apply(expression_data_UQ_batch_normalized, 1, function(i) i*100/sum(i))
input_mat = t(as.matrix(expression_data_gene_norm))

## Part 2, WGCNA:

allowWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

# Plots:
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point(size=3) +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.05) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  theme(text = element_text(size=18))
  # This adds some nicer aesthetics to our plot
  # theme_classic()

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# WGCNA:
picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
sub<-input_mat[,1:10000]
netwk <- blockwiseModules(input_mat,                # <= input here
                          randomSeed = 1234,
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree Options ==
                          deepSplit = 2,
                          pamStage = T,
                          pamRespectsDendro = T,
                          #detectCutHeight = 0.9,
                          #minModuleSize = 30,
                          
                          # == Block Options ==
                          maxBlockSize = 5000,
                          corType = "pearson",
                          
                          # == Module Adjustments ==
                          reassignThreshold = 1e-6, ###
                          minKMEtoStay = 0.01,
                          
                          # == Module Merging ==
                          mergeCutHeight = 0.3, ###
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

readr::write_rds(netwk,file = "Output/WGCNA_netwk_TPM_v4.RDS")


cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[3]],
  mergedColors[netwk$blockGenes[[3]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  cex.colorLabels = 1.2)

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

write_delim(module_df,
            file = "Output/gene_modules_TPM_UQ_batch_v4.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

write_delim(MEs0,
              file = "Output/eigengene_modules_TPM_UQ_batch_v4.txt",
            delim = "\t")

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

# mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
#   geom_tile() +
#   theme_bw() +
#   scale_fill_gradient2(
#     low = "blue",
#     high = "red",
#     mid = "white",
#     midpoint = 0,
#     limit = c(-1,1)) +
#   theme(axis.text.x = element_text(angle=90)) +
#   labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

modules_eigengene<-t(MEs0)

# Plots: 
pheatmap(MEs0, cluster_cols=TRUE, cluster_rows=FALSE, show_colnames=TRUE, kmeans_k=18,fontsize=13)

module_df<-read.table("Output/gene_modules_TPM_UQ_batch.txt", sep = "\t", header=TRUE)
correspondence_table<-read.table("Data/Traes_2_EGI.csv", sep = ",", header=TRUE, row.names = 1)
module_df<-merge(module_df,correspondence_table, by.x="gene_id", by.y="EGI_LOC",all.x = TRUE)
module_df<-module_df %>% drop_na()
module_df$Traes<-unlist(lapply(module_df$Traes_long,function(i){str_sub(i,end=-3)}))
df<-module_df %>%
  group_by(colors) %>%
  summarize(Traes = list(Traes))
l<-list(df$Traes)[[1]]
names(l)<-df$colors

## Enrichment Analysis:
ck <- compareCluster(geneCluster = l, fun = enrichGO, OrgDb = "org.Taestivum.eg.db", keyType = "GID", ont="BP")
ck <- setReadable(ck, OrgDb = "org.Taestivum.eg.db", keyType="GID")
res<-ck@compareClusterResult
res<-res %>% # take the dataframe
  group_by(Cluster) %>% # group it by the grouping variable
  slice(1:3) # and pick rows 1 to 3 per group
mat<-as.data.frame.matrix(xtabs(-log10(p.adjust)~ID+Cluster,res,addNA=TRUE))
mat[mat>50]<-50
mat_filter<-mat[rowSums(mat > 10) >= 1, ]
pheatmap(mat, cluster_cols=TRUE, cluster_rows=TRUE, show_colnames=TRUE, show_rownames=TRUE,fontsize=13,fontsize_row=9)

dotplot(ck)

cnetplot(ck,node_label="category",cex_category=10, showCategory=2)
