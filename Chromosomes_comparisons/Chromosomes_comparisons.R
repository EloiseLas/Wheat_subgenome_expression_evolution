library(stringr)
library(gt)
library(gtExtras)
library(dplyr)
library(palettes)

## Load the data:
homeologs<-read.table("Data/homoeolog_Traes_EGI.csv", sep = ",", header=TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("X")]

## Identify the location on chromosomes:
homeologs$A_chr <- unlist(lapply(homeologs[,3], function(i){strsplit(i, '')[[1]][8]}))
homeologs$B_chr <- unlist(lapply(homeologs[,2], function(i){strsplit(i, '')[[1]][8]}))
homeologs$D_chr <- unlist(lapply(homeologs[,1], function(i){strsplit(i, '')[[1]][8]}))

## Look at the different configurations:
chr_group_A <- homeologs %>% group_by(A_chr) %>% tally()
chr_group_B <- homeologs %>% group_by(B_chr) %>% tally()
chr_group_D <- homeologs %>% group_by(D_chr) %>% tally()
chr_group <- homeologs %>% group_by(A_chr,B_chr,D_chr) %>% tally()
chr_group <- chr_group[order(chr_group$n,decreasing = TRUE),]
chr_group <- data.frame(chr_group)
chr_group$A_chr <- as.numeric(chr_group$A_chr)
chr_group$B_chr <- as.numeric(chr_group$B_chr)
chr_group$D_chr <- as.numeric(chr_group$D_chr)

color_palette <- rev(met_palettes$Cassatt1)
character_vector <- as.character(color_palette)

chr_group[1:13,] %>%
  gt() %>%
  cols_label(
    A_chr = 'A',
    B_chr = 'B',
    D_chr = 'D',
    n = 'Number of Triads',
  ) %>%
  gt_theme_538() %>% 
  data_color(
    columns = c(A_chr, B_chr, D_chr),
    colors  = character_vector
  ) %>% 
  tab_spanner(
    label = md('**Chromosome**'),
    columns = c(A_chr, B_chr, D_chr)
  ) %>% 
  tab_header(
    title = 'Homeologs locations on chromosomes',
  ) 

## Focus on the possible rearrangement between chr 4 and 5 in subgenome A:
homeologs <- homeologs %>%
  mutate(special=case_when(
    A_chr == 5 & B_chr == 4 & D_chr == 4 ~ "A_5",
    A_chr == 4 & B_chr == 5 & D_chr == 5 ~ "A_4",
    TRUE ~ "N"))

expression_level<-read.table("Output/homeologs_pvalue.csv", sep = ",", header=TRUE, row.names = 1)
merged<-merge(expression_level,homeologs,by="ID")
ggplot(merged, aes(x=special, y=means_A, fill = special)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#FF5733", "#24D2F9", "#CAA0EE")) +
  xlab('Expression Profile Categories') + ylab('Mean Expression Level of the Triad') +
  guides(fill="none") +
  theme(text = element_text(size=15))

pairwise.wilcox.test(merged$means_A, merged$special)

# homeologs$category <- ifelse(homeologs$A_chr == homeologs$B_chr & homeologs$B_chr == homeologs$D_chr, "S", NA)
homeologs <- homeologs %>%
     mutate(category=case_when(
       A_chr == B_chr & B_chr == D_chr ~ "S",
       A_chr == B_chr & B_chr != D_chr ~ "D",
       A_chr != B_chr & B_chr == D_chr ~ "A",
       A_chr == D_chr & B_chr != D_chr ~ "B",
       A_chr != B_chr & B_chr != D_chr ~ "F",
       TRUE ~ "other"))

write.csv(homeologs,file = "Output/homeo_chr_category.csv")

