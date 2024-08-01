library(dplyr)
library(plyr)
library(tidyr)

correspondence_table<-read.table("Data/Traes_2_EGI.csv", sep = ",", header=TRUE, row.names = 1)
homeologs<-read.table("Data/homoeolog_CSv2_ABD.list", sep = "\t", header=FALSE,col.names=c("A","B","D"))
homeologs<-merge(homeologs,correspondence_table, by.x="A", by.y="Traes_long",all.x = TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("RefSeqID", "EGI")]
names(homeologs)[names(homeologs) == 'EGI_LOC'] <- 'EGI_LOC_A'
homeologs<-merge(homeologs,correspondence_table, by.x="B", by.y="Traes_long", all.x = TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("RefSeqID", "EGI")]
names(homeologs)[names(homeologs) == 'EGI_LOC'] <- 'EGI_LOC_B'
homeologs<-merge(homeologs,correspondence_table, by.x="D", by.y="Traes_long", all.x = TRUE)
homeologs<-homeologs[ ,  !names(homeologs) %in% 
                        c("RefSeqID", "EGI")]
names(homeologs)[names(homeologs) == 'EGI_LOC'] <- 'EGI_LOC_D'

# Add ID:
homeologs$ID<-c(1:24990)

homeologs_same_name<-homeologs

homeologs_same_name<-homeologs_same_name %>%
  mutate(logical=case_when(
    EGI_LOC_A == EGI_LOC_B | EGI_LOC_B == EGI_LOC_D | EGI_LOC_A == EGI_LOC_D ~ "Same_names",
    TRUE ~ "Different_names"))

group <- homeologs_same_name %>% group_by(logical) %>% tally()

homeologs_clean<-subset(homeologs_same_name,logical=="Different_names")
homeologs_clean<-homeologs_clean[ ,  !names(homeologs_clean) %in% 
                        c("logical")]

## Count what was lost in translation:

homeologs_count<-homeologs_clean

homeologs_count[is.na(homeologs_count)] = "NA"

homeologs_count<-homeologs_count %>%
  mutate(logical=case_when(
    EGI_LOC_A == "NA" & EGI_LOC_B == "NA" & EGI_LOC_D == "NA" ~ "All_missing",
    EGI_LOC_A == "NA" & EGI_LOC_B != "NA" & EGI_LOC_D != "NA" ~ "A_missing",
    EGI_LOC_A != "NA" & EGI_LOC_B == "NA" & EGI_LOC_D != "NA" ~ "B_missing",
    EGI_LOC_A != "NA" & EGI_LOC_B != "NA" & EGI_LOC_D == "NA" ~ "D_missing",
    EGI_LOC_A == "NA" & EGI_LOC_B == "NA" & EGI_LOC_D != "NA" ~ "AB_missing",
    EGI_LOC_A == "NA" & EGI_LOC_B != "NA" & EGI_LOC_D == "NA" ~ "AD_missing",
    EGI_LOC_A != "NA" & EGI_LOC_B == "NA" & EGI_LOC_D == "NA" ~ "BD_missing",
    TRUE ~ "Complete"))

group <- homeologs_count %>% group_by(logical) %>% tally()

sum(is.na(homeologs_clean$EGI_LOC_A))
sum(is.na(homeologs_clean$EGI_LOC_B))
sum(is.na(homeologs_clean$EGI_LOC_D))

sum(is.na(homeologs_clean))

write.csv(homeologs_clean,
          file = "Data/homoeolog_Traes_EGI.csv")
