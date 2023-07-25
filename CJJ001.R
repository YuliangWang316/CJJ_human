library(Seurat)
library(cowplot)
NB <- read.csv("e:/P23042711/NB-B/20230630/Matrix/NB-B_matrix.csv", 
               stringsAsFactors = FALSE)
DN <- read.csv("e:/P23042711/DN-B/20230630/Matrix/DN-B_matrix.csv", 
               stringsAsFactors = FALSE)
DP <- read.csv("e:/P23042711/DP-B/20230630/Matrix/DP-B_matrix.csv", 
               stringsAsFactors = FALSE)
PC <- read.csv("e:/P23042711/PC-B/20230630/Matrix/PC-B_matrix.csv", 
               stringsAsFactors = FALSE)
NB_igh <- NB[NB$chain == "IGH", ]
DN_igh <- DN[DN$chain == "IGH", ]
DP_igh <- DP[DP$chain == "IGH", ]
PC_igh <- PC[PC$chain == "IGH", ]

NB_igh_prod <- NB_igh[NB_igh$productive == "True",]
DN_igh_prod <- DN_igh[DN_igh$productive == "True",]
DP_igh_prod <- DP_igh[DP_igh$productive == "True",]
PC_igh_prod <- PC_igh[PC_igh$productive == "True",]


NB_igh_prod_vgene <- NB_igh_prod[NB_igh_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
NB_igh_prod_vgene_dedup <- NB_igh_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(NB_igh_prod_vgene_dedup$contig_id, 
            "NB_igh_valid_contig_id.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

DN_igh_prod_vgene <- DN_igh_prod[DN_igh_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
DN_igh_prod_vgene_dedup <- DN_igh_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(DN_igh_prod_vgene_dedup$contig_id, 
            "DN_igh_valid_contig_id.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
DP_igh_prod_vgene <- DP_igh_prod[DP_igh_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
DP_igh_prod_vgene_dedup <- DP_igh_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(DP_igh_prod_vgene_dedup$contig_id, 
            "DP_igh_valid_contig_id.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
PC_igh_prod_vgene <- PC_igh_prod[PC_igh_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
PC_igh_prod_vgene_dedup <- PC_igh_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(PC_igh_prod_vgene_dedup$contig_id, 
            "PC_igh_valid_contig_id.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

DN_igh_prod_vgene_dedup$group<-rep("DN",length(rownames(DN_igh_prod_vgene_dedup)))
DP_igh_prod_vgene_dedup$group<-rep("DP",length(rownames(DP_igh_prod_vgene_dedup)))
NB_igh_prod_vgene_dedup$group<-rep("NB",length(rownames(NB_igh_prod_vgene_dedup)))
PC_igh_prod_vgene_dedup$group<-rep("PC",length(rownames(PC_igh_prod_vgene_dedup)))


Total_igh_prod_vgene_dedup<-rbind(DN_igh_prod_vgene_dedup,DP_igh_prod_vgene_dedup,NB_igh_prod_vgene_dedup,PC_igh_prod_vgene_dedup)

NB_prod <- NB[NB$productive == "True",]
DN_prod <- DN[DN$productive == "True",]
DP_prod <- DP[DP$productive == "True",]
PC_prod <- PC[PC$productive == "True",]


NB_prod_vgene <- NB_prod[NB_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
NB_prod_vgene_dedup <- NB_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(NB_prod_vgene_dedup$contig_id, 
            "NB_valid_contig_id_with_igk.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

DN_prod_vgene <- DN_prod[DN_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
DN_prod_vgene_dedup <- DN_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(DN_prod_vgene_dedup$contig_id, 
            "DN_valid_contig_id_with_igk.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
DP_prod_vgene <- DP_prod[DP_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
DP_prod_vgene_dedup <- DP_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(DP_prod_vgene_dedup$contig_id, 
            "DP_valid_contig_id_with_igk.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
PC_prod_vgene <- PC_prod[PC_prod$v_gene != "None", ]
## Remove cell barcode with more than 1 heavy chain assembled
library(dplyr)
PC_prod_vgene_dedup <- PC_prod_vgene %>%
  group_by(barcode) %>% 
  filter(n() == 1)
## Export valid contig ID
write.table(PC_prod_vgene_dedup$contig_id, 
            "PC_valid_contig_id_with_igk.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

DN_prod_vgene_dedup$group<-rep("DN",length(rownames(DN_prod_vgene_dedup)))
DP_prod_vgene_dedup$group<-rep("DP",length(rownames(DP_prod_vgene_dedup)))
NB_prod_vgene_dedup$group<-rep("NB",length(rownames(NB_prod_vgene_dedup)))
PC_prod_vgene_dedup$group<-rep("PC",length(rownames(PC_prod_vgene_dedup)))


Total_prod_vgene_dedup<-rbind(DN_prod_vgene_dedup,DP_prod_vgene_dedup,NB_prod_vgene_dedup,PC_prod_vgene_dedup)

library(dplyr)
library(tidyverse)
DN_igblast <- read.delim("e:/P23042711/DN-B/DN_blast_01_db-pass.tsv", 
                         stringsAsFactors = FALSE)
DP_igblast <- read.delim("e:/P23042711/DP-B/DP_blast_01_db-pass.tsv", 
                         stringsAsFactors = FALSE)
NB_igblast <- read.delim("e:/P23042711/NB-B/NB_blast_01_db-pass.tsv", 
                         stringsAsFactors = FALSE)
PC_igblast <- read.delim("e:/P23042711/PC-B/PC_blast_01_db-pass.tsv", 
                         stringsAsFactors = FALSE)
DN_igblast$group<-rep("DN",length(rownames(DN_igblast)))
DP_igblast$group<-rep("DP",length(rownames(DP_igblast)))
NB_igblast$group<-rep("NB",length(rownames(NB_igblast)))
PC_igblast$group<-rep("PC",length(rownames(PC_igblast)))
Total_igblast<-rbind(DN_igblast,DP_igblast,NB_igblast,PC_igblast)
write.table(Total_igblast,file = "e:/P23042711/Total.txt",sep = "\t")
# library(tigger)
# novel <- findNovelAlleles(Total_igblast, SampleGermlineIGHV, nproc=40)
# novel_rows <- selectNovel(novel)
# novel_row <- which(!is.na(novel$polymorphism_call))[1]
# plotNovel(AIRRDb, novel[novel_row, ])
#whether choose c gene?
library(alakazam)
library(dowser)
library(dplyr)
library(scales)
library(shazam)
Total_dist_ham <- distToNearest(Total_igblast , 
                                sequenceColumn="junction", 
                                vCallColumn="v_call", jCallColumn="j_call",
                                model="ham", normalize="len", nproc=40)
Total_dist_s5f <- distToNearest(Total_igblast, 
                                sequenceColumn="junction", 
                                vCallColumn="v_call", jCallColumn="j_call",
                                model="hh_s5f", normalize="none", nproc=40)
library(ggplot2)
p1 <- ggplot(subset(Total_dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)
p2 <- ggplot(subset(Total_dist_s5f, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("HH_S5F distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 50, 5)) +
  geom_histogram(color="white", binwidth=1) +
  geom_vline(xintercept=7, color="firebrick", linetype=2)
plot(p2)
output <- findThreshold(Total_dist_ham$dist_nearest, method="density")
threshold <- output@threshold
plot(output, title="Density Method")
print(output)
output <- findThreshold(Total_dist_ham$dist_nearest, method="gmm", model="gamma-gamma")
plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
print(output)

library(tidyverse)

Total_new_igblast <- add_column(Total_igblast, 
                                c_gene = Total_igh_prod_vgene_dedup$c_gene[match(Total_igblast$sequence_id,Total_igh_prod_vgene_dedup$contig_id)], 
                                .after = "j_support")
Total_new_igblast <- add_column(Total_new_igblast, 
                     length = sapply(strsplit(Total_new_igblast$sequence, ""), length), 
                     .after = "locus")

Total_new_igblast_new<-Total_new_igblast[Total_new_igblast$c_gene != "",]
Total_new_igblast_new$group<-factor(Total_new_igblast_new$group,levels = c("NB","DN","DP","PC"))
Total_new_igblast_new<-Total_new_igblast_new[which(Total_new_igblast_new$c_gene == "IGHA1" | Total_new_igblast_new$c_gene == "IGHA2" | Total_new_igblast_new$c_gene == "IGHG1" | Total_new_igblast_new$c_gene == "IGHG2" | Total_new_igblast_new$c_gene == "IGHG3" | Total_new_igblast_new$c_gene == "IGHG4" | Total_new_igblast_new$c_gene == "IGHM" ),]
Total_new_igblast_new$c_gene<-factor(Total_new_igblast_new$c_gene,levels = c("IGHM","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4"))
write.table(Total_new_igblast_new,file = "e:/P23042711/Total_new.txt",sep = "\t")


Total_dist_ham <- distToNearest(Total_new_igblast_new , 
                                sequenceColumn="junction", 
                                vCallColumn="v_call", jCallColumn="j_call",
                                model="ham", normalize="len", nproc=40)
Total_dist_s5f <- distToNearest(Total_new_igblast_new, 
                                sequenceColumn="junction", 
                                vCallColumn="v_call", jCallColumn="j_call",
                                model="hh_s5f", normalize="none", nproc=40)
library(ggplot2)
p1 <- ggplot(subset(Total_dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("Hamming distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 1, 0.1)) +
  geom_histogram(color="white", binwidth=0.02) +
  geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)
p2 <- ggplot(subset(Total_dist_s5f, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
  theme_bw() + 
  xlab("HH_S5F distance") + 
  ylab("Count") +
  scale_x_continuous(breaks=seq(0, 50, 5)) +
  geom_histogram(color="white", binwidth=1) +
  geom_vline(xintercept=7, color="firebrick", linetype=2)
plot(p2)
output <- findThreshold(Total_dist_ham$dist_nearest, method="density")
threshold <- output@threshold
plot(output, title="Density Method")
print(output)
output <- findThreshold(Total_dist_ham$dist_nearest, method="gmm", model="gamma-gamma")
plot(output, binwidth=0.02, title="GMM Method: gamma-gamma")
print(output)

Total_new_igblast_new<-read.delim("e:/P23042711/Total_new_clone-pass.tsv", 
                                        stringsAsFactors = FALSE)
clone_number<-as.data.frame(table(Total_new_igblast_new$clone_id))
clone_use<-clone_number[which(clone_number$Freq >2 ),]
clone_use<-clone_use$Var1
Total_new_igblast_new$Hits<-"No_Hit"
for (i in 1:length(clone_use)) {
  for (j in 1:length(rownames(Total_new_igblast_new))) {
    
    if(clone_use[i] == Total_new_igblast_new$clone_id[j]){
      Total_new_igblast_new$Hits[j]<-"Hit"      
    }
  }
  
}
Total_new_igblast_new_clone_use<-Total_new_igblast_new[which(Total_new_igblast_new$Hits == "Hit"),]
write.table(Total_new_igblast_new_clone_use,file = "e:/P23042711/Total_new_clone_use.txt",sep = "\t")

Total_igblast_new<-read.delim("e:/P23042711/Total_clone-pass.tsv", 
                                  stringsAsFactors = FALSE)
clone_number<-as.data.frame(table(Total_igblast_new$clone_id))
clone_use<-clone_number[which(clone_number$Freq >2 ),]
clone_use<-clone_use$Var1
Total_igblast_new$Hits<-"No_Hit"
for (i in 1:length(clone_use)) {
  for (j in 1:length(rownames(Total_igblast_new))) {
    
    if(clone_use[i] == Total_igblast_new$clone_id[j]){
      Total_igblast_new$Hits[j]<-"Hit"      
    }
  }
  
}
Total_igblast_new_clone_use<-Total_igblast_new[which(Total_igblast_new$Hits == "Hit"),]
write.table(Total_igblast_new_clone_use,file = "e:/P23042711/Total_clone_use.txt",sep = "\t")
