library(Seurat)
library(cowplot)
pbmc<-readRDS("d:/CJJ_humanTonsil_DP_DN_NB_PC_scRNA.rds")
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
DN_igblast <- read.delim("e:/P23042711/DN-B/DN_igh_valid_contig_igblast_annotation.tsv", 
                      stringsAsFactors = FALSE)
DP_igblast <- read.delim("e:/P23042711/DP-B/DP_igh_valid_contig_igblast_annotation.tsv", 
                         stringsAsFactors = FALSE)
NB_igblast <- read.delim("e:/P23042711/NB-B/NB_igh_valid_contig_igblast_annotation.tsv", 
                         stringsAsFactors = FALSE)
PC_igblast <- read.delim("e:/P23042711/PC-B/PC_igh_valid_contig_igblast_annotation.tsv", 
                         stringsAsFactors = FALSE)
DN_result <- DN_igblast %>% 
  dplyr::select(barcode = sequence_id,
                chain = locus,
                stop_codon,
                productive,
                v_gene = v_call,
                v_score,
                v_evalue = v_support,
                d_gene = d_call,
                d_score,
                d_evalue = d_support,
                j_gene = j_call,
                j_score,
                j_evalue = j_support,
                v_germline_start)
DP_result <- DP_igblast %>% 
  dplyr::select(barcode = sequence_id,
                chain = locus,
                stop_codon,
                productive,
                v_gene = v_call,
                v_score,
                v_evalue = v_support,
                d_gene = d_call,
                d_score,
                d_evalue = d_support,
                j_gene = j_call,
                j_score,
                j_evalue = j_support,
                v_germline_start)
NB_result <- NB_igblast %>% 
  dplyr::select(barcode = sequence_id,
                chain = locus,
                stop_codon,
                productive,
                v_gene = v_call,
                v_score,
                v_evalue = v_support,
                d_gene = d_call,
                d_score,
                d_evalue = d_support,
                j_gene = j_call,
                j_score,
                j_evalue = j_support,
                v_germline_start)
PC_result <- PC_igblast %>% 
  dplyr::select(barcode = sequence_id,
                chain = locus,
                stop_codon,
                productive,
                v_gene = v_call,
                v_score,
                v_evalue = v_support,
                d_gene = d_call,
                d_score,
                d_evalue = d_support,
                j_gene = j_call,
                j_score,
                j_evalue = j_support,
                v_germline_start)
library(tidyverse)
DN_result <- add_column(DN_result, 
                     c_gene = DN_igh_prod_vgene_dedup$c_gene[match(DN_igblast$sequence_id, 
                                                                DN_igh_prod_vgene_dedup$contig_id)], 
                     .after = "j_evalue")
DP_result <- add_column(DP_result, 
                        c_gene = DP_igh_prod_vgene_dedup$c_gene[match(DP_igblast$sequence_id, 
                                                                      DP_igh_prod_vgene_dedup$contig_id)], 
                        .after = "j_evalue")
NB_result <- add_column(NB_result, 
                        c_gene = NB_igh_prod_vgene_dedup$c_gene[match(NB_igblast$sequence_id, 
                                                                      NB_igh_prod_vgene_dedup$contig_id)], 
                        .after = "j_evalue")
PC_result <- add_column(PC_result, 
                        c_gene = PC_igh_prod_vgene_dedup$c_gene[match(PC_igblast$sequence_id, 
                                                                      PC_igh_prod_vgene_dedup$contig_id)], 
                        .after = "j_evalue")
DN_result <- add_column(DN_result, 
                     length = sapply(strsplit(DN_igblast$sequence, ""), length), 
                     .after = "chain")
DP_result <- add_column(DP_result, 
                        length = sapply(strsplit(DP_igblast$sequence, ""), length), 
                        .after = "chain")
NB_result <- add_column(NB_result, 
                        length = sapply(strsplit(NB_igblast$sequence, ""), length), 
                        .after = "chain")
PC_result <- add_column(PC_result, 
                        length = sapply(strsplit(PC_igblast$sequence, ""), length), 
                        .after = "chain")
DN_result$Total_length <- nchar(DN_igblast$v_sequence_alignment)
DP_result$Total_length <- nchar(DP_igblast$v_sequence_alignment)
NB_result$Total_length <- nchar(NB_igblast$v_sequence_alignment)
PC_result$Total_length <- nchar(PC_igblast$v_sequence_alignment)
DN_vgene_nt_match <- sapply(1:nrow(DN_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(DN_igblast[x,]$v_sequence_alignment, DN_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
DN_result$Total_match <- sapply(DN_vgene_nt_match, sum)
DP_vgene_nt_match <- sapply(1:nrow(DP_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(DP_igblast[x,]$v_sequence_alignment, DP_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
DP_result$Total_match <- sapply(DP_vgene_nt_match, sum)
NB_vgene_nt_match <- sapply(1:nrow(NB_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(NB_igblast[x,]$v_sequence_alignment, NB_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
NB_result$Total_match <- sapply(NB_vgene_nt_match, sum)
PC_vgene_nt_match <- sapply(1:nrow(PC_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(PC_igblast[x,]$v_sequence_alignment, PC_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
PC_result$Total_match <- sapply(PC_vgene_nt_match, sum)
DN_vgene_nt_mismatch <- sapply(1:nrow(DN_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(DN_igblast[x,]$v_sequence_alignment, DN_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] != i[2] & i[1] != "-" & i[2] != "-"
        }
  )
}
)
DN_result$Total_mismatch <- sapply(DN_vgene_nt_mismatch, sum)
DP_vgene_nt_mismatch <- sapply(1:nrow(DP_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(DP_igblast[x,]$v_sequence_alignment, DP_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] != i[2] & i[1] != "-" & i[2] != "-"
        }
  )
}
)
DP_result$Total_mismatch <- sapply(DP_vgene_nt_mismatch, sum)
NB_vgene_nt_mismatch <- sapply(1:nrow(NB_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(NB_igblast[x,]$v_sequence_alignment, NB_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] != i[2] & i[1] != "-" & i[2] != "-"
        }
  )
}
)
NB_result$Total_mismatch <- sapply(NB_vgene_nt_mismatch, sum)
PC_vgene_nt_mismatch <- sapply(1:nrow(PC_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(PC_igblast[x,]$v_sequence_alignment, PC_igblast[x,]$v_germline_alignment), "")), 
        2, 
        function(i){
          i[1] != i[2] & i[1] != "-" & i[2] != "-"
        }
  )
}
)
PC_result$Total_mismatch <- sapply(PC_vgene_nt_mismatch, sum)
DN_gap_query <- str_count(DN_igblast$v_sequence_alignment, "-")
DN_gap_ref <- str_count(DN_igblast$v_germline_alignment, "-")
DN_result <- DN_result %>%
  mutate(Total_diff = Total_length-Total_match,
         Total_gap = DN_gap_query + DN_gap_ref,
         Total_identity = Total_match/Total_length)
DP_gap_query <- str_count(DP_igblast$v_sequence_alignment, "-")
DP_gap_ref <- str_count(DP_igblast$v_germline_alignment, "-")
DP_result <- DP_result %>%
  mutate(Total_diff = Total_length-Total_match,
         Total_gap = DP_gap_query + DP_gap_ref,
         Total_identity = Total_match/Total_length)
NB_gap_query <- str_count(NB_igblast$v_sequence_alignment, "-")
NB_gap_ref <- str_count(NB_igblast$v_germline_alignment, "-")
NB_result <- NB_result %>%
  mutate(Total_diff = Total_length-Total_match,
         Total_gap = NB_gap_query + NB_gap_ref,
         Total_identity = Total_match/Total_length)
PC_gap_query <- str_count(PC_igblast$v_sequence_alignment, "-")
PC_gap_ref <- str_count(PC_igblast$v_germline_alignment, "-")
PC_result <- PC_result %>%
  mutate(Total_diff = Total_length-Total_match,
         Total_gap = PC_gap_query + PC_gap_ref,
         Total_identity = Total_match/Total_length)
DN_result$indel <- ifelse(nchar(DN_igblast$v_sequence_alignment_aa)==nchar(DN_igblast$v_germline_alignment_aa), 
                       0, 
                       1)
NB_result$indel <- ifelse(nchar(NB_igblast$v_sequence_alignment_aa)==nchar(NB_igblast$v_germline_alignment_aa), 
                          0, 
                          1)
DP_result$indel <- ifelse(nchar(DP_igblast$v_sequence_alignment_aa)==nchar(DP_igblast$v_germline_alignment_aa), 
                          0, 
                          1)
PC_result$indel <- ifelse(nchar(PC_igblast$v_sequence_alignment_aa)==nchar(PC_igblast$v_germline_alignment_aa), 
                          0, 
                          1)
DN_vgene_mutation_offset <- (DN_igblast$v_germline_start + 1) %/% 3
DP_vgene_mutation_offset <- (DP_igblast$v_germline_start + 1) %/% 3
NB_vgene_mutation_offset <- (NB_igblast$v_germline_start + 1) %/% 3
PC_vgene_mutation_offset <- (PC_igblast$v_germline_start + 1) %/% 3
DN_vgene_aa_match <- sapply(1:nrow(DN_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(DN_igblast[x,]$v_sequence_alignment_aa, 
                                  DN_igblast[x,]$v_germline_alignment_aa), 
                                "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
DN_vgene_mutation_pos <- lapply(DN_vgene_aa_match, function(x) which(!x))
DP_vgene_aa_match <- sapply(1:nrow(DP_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(DP_igblast[x,]$v_sequence_alignment_aa, 
                                  DP_igblast[x,]$v_germline_alignment_aa), 
                                "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
DP_vgene_mutation_pos <- lapply(DP_vgene_aa_match, function(x) which(!x))
NB_vgene_aa_match <- sapply(1:nrow(NB_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(NB_igblast[x,]$v_sequence_alignment_aa, 
                                  NB_igblast[x,]$v_germline_alignment_aa), 
                                "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
NB_vgene_mutation_pos <- lapply(NB_vgene_aa_match, function(x) which(!x))
PC_vgene_aa_match <- sapply(1:nrow(PC_igblast), function(x) {
  apply(do.call(rbind, strsplit(c(PC_igblast[x,]$v_sequence_alignment_aa, 
                                  PC_igblast[x,]$v_germline_alignment_aa), 
                                "")), 
        2, 
        function(i){
          i[1] == i[2]
        }
  )
}
)
PC_vgene_mutation_pos <- lapply(PC_vgene_aa_match, function(x) which(!x))
DN_vgene_mutation_pos <- lapply(1:nrow(DN_igblast), 
                             function(i) DN_vgene_mutation_pos[[i]] + 
                               DN_vgene_mutation_offset[i])
DP_vgene_mutation_pos <- lapply(1:nrow(DP_igblast), 
                                function(i) DP_vgene_mutation_pos[[i]] + 
                                  DP_vgene_mutation_offset[i])
NB_vgene_mutation_pos <- lapply(1:nrow(NB_igblast), 
                                function(i) NB_vgene_mutation_pos[[i]] + 
                                  NB_vgene_mutation_offset[i])
PC_vgene_mutation_pos <- lapply(1:nrow(PC_igblast), 
                                function(i) PC_vgene_mutation_pos[[i]] + 
                                  PC_vgene_mutation_offset[i])
library(Seurat)
library(tidyverse)
library(RColorBrewer)
DN_igh_prod_vgene_dedup_IgG1<-DN_igh_prod_vgene_dedup[DN_igh_prod_vgene_dedup$c_gene == "IGHG1",]
DN_igh_prod_vgene_dedup_IgG1 %>%
  group_by(v_gene) %>%
  summarise(count = n()) %>%
  mutate(fraction = count/sum(count),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1))) %>%
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=v_gene)) +
  geom_rect(color="black") +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(61)) +
  theme_void()

DP_igh_prod_vgene_dedup_IgG1<-DP_igh_prod_vgene_dedup[DP_igh_prod_vgene_dedup$c_gene == "IGHG1",]
DP_igh_prod_vgene_dedup_IgG1 %>%
  group_by(v_gene) %>%
  summarise(count = n()) %>%
  mutate(fraction = count/sum(count),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1))) %>%
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=v_gene)) +
  geom_rect(color="black") +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(61)) +
  theme_void()
NB_igh_prod_vgene_dedup_IgG1<-NB_igh_prod_vgene_dedup[NB_igh_prod_vgene_dedup$c_gene == "IGHG1",]
NB_igh_prod_vgene_dedup_IgG1 %>%
  group_by(v_gene) %>%
  summarise(count = n()) %>%
  mutate(fraction = count/sum(count),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1))) %>%
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=v_gene)) +
  geom_rect(color="black") +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(61)) +
  theme_void()
PC_igh_prod_vgene_dedup_IgG1<-PC_igh_prod_vgene_dedup[PC_igh_prod_vgene_dedup$c_gene == "IGHG1",]
PC_igh_prod_vgene_dedup_IgG1 %>%
  group_by(v_gene) %>%
  summarise(count = n()) %>%
  mutate(fraction = count/sum(count),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n=-1))) %>%
  ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=v_gene)) +
  geom_rect(color="black") +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(61)) +
  theme_void()
DP_IGHG1<-DP_result[DP_result$c_gene == "IGHG1",]
NB_IGM<-NB_result[NB_result$c_gene == "IGHM",]
DN_IGHG1<-DN_result[DN_result$c_gene == "IGHG1",]
DP_IGHG1$Total<-rep(1,length(rownames(DP_IGHG1)))
NB_IGM$Total<-rep(1,length(rownames(NB_IGM)))
DN_IGHG1$Total<-rep(1,length(rownames(DN_IGHG1)))
DP_IGHG1$Total_unidentify<-DP_IGHG1$Total - DP_IGHG1$Total_identity
NB_IGM$Total_unidentify<-NB_IGM$Total - NB_IGM$Total_identity
DN_IGHG1$Total_unidentify<-DN_IGHG1$Total - DN_IGHG1$Total_identity
library(beeswarm)
library(ggbeeswarm)
ggplot(DP_IGHG1,aes(x= c_gene,y=Total_unidentify))+ geom_quasirandom(size=0.5,alpha=0.8) + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +geom_boxplot(alpha=0) + scale_color_manual(values = c("Gray"))
ggplot(NB_IGM,aes(x= c_gene,y=Total_unidentify))+ geom_quasirandom(size=0.5,alpha=0.8) + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +geom_boxplot(alpha=0) + scale_color_manual(values = c("Gray"))
ggplot(DN_IGHG1,aes(x= c_gene,y=Total_unidentify))+ geom_quasirandom(size=0.5,alpha=0.8) + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +geom_boxplot(alpha=0) + scale_color_manual(values = c("Gray"))
write.table(DP_IGHG1,file = "DP_IGHG1mutation.txt",sep = "\t")
write.table(DN_IGHG1,file = "DN_IGHG1mutation.txt",sep = "\t")
write.table(NB_IGM,file = "NB_IGMmutation.txt",sep = "\t")
DN_result$group<-rep("DN",length(rownames(DN_result)))
DP_result$group<-rep("DP",length(rownames(DP_result)))
NB_result$group<-rep("NB",length(rownames(NB_result)))
PC_result$group<-rep("PC",length(rownames(PC_result)))
DN_result$Total<-rep(1,length(rownames(DN_result)))
DP_result$Total<-rep(1,length(rownames(DP_result)))
NB_result$Total<-rep(1,length(rownames(NB_result)))
PC_result$Total<-rep(1,length(rownames(PC_result)))
DN_result$Total_unidentify<-DN_result$Total - DN_result$Total_identity
DP_result$Total_unidentify<-DP_result$Total - DP_result$Total_identity
NB_result$Total_unidentify<-NB_result$Total - NB_result$Total_identity
PC_result$Total_unidentify<-PC_result$Total - PC_result$Total_identity
Total<-rbind(DN_result,DP_result,NB_result,PC_result)
Total<-Total[Total$c_gene != "",]
Total$group<-factor(Total$group,levels = c("NB","DN","DP","PC"))

library(ggplot2)
library(ggprism)
ggplot(data = Total, aes(x = c_gene, y = Total_unidentify, fill = group)) +
  geom_violin()+
  geom_boxplot()+
  theme_prism()+
  scale_color_manual(values = brewer.pal(9,"Set1")[3:8])
DP_IGHG1<-DP_result[DP_result$c_gene == "IGHG1",]
DN_IGHG1<-DN_result[DN_result$c_gene == "IGHG1",]
PC_IGHG1<-PC_result[PC_result$c_gene == "IGHG1",]
Total_IGHG1<-Total[Total$c_gene == "IGHG1",]
IGHG1_Cellratio <- as.data.frame(prop.table(table(Total_IGHG1$v_gene, Total_IGHG1$group), margin = 2))
colourCount = length(unique(IGHG1_Cellratio$Var1))
library(ggplot2)
ggplot(IGHG1_Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
rownames(DN_igblast)<-DN_igblast[,1]
rownames(DP_igblast)<-DP_igblast[,1]
rownames(PC_igblast)<-PC_igblast[,1]
rownames(DN_result)<-DN_result[,1]
rownames(DP_result)<-DP_result[,1]
rownames(PC_result)<-PC_result[,1]
DN_result<-DN_result[rownames(DN_igblast),]
DP_result<-DP_result[rownames(DP_igblast),]
PC_result<-PC_result[rownames(PC_igblast),]
DN_result_igblast<-cbind(DN_result,DN_igblast)
DP_result_igblast<-cbind(DP_result,DP_igblast)
PC_result_igblast<-cbind(PC_result,PC_igblast)
Total_result_igblast<-rbind(DN_result_igblast,DP_result_igblast,PC_result_igblast)
Total_result_igblast_IGHG1<-Total_result_igblast[which(Total_result_igblast$c_gene == "IGHG1"),]
DN_igblast_IGHG1_vgene<-Total_result_igblast_IGHG1[which(Total_result_igblast_IGHG1$group == "DN" ),]
DP_igblast_IGHG1_vgene<-Total_result_igblast_IGHG1[which(Total_result_igblast_IGHG1$group == "DP" ),]
PC_igblast_IGHG1_vgene<-Total_result_igblast_IGHG1[which(Total_result_igblast_IGHG1$group == "PC" ),]

PC_DP_IGHG1_overlapVgene<-intersect(PC_igblast_IGHG1_vgene$v_sequence_alignment,DP_igblast_IGHG1_vgene$v_sequence_alignment)
PC_DN_IGHG1_overlapVgene<-intersect(PC_igblast_IGHG1_vgene$v_sequence_alignment,DN_igblast_IGHG1_vgene$v_sequence_alignment)

PC_DP_result_igblast_IGHG1<-Total_result_igblast_IGHG1[which(Total_result_igblast_IGHG1$group == "DP" | Total_result_igblast_IGHG1$group == "PC"),]
PC_DN_result_igblast_IGHG1<-Total_result_igblast_IGHG1[which(Total_result_igblast_IGHG1$group == "DN" | Total_result_igblast_IGHG1$group == "PC"),]

PC_DP_result_igblast_IGHG1$Hits<-"No_Hit"
PC_DN_result_igblast_IGHG1$Hits<-"No_Hit"

for (i in 1:length(PC_DP_IGHG1_overlapVgene)) {
  for (j in 1:length(rownames(PC_DP_result_igblast_IGHG1))) {
    
      if(PC_DP_IGHG1_overlapVgene[i] ==  PC_DP_result_igblast_IGHG1$v_sequence_alignment[j] ){
        PC_DP_result_igblast_IGHG1$Hits[j]<-"Hit"      
      }
    }
    
  }


for (i in length(PC_DN_IGHG1_overlapVgene)) {
  for (j in 1:length(rownames(PC_DN_result_igblast_IGHG1))) {
    
      if(PC_DN_IGHG1_overlapVgene[i] == PC_DN_result_igblast_IGHG1$v_sequence_alignment[j]){
        PC_DN_result_igblast_IGHG1$Hits[j]<-"Hit"      
      }
    }
    
  }
PC_DP_same_IGHG1_Cellratio <- as.data.frame(prop.table(table(PC_DP_result_igblast_IGHG1$Hits, PC_DP_result_igblast_IGHG1$group), margin = 2))
colourCount = length(unique(PC_DP_same_IGHG1_Cellratio$Var1))
library(ggplot2)
ggplot(PC_DP_same_IGHG1_Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
PC_DN_same_IGHG1_Cellratio <- as.data.frame(prop.table(table(PC_DN_result_igblast_IGHG1$Hits, PC_DN_result_igblast_IGHG1$group), margin = 2))
colourCount = length(unique(PC_DN_same_IGHG1_Cellratio$Var1))
library(ggplot2)
ggplot(PC_DN_same_IGHG1_Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

