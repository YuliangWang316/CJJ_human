library(dplyr)
library(tidyverse)
DN_igblast <- read.delim("e:/P23042711/DN-B/DN_igh_valid_contig_igblast_annotation_new.tsv", 
                         stringsAsFactors = FALSE)
DP_igblast <- read.delim("e:/P23042711/DP-B/DP_igh_valid_contig_igblast_annotation_new.tsv", 
                         stringsAsFactors = FALSE)
NB_igblast <- read.delim("e:/P23042711/NB-B/NB_igh_valid_contig_igblast_annotation_new.tsv", 
                         stringsAsFactors = FALSE)
PC_igblast <- read.delim("e:/P23042711/PC-B/PC_igh_valid_contig_igblast_annotation_new.tsv", 
                         stringsAsFactors = FALSE)
DN_igblast$group<-rep("DN",length(rownames(DN_igblast)))
DP_igblast$group<-rep("DP",length(rownames(DP_igblast)))
NB_igblast$group<-rep("NB",length(rownames(NB_igblast)))
PC_igblast$group<-rep("PC",length(rownames(PC_igblast)))
Total_igblast<-rbind(DN_igblast,DP_igblast,NB_igblast,PC_igblast)
write.table(Total_igblast,file = "e:/P23042711/Total_new.txt",sep = "\t")
# library(tigger)
# novel <- findNovelAlleles(Total_igblast, SampleGermlineIGHV, nproc=40)
# novel_rows <- selectNovel(novel)
# novel_row <- which(!is.na(novel$polymorphism_call))[1]
# plotNovel(AIRRDb, novel[novel_row, ])
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
