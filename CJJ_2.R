library(ggplot2)
library(ggpubr)

NB<-read.csv("e:/P23042711/NB-B/20230630/Matrix/NB-B_matrix.csv",sep = ",",header = TRUE)
DN<-read.csv("e:/P23042711/DN-B/20230630/Matrix/DN-B_matrix.csv",sep = ",",header = TRUE)
DP<-read.csv("e:/P23042711/DP-B/20230630/Matrix/DP-B_matrix.csv",sep = ",",header = TRUE)
PC<-read.csv("e:/P23042711/PC-B/20230630/Matrix/PC-B_matrix.csv",sep = ",",header = TRUE)

NB<-NB[which(NB$chain == "IGH"),]
NB<-NB[which(NB$c_gene != ""),]
DN<-DN[which(DN$chain == "IGH"),]
DN<-DN[which(DN$c_gene != ""),]
DP<-DP[which(DP$chain == "IGH"),]
DP<-DP[which(DP$c_gene != ""),]
PC<-PC[which(PC$chain == "IGH"),]
PC<-PC[which(PC$c_gene != ""),]

NB$group<-rep("NB",length(rownames(NB)))
DN$group<-rep("DN",length(rownames(DN)))
DP$group<-rep("DP",length(rownames(DP)))
PC$group<-rep("PC",length(rownames(PC)))

Total<-rbind(NB,DN,DP,PC)
Cellratio <- prop.table(table(Total$c_gene, Total$group), margin = 2)

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
