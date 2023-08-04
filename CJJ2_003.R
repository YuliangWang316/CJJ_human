Total <- read.delim("e:/P23042711/Total_clone-pass_germ-pass.tsv", 
                    stringsAsFactors = FALSE)
Total$Total_length <- nchar(Total$sequence_alignment)
Total$germ_length <- nchar(Total$germline_alignment)
Total$Input_ger<-Total$Total_length- Total$germ_length

PC<-Total[which(Total$group == "PC"),]
DP<-Total[which(Total$group == "DP"),]
DN<-Total[which(Total$group == "DN"),]
NB<-Total[which(Total$group == "NB"),]
PC_clone<-as.data.frame(table(PC$clone_id))
DP_clone<-as.data.frame(table(DP$clone_id))
DN_clone<-as.data.frame(table(DN$clone_id))
NB_clone<-as.data.frame(table(NB$clone_id))

remove(DN,DN_clone,DP,DP_clone,NB,NB_clone,PC,PC_clone)
SHARE<-read.table("E:/P23042711/2_NEW.txt",sep = "\t",header = TRUE)

Total$two<-"No_Hit"

for (i in 1:93) {
  for (j in 1:length(rownames(Total))) {
    
    if(SHARE$DN_PC[i] ==  Total$clone_id[j] ){
      Total$two[j]<-"Hit"      
    }
  }
  
}


remove(SHARE,i,j)
Total_two<-Total[which(Total$two == "Hit"),]



library(dowser)
clones_two = formatClones(Total_two, traits=c("group","c_call"))
clones_two_phylip = getTrees(clones_two, build="dnapars", exec="e:/P23042711/imgt/phylip-3.698/exe/dnapars.exe")
clones_two_1 = getTrees(clones_two)
clones_two_pml = getTrees(clones_two, build="pml")
clones_two_dnaml = getTrees(clones_two, build="dnaml", exec="e:/P23042711/imgt/phylip-3.698/exe/dnaml.exe")
plots_two_1=plotTrees(clones_two_1, tips="group")
plots_two_phylip=plotTrees(clones_two_phylip, tips="group")
plots_two_pml=plotTrees(clones_two_pml, tips="group")
plots_two_dnaml=plotTrees(clones_two_dnaml, tips="group")


for (i in 1:93) {
  a<-plots_two_1[[i]][["labels"]][["title"]]
  pdf(paste0("Total_two_1_",i,"_",a,".pdf"))
  remove(a)
  print(plots_two_1[[i]])
  dev.off()
}
for (i in 1:93) {
  a<-plots_two_phylip[[i]][["labels"]][["title"]]
  pdf(paste0("Total_two_phylip_",i,"_",a,".pdf"))
  remove(a)
  print(plots_two_phylip[[i]])
  dev.off()
}
for (i in 1:91) {
  a<-plots_two_pml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_two_pml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_two_pml[[i]])
  dev.off()
}
for (i in 1:93) {
  a<-plots_two_dnaml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_two_dnaml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_two_dnaml[[i]])
  dev.off()
}

write.table(Total_two,file = "e:/P23042711/Total_two_DN_PC.txt",sep = "\t")

library(alakazam)
library(igraph)
library(dplyr)
phylip_exec <- "e:/P23042711/imgt/phylip-3.698/exe/dnapars.exe"
Two_clone<-unique(Total_two$clone_id)
for (i in Two_clone) {
  sub_db <- subset(Total_two, clone_id == i)
  clone <- makeChangeoClone(sub_db, text_fields=c("group", "c_call") 
  )
  if(length(clone@data$sequence_id)>1){
    if(nchar(clone@data$sequence)[1] == nchar(clone@germline)){
      pdf(paste0("Alakazam_two_",i,".pdf"))
      graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
      
      # Set node colors
      V(graph)$color[V(graph)$group == "DN"] <- "seagreen"
      V(graph)$color[V(graph)$group == "PC"] <- "steelblue"
      V(graph)$color[V(graph)$name == "Germline"] <- "black"
      V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
      
      # Set node labels
      V(graph)$label <- paste(V(graph)$group, V(graph)$c_call, sep=", ")
      V(graph)$label[V(graph)$name == "Germline"] <- ""
      V(graph)$label[grepl("Inferred", V(graph)$name)] <- ""
      
      # Set node shapes
      V(graph)$shape <- "circle"
      V(graph)$shape[V(graph)$name == "Germline"] <- "circle"
      V(graph)$shape[grepl("Inferred", V(graph)$name)] <- "circle"
      
      # Set node sizes
      V(graph)$size <- 20
      V(graph)$size[V(graph)$name == "Germline"] <- 10
      V(graph)$size[grepl("Inferred", V(graph)$name)] <- 5 
      
      # Remove large default margins
      par(mar=c(0, 0, 0, 0) + 0.05)
      
      # Plot the example tree
      plot(graph, layout=layout_as_tree, vertex.frame.color="grey", 
           vertex.label.color="black", edge.label.color="black", 
           edge.arrow.mode=0) 
      # legend("topleft", c("Germline", "Inferred", "DP", "PC"), 
      #        fill=c("black", "white", "seagreen", "steelblue"), cex=0.75)
      dev.off()
    }
    
  }
  
}
remove(Two_clone,i,sub_db,graph,clone)

remove(phylip_exec)

# clones_two_igphyml = readRDS("e:/clone_two_igphyml.rds")
# plots_two_igphyml=plotTrees(clones_two_igphyml, tips="group")
# for (i in 1:86) {
#   a<-plots_two_igphyml[[i]][["labels"]][["title"]]
#   pdf(paste0("Total_two_igphyml_",i,"_",a,".pdf"))
#   remove(a)
#   print(plots_two_igphyml[[i]])
#   dev.off()
# }
