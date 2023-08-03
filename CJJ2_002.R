Total <- read.delim("e:/P23042711/Total_DN_DP_PC_clone-pass_germ-pass.tsv", 
                    stringsAsFactors = FALSE)
Total$Total_length <- nchar(Total$sequence_alignment)
Total$germ_length <- nchar(Total$germline_alignment)
Total$Input_ger<-Total$Total_length- Total$germ_length

PC<-Total[which(Total$group == "PC"),]
DP<-Total[which(Total$group == "DP"),]
DN<-Total[which(Total$group == "DN"),]
PC_clone<-as.data.frame(table(PC$clone_id))
DP_clone<-as.data.frame(table(DP$clone_id))
DN_clone<-as.data.frame(table(DN$clone_id))

remove(DN,DN_clone,DP,DP_clone,PC,PC_clone)
SHARE<-read.table("E:/P23042711/3_new_new.txt",sep = "\t",header = TRUE)
Total$DN_PC<-"No_Hit"
Total$DN_DP<-"No_Hit"

for (i in 1:97) {
  for (j in 1:length(rownames(Total))) {
    
    if(SHARE$DN_PC[i] ==  Total$clone_id[j] ){
      Total$DN_PC[j]<-"Hit"      
    }
  }
  
}

for (i in 1:length(SHARE$DN_DP)) {
  for (j in 1:length(rownames(Total))) {
    
    if(SHARE$DN_DP[i] ==  Total$clone_id[j] ){
      Total$DN_DP[j]<-"Hit"      
    }
  }
  
}
remove(SHARE,i,j)
Total_DN_PC<-Total[which(Total$DN_PC == "Hit"),]
Total_DN_DP<-Total[which(Total$DN_DP == "Hit"),]


library(dowser)
clones_DN_PC = formatClones(Total_DN_PC, traits=c("group","c_call"))
clones_DN_PC_phylip = getTrees(clones_DN_PC, build="dnapars", exec="e:/P23042711/imgt/phylip-3.698/exe/dnapars.exe")
clones_DN_PC_1 = getTrees(clones_DN_PC)
clones_DN_PC_pml = getTrees(clones_DN_PC, build="pml")
clones_DN_PC_dnaml = getTrees(clones_DN_PC, build="dnaml", exec="e:/P23042711/imgt/phylip-3.698/exe/dnaml.exe")
plots_DN_PC_1=plotTrees(clones_DN_PC_1, tips="group")
plots_DN_PC_phylip=plotTrees(clones_DN_PC_phylip, tips="group")
plots_DN_PC_pml=plotTrees(clones_DN_PC_pml, tips="group")
plots_DN_PC_dnaml=plotTrees(clones_DN_PC_dnaml, tips="group")
# plots_DN_PC_1_c=plotTrees(clones_DN_PC_1, tips="c_call")
# plots_DN_PC_phylip_c=plotTrees(clones_DN_PC_phylip, tips="c_call")
# plots_DN_PC_pml_c=plotTrees(clones_DN_PC_pml, tips="c_call")
# plots_DN_PC_dnaml_c=plotTrees(clones_DN_PC_dnaml, tips="c_call")

clones_DN_DP = formatClones(Total_DN_DP, traits=c("group","c_call"))
clones_DN_DP_phylip = getTrees(clones_DN_DP, build="dnapars", exec="e:/P23042711/imgt/phylip-3.698/exe/dnapars.exe")
clones_DN_DP_1 = getTrees(clones_DN_DP)
clones_DN_DP_pml = getTrees(clones_DN_DP, build="pml")
clones_DN_DP_dnaml = getTrees(clones_DN_DP, build="dnaml", exec="e:/P23042711/imgt/phylip-3.698/exe/dnaml.exe")
plots_DN_DP_1=plotTrees(clones_DN_DP_1, tips="group")
plots_DN_DP_phylip=plotTrees(clones_DN_DP_phylip, tips="group")
plots_DN_DP_pml=plotTrees(clones_DN_DP_pml, tips="group")
plots_DN_DP_dnaml=plotTrees(clones_DN_DP_dnaml, tips="group")
plots_DN_DP_1_c=plotTrees(clones_DN_DP_1, tips="c_call")
plots_DN_DP_phylip_c=plotTrees(clones_DN_DP_phylip, tips="c_call")
plots_DN_DP_pml_c=plotTrees(clones_DN_DP_pml, tips="c_call")
plots_DN_DP_dnaml_c=plotTrees(clones_DN_DP_dnaml, tips="c_call")

for (i in 1:97) {
  a<-plots_DN_PC_1[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_PC_1_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_PC_1[[i]])
  dev.off()
}
for (i in 1:97) {
  a<-plots_DN_PC_phylip[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_PC_phylip_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_PC_phylip[[i]])
  dev.off()
}
for (i in 1:95) {
  a<-plots_DN_PC_pml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_PC_pml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_PC_pml[[i]])
  dev.off()
}
for (i in 1:97) {
  a<-plots_DN_PC_dnaml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_PC_dnaml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_PC_dnaml[[i]])
  dev.off()
}

for (i in 1:219) {
  a<-plots_DN_DP_1[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_1_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_1[[i]])
  dev.off()
}
for (i in 1:219) {
  a<-plots_DN_DP_phylip[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_phylip_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_phylip[[i]])
  dev.off()
}
for (i in 1:218) {
  a<-plots_DN_DP_pml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_pml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_pml[[i]])
  dev.off()
}
for (i in 1:219) {
  a<-plots_DN_DP_dnaml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_dnaml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_dnaml[[i]])
  dev.off()
}

# for (i in 1:88) {
#   a<-plots_DN_PC_1_c[[i]][["labels"]][["title"]]
#   pdf(paste0("Total_DN_PC_c_1_",i,"_",a,".pdf"))
#   remove(a)
#   print(plots_DN_PC_1_c[[i]])
#   dev.off()
# }
# for (i in 1:88) {
#   a<-plots_DN_PC_phylip_c[[i]][["labels"]][["title"]]
#   pdf(paste0("Total_DN_PC_phylip_c_",i,"_",a,".pdf"))
#   remove(a)
#   print(plots_DN_PC_phylip_c[[i]])
#   dev.off()
# }
# for (i in 1:86) {
#   a<-plots_DN_PC_pml_c[[i]][["labels"]][["title"]]
#   pdf(paste0("Total_DN_PC_pml_c_",i,"_",a,".pdf"))
#   remove(a)
#   print(plots_DN_PC_pml_c[[i]])
#   dev.off()
# }
# for (i in 1:88) {
#   a<-plots_DN_PC_dnaml_c[[i]][["labels"]][["title"]]
#   pdf(paste0("Total_DN_PC_dnaml_c_",i,"_",a,".pdf"))
#   remove(a)
#   print(plots_DN_PC_dnaml_c[[i]])
#   dev.off()
# }

for (i in 1:219) {
  a<-plots_DN_DP_1_c[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_c_1_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_1_c[[i]])
  dev.off()
}
for (i in 1:219) {
  a<-plots_DN_DP_phylip_c[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_phylip_c_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_phylip_c[[i]])
  dev.off()
}
for (i in 1:218) {
  a<-plots_DN_DP_pml_c[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_pml_c_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_pml_c[[i]])
  dev.off()
}
for (i in 1:219) {
  a<-plots_DN_DP_dnaml_c[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_dnaml_c_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_dnaml_c[[i]])
  dev.off()
}
write.table(Total_DN_DP,file = "e:/P23042711/Total_DN_DP_PC_DN_DP.txt",sep = "\t")
write.table(Total_DN_PC,file = "e:/P23042711/Total_DN_DP_PC_DN_PC.txt",sep = "\t")

library(alakazam)
library(igraph)
library(dplyr)
phylip_exec <- "e:/P23042711/imgt/phylip-3.698/exe/dnapars.exe"
DN_PC_clone<-unique(Total_DN_PC$clone_id)
for (i in DN_PC_clone) {
  sub_db <- subset(Total_DN_PC, clone_id == i)
  clone <- makeChangeoClone(sub_db, text_fields=c("group", "c_call") 
  )
  if(length(clone@data$sequence_id)>1){
    if(nchar(clone@data$sequence)[1] == nchar(clone@germline)){
      pdf(paste0("Alakazam_DN_PC_",i,".pdf"))
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
remove(DN_PC_clone,i,sub_db,graph,clone)

DN_DP_clone<-unique(Total_DN_DP$clone_id)
for (i in DN_DP_clone) {
  sub_db <- subset(Total_DN_DP, clone_id == i)
  clone <- makeChangeoClone(sub_db, text_fields=c("group", "c_call") 
  )
  if(length(clone@data$sequence_id)>1){
    if(nchar(clone@data$sequence)[1] == nchar(clone@germline)){
      pdf(paste0("Alakazam_DN_DP_",i,".pdf"))
      graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
      
      # Set node colors
      V(graph)$color[V(graph)$group == "DN"] <- "yellow"
      V(graph)$color[V(graph)$group == "DP"] <- "seagreen"
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
remove(DN_DP_clone,i,sub_db,graph,clone)
remove(phylip_exec)

clones_DN_PC_igphyml = readRDS("e:/clone_DN_DP_PC_DN_PC.rds")
plots_DN_PC_igphyml=plotTrees(clones_DN_PC_igphyml, tips="group")
for (i in 1:97) {
  a<-plots_DN_PC_igphyml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_PC_igphyml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_PC_igphyml[[i]])
  dev.off()
}
# plots_DN_PC_igphyml_c=plotTrees(clones_DN_PC_igphyml, tips="c_call")
# for (i in 1:88) {
#   a<-plots_DN_PC_igphyml_c[[i]][["labels"]][["title"]]
#   pdf(paste0("Total_DN_PC_igphyml_c_",i,"_",a,".pdf"))
#   remove(a)
#   print(plots_DN_PC_igphyml_c[[i]])
#   dev.off()
# }
clones_DN_DP_igphyml = readRDS("e:/clone_DN_DP_PC_DN_DP_partial.rds")
plots_DN_DP_igphyml=plotTrees(clones_DN_DP_igphyml, tips="group")
for (i in 1:219) {
  a<-plots_DN_DP_igphyml[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_igphyml_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_igphyml[[i]])
  dev.off()
}
plots_DN_DP_igphyml_c=plotTrees(clones_DN_DP_igphyml, tips="c_call")
for (i in 1:219) {
  a<-plots_DN_DP_igphyml_c[[i]][["labels"]][["title"]]
  pdf(paste0("Total_DN_DP_igphyml_c_",i,"_",a,".pdf"))
  remove(a)
  print(plots_DN_DP_igphyml_c[[i]])
  dev.off()
}
