library(Seurat)
library(cowplot)
PC.data <- Read10X(data.dir = "E:/P23042711/PC/20230630/Matrix/PC_matrix/")
DP.data <- Read10X(data.dir = "E:/P23042711/DP/20230630/Matrix/DP_matrix/")
DN.data <- Read10X(data.dir = "E:/P23042711/DN/20230630/Matrix/DN_matrix/")
NB.data <- Read10X(data.dir = "E:/P23042711/NB/20230630/Matrix/NB_matrix/")

PC.data <- as.data.frame(PC.data)
DP.data <- as.data.frame(DP.data)
DN.data <- as.data.frame(DN.data)
NB.data <- as.data.frame(NB.data)

for (i in 1:length(colnames(PC.data))) {
  colnames(PC.data)[i] <- paste(colnames(PC.data)[i],"PC",i,sep = "-")  
}
for (i in 1:length(colnames(DP.data))) {
  colnames(DP.data)[i] <- paste(colnames(DP.data)[i],"DP",i,sep = "-")  
}
for (i in 1:length(colnames(DN.data))) {
  colnames(DN.data)[i] <- paste(colnames(DN.data)[i],"DN",i,sep = "-")  
}
for (i in 1:length(colnames(NB.data))) {
  colnames(NB.data)[i] <- paste(colnames(NB.data)[i],"NB",i,sep = "-")  
}


PC.metadata<-data.frame(colnames(PC.data),rep("PC",length(colnames(PC.data))))
colnames(PC.metadata)<-c("barcode","group")
DP.metadata<-data.frame(colnames(DP.data),rep("DP",length(colnames(DP.data))))
colnames(DP.metadata)<-c("barcode","group")
DN.metadata<-data.frame(colnames(DN.data),rep("DN",length(colnames(DN.data))))
colnames(DN.metadata)<-c("barcode","group")
NB.metadata<-data.frame(colnames(NB.data),rep("NB",length(colnames(NB.data))))
colnames(NB.metadata)<-c("barcode","group")

pbmc.metadata<-rbind(PC.metadata,DP.metadata,DN.metadata,NB.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(PC.data,DP.data,DN.data,NB.data)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata,min.cells = 3, min.features = 200)
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 1.2)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")
DotPlot(pbmc,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"))
VlnPlot(pbmc,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"),pt.size = 0,sort = TRUE)
DimPlot(pbmc,reduction = "pca")
DimPlot(pbmc,reduction = "pca",split.by = "group")
pbmc <- RunLDA(pbmc, labels = pbmc$group)
pbmc <- RunUMAP(pbmc ,reduction = "lda",dims = 1:3,reduction.name = "lda_umap")
pbmc <- RunTSNE(pbmc, reduction = "lda",dims = 1:3,reduction.name = "lda_tsne" )
DimPlot(pbmc,reduction = "lda")
DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")
saveRDS(pbmc,file = "d:/CJJ_humanTonsil_DP_DN_NB_PC_scRNA.rds")