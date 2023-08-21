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
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
remove(DN.data,DP.data,NB.data,PC.data)
remove(DN.metadata,DP.metadata,NB.metadata,PC.metadata)
remove(pbmc.data,pbmc.metadata,i)
Idents(pbmc)<-pbmc$group
NB<-subset(pbmc,idents = c("NB"))
DN<-subset(pbmc,idents = c("DN"))
DP<-subset(pbmc,idents = c("DP"))
PC<-subset(pbmc,idents = c("PC"))

NB_new<-NB[,sample(1:ncol(NB),2000)]
DN_new<-DN[,sample(1:ncol(DN),2000)]
DP_new<-DP[,sample(1:ncol(DP),2000)]
PC_new<-PC[,sample(1:ncol(PC),2000)]
pbmc_new<-pbmc.big <- merge(NB_new, y = c(DN_new, DP_new,PC_new), add.cell.ids = c("NB", "DN", "DP","PC"), project = "Total")
pbmc_new <- NormalizeData(pbmc_new)
pbmc_new <- FindVariableFeatures(pbmc_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_new)
pbmc_new <- ScaleData(pbmc_new, features = all.genes)
pbmc_new <- RunPCA(pbmc_new, features = VariableFeatures(object = pbmc_new))
ElbowPlot(pbmc_new)
pbmc_new <- FindNeighbors(pbmc_new, dims = 1:20)
pbmc_new <- FindClusters(pbmc_new, resolution = 1.2)
pbmc_new <- RunUMAP(pbmc_new, dims = 1:20)
pbmc_new <- RunTSNE(pbmc_new, dims = 1:20)
DimPlot(pbmc_new, reduction = "umap")
DimPlot(pbmc_new, reduction = "umap",split.by = "group")
DimPlot(pbmc_new, reduction = "tsne")
DimPlot(pbmc_new, reduction = "tsne",split.by = "group")
DotPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"))
VlnPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"),pt.size = 0,sort = TRUE)
DimPlot(pbmc_new,reduction = "pca")
DimPlot(pbmc_new,reduction = "pca",split.by = "group")
pbmc_new <- RunLDA(pbmc_new, labels = pbmc_new$group)
pbmc_new <- RunUMAP(pbmc_new ,reduction = "lda",dims = 1:3,reduction.name = "lda_umap")
pbmc_new <- RunTSNE(pbmc_new, reduction = "lda",dims = 1:3,reduction.name = "lda_tsne" )
Idents(pbmc_new)<-pbmc_new$group
DimPlot(pbmc_new,reduction = "lda")
DimPlot(pbmc_new,reduction = "lda_umap")
DimPlot(pbmc_new,reduction = "lda_tsne")
pbmc_new_new<-subset(pbmc_new,idents =c("DN","DP","PC"))
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc_new_new@assays$RNA@counts)
pd <-pbmc_new_new@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
                              expressionFamily = VGAM::negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~group",cores = 20)

pbmcmarkers_new<-pbmc.marker[which(pbmc.marker$p_val_adj < 1e-10 ),]
# pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
ordering_genes <- pbmcmarkers_new$gene
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-50))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2,norm_method = "log")

monocle_cds <-orderCells(monocle_cds,reverse = TRUE)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75,)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75,)+ facet_wrap(~group, nrow = 1)

