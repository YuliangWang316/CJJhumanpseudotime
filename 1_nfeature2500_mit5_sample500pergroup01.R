library(Seurat)
library(cowplot)
set.seed(123)
PC.data <- Read10X(data.dir = "f:/Newfile3/P23042711/PC/20230630/Matrix/PC_matrix/")
DP.data <- Read10X(data.dir = "f:/Newfile3/P23042711/DP/20230630/Matrix/DP_matrix/")
DN.data <- Read10X(data.dir = "f:/Newfile3/P23042711/DN/20230630/Matrix/DN_matrix/")
NB.data <- Read10X(data.dir = "f:/Newfile3/P23042711/NB/20230630/Matrix/NB_matrix/")

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

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
remove(DN.data,DP.data,NB.data,PC.data)
remove(DN.metadata,DP.metadata,NB.metadata,PC.metadata)
remove(pbmc.data,pbmc.metadata,i)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
Idents(pbmc)<-pbmc$group
gc()
NB<-subset(pbmc,idents = c("NB"))
DN<-subset(pbmc,idents = c("DN"))
DP<-subset(pbmc,idents = c("DP"))
PC<-subset(pbmc,idents = c("PC"))

NB_new<-NB[,sample(1:ncol(NB),500)]
DN_new<-DN[,sample(1:ncol(DN),500)]
DP_new<-DP[,sample(1:ncol(DP),500)]
PC_new<-PC[,sample(1:ncol(PC),500)]
pbmc_new <- merge(NB_new, y = c(DN_new, DP_new,PC_new), add.cell.ids = c("NB", "DN", "DP","PC"), project = "Total")
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
remove(all.genes,DN,DN_new,DP,DP_new,NB,NB_new,pbmc,pbmc_Raw,PC,PC_new,pbmc.big)
remove(g2m.genes,s.genes)
gc()
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
remove(data,fd,fData,pbmc_new,pbmc_new_new,pd)
gc()
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~group",cores = 20)

pbmcmarkers_new<-pbmc.marker[which(pbmc.marker$p_val_adj < 0.05 ),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# ordering_genes <- diff_test_res$gene_short_name
# ordering_genes <- row.names (subset(diff_test_res, qval < 1e-10))
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-200))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)

monocle_cds <-orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)

plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75,)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75,)+ facet_wrap(~group, nrow = 1)

my_genes <- row.names(subset(fData(monocle_cds),
                             gene_short_name %in% c("Irf4","Ly75","Kdm6b","Prdm1")))
cds_subset <- monocle_cds[my_genes,]
# plot_genes_in_pseudotime(cds_subset, color_by = "group")
plotdf=pData(monocle_cds)
library(ggridges)
mycolor<-c("#619CFF","#00BA38","#F8766D")
ggplot(plotdf, aes(x=Pseudotime,y=group,fill=group))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = mycolor)

pbmc_new_new_new<-as.SingleCellExperiment(pbmc_new_new)

library(slingshot)

sce_6<-slingshot(pbmc_new_new_new,clusterLabels = "group",reducedDim = 'LDA_TSNE',start.clus = "DN",end.clus = c("PC"),dist.method= "slingshot") 
library(TrajectoryUtils)
a<-data.frame(sce_6$slingPseudotime_1,sce_6$slingPseudotime_2)

b<-averagePseudotime(a)
remove(a)
sce_6$slingPseudotime_3<-b
remove(b)
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_6$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce_6)$LDA_TSNE, col = plotcol, pch=16, asp = 1) 
lines(SlingshotDataSet(sce_6), lwd=2, col='black')
data<-sce_6@int_colData$reducedDims@listData$LDA_TSNE
color<-data.frame(sce_6$slingPseudotime_1,plotcol)
rownames(color)<-rownames(data)
data<-cbind(data,color)
colnames(data)[3]<-"Pseudotime"
colnames(data)[4]<-"color"
remove(color)
library(ggplot2)
library(scales)
library(ggsn)
ggplot(data, aes(x = ldatsne_1, y = ldatsne_2,fill=Pseudotime)) +
  geom_point(col = plotcol ) + 
  theme_classic() +
  scale_fill_gradientn(colors = colors)
DimPlot(pbmc_new_new,reduction = "lda_tsne")