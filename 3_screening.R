library(Seurat)
library(cowplot)
set.seed(123)
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))

library(ggridges)
library(ggplot2)

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
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
remove(DN.data,DP.data,NB.data,PC.data)
remove(DN.metadata,DP.metadata,NB.metadata,PC.metadata)
remove(pbmc.data,pbmc.metadata,i)
gc()
setwd("e:/CJJscreen/")
pbmc_raw<-pbmc
remove(pbmc)
gc()
for (i in c(2500,3000,3500,4000,4500,5000,5500,6000)) {
  for (j in c(5,10,15,20,25,30,35,40)) {
    for (p in c(1e-5,1e-10,1e-20,1e-30,1e-40)) {
      
        pbmc <- subset(pbmc_raw, subset = nFeature_RNA > 200 & nFeature_RNA < i & percent.mt < j)
        Idents(pbmc)<-pbmc$group
        
        
        DN<-subset(pbmc,idents = c("DN"))
        DP<-subset(pbmc,idents = c("DP"))
        PC<-subset(pbmc,idents = c("PC"))
      
        for (q in seq(500,min(ncol(DN),ncol(DP),ncol(PC))%/%100*100,100)) {
        DN_new<-DN[,sample(1:ncol(DN),q)]
        DP_new<-DP[,sample(1:ncol(DP),q)]
        PC_new<-PC[,sample(1:ncol(PC),q)]
        pbmc_new <- merge( DN_new,y = c( DP_new,PC_new), add.cell.ids = c( "DN", "DP","PC"), project = "Total")
        remove(DN_new,DP_new,PC_new)
        gc()
        Idents(pbmc_new)<-pbmc_new$group
        pbmc_new_new<-subset(pbmc_new,idents =c("DN","DP","PC"))
        
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
        diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~group",cores = 20)
        ordering_genes <- row.names (subset(diff_test_res, qval < p))
        monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
        monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)
        
        monocle_cds <-orderCells(monocle_cds)
        plotdf=pData(monocle_cds)
        mycolor<-c("#619CFF","#00BA38","#F8766D")
        pdf(paste("vocano",i,j,p,q,".pdf",sep = "_"))
        g=ggplot(plotdf, aes(x=Pseudotime,y=group,fill=group))+
          geom_density_ridges(scale=1) +
          geom_vline(xintercept = c(5,10),linetype=2)+
          scale_y_discrete("")+
          theme_minimal()+
          theme(
            panel.grid = element_blank()
          )+scale_fill_manual(values = mycolor)
        print(g)
        dev.off()
        }
        remove(DM,DP,PC)
        gc()
    }
  }
}
