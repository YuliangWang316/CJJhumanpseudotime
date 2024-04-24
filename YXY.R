library(tidyverse)
library(DESeq2)
setwd("D:/kdm6a")
mycounts_lz<-read.table("LZ_counts2.txt",header = TRUE,row.names = 1,sep ="\t" )
condition_lz<-factor(c(rep("WT",3),rep("KO",2)),levels = c("WT","KO"))
colData_lz<-data.frame(row.names = colnames(mycounts_lz),condition_lz)

dds_lz <- DESeqDataSetFromMatrix(mycounts_lz, colData_lz, design= ~ condition_lz)
dds_lz <- DESeq(dds_lz)

res_lz= results(dds_lz)
res_lz = res_lz[order(res_lz$pvalue),]
head(res_lz)
summary(res_lz)
write.csv(res_lz,file="LZ2.csv")
setwd("D:/kdm6a")
mycounts_dz<-read.table("DZ_counts.txt",header = TRUE,row.names = 1,sep ="\t" )
condition_dz<-factor(c(rep("WT",3),rep("KO",4)),levels = c("WT","KO"))
colData_dz<-data.frame(row.names = colnames(mycounts_dz),condition_dz)

dds_dz <- DESeqDataSetFromMatrix(mycounts_dz, colData_dz, design= ~ condition_dz)
dds_dz <- DESeq(dds_dz)

res_dz= results(dds_dz)
res_dz = res_dz[order(res_dz$pvalue),]
head(res_dz)
summary(res_dz)
write.csv(res_dz,file="DZ.csv")
