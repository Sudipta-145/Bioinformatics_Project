
library(BiocManager)
library(Seurat)
library(dplyr)
library(ggplot2)

setwd("C:\\Users\\WINDOWS 10\\Documents\\TTT\\Sudipta")
system("GSE176078_RAW.TAR")

#counts<-read.csv("C:\\Users\\WINDOWS 10\\Documents\\TTT\\Sudipta\\metadata.csv",header=TRUE,row.names = 1)
#seurat_object<-CreateSeuratObject(counts = counts,project = "GSE176078")

counts <- Read10X(data.dir = "C:/Users/WINDOWS 10/Documents/TTT/Sudipta/",gene.column = 1)
seurat_object <- CreateSeuratObject(counts = counts,project = "GSE176078")

seurat_object[["percent.mito"]]<-PercentageFeatureSet(seurat_object,pattern = "^MT-")
seurat_object
VlnPlot(seurat_object,features = c("nFeature_RNA","nCount_RNA","percent.mito"))
FeatureScatter(seurat_object,feature1 ="nCount_RNA",feature2 = "percent.mito" )

seurat_object<-subset(seurat_object,subset =nFeature_RNA>200 & nFeature_RNA<6000 & percent.mito<10 )
seurat_object<-NormalizeData(seurat_object,normalization.method = "LogNormalize")

seurat_object<-FindVariableFeatures(seurat_object,selection.method = "vst",nfeatures = 2000)
seurat_object<-ScaleData(seurat_object)
VlnPlot(seurat_object,features = c("nFeature_RNA","nCount_RNA","percent.mito"))

seurat_object<-RunPCA(seurat_object,features = VariableFeatures(seurat_object))
ElbowPlot(seurat_object)

seurat_object<-FindNeighbors(seurat_object,dims = 1:15)
seurat_object<-FindClusters(seurat_object,resolution = 0.5)

seurat_object<-RunUMAP(seurat_object,dims = 1:15)
DimPlot(seurat_object,reduction = "umap",label = TRUE)

markers<-FindAllMarkers(seurat_object,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
FeaturePlot(seurat_object,features = "NADK")

DimPlot(seurat_object, group.by = "celltype_major")

new_cluster_ids<-c("cells")
names(new_cluster_ids)<-levels(seurat_object)
seurat_object<-RenameIdents(seurat_object,new.cluster.ids)
DimPlot(seurat_object,reduction = "umap",label=TRUE)
