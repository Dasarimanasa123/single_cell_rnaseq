# load lybrarie
BiocManager::install("Seurat")
library(Seurat)
library(tidyverse)
BiocManager::install("hdf5r")

library(hdf5r)
# load the NSCLC DATASET
nsclc.sparse.m <- Read10X_h5(filename = 'data/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$`Gene Expression`

#Initialize the seurat object with the raw 9normalized data)

nsclc.seurat.obj <- CreateSeuratObject(counts = cts ,project="NSCLC",min.cell=3,min.features = 200)
str(nsclc.seurat.obj)

#1. qc(Quality control)
View(nsclc.seurat.obj@meta.data)
library(ggplot2)
# % mt reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj , pattern = "^MT-")

View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj,features =c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

FeatureScatter(nsclc.seurat.obj,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")+ geom_smooth(method = "lm")

# 2. filtering
nsclc.seurat.obj <- subset(nsclc.seurat.obj ,subset = nFeature_RNA >200 & nFeature_RNA < 2500 & percent.mt<5)

# 3. Normaliz data
# nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj ,normalization.method="logNormalize",scale.factor=1000)
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)
# Identification of highly variable features
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj ,selection.method = "vst",nfeatures = 2000)

# identification of most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj),10)

#plot variable features with and withoutlabels

plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot= plot1,points = top10 ,repel =TRUE)

#?LabelPoints
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj ,features = all.genes)


str(nsclc.seurat.obj)

# perform liniear dimentionlity reduction

nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj,features = VariableFeatures(object = nsclc.seurat.obj))

# visualize pca results
print(nsclc.seurat.obj[["pca"]], dims = 1:5,nfeatures =5)
DimHeatmap(nsclc.seurat.obj , dim = 1 ,cells = 500 , balanced = TRUE)

#determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)

#clustaring
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj , dims = 1:15)

#understanding clusters

nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj , resolution = c(0.1,0.3,0.5,0.7,1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj , group.by = "RNA_snn_res.0.5",label =TRUE)

# scaling identi of clusters
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

#non-linear dimensionality reduction

?reticulate
install.packages("reticulate")
reticulate::py_install(packages =  "umap-learn")

library(umap-learn)
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj ,dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")


saveRDS(nsclc.seurat.obj@meta.data, file = "/data/nsclc.rds")


