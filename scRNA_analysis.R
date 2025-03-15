install.packages("Seurat")
install.packages("tidyverse")
install.packages("hdf5r")

library(Seurat)
library(tidyverse)
library(hdf5r)
library(dplyr)

#loading data
nsclc <- Read10X_h5(filename= "C:/Users/shwet/Desktop/scRNA/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc)
cts <- nsclc$'Gene Expression'
str(cts)

#Initialize the Seurat object with the raw (non-normalized data)

nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)

View(nsclc.seurat.obj@meta.data)

#Quality Control
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

#Filtering
nsclc.seurat.obj <- subset(nsclc.seurat.obj , subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalize
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor =  10000)
# Identify variable features
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

#Id the 10 most variable genes
top10_genes <- head(VariableFeatures(nsclc.seurat.obj), 10)
str(top10_genes)

#plot variable features with and without labels
top_feat_plot <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = top_feat_plot, points = top10_genes, repel = TRUE)

#Scaling
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

#Perform Linear Dimensionality reduction
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))
print(nsclc.seurat.obj[["pca"]], dims = 1:15, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

#determine data dimensionality
ElbowPlot(nsclc.seurat.obj)

#Clustering
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:13)
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = 0.5)
#identify clusters
nsclc.seurat.obj <- SetIdent(nsclc.seurat.obj, value = "seurat_clusters")
levels(nsclc.seurat.obj)
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:13)

DimPlot(nsclc.seurat.obj, reduction = "umap", label = TRUE)
levels(nsclc.seurat.obj)
#find markers

nsclc.markers <- FindAllMarkers( nsclc.seurat.obj , only.pos = TRUE)

#Visualize first 10 and features per cluster in a heatmap
top10 <- nsclc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()
DoHeatmap(nsclc.seurat.obj, features = top10$gene) + NoLegend()

#top 2 features per cluster
top_markers <- nsclc.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 2)

print(top_markers, n =30)

#Identifying Clusters

new.cluster.ids <- c("Lung Epithelial", "CAFs", "Neuronal", "NK", "Microtubular", "Macrophage","Cancer1", "B", "Cancer2","Cancer3","Lung Epithelial","Cancer4","Endothelial","Stromal","Mast")
names(new.cluster.ids) <- levels(nsclc.seurat.obj)
nsclc.seurat.obj <- RenameIdents(nsclc.seurat.obj, new.cluster.ids)
DimPlot(nsclc.seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()





