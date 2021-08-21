library(dplyr)
library(Seurat)

#example data
pmbc.data2 <- Read10X(data.dir = "C:/Users/james/Desktop/Projects/scRNA-seq/data/filtered_gene_bc_matrices/hg19")
pbmc2 <- CreateSeuratObject(counts = pmbc.data2, project = "pbmc_B6", min.cells = 3, min.features = 200)
pbmc2

# Mouse peripheral blood mononuclear cells (PBMCs) from a C57BL/6 strain (age 8 weeks, male)
pbmc.data <- Read10X(data.dir = "C:/Users/james/Desktop/Projects/scRNA-seq/data/Mouse_PBMC_10K_Multiplex_count_raw_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = pbmc.data[["Gene Expression"]], project = "pbmc_B6", min.cells = 4, min.features = 200)
pbmc

#mouse 1 
pbmc.data3 <- Read10X(data.dir = "C:/Users/james/Desktop/Projects/scRNA-seq/data/sample_feature_bc_matrix")
pbmc3 <- CreateSeuratObject(counts = pbmc.data3[["Gene Expression"]], project = "pbmc_B6", min.cells = 3, min.features = 200)
pbmc3


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^UBR-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")



plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))




###  Noramlizing Data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


### Identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))








### Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


### PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:4, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

ElbowPlot(pbmc)

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)





### Cell Type Identification
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)


pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")







### Finding differentially expressed features (cluster biomarkers)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25)
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25)

marker_list1 <- data.frame(row.names(cluster1.markers))


head(cluster1.markers, n = 5)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


VlnPlot(pbmc, features = c("Cd74", "Lamp1"))


head(cluster1.markers, n = 15)
FeaturePlot(pbmc, features = c("Cd3e", "Cd4", "Cd8a", "Foxp3", "Cd40", "Cd44", "Cd19", "Cd74", "Cd9", "Lamp1", "Trim7"))
# Cd3epsilon, Cd4, T Cell
# Cd8a T Cell
# Foxp3 Regualtory T cell
#
# Cd40 activated B cell
# Cd19 (activated B cells), Cd74 for B Cells
# Cd74 also shows up in some T cells and Macrophages (monocytes on the right)
# Cd9 on T, B, Monocytes. Platelet also has Cd9 (the small island on the bottom right side)
# Lamp1 -> on most cells but mostly on monocyte. also on RBC
#
# 


top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()



new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
