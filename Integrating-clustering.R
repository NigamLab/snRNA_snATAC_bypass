library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(reticulate)
library(Matrix)
library(harmony)
library(future)
use_python("~/miniconda3/bin/")
source("/home/wli5/07_patient3_patient4/custom_seurat_functions.R")

alldata <- readRDS("data/Dim.rds")
plan("multicore", workers = 15)
options(future.globals.maxSize = 256000 * 1024^2)

alldata.list <- SplitObject(alldata, split.by = "orig.ident")
for (i in 1:length(alldata.list)) {
    alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst",
        nfeatures = 2000, verbose = FALSE)
}
#alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,reduction = "cca")
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,reduction = "rpca")
alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "rpca")
#alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)
alldata.int <- RunUMAP(alldata.int,reduction = "pca", dims = 1:30, verbose = F)
alldata.int <- RunTSNE(alldata.int,reduction = "pca", dims = 1:30,verbose = F)
saveRDS(alldata.int, "data/Integrating.rds")
alldata.int <- readRDS("data/Integrating.rds")
alldata.int <- FindNeighbors(alldata.int, dims = 1:30, k.param = 10, verbose = F)
alldata.int <- FindClusters(alldata.int, resolution = 0.2, verbose = F)
DimPlot(alldata.int,reduction = "umap") 
plot_integrated_clusters(alldata.int)
DimPlot(alldata.int, reduction = "umap", split.by = "orig.ident") + NoLegend()
saveRDS(alldata.int, "data/Clustering.rds")

#### Run harmony.####
alldata.harmony <- RunHarmony(alldata, group.by.vars = "orig.ident", reduction = "pca",
    dims.use = 1:30, assay.use = "RNA")
# Here we use all PCs computed from Harmony for UMAP calculation
alldata.int[["harmony"]] <- alldata.harmony[["harmony"]]
harmony_embeddings <- Embeddings(alldata.harmony, 'harmony')
harmony_embeddings[1:5, 1:5]

alldata.harmony <- alldata.harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>%
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>%
  FindClusters() %>%
  identity()
alldata.harmony <- SetIdent(alldata.harmony,value = "orig.ident")
DimPlot(alldata.harmony,reduction = "umap") + plot_annotation(title = "PBMC cells, after integration (Harmony)")
DimPlot(alldata.harmony, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()
alldata.harmony <- SetIdent(alldata.harmony,value = "seurat_clusters")
DimPlot(alldata.harmony,label = T) + NoLegend()
plot_integrated_clusters(alldata.harmony)

saveRDS(alldata.harmony, "data/Harmony-clustering.rds")

