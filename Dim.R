#!/usr/bin/Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(reticulate)
alldata <- readRDS("data/scRNA_qc_remove_doublets.rds")

suppressWarnings(suppressMessages(alldata <- FindVariableFeatures(alldata, selection.method = "vst",
    nfeatures = 2000, verbose = FALSE, assay = "RNA")))


alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA", "S.Score", "G2M.Score"),
    assay = "RNA")
alldata <- RunPCA(alldata, npcs = 50, verbose = F)
alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30, perplexity = 30, max_iter = 1000,
    theta = 0.5, eta = 200, num_threads = 10)
alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
    n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)

alldata <- RunUMAP(alldata, reduction.name = "UMAP10_on_PCA", reduction = "pca",
    dims = 1:30, n.components = 10, n.neighbors = 30, n.epochs = 200, min.dist = 0.3,
    learning.rate = 1, spread = 1)

alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_ScaleData", features = alldata@assays$RNA@var.features,
    assay = "RNA", n.components = 2, n.neighbors = 30, n.epochs = 200, min.dist = 0.3,
    learning.rate = 1, spread = 1)
# Build Graph
alldata <- FindNeighbors(alldata, reduction = "pca", graph.name = "SNN", assay = "RNA",
    k.param = 20, features = alldata@assays$RNA@var.features)

# Run UMAP on a graph
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph", graph = "SNN", assay = "RNA")

saveRDS(alldata, "data/Dim.rds")
