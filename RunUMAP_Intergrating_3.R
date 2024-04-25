#!/usr/bin/Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(reticulate)
library(Matrix)
library(DoubletFinder)

source("/home/wli5/05_patient1_patient2/custom_seurat_functions.R")
#use_python("/Users/wli5/miniconda3/bin/")
#use_python("/Users/wli5/miniconda3/lib/python3.9/")
#py_install("umap-learn")

## Patient 1 is multiomic data including gene expression and peaks. It needs to specify the matrix 
dat_01_pre <- Read10X_h5("/home/wli5/02_Seurat/02_cellranger_ARC/Run_01_pre/outs/filtered_feature_bc_matrix.h5")
dat_01_post <- Read10X_h5("/home/wli5/02_Seurat/02_cellranger_ARC/Run_01_post/outs/filtered_feature_bc_matrix.h5")
dat_01_8h <- Read10X_h5("/home/wli5/02_Seurat/02_cellranger_ARC/Run_01_8h/outs/filtered_feature_bc_matrix.h5")
dat_01_24h <- Read10X_h5("/home/wli5/02_Seurat/02_cellranger_ARC/Run_01_24h/outs/filtered_feature_bc_matrix.h5")
## Patient 2 is snRNA data. It dosen't need to specify.
dat_02_pre <- Read10X_h5("/home/wli5/04_patient2_3_GEX_data/20220420_Nigam_GEX/Run_02_pre_scRNA/outs/filtered_feature_bc_matrix.h5")
dat_02_post <- Read10X_h5("/home/wli5/04_patient2_3_GEX_data/20220420_Nigam_GEX/Run_02_post_scRNA/outs/filtered_feature_bc_matrix.h5")
dat_02_8h <- Read10X_h5("/home/wli5/04_patient2_3_GEX_data/20220420_Nigam_GEX/Run_02_6h_scRNA/outs/filtered_feature_bc_matrix.h5")
dat_02_24h <- Read10X_h5("/home/wli5/04_patient2_3_GEX_data/20220420_Nigam_GEX/Run_02_24h_scRNA/outs/filtered_feature_bc_matrix.h5")
### Patient 3 is mutiomic data including gene expression and peask. It needs to specify the matrix.
dat_03_pre <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_03_pre_arc/outs/filtered_feature_bc_matrix.h5")
dat_03_post <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_03_post_arc/outs/filtered_feature_bc_matrix.h5")
dat_03_8h <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_03_08h_arc/outs/filtered_feature_bc_matrix.h5")
dat_03_24h <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_03_24h_arc/outs/filtered_feature_bc_matrix.h5")
### Patient 4 is mutiomic data including gene expression and peask. It needs to specify the matrix.
dat_04_pre <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_04_pre_arc/outs/filtered_feature_bc_matrix.h5")
dat_04_post <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_04_post_arc/outs/filtered_feature_bc_matrix.h5")
dat_04_8h <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_04_08h_arc/outs/filtered_feature_bc_matrix.h5")
dat_04_24h <- Read10X_h5("/home/wli5/07_patient3_patient4/05_cell_ranger/Run_04_24h_arc/outs/filtered_feature_bc_matrix.h5")
### Creating Seurat objects. ####
SPre01 <- CreateSeuratObject(counts = dat_01_pre$`Gene Expression`, project = "01_pre", min.cells = 3, min.features = 200)
SPost01 <- CreateSeuratObject(counts = dat_01_post$`Gene Expression`, project = "01_post", min.cells = 3, min.features = 200)
S8h01 <- CreateSeuratObject(counts = dat_01_8h$`Gene Expression`, project = "01_8h", min.cells = 3, min.features = 200)
S24h01 <- CreateSeuratObject(counts = dat_01_24h$`Gene Expression`, project = "01_24h", min.cells = 3, min.features = 200)
SPre02 <- CreateSeuratObject(counts = dat_02_pre, project = "02_pre", min.cells = 3, min.features = 200)
SPost02 <- CreateSeuratObject(counts = dat_02_post, project = "02_post", min.cells = 3, min.features = 200)
S8h02 <- CreateSeuratObject(counts = dat_02_8h, project = "02_8h", min.cells = 3, min.features = 200)
S24h02 <- CreateSeuratObject(counts = dat_02_24h, project = "02_24h", min.cells = 3, min.features = 200)
SPre03 <- CreateSeuratObject(counts = dat_03_pre$`Gene Expression`, project = "03_pre", min.cells = 3, min.features = 200)
SPost03 <- CreateSeuratObject(counts = dat_03_post$`Gene Expression`, project = "03_post", min.cells = 3, min.features = 200)
S8h03 <- CreateSeuratObject(counts = dat_03_8h$`Gene Expression`, project = "03_8h", min.cells = 3, min.features = 200)
S24h03 <- CreateSeuratObject(counts = dat_03_24h$`Gene Expression`, project = "03_24h", min.cells = 3, min.features = 200)
SPre04 <- CreateSeuratObject(counts = dat_04_pre$`Gene Expression`, project = "04_pre", min.cells = 3, min.features = 200)
SPost04 <- CreateSeuratObject(counts = dat_04_post$`Gene Expression`, project = "04_post", min.cells = 3, min.features = 200)
S8h04 <- CreateSeuratObject(counts = dat_04_8h$`Gene Expression`, project = "04_8h", min.cells = 3, min.features = 200)
S24h04 <- CreateSeuratObject(counts = dat_04_24h$`Gene Expression`, project = "04_24h", min.cells = 3, min.features = 200)

SPre01$time = "pre01"
SPost01$time = "post01"
S8h01$time = "8h01"
S24h01$time = "24h01"
SPre02$time = "pre02"
SPost02$time = "post02"
S8h02$time = "8h02"
S24h02$time = "24h02"
SPre03$time = "pre03"
SPost03$time = "post03"
S8h03$time = "8h03"
S24h03$time = "24h03"
SPre04$time = "pre04"
SPost04$time = "post04"
S8h04$time = "8h04"
S24h04$time = "24h04"
alldata <- merge(SPre01, c(SPost01, S8h01, S24h01,SPre02,SPost02, S8h02, S24h02,SPre03,SPost03, S8h03, S24h03,SPre04,SPost04, S8h04, S24h04), add.cell.ids = c("pre01", "post01", "8h01", "24h01","pre02", "post02", "8h02", "24h02","pre03", "post03", "8h03", "24h03","pre04", "post04", "8h04", "24h04"))
#rm(dat_01_pre,dat_01_post,dat_01_8h,dat_01_24h,dat_02_pre,dat_02_post,dat_02_8h,dat_02_24h,SPre01,SPost01,S8h01,S24h01,SPre02,SPost02,S8h02,S24h02)
#gc()

total_counts_per_cell <- colSums(alldata@assays$RNA@counts)
mito_genes <- rownames(alldata)[grep("^MT-", rownames(alldata))] # Using mito genes not tRNA or rNRA
alldata$percent_mito <- colSums(alldata@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
ribo_genes <- rownames(alldata)[grep("^RP[SL]", rownames(alldata))]
head(ribo_genes, 10)
alldata$percent_ribo <- colSums(alldata@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
pdf("Feats.pdf")
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()
FeatureScatter(alldata, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()

selected_c <- WhichCells(alldata, expression = nFeature_RNA > 200)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]
data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)
selected_mito <- WhichCells(alldata, expression = percent_mito < 0.2)
#selected_ribo <- WhichCells(alldata, expression = percent_ribo > 0.05)
#selected_hb <- WhichCells(alldata, expression = percent_hb < 0.1)

# and subset the object to only keep those cells
data.filt <- subset(alldata, features = selected_f, cells = selected_c)
data.filt <- subset(data.filt, cells = selected_mito)
#data.filt <- subset(data.filt, cells = selected_ribo)
#data.filt <- subset(data.filt, cells = selected_hb)
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ] ## Remove all MT- genes including tRNA and rRNA
#data.filt <- data.filt[!grepl("^HBB$|^HBB-2$|^HBE$|^LOC117884394$|^HBAD$|^HBAA$", rownames(data.filt)), ]

# Before running CellCycleScoring the data need to be normalized and logtransformed.
data.filt = NormalizeData(data.filt)
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes, match=rownames(data.filt))
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(data.filt))
data.filt <- CellCycleScoring(data.filt, g2m.features=g2m_genes, s.features=s_genes)

# Remove doublets
suppressMessages(require(DoubletFinder))
data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito", "S.Score", "G2M.Score"), verbose = F)
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)
sweep.res <- paramSweep_v3(data.filt) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp <- round(ncol(data.filt) * 0.04)  # expect 4% doublets
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
saveRDS(data.filt, "/home/wli5/07_patient3_patient4/data/scRNA_qc_remove_doublets.rds")
