## scRNA-seq analysis of human kidney lymphatic vasculature for Jafree et al. 2023
## Merged dataset derived from multiple scRNA-seq studies of human kidney (see Methods for details)
## Authors: Daniyal Jafree (University College London), Dr Benjamin Stewart (University of Cambridge)
## Version 2: 07/01/2023 - for submission to Cell Reports Medicine

#----------------------------------------------------------------------------------------------------------------#

## Load packages and set working directory. Please change working directory as required.
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)
library(ggrepel)
library(patchwork)
library(tidyselect)
library(ktplots)
library(harmony)

#----------------------------------------------------------------------------------------------------------------#

## Load dataset and set cell types as active.ident, before creating UMAP of merged dataset
load("/Users/daniyaljafree/Ben_data/Sobj/combined_data.RData")
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$celltype)
DimPlot(seurat.object, pt.size = 0.5, raster= F, label = F)  # Figure 2A
VlnPlot(seurat.object, features = c("FLT4"), raster = F, pt.size = 0, split.by = "disease_state")
table(seurat.object@meta.data$donor_ID)

## Derive number of cells by condition
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$pathology)
pathology.grouping_seurat.object <- c("Ctrl", "Ctrl", "CKD", "Rej_nonimmune", "Rej_immune", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Rej_nonimmune", "Rej_immune", "Ctrl")
names(pathology.grouping_seurat.object) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, pathology.grouping_seurat.object)
seurat.object@meta.data$pathology.grouping_seurat.object <- Idents(seurat.object)
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$celltype)
table(seurat.object@meta.data$pathology.grouping_seurat.object, seurat.object@meta.data$celltype) # Input data for Figure S1B
table(seurat.object@meta.data$pathology) # Input data for Figure S1B

## Split UMAP according to individual datasets
DimPlot(seurat.object, pt.size = 0.5, raster= F, label = F, split.by = "pathology", ncol = 6) # Figure S1A
table(seurat.object@meta.data$pathology)

## Heatmap to visualise all cell types by top two markers
my_levels <- c("lymphatic_endothelium",
               "arterial_endothelium", "GEC", "PCE", "VRE", "venular_endothelium", # ENDOTHELIUM
               "podocyte", "GPEC", "PTEC", "VCAM1+_PTEC", "DTL_LOH", "TAL_LOH", "DCT", "CNT", "collecting_duct", "PC", "IC-A", "IC-B", "pelvic_epithelium", # EPITHELIUM
               "mesangial_cell", "fibroblast", "myofibroblast", # STROMAL
               "classical_monocyte", "nonclassical_monocyte", "macrophage", "TREM2_macrophage", "cDC2", "pDC", "mast_cell", "NK1", "NK2", "NKT_cell", "CD8_T_cell", "naive_CD4_T_cell", "effector_CD4_T_cell", "Treg", "B_cell", "plasma_cell"  # IMMUNE
)
seurat.object@active.ident <- factor(x = seurat.object@active.ident, levels = my_levels)

## Subset healthy cells and compute differential expression relative to other cell types for GO term analysis
seurat.object.control <- subset(seurat.object, subset = tissue_state == "healthy")
LEC_cells_healthy <- WhichCells(seurat.object.control, ident = "lymphatic_endothelium")
LEC.DE <- FindMarkers(seurat.object.control, ident.1 = LEC_cells_healthy, only.pos = TRUE)
write.csv(LEC.DE, file = "LEC.DE.csv") # Output file used for Figure 2B

## DotPlot for visualisation of cell type specificity of markers in Figure 2C
LEC_DE_list <- c("CCL21", "FABP4", "TFF3", "MMRN1", "TFPI", "APOD", "PROX1", "CLDN5", "NRP2", "GNG11", "SNCG", "PDPN", "AKAP12", "LAMA4", "CAV1", "FABP5", "NNMT", "VIM", "ANGPT2", "EFEMP1")
DotPlot(seurat.object.control, features = LEC_DE_list)

## Export LECs from 'healthy' kidneys for comparative analysis with LECs from other organs
LEC_control <- subset(seurat.object.control, idents = "lymphatic_endothelium")
LEC_raw <- LEC_control@assays$RNA@counts
kidney_LEC <- CreateSeuratObject(LEC_raw, project = "kidney_LEC")
saveRDS(kidney_LEC, file = "/Users/daniyaljafree/Ben_data/Sobj/Crossorgancomp/Comparison/kidney_LEC.rds") # Output required for Figure 2E

## Export LECs, group by conditions by pathology and assess expression of MHC Class II-encoding transcripts
LEC <- subset(seurat.object, idents = "lymphatic_endothelium") # Use this subset for downstream analyses
LEC <- SetIdent(LEC, value = LEC@meta.data$pathology)
pathology.grouping_LEC <- c("Ctrl", "Ctrl", "CKD", "Rej_immune", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Rej_nonimmune", "Rej_immune", "Ctrl")
names(pathology.grouping_LEC) <- levels(LEC)
LEC <- RenameIdents(LEC, pathology.grouping_LEC)
LEC@meta.data$pathology.grouping_LEC <- Idents(LEC)
LEC <- SetIdent(LEC, value = LEC@meta.data$pathology)
my_levels <- c("no_rejn", "live_donor_kidney_biopsy", "3_month_transplant_bx", "12_month_transplant_bx", "healthy", "CKD", "mycotic pseudoaneurysm", "renal vein thrombus", "recurrent FSGS", "rejn", "chronic rejection")
LEC@meta.data$pathology <- factor(x = LEC@meta.data$pathology, levels = my_levels)
LEC <- SetIdent(LEC, value = LEC@meta.data$pathology.grouping_LEC)
DotPlot(LEC, features = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5"), dot.scale = 10, group.by = "pathology.grouping_LEC", col.min = 0, scale.min = 0, scale = F)
HealthyvsCKD <- FindMarkers(LEC, ident.1 = "Ctrl", ident.2 = "CKD")
HealthyvsRej_immune <- FindMarkers(LEC, ident.1 = "Ctrl", ident.2 = "Rej_immune")
HealthyvsRej_nonimmune <- FindMarkers(LEC, ident.1 = "Ctrl", ident.2 = "Rej_nonimmune")
CKDvsRej_immune <- FindMarkers(LEC, ident.1 = "CKD", ident.2 = "Rej_immune")
write.csv(HealthyvsCKD, file = "HealthyvsCKD.csv")
write.csv(HealthyvsRej_immune, file = "HealthyvsRej_immune.csv")
write.csv(HealthyvsRej_nonimmune, file = "HealthyvsRej_nonimmune.csv")
write.csv(CKDvsRej_immune, file = "CKDvsRej_immune.csv")

## Information on donor or recipient origin
table(LEC@meta.data$genotype, LEC@meta.data$pathology)

## Set-up for CellPhoneDB
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$pathology)
pathology.grouping <- c("Ctrl_comb", "Ctrl_comb", "CKD", "Rej_non", "Rej", "Ctrl_comb", "Ctrl_comb", "Rej_non", "Ctrl_comb", "Rej_non", "Rej", "Rej_non")
names(pathology.grouping) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, pathology.grouping)
seurat.object@meta.data$pathology_grouping <- Idents(seurat.object)
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$celltype)
seurat.object.immune = subset(seurat.object, idents = c("lymphatic_endothelium", "CD8_T_cell", "B_cell", "Treg", "effector_CD4_T_cell", "naive_CD4_T_cell", "plasma_cell"))
seurat.object.immune@meta.data$celltype_aggregate = paste(seurat.object.immune@active.ident, seurat.object.immune@meta.data$pathology_grouping, sep = "_")
seurat.object.immune <- SetIdent(seurat.object.immune, value = "celltype_aggregate")
## Generate output files for CellPhoneDB
count_raw <- seurat.object.immune@assays$RNA@counts
count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
write.table(count_norm, "cellphonedb_count.txt", sep="\t", quote=F)
meta_data <- cbind(rownames(seurat.object.immune@meta.data), as.vector(as.data.frame(seurat.object.immune@active.ident)[,1]))
write.table(meta_data, "cellphonedb_meta.txt", sep="\t", quote=F, row.names = F)

## Kelvin scripts
pvals <- read.delim("/Users/daniyaljafree/Ben_data/Sobj/CellPhoneDBVERSION2/out/pvalues.txt", check.names = FALSE)
regmeans <- read.delim("/Users/daniyaljafree/Ben_data/Sobj/CellPhoneDBVERSION2/out/Filtered_means/regulatory_means.txt", check.names = FALSE)
decon <- read.delim("/Users/daniyaljafree/Ben_data/Sobj/CellPhoneDBVERSION2/out/deconvoluted.txt", check.names = FALSE)
interaction_annotation <- read.delim("/Users/daniyaljafree/Ben_data/Sobj/CellPhoneDBVERSION2/out/interaction_annotation.txt", check.names = FALSE)
sce <- Seurat::as.SingleCellExperiment(seurat.object.immune)

## ABMR rejection
plot_cpdb2(cell_type1 = 'lymphatic_endothelium_Rej', cell_type2 = 'effector_CD4_T_cell_Rej|naive_CD4_T_cell_Rej|Treg_Rej',
           scdata = sce,
           idents = 'celltype_aggregate', # column name where the cell ids are located in the metadata
           means = regmeans,
           pvals = pvals,
           deconvoluted = decon,
           desiredInteractions = list(
             c('lymphatic_endothelium_Rej', 'effector_CD4_T_cell_Rej'),
             c('lymphatic_endothelium_Rej', 'naive_CD4_T_cell_Rej'),
             c('lymphatic_endothelium_Rej', 'Treg_Rej')),
           interaction_grouping = interaction_annotation,
           edge_group_colors = c(
             "Stimulatory" = "#e15759",
             "Inhibitory" = "#4e79a7"
           ),
           node_group_colors = c(
             "lymphatic_endothelium_Rej" = "black",
             "effector_CD4_T_cell_Rej" = "grey",
             'naive_CD4_T_cell_Rej' = "grey",
             "Treg_Rej" = "grey"),
           frac = 0.1,
           keep_significant_only = F,
           standard_scale = T,
           remove_self = TRUE,
           return_df = F,
)

## ABMR rejection non-CD4+ T cell types
plot_cpdb2(cell_type1 = 'lymphatic_endothelium_Rej', cell_type2 = 'CD8_T_cell_Rej|B_cell_Rej|plasma_cell_Rej',
           scdata = sce,
           idents = 'celltype_aggregate', # column name where the cell ids are located in the metadata
           means = regmeans,
           pvals = pvals,
           deconvoluted = decon,
           desiredInteractions = list(
             c('lymphatic_endothelium_Rej', 'CD8_T_cell_Rej'),
             c('lymphatic_endothelium_Rej', 'B_cell_Rej'),
             c('lymphatic_endothelium_Rej', 'plasma_cell_Rej')),
           interaction_grouping = interaction_annotation,
           edge_group_colors = c(
             "Stimulatory" = "#e15759",
             "Inhibitory" = "#4e79a7"
           ),
           node_group_colors = c(
             "lymphatic_endothelium_Rej" = "black",
             "CD8_T_cell_Rej" = "grey",
             'B_cell_Rej' = "green",
             "plasma_cell_Rej" = "grey"),
           frac = 0.1,
           keep_significant_only = F,
           standard_scale = T,
           remove_self = TRUE,
           return_df = F,
)

## CKD
plot_cpdb2(cell_type1 = 'lymphatic_endothelium_CKD', cell_type2 = 'effector_CD4_T_cell_CKD|naive_CD4_T_cell_CKD|Treg_CKD',
           scdata = sce,
           idents = 'celltype_aggregate', # column name where the cell ids are located in the metadata
           means = regmeans,
           pvals = pvals,
           deconvoluted = decon,
           desiredInteractions = list(
             c('lymphatic_endothelium_CKD', 'effector_CD4_T_cell_CKD'),
             c('lymphatic_endothelium_CKD', 'naive_CD4_T_cell_CKD'),
             c('lymphatic_endothelium_CKD', 'Treg_CKD')),
           interaction_grouping = interaction_annotation,
           edge_group_colors = c(
             "Stimulatory" = "#e15759",
             "Inhibitory" = "#4e79a7"
           ),
           node_group_colors = c(
             "lymphatic_endothelium_CKD" = "black",
             "effector_CD4_T_cell_CKD" = "grey",
             'naive_CD4_T_cell_CKD' = "grey",
             "Treg_CKD" = "grey"),
           frac = 0.1,
           keep_significant_only = F,
           standard_scale = T,
           remove_self = TRUE,
           return_df = F,
)

## Non-ABMR
plot_cpdb2(cell_type1 = 'lymphatic_endothelium_Rej_non', cell_type2 = 'effector_CD4_T_cell_Rej_non|naive_CD4_T_cell_Rej_non|Treg_Rej_non',
           scdata = sce,
           idents = 'celltype_aggregate', # column name where the cell ids are located in the metadata
           means = regmeans,
           pvals = pvals,
           deconvoluted = decon,
           desiredInteractions = list(
             c('lymphatic_endothelium_Rej_non', 'effector_CD4_T_cell_Rej_non'),
             c('lymphatic_endothelium_Rej_non', 'naive_CD4_T_cell_Rej_non'),
             c('lymphatic_endothelium_Rej_non', 'Treg_Rej_non')),
           interaction_grouping = interaction_annotation,
           edge_group_colors = c(
             "Stimulatory" = "#e15759",
             "Inhibitory" = "#4e79a7"
           ),
           node_group_colors = c(
             "lymphatic_endothelium_Rej_non" = "black",
             "effector_CD4_T_cell_Rej_non" = "grey",
             'naive_CD4_T_cell_Rej_non' = "grey",
             "Treg_Rej_non" = "grey"),
           frac = 0.1,
           keep_significant_only = F,
           standard_scale = T,
           remove_self = TRUE,
           return_df = F,
)

## Auxiliary function for plotting Dot Plots with cell-cell interaction strength
plot_cpdb(cell_type1 = 'lymphatic_endothelium_Rej', cell_type2 = 'effector_CD4_T_cell_Rej', scdata = seurat.object.immune,
          idents = 'celltype_aggregate', # column name where the cell ids are located in the metadata
          means = regmeans, pvals = pvals,
          genes = c("PVR", "LGALS9")) +
  small_axis(fontsize = 10) + small_grid() + small_guide() + small_legend(fontsize = 5)

## Subset rejecting endothelium to identify putative HEVs
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$pathology_grouping)
seurat.object.rej <- subset(seurat.object, idents = c("Rej", "Rej_non"))
seurat.object.rej <- SetIdent(seurat.object.rej, value = seurat.object.rej@meta.data$celltype)
VlnPlot(seurat.object.rej, features = c("NTAN1", "CCL21"), split.by = "pathology_grouping")
HEV <- subset(seurat.object.rej, subset = CCL21 > 1 & NTAN1 > 1, idents = c("PCE", "arterial_endothelium"))
LEC_rej <- subset(seurat.object.rej, idents = "lymphatic_endothelium")
HEV_LEC <- merge(LEC_rej, y = HEV)
HEV_LEC <- NormalizeData(HEV_LEC, normalization.method = "LogNormalize", scale.factor = 10000)
HEV_LEC <- FindVariableFeatures(HEV_LEC, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(HEV_LEC), 20)
plot1 <- VariableFeaturePlot(HEV_LEC)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(HEV_LEC)
HEV_LEC <- ScaleData(HEV_LEC, features = all.genes)
HEV_LEC <- RunPCA(HEV_LEC, features = VariableFeatures(object = HEV_LEC))
ElbowPlot(HEV_LEC)
HEV_LEC <- HEV_LEC %>%
  RunHarmony("donor_ID", plot_convergence = TRUE)
HEV_LEC <- FindNeighbors(HEV_LEC, reduction = 'harmony', dims = 1:6)
HEV_LEC <- FindClusters(HEV_LEC, resolution = 0.1)
HEV_LEC <- RunUMAP(HEV_LEC, reduction = 'harmony', dims = 1:6)
DimPlot(HEV_LEC, reduction = "umap", label = TRUE, pt.size = 2, cols = c('0' = 'green3', '1' = 'orange2'), label.size = 0)
FeaturePlot(HEV_LEC, features = c("UBD", "ACKR1", "PLVAP", "CCL21", "NTAN1"))
VlnPlot(HEV_LEC, features = c("UBD", "ACKR1", "PLVAP", "CCL21", "NTAN1"))
HEV_LEC.markers <- FindAllMarkers(HEV_LEC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HEV_LEC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(HEV_LEC.markers, file = "HEV_LEC.markers_annotated.csv")
VlnPlot(HEV_LEC, features = c("PROX1", "PDPN", "ANGPT2", "EMCN", "TEK", "FLT1"), pt.size = 0, cols = c('0' = 'green3', '1' = 'orange2')) # marker genes
VlnPlot(HEV_LEC, features = c("LYVE1", "CLDN5", "SELE", "ICAM1", "MCAM", "PLVAP"),  pt.size = 0, cols = c('0' = 'green3', '1' = 'orange2')) # adhesion molecules
VlnPlot(HEV_LEC, features = c("TFF3", "RELN", "CX3CL1", "CXCL16", "CD40", "IL32"),  pt.size = 0, cols = c('0' = 'green3', '1' = 'orange2')) # inflammatory molecules

#----------------------------------------------------------------------------------------------------#

## CROSS-ORGAN LYMPHATIC COMPARISON
# At this point, if necessary, delete irrelevant variables for downstream analysis
## Load the datasets and merge into a single Seurat object
kidney_LEC <- readRDS(file = "/Users/daniyaljafree/Ben_data/Sobj/Crossorgancomp/Comparison/kidney_LEC.rds")
heart_LEC <- readRDS(file = "/Users/daniyaljafree/Ben_data/Sobj/Crossorgancomp/Comparison/heart_LEC.rds")
skin_LEC <- readRDS(file = "/Users/daniyaljafree/Ben_data/Sobj/Crossorgancomp/Comparison/skin_LEC.rds")
lung_LEC <- readRDS(file = "/Users/daniyaljafree/Ben_data/Sobj/Crossorgancomp/Comparison/lung_LEC.rds")
sobj <- merge(kidney_LEC, y = c(heart_LEC, skin_LEC, lung_LEC))
sobj

## Normalization, scaling and PCA before Harmony batch integration
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(sobj), 20)
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
ElbowPlot(sobj)
pct <- sobj[["pca"]]@stdev / sum(sobj[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
p1 <- DimPlot(object = sobj, reduction = "pca", pt.size = .1)
p2 <- VlnPlot(object = sobj, features = "PC_1", pt.size = .1)
p1
p2
sobj$organ <- Idents(sobj)
sobj <- sobj %>% 
  RunHarmony("organ", plot_convergence = TRUE)
sobj_harmony <- Embeddings(sobj, 'harmony')
sobj_harmony[1:5, 1:5]
p1 <- DimPlot(object = sobj, reduction = "harmony", pt.size = .1, group.by = "organ")
p2 <- VlnPlot(object = sobj, features = "harmony_1", group.by = "organ", pt.size = .1)
p1
p2

## Clustering and further analysis
sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:12)
sobj <- FindClusters(sobj, resolution = 0.4)
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:12)
DimPlot(sobj, reduction = "umap", label = F, pt.size = 1)
DimPlot(sobj, reduction = "umap", label = TRUE, pt.size = 2, group.by = "organ")
Idents(sobj) <- sobj@meta.data$organ
comparison.markers <- FindAllMarkers(sobj, only.pos = F, logfc.threshold = 0.1)
comparison.markers %>% group_by(cluster) %>% top_n(n = 20)
top20 <- comparison.markers %>% group_by(cluster) %>% top_n(n = 20)
write.csv(comparison.markers,"comparison.markers.csv", row.names = FALSE)
my_levels <- c("kidney_LEC", "heart_LEC", "lung_LEC", "skin_LEC")
sobj@meta.data$organ <- factor(x = sobj@meta.data$organ, levels = my_levels)
VlnPlot(sobj, group.by = "organ", features = "LYVE1", pt.size = 0) # Violin plot for Figure 2E
VlnPlot(sobj, group.by = "organ", features = "DNASE1L3", pt.size = 0)
VlnPlot(sobj, group.by = "organ", features = "CCL14", pt.size = 0)
VlnPlot(sobj, group.by = "organ", features = c("LYVE1", "DNASE1L3", "CCL14"), pt.size = 0, stack = F)

#----------------------------------------------------------------------------------------------------#
#----------------------------------------------### END ###-------------------------------------------#
#----------------------------------------------------------------------------------------------------#



