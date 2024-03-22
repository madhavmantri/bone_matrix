# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(fitdistrplus); library(ggplot2)
library(nichenetr)
library(EnhancedVolcano)
library(clusterProfiler)

setwd("/fs/cbsuvlaminck3/workdir/mm2937/FishbachLab/")

cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")

load("Robjs_ImmunePos/seurat.object.immunePos.Robj")

####################################  Monocyte cell analysis starts here! ################################################
DefaultAssay(seurat.object) <- "RNA"
fibroblasts <- subset(seurat.object, subset = RNA_Celltype == "Fibroblasts")
fibroblasts <- NormalizeData(fibroblasts) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
fibroblasts <- CellCycleScoring(fibroblasts, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
ElbowPlot(object = fibroblasts, ndims = 50)
DimPlot(fibroblasts, reduction = "pca", group.by = "Phase")
DimPlot(fibroblasts, reduction = "pca", group.by = "Sample")
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:20, force.recalc = TRUE)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.5)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.3)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.2)

table(Idents(fibroblasts))
fibroblasts <- RunTSNE(object = fibroblasts, reduction = "pca", dims = 1:20)

DimPlot(fibroblasts, reduction = "tsne", group.by = "Phase")
DimPlot(fibroblasts, reduction = "tsne", group.by = "RNA_snn_res.0.2") 
DimPlot(fibroblasts, reduction = "tsne", group.by = "Sample")
DimPlot(fibroblasts, reduction = "tsne", group.by = "RNA_snn_res.0.2") | DimPlot(fibroblasts, reduction = "tsne", group.by = "Sample") |  DimPlot(fibroblasts, reduction = "tsne", group.by = "Phase")
FeaturePlot(fibroblasts, features = "Assignment_Probability", ncol = 1)
FeaturePlot(fibroblasts, features = c("Dcn", "Postn", "Tcf21", "Col3a1", "Acta2"), ncol = 2)

fibroblasts$Assignment_Probability
fibroblasts <- subset(fibroblasts, subset = RNA_snn_res.0.2 %in% c("0", "1", "3"))

# save(fibroblasts, file="./Robjs_ImmunePos/fibroblasts.Robj")
load("./Robjs_ImmunePos/fibroblasts.Robj")

Idents(fibroblasts) <- fibroblasts$RNA_snn_res.0.2
fibroblasts.markers <- FindAllMarkers(fibroblasts, logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T)
markers.top10.fibroblasts = fibroblasts.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.fibroblasts = fibroblasts.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DotPlot(fibroblasts, features = markers.top10.fibroblasts$gene) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0))

summary <- as.data.frame.matrix(table(fibroblasts@meta.data[c("Sample", "RNA_snn_res.0.2")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity") + 
  labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic() + theme(text = element_text(size = 14))


CAF_markers <- FindMarkers(fibroblasts, assay = "RNA", group.by = "RNA_snn_res.0.2", ident.1 = "1")
EnhancedVolcano(CAF_markers,
                lab = rownames(CAF_markers),
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, labSize = 3,
                pointSize = 2.0, FCcutoff = 1.0, colAlpha = 0.8) + xlim(0.8,3.6)



DimPlot(monocytes, reduction = "umap", group.by = "RNA_snn_res.0.2")
FeaturePlot(fibroblasts, c( "Acta2", "Tagln", "Myl9", "Ctgf"))


####################################  Monocyte cell analysis starts here! ################################################
DefaultAssay(monocytes) <- "RNA"
tams <- subset(monocytes, subset = RNA_snn_res.0.2 == "1")

tams <- subset(tams, subset = Sample %in% c("Tumor_Bone", "Tumor_DMB"))


tams <- NormalizeData(tams) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
tams <- CellCycleScoring(tams, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
ElbowPlot(object = tams, ndims = 50)
DimPlot(tams, reduction = "pca", group.by = "Phase")
DimPlot(tams, reduction = "pca", group.by = "Sample")
tams <- FindNeighbors(tams, dims = 1:20, force.recalc = TRUE)
tams <- FindClusters(tams, resolution = 0.5)
tams <- FindClusters(tams, resolution = 0.3)
tams <- FindClusters(tams, resolution = 0.2)

table(Idents(tams))
tams <- RunUMAP(object = tams, reduction = "pca", dims = 1:20)
DimPlot(tams, reduction = "umap", group.by = "Sample")
DimPlot(tams, reduction = "umap", group.by = "Phase")
DimPlot(tams, reduction = "umap", group.by = "RNA_snn_res.0.2") / DimPlot(tams, reduction = "umap", group.by = "Sample")
FeaturePlot(tams, features = c("Cxcl1", "Cxcl2", "Cxcl3", "Hmox1", "Ptgs2", "Tgm2"))

Idents(tams) <- tams$RNA_snn_res.0.2
tams.markers <- FindAllMarkers(tams, logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T)
markers.top10.tams = tams.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.tams = tams.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(tams, features = markers.top10.tams$gene)
write.csv(tams.markers, file = "./csvs.round3/tams.markers.RNA_snn_res.0.2.csv", row.names = TRUE)


summary <- as.data.frame.matrix(table(tams@meta.data[c("Sample", "RNA_snn_res.0.2")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
write.csv(summary, "counts_celltypes_tams.csv")
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)

ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity") + 
  labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic()  + theme(text = element_text(size = 20))



