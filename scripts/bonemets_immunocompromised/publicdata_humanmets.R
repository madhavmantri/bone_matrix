library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(here); library(fitdistrplus); library(ggplot2)
library(nichenetr)
library(future)
library(harmony)
library(EnhancedVolcano)
library(pals)
setwd("/workdir/mm2937/FishbachLab/")


bonemets_public <- read.csv(file = "./PublicData/BoneMets/GSE190772_BoM_logCounts.txt", sep = "\t")
rownames(bonemets_public) <- bonemets_public$gene
bonemets_public$gene <- NULL
max(bonemets_public)

bonemets_public_metadata <- read.csv(file = "./PublicData/BoneMets/GSE190772_BoM_MetaData.txt", sep = "\t") 
rownames(bonemets_public_metadata) <- bonemets_public_metadata$Cell
dim(bonemets_public_metadata)
colnames(bonemets_public_metadata)


public_obj <- CreateSeuratObject(counts = bonemets_public, meta.data = bonemets_public_metadata)
public_obj$Sample <- str_split_fixed(colnames(public_obj), pattern = "_", 2)[,1]
dim(public_obj)

public_obj <- FindVariableFeatures(public_obj) %>% ScaleData() %>% RunPCA()
ElbowPlot(object = public_obj, ndims = 50)
DimPlot(public_obj, reduction = "pca", group.by = "cellType")

public_obj <- FindNeighbors(public_obj, dims = 1:20, force.recalc = TRUE)
public_obj <- FindClusters(public_obj, resolution = 0.5)
public_obj <- FindClusters(public_obj, resolution = 0.3)
public_obj <- FindClusters(public_obj, resolution = 0.2)

table(Idents(public_obj))
public_obj <- RunTSNE(object = public_obj, reduction = "pca", dims = 1:30)
DimPlot(public_obj, reduction = "tsne", group.by = "Sample") | DimPlot(public_obj, reduction = "tsne", group.by = "seurat_clusters")  | DimPlot(public_obj, reduction = "tsne", group.by = "cellType") 
DimPlot(public_obj, reduction = "tsne", group.by = "cellType") 


DimPlot(public_obj[,public_obj$nCount_RNA > 1000], reduction = "tsne", group.by = "cellType")
dim(public_obj)


library(pals)
pdf(file = "./plots/publioc_bone_met_umap_without_controls.pdf", width = 5.0, height = 1.5)
(DimPlot(public_obj, reduction = "tsne", group.by = "cellType", pt.size = 0.001) + theme_void(base_size = 6)  + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), 
        legend.position = 'right', 
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Cell types")) | (DimPlot(public_obj, reduction = "tsne", group.by = "Sample", pt.size = 0.001) + theme_void(base_size = 6) + theme_void(base_size = 6) + 
                                                                                             theme(text = element_text(size = 6), 
                                                                                                   legend.position = 'right', 
                                                                                                   legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Sample"))
dev.off()


dim(public_obj)
summary <- as.data.frame.matrix(table(public_obj@meta.data[c("Sample", "cellType")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$sample)
# summary_prop$sample <- factor(summary_prop$sample, levels = rev(c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB")))
pdf(file = "./plots/immune_neg_fibroblast_composition.pdf", width = 1.9, height = 1.6)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity", position = "dodge") + labs(x = "Condition", y = "Composition", fill = "subtypes") + theme_classic() + theme(plot.margin = margin(1,1,1,1), text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.text.x = element_text(size = 6, angle = 0), axis.line = element_line(size = 0.4), legend.position = "none", legend.text = element_text(size = 6), legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0)) + scale_fill_manual("Fibroblast clusters", values = c(as.vector(alphabet()[17:18]), as.vector(alphabet()[22:23]))) + coord_flip()
dev.off()


FeaturePlot(public_obj, features = c("nFeature_RNA", "nCount_RNA"))
min(public_obj$nCount_RNA)
rownames(public_obj)
colnames(public_obj)

pdf(file = "./plots/publioc_bone_met_umap_featureplots1.pdf", width = 6.0, height = 1.5)
FeaturePlot(public_obj, c("COL8A1", "COL11A1", "THBS1", "MGP"), cols = c("grey", "black"), ncol = 4, pt.size = 0.00001, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

pdf(file = "./plots/publioc_bone_met_umap_featureplots3.pdf", width = 3.0, height = 1.5)
FeaturePlot(public_obj, c("MFAP4", "COL8A1"), cols = c("grey", "black"), ncol = 2, pt.size = 0.00001, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

pdf(file = "./plots/publioc_bone_met_umap_featureplots4.pdf", width = 4.5, height = 1.5)
FeaturePlot(public_obj, c("COL11A1", "THBS1", "MGP"), cols = c("grey", "black"), ncol = 3, pt.size = 0.00001, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()


pdf(file = "./plots/publioc_bone_met_umap_featureplots2.pdf", width = 3.0, height = 1.5)
FeaturePlot(public_obj, c("ACTA2", "MYH11"), cols = c("grey", "black"), ncol = 2, pt.size = 0.00001, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

dim(public_obj)

FeaturePlot(public_obj, c("ACTA2", "TAGLN"))

save(public_obj, file="Robjs_ImmuneNeg/public_obj.Robj")
load("Robjs_ImmuneNeg/public_obj.Robj")
dim(public_obj)


cafs <- subset(public_obj, cellType %in% c("Fibroblast"))
cafs <- FindVariableFeatures(cafs) %>% ScaleData() %>% RunPCA()
ElbowPlot(object = cafs, ndims = 50)
DimPlot(cafs, reduction = "pca", group.by = "cellType")

cafs <- FindNeighbors(cafs, dims = 1:20, force.recalc = TRUE)
cafs <- FindClusters(cafs, resolution = 0.5)
cafs <- FindClusters(cafs, resolution = 0.3)
cafs <- FindClusters(cafs, resolution = 0.2)

table(Idents(cafs))
cafs <- RunTSNE(object = cafs, reduction = "pca", dims = 1:30)
DimPlot(cafs, reduction = "tsne", group.by = "Sample") | DimPlot(cafs, reduction = "tsne", group.by = "seurat_clusters")  | DimPlot(cafs, reduction = "tsne", group.by = "cellType") 

cafs$RNA_snn_res.0.2
DimPlot(cafs, reduction = "tsne", group.by = "RNA_snn_res.0.3") 
FeaturePlot(cafs, c("PDPN", "S100A4", "TAGLN", "ACTA2"))
FeaturePlot(cafs, c("nCount_RNA", "nFeature_RNA"))

Idents(cafs) cafs$RNA_snn_res.0.2
markers.all.RNA = FindAllMarkers(object = cafs, assay = "RNA", logfc.threshold = 0.5)
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DotPlot(cafs, features = unique(markers.top10.RNA$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0))

            