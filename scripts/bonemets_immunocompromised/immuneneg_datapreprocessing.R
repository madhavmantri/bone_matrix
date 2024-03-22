# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(here); library(fitdistrplus); library(ggplot2)
library(nichenetr)
library(future)

setwd("/workdir/mm2937/FishbachLab/")

load(file = "cc.genes.rda") 
cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
samples_round1= c("Tumor_Bone_111425_HNYY2BGXC_1M", "Tumor_DMB_111426_HNYY2BGXC_2C")
samples_round2 = c("Tumor_DMB_116340_H7HWFBGXF_1A", "Tumor_Bone_116341_H7HWFBGXF_2B", "Control_DMB_116342_H7HWFBGXF_3C", "Control_Bone_116343_H7HWFBGXF_4D")
cc.genes.humans <- list("s.genes" = c(paste("GRCh38-", cc.genes.updated.2019$s.genes, sep = "")), "g2m.genes" = c(paste("GRCh38-", cc.genes.updated.2019$g2m.genes, sep = ""))) 
cc.genes.mouse <- list("s.genes" = c(paste("mm10---", cc.genes.mouse$s.genes, sep = "")), "g2m.genes" = c(paste("mm10---", cc.genes.mouse$g2m.genes, sep = ""))) 

prepare_datasets <- function(counts_path, project_name, sample_id, sample_name){
  # Read data
  data <- Read10X_h5(filename = paste0(counts_path, sample_id, "/outs/raw_feature_bc_matrix.h5", sep = ""))
  seurat.object <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 1, project = project_name)
  
  # Human v/s mouse
  seurat.object$percent_HUMAN_RNA <- PercentageFeatureSet(object = seurat.object, pattern = "^GRCh38-")
  seurat.object$percent_MOUSE_RNA <- PercentageFeatureSet(object = seurat.object, pattern = "^mm10---")
  seurat.object$Sample = sample_name
  
  print(sample_name)
  print(dim(seurat.object))
  return(seurat.object)
  # return(data)
}

# tumor_bone_nc = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg_NC", samples_round1[1], "Tumor_Bone")
# tumor_dmb_nc = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg_NC", samples_round1[2], "Tumor_DMB")

tumor_dmb = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[1], "Tumor_DMB")
tumor_bone = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[2], "Tumor_Bone")
control_dmb = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[3], "Control_DMB")
control_bone = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[4], "Control_Bone")

dim(control_bone)
dim(control_dmb)
dim(tumor_bone)
dim(tumor_dmb)

min(control_bone$nCount_RNA)
min(control_dmb$nCount_RNA)
min(tumor_bone$nCount_RNA)
min(tumor_dmb$nCount_RNA)

# save.image("Robjs_ImmuneNeg/all_data.Robj")
load("Robjs_ImmuneNeg/all_data.Robj")

seurat.object <- merge(tumor_dmb, c(tumor_bone, control_dmb, control_bone), add.cell.ids = samples_round2)

human.mito.genes <- grep(pattern = "GRCh38-MT-", x = rownames(seurat.object), value = TRUE)
mouse.mito.genes <- grep(pattern = "mm10---mt-", x = rownames(seurat.object), value = TRUE)
mito_genes <- c(human.mito.genes, mouse.mito.genes)
seurat.object$percent.mito <- PercentageFeatureSet(seurat.object, features = mito_genes)

# VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "Sample", pt.size = 0.00)
seurat.object <- subset(seurat.object, subset = nFeature_RNA >= 200 & nCount_RNA >= 1000 & percent.mito <= 30)
# VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "Sample", pt.size = 0.00)

table(seurat.object$Sample)
min(control_bone$nCount_RNA)
min(control_dmb$nCount_RNA)
min(tumor_bone$nCount_RNA)
min(tumor_dmb$nCount_RNA)

human_genes <- grep("GRCh38-", x = rownames(seurat.object), value = TRUE)
mouse_genes <- grep("mm10---", x = rownames(seurat.object), value = TRUE)

# save(human_genes, mouse_genes, file = "Robjs_ImmuneNeg/gene_names.Robj")
load("Robjs_ImmuneNeg/gene_names.Robj")

seurat.object$nCount_HUMAN_RNA <- colSums(GetAssayData(seurat.object)[human_genes,])
seurat.object$nCount_MOUSE_RNA <- colSums(GetAssayData(seurat.object)[mouse_genes,])

VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA"), group.by = "Sample",)
FeatureScatter(object = seurat.object, feature1 = "nCount_HUMAN_RNA", feature2 = "nCount_MOUSE_RNA", group.by = "orig.ident")
FeatureScatter(object = seurat.object, feature1 = "percent_HUMAN_RNA", feature2 = "percent_MOUSE_RNA", group.by = "orig.ident") + 
  geom_hline(yintercept=90, linetype="dashed") +  geom_vline(xintercept = 90, linetype="dashed")

seurat.object$organism <- NA
seurat.object$organism[seurat.object$percent_HUMAN_RNA >= 90] <- "HUMAN"
seurat.object$organism[seurat.object$percent_MOUSE_RNA >= 90] <- "MOUSE"
seurat.object$organism[is.na(seurat.object$organism)] <- "DOUBLET"
table(seurat.object$organism)

FeatureScatter(object = seurat.object, feature1 = "nCount_HUMAN_RNA", feature2 = "nCount_MOUSE_RNA", group.by = "organism")
FeatureScatter(object = seurat.object, feature1 = "percent_HUMAN_RNA", feature2 = "percent_MOUSE_RNA", group.by = "organism")

# Normalize the data, then center and scale.
seurat.object <- NormalizeData(object = seurat.object)
seurat.object <- CellCycleScoring(seurat.object, s.features = c(cc.genes.mouse$s.genes, cc.genes.humans$s.genes), g2m.features = c(cc.genes.mouse$g2m.genes, cc.genes.humans$g2m.genes))
seurat.object <- FindVariableFeatures(object = seurat.object, selection.method = "vst", nfeatures = 2000)
sum(VariableFeatures(seurat.object) %in% human_genes)
sum(VariableFeatures(seurat.object) %in% mouse_genes)

colnames(seurat.object@meta.data)
# Run Principal Component Analysis.
seurat.object <- ScaleData(object = seurat.object)
dim(GetAssayData(seurat.object, slot = "scale.data"))
seurat.object <- RunPCA(object = seurat.object, npcs = 50)
ElbowPlot(object = seurat.object, ndims = 50)
n.pcs = 20
seurat.object <- FindNeighbors(seurat.object, dims = 1:n.pcs, force.recalc = TRUE)
seurat.object <- FindClusters(seurat.object, resolution = 0.10)
seurat.object <- FindClusters(seurat.object, resolution = 0.50)

table(Idents(seurat.object))
# To visualize

seurat.object <- RunUMAP(object = seurat.object, reduction = "pca", dims = 1:n.pcs)
seurat.object <- RunTSNE(object = seurat.object, reduction = "pca", dims = 1:n.pcs)

DimPlot(seurat.object, reduction = "tsne", group.by = "Phase")
DimPlot(seurat.object, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE) | DimPlot(seurat.object, reduction = "umap", group.by = "Sample")
DimPlot(seurat.object, reduction = "tsne", group.by = "RNA_snn_res.0.1", label = TRUE) | DimPlot(seurat.object, reduction = "tsne", group.by = "Sample")

FeaturePlot(seurat.object, reduction = "tsne", features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))
VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "RNA_snn_res.0.1", pt.size = 0.01)

DefaultAssay(seurat.object) <- "RNA"
seurat.object$Condition <- str_split_fixed(seurat.object$Sample, pattern = "_", n = 2)[,2] 


library(ggrepel)
Idents(seurat.object) <- seurat.object$RNA_snn_res.0.1
markers.all.RNA = FindAllMarkers(object = seurat.object, assay = "RNA", logfc.threshold = 2.0, min.pct = 0.5, min.diff.pct = 0.5)
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DotPlot(seurat.object, features = unique(markers.top5.RNA$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0)) 

Idents(seurat.object) <- seurat.object$RNA_snn_res.0.1
seurat.object <- RenameIdents(object = seurat.object, `0` = "Monocytes", `1` = "Fibroblasts-1", 
                              `2` = "Fibroblasts-2", `3` = "Smooth muscle cells", `4` = "Cancer cells", `6` = "Endothelial cells",
                              `5` = "Fibroblasts-3", `7` = "NK Cells/ T cells", `8` = "Monocytes")
seurat.object$RNA_Celltype <- Idents(seurat.object)

DimPlot(seurat.object, reduction = "tsne", group.by = "RNA_Celltype") | DimPlot(seurat.object, reduction = "tsne", group.by = "Sample")

Idents(seurat.object) <- seurat.object$RNA_snn_res.0.1
seurat.object <- RenameIdents(object = seurat.object, `0` = "Monocytes", `1` = "Fibroblasts", 
                              `2` = "Fibroblasts", `3` = "Smooth muscle", `4` = "Cancer cells", `6` = "Endothelial cells",
                              `5` = "Fibroblasts", `7` = "NK Cells/ T cells", `8` = "Monocytes")
seurat.object$Celltype_groups <- Idents(seurat.object)
table(seurat.object$Celltype_groups)

library(pals)
pdf(file = "./plots/immune_neg_umap.pdf", width = 2.5, height = 1.5)
DimPlot(seurat.object, reduction = "tsne", group.by = "Celltype_groups", cols = "Set1", pt.size = 0.001) + labs(color = "Celltypes") + theme_void(base_size = 6) + 
        theme(text = element_text(size = 6), 
        legend.position = 'right', 
        legend.key.size = unit(0.3,"cm"), legend.text = element_text(size = 6)) + ggtitle(NULL)
dev.off()

library(pals)
pdf(file = "./plots/immune_neg_umap_without_controls.pdf", width = 2.5, height = 1.5)
DimPlot(subset(seurat.object, Sample  %in% c("Tumor_Bone", "Tumor_DMB")), reduction = "tsne", group.by = "Celltype_groups", cols = "Set1", pt.size = 0.001) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), 
        legend.position = 'right', 
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL)
dev.off()



pdf(file = "./plots/immune_neg_umap_split.pdf", width = 1.5, height = 1.5)
DimPlot(seurat.object, reduction = "tsne", split.by = "Sample",cols = "Set1", pt.size = 0.001, ncol = 2) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), 
        legend.position = 'none', 
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL)
dev.off()


FeaturePlot(seurat.object, reduction = "tsne", features = c("mm10---Ms4a1", "mm10---Cd72"))


library(ggrepel)
Idents(seurat.object) <- seurat.object$RNA_Celltype
markers.all.RNA = FindAllMarkers(object = seurat.object, assay = "RNA", logfc.threshold = 2.0, min.pct = 0.5, min.diff.pct = 0.5)
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)

pdf(file = "./plots/immune_neg_markers.pdf", width = 6.5, height = 2.5)
DotPlot(seurat.object, features = unique(markers.top5.RNA$gene), cols = c("lightgray", "brown"), dot.min = .10, dot.scale = 3.0) + xlab("Genes") + ylab("Cell types") + theme(plot.margin = margin(1,1,1,1), axis.text.x = element_text(size = 6, angle = 90, hjust = 1.0), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), legend.position = "bottom", legend.key.size = unit(0.2,"cm"), legend.title = element_text(size = 6), legend.text = element_text(size = 6), legend.justification = "center", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10, 0, 0, 0)) + scale_color_gradient(low = "lightgray", high = "black", trans = "exp")
dev.off()

summary <- as.data.frame.matrix(table(seurat.object@meta.data[c("Sample", "Celltype_groups")]))
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)

pdf(file = "./plots/immune_neg_composition.pdf", width = 4.0, height = 1.5)
ggplot(data = summary_prop) + geom_bar(mapping = aes(fill = sample, y = value, x = variable), stat = "identity", width = 0.7, position = "dodge") + labs(y = "Composition", fill = "Samples", x = "Celltypes") + theme_classic() + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "top", legend.key.size = unit(0.2,"cm"), legend.text = element_text(size = 6), legend.justification = "center", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,-5,-5,-5)) + scale_fill_brewer(palette = "Set2") 
dev.off()
  
pdf(file = "./plots/immune_neg_dataqc.pdf", width = 6.5, height = 2.0)
VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "Sample", pt.size = 0.00) & theme_classic() & theme(plot.margin = margin(0,0,0,10),text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.text = element_text(size = 6), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-5,-5,-5,-5)) & scale_fill_brewer(palette = "Set2") & xlab("Condition")
dev.off()

# save(seurat.object, file="Robjs_ImmuneNeg/seurat.object.immuneNeg.exp1.Robj")
load("Robjs_ImmuneNeg/seurat.object.immuneNeg.exp1.Robj")
table(seurat.object$Sample)


dim(seurat.object)

####################################################################
seurat.object.nc <- merge(tumor_bone_nc, c(tumor_dmb_nc), add.cell.ids = samples_round1)

human.mito.genes <- grep(pattern = "GRCh38-MT-", x = rownames(seurat.object.nc), value = TRUE)
mouse.mito.genes <- grep(pattern = "mm10---mt-", x = rownames(seurat.object.nc), value = TRUE)
mito_genes <- c(human.mito.genes, mouse.mito.genes)
seurat.object.nc$percent.mito <- PercentageFeatureSet(seurat.object.nc, features = mito_genes)

VlnPlot(seurat.object.nc, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "Sample", pt.size = 0.00)
seurat.object.nc <- subset(seurat.object.nc, subset = nFeature_RNA >= 200 & nCount_RNA >= 1000 & percent.mito <= 30)
VlnPlot(seurat.object.nc, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "Sample", pt.size = 0.00)

load("Robjs_ImmuneNeg/gene_names.Robj")

seurat.object.nc$nCount_HUMAN_RNA <- colSums(GetAssayData(seurat.object.nc)[human_genes,])
seurat.object.nc$nCount_MOUSE_RNA <- colSums(GetAssayData(seurat.object.nc)[mouse_genes,])

VlnPlot(seurat.object.nc, features = c("nCount_RNA", "nFeature_RNA"), group.by = "Sample",)
FeatureScatter(object = seurat.object.nc, feature1 = "nCount_HUMAN_RNA", feature2 = "nCount_MOUSE_RNA", group.by = "Sample")
FeatureScatter(object = seurat.object.nc, feature1 = "percent_HUMAN_RNA", feature2 = "percent_MOUSE_RNA", group.by = "Sample") + 
  geom_hline(yintercept=90, linetype="dashed") +  geom_vline(xintercept = 90, linetype="dashed")

seurat.object.nc$organism <- NA
seurat.object.nc$organism[seurat.object.nc$percent_HUMAN_RNA >= 90] <- "HUMAN"
seurat.object.nc$organism[seurat.object.nc$percent_MOUSE_RNA >= 90] <- "MOUSE"
seurat.object.nc$organism[is.na(seurat.object.nc$organism)] <- "DOUBLET"
table(seurat.object.nc$organism)

FeatureScatter(object = seurat.object.nc, feature1 = "nCount_HUMAN_RNA", feature2 = "nCount_MOUSE_RNA", group.by = "organism")
FeatureScatter(object = seurat.object.nc, feature1 = "percent_HUMAN_RNA", feature2 = "percent_MOUSE_RNA", group.by = "organism")


##############################################################

seurat.object.immuneNeg = merge(seurat.object, c(seurat.object.nc))
seurat.object.immuneNeg

# Normalize the data, then center and scale.
seurat.object.immuneNeg <- NormalizeData(object = seurat.object.immuneNeg)
seurat.object.immuneNeg <- CellCycleScoring(seurat.object.immuneNeg, s.features = c(cc.genes.mouse$s.genes, cc.genes.humans$s.genes), g2m.features = c(cc.genes.mouse$g2m.genes, cc.genes.humans$g2m.genes))
seurat.object.immuneNeg <- FindVariableFeatures(object = seurat.object.immuneNeg, selection.method = "vst", nfeatures = 2000)
sum(VariableFeatures(seurat.object.immuneNeg) %in% human_genes)
sum(VariableFeatures(seurat.object.immuneNeg) %in% mouse_genes)

colnames(seurat.object.immuneNeg@meta.data)
# Run Principal Component Analysis.
seurat.object.immuneNeg <- ScaleData(object = seurat.object.immuneNeg)
dim(GetAssayData(seurat.object.immuneNeg, slot = "scale.data"))
seurat.object.immuneNeg <- RunPCA(object = seurat.object.immuneNeg, npcs = 50)
ElbowPlot(object = seurat.object.immuneNeg, ndims = 50)
n.pcs = 20
seurat.object.immuneNeg <- FindNeighbors(seurat.object.immuneNeg, dims = 1:n.pcs, force.recalc = TRUE)
seurat.object.immuneNeg <- FindClusters(seurat.object.immuneNeg, resolution = 0.10)
seurat.object.immuneNeg <- FindClusters(seurat.object.immuneNeg, resolution = 0.50)

table(Idents(seurat.object.immuneNeg))
# To visualize

seurat.object.immuneNeg <- RunUMAP(object = seurat.object.immuneNeg, reduction = "pca", dims = 1:n.pcs)
seurat.object.immuneNeg <- RunTSNE(object = seurat.object.immuneNeg, reduction = "pca", dims = 1:n.pcs)

DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "orig.ident")
DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "organism")

DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "RNA_snn_res.0.1", label = TRUE) | DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "Sample")

FeaturePlot(seurat.object.immuneNeg, reduction = "tsne", features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))
VlnPlot(seurat.object.immuneNeg, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "RNA_snn_res.0.1", pt.size = 0.01)

library(ggrepel)
Idents(seurat.object.immuneNeg) <- seurat.object.immuneNeg$RNA_snn_res.0.1
markers.all.RNA = FindAllMarkers(object = seurat.object.immuneNeg, assay = "RNA", logfc.threshold = 2.0, min.pct = 0.5, min.diff.pct = 0.5)
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DotPlot(seurat.object.immuneNeg, features = unique(markers.top5.RNA$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0)) 

Idents(seurat.object.immuneNeg) <- seurat.object.immuneNeg$RNA_snn_res.0.1
seurat.object.immuneNeg <- RenameIdents(object = seurat.object.immuneNeg, `1` = "Monocytes", `0` = "Fibroblasts-1", 
                              `2` = "Fibroblasts-2", `3` = "Smooth muscle cells", `5` = "Cancer cells", `4` = "Endothelial cells",
                              `6` = "Fibroblasts-3", `7` = "NK Cells/ T cells", `8` = "Basal cells", `9` = "Monocytes")
seurat.object.immuneNeg$RNA_Celltype <- Idents(seurat.object.immuneNeg)
DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "RNA_Celltype") | DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "Sample")

unique(seurat.object.immuneNeg$orig.ident)
seurat.object.immuneNeg$Experiment <- NA
seurat.object.immuneNeg$Experiment[seurat.object.immuneNeg$orig.ident == "Matrix_ImmuneNeg"] <- "Experiment1"
seurat.object.immuneNeg$Experiment[seurat.object.immuneNeg$orig.ident == "Matrix_ImmuneNeg_NC"] <- "Experiment2"
DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "Experiment")

FeaturePlot(seurat.object.immuneNeg, reduction = "tsne", features = c("mm10---C1qa", "mm10---Cdh5", "mm10---Postn", "mm10---Dcn", "mm10---Acta2"))

library(ggrepel)
Idents(seurat.object.immuneNeg) <- seurat.object.immuneNeg$RNA_Celltype
markers.all.RNA = FindAllMarkers(object = seurat.object.immuneNeg, assay = "RNA", logfc.threshold = 2.0, min.pct = 0.5, min.diff.pct = 0.5)
write.csv(markers.all.RNA, file = "csvs/celltype_markers_matrix_immune_neg.csv")
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DotPlot(seurat.object.immuneNeg, features = unique(markers.top5.RNA$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0))

# save(seurat.object.immuneNeg, file="Robjs_ImmuneNeg/seurat.object.immuneNeg.allexp.Robj")
load("Robjs_ImmuneNeg/seurat.object.immuneNeg.allexp.Robj")

Idents(seurat.object.immuneNeg) <- seurat.object.immuneNeg$RNA_snn_res.0.1
seurat.object.immuneNeg <- RenameIdents(object = seurat.object.immuneNeg, `1` = "Monocytes", `0` = "Fibroblasts", 
                                        `2` = "Fibroblasts", `3` = "Smooth muscle cells", `5` = "Cancer cells", `4` = "Endothelial cells",
                                        `6` = "Fibroblasts", `7` = "NK Cells/ T cells", `8` = "Basal cells", `9` = "Monocytes")
seurat.object.immuneNeg$Celltype_groups <- Idents(seurat.object.immuneNeg)
DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "Celltype_groups") / DimPlot(seurat.object.immuneNeg, reduction = "tsne", group.by = "Sample")


seurat.object.immuneNeg$Condition <- paste(seurat.object.immuneNeg$Sample, seurat.object.immuneNeg$Experiment, sep = ":")
head(seurat.object.immuneNeg@meta.data)
summary <- as.data.frame.matrix(table(seurat.object.immuneNeg@meta.data[c("Condition", "Celltype_groups")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)
ggplot(data = summary_prop) + geom_bar(mapping = aes(fill = sample, y = value, x = variable), stat = "identity", width = 0.6, position = "dodge") + labs(y = "Composition", fill = "Sample", x = "Celltypes") + theme_classic() + theme(plot.margin = margin(0,0,0,0), legend.position = "top", text = element_text(color = "black", size = 14), axis.text.x = element_text(angle = 0))

table(seurat.object.immuneNeg$orig.ident)
