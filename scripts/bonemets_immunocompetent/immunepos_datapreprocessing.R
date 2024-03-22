# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(fitdistrplus); library(ggplot2)
library(nichenetr)
library(EnhancedVolcano)
library(clusterProfiler)

setwd("/workdir/mm2937/FishbachLab/")

cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
# cc.genes.humans <- list("s.genes" = c(paste("HUMAN-", cc.genes.updated.2019$s.genes, sep = "")), "g2m.genes" = c(paste("HUMAN-", cc.genes.updated.2019$g2m.genes, sep = ""))) 
# cc.genes.mouse <- list("s.genes" = c(paste("MOUSE-", cc.genes.mouse$s.genes, sep = "")), "g2m.genes" = c(paste("MOUSE-", cc.genes.mouse$g2m.genes, sep = ""))) 
# samples = c("111425_HNYY2BGXC_1M_22Nov19", "111426_HNYY2BGXC_2C_22Nov19")
# samples = c("Tumor_DMB_116340_H7HWFBGXF_1A", "Tumor_Bone_116341_H7HWFBGXF_2B", "Control_DMB_116342_H7HWFBGXF_3C", "Control_Bone_116343_H7HWFBGXF_4D")
# samples = c("DT", "BT", "DC", "BC")

cmo_list = c("CMO301", "CMO302", "CMO303", "CMO304")
sample_list = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB")
  
# sample_id,cmo_ids
# BC,CMO301
# DC,CMO302
# BT,CMO303
# DT,CMO304

prepare_datasets <- function(counts_path, project_name, sample_id, sample_name){
  # Read data
  data <- Read10X(data.dir = paste0(counts_path, "Cellplex_AAAMNMYM5_10444254/outs/per_sample_outs/", 
                                       sample_id, "/count/sample_feature_bc_matrix/", sep = ""))
  seurat.object <- CreateSeuratObject(counts = data$`Gene Expression`, project = project_name)
  seurat.object[["CMO"]] <- CreateAssayObject(counts = data$`Multiplexing Capture`)
  seurat.object$Sample = sample_name
  
  # Mitrocondria
  seurat.object = PercentageFeatureSet(object = seurat.object, pattern = "^mt-", col.name = "percent.mito")
  
  # Human v/s mouse
  # seurat.object$percent_HUMAN_RNA <- PercentageFeatureSet(object = seurat.object, pattern = "^HUMAN-")
  # seurat.object$percent_MOUSE_RNA <- PercentageFeatureSet(object = seurat.object, pattern = "^MOUSE-")
  
  print(sample_name)
  print(dim(seurat.object))
  return(seurat.object)
}

# tumor_bone_nc = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg_NC", samples_round1[1], "Tumor_Bone")
# tumor_dmb_nc = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg_NC", samples_round1[2], "Tumor_DMB")
# 
# tumor_dmb = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[1], "Tumor_DMB")
# tumor_bone = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[2], "Tumor_Bone")
# control_dmb = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[3], "Control_DMB")
# control_bone = prepare_datasets("./10x_counts/", "Matrix_ImmuneNeg", samples_round2[4], "Control_Bone")

tumor_dmb = prepare_datasets("10x_counts/", "Matrix_ImmunePos", sample_list[4], "Tumor_DMB")
tumor_bone = prepare_datasets("10x_counts/", "Matrix_ImmunePos", sample_list[3], "Tumor_Bone")
control_dmb = prepare_datasets("10x_counts/", "Matrix_ImmunePos", sample_list[2], "Control_DMB")
control_bone = prepare_datasets("10x_counts/", "Matrix_ImmunePos", sample_list[1], "Control_Bone")

# save.image("Robjs_ImmunePos/all_data.Robj")
load("Robjs_ImmunePos/all_data.Robj")

seurat.object <- merge(tumor_dmb, c(tumor_bone, control_dmb, control_bone), add.cell.ids = sample_list)
dim(seurat.object)
table(seurat.object$Sample)

cmo_data = read.csv("10x_counts/Cellplex_AAAMNMYM5_10444254/outs/multi/multiplexing_analysis/assignment_confidence_table.csv", row.names = 1)
rownames(cmo_data) <- cmo_data$Barcodes
dim(cmo_data)

ge_data <- Read10X_h5(filename = paste0("10x_counts/Cellplex_AAAMNMYM5_10444254/outs/multi/count/raw_feature_bc_matrix.h5"))
seurat.object.ge_data <- CreateSeuratObject(counts = ge_data$`Gene Expression`, assay = "RNA", project = "Matrix_ImmunePos")
seurat.object.ge_data[["CMO"]] <- CreateAssayObject(counts = ge_data$`Multiplexing Capture`[cmo_list,])
dim(seurat.object.ge_data)

seurat.object.immunePos <- seurat.object.ge_data[,cmo_data$Barcodes]
seurat.object.immunePos <- AddMetaData(seurat.object.immunePos, metadata = cmo_data)
seurat.object.immunePos = PercentageFeatureSet(object = seurat.object.immunePos, pattern = "^mt-", col.name = "percent.mito")
DefaultAssay(seurat.object.immunePos)
dim(seurat.object.immunePos)

VlnPlot(seurat.object.immunePos, features = c("nCount_RNA", "nFeature_RNA", "nCount_CMO", "nFeature_CMO", "eGFP"), pt.size = 0.0, ncol = 4, group.by = "Assignment")
dim(seurat.object.immunePos)
seurat.object.immunePos <- subset(seurat.object.immunePos, subset = nCount_RNA > 0 & nCount_CMO > 0)
dim(seurat.object.immunePos)
# seurat.object <- subset(seurat.object, subset = nCount_RNA >= 1000 & percent.mito <= 30)
VlnPlot(seurat.object.immunePos, features = c("nCount_RNA", "nFeature_RNA", "nCount_CMO", "nFeature_CMO", "eGFP"), group.by = "Assignment", pt.size = 0)
dim(seurat.object.immunePos)

head(seurat.object.immunePos@meta.data)
FeatureScatter(seurat.object.immunePos, feature1 = "nCount_RNA", feature2 = "nCount_CMO", group.by = "Assignment")
FeatureScatter(seurat.object.immunePos, feature1 = "nCount_RNA", feature2 = "Assignment_Probability", group.by = "Assignment")
FeatureScatter(seurat.object.immunePos, feature1 = "nCount_CMO", feature2 = "Assignment_Probability", group.by = "Assignment")

rowSums(seurat.object.immunePos@meta.data[,c(cmo_list, "Multiplet", "Blanks")])
seurat.object.immunePos@meta.data["Max_CMO"] <- colnames(seurat.object.immunePos@meta.data[,cmo_list])[max.col(seurat.object.immunePos@meta.data[,cmo_list], ties.method="first")]
seurat.object.immunePos@meta.data["Max_Assignment"] <- colnames(seurat.object.immunePos@meta.data[,c(cmo_list, "Multiplet", "Blanks")])[max.col(seurat.object.immunePos@meta.data[,c(cmo_list, "Multiplet", "Blanks")], ties.method="first")]

table(seurat.object.immunePos$Assignment)
table(seurat.object.immunePos$Max_CMO)
table(seurat.object.immunePos$Max_Assignment)

FeatureScatter(seurat.object.immunePos, feature1 = "nCount_RNA", feature2 = "nCount_CMO", group.by = "Max_Assignment")
FeatureScatter(seurat.object.immunePos, feature1 = "nCount_RNA", feature2 = "Assignment_Probability", group.by = "Max_Assignment")
FeatureScatter(seurat.object.immunePos, feature1 = "nCount_CMO", feature2 = "Assignment_Probability", group.by = "Max_Assignment")

VlnPlot(seurat.object.immunePos, features = c("eGFP"), pt.size = 0.0001, ncol = 1, group.by = "Max_Assignment")
VlnPlot(subset(seurat.object.immunePos, subset = Assignment_Probability > 0.90), features = c("eGFP"), pt.size = 0.001, ncol = 1, group.by = "Max_Assignment")
VlnPlot(seurat.object.immunePos, features = c("Assignment_Probability"), pt.size = 0.001, ncol = 1, group.by = "Max_Assignment")

DefaultAssay(seurat.object.immunePos) <- "CMO"
# Normalize the data, then center and scale.
# seurat.object <- NormalizeData(object = seurat.object)
seurat.object.immunePos<-
  NormalizeData(seurat.object.immunePos,
                assay = "CMO",
                normalization.method = "CLR")
dim(GetAssayData(seurat.object.immunePos))

VariableFeatures(seurat.object.immunePos) <- cmo_list

# Run Principal Component Analysis.
seurat.object.immunePos <- ScaleData(object = seurat.object.immunePos)
dim(GetAssayData(seurat.object.immunePos, slot = "scale.data"))
seurat.object.immunePos <- RunPCA(object = seurat.object.immunePos)
ElbowPlot(object = seurat.object.immunePos, ndims = 50)
n.pcs = 3
seurat.object.immunePos <- FindNeighbors(seurat.object.immunePos, dims = 1:n.pcs, force.recalc = TRUE)
# seurat.object <- FindClusters(seurat.object, resolution = 0.10)
# seurat.object <- FindClusters(seurat.object, resolution = 0.30)
# seurat.object <- FindClusters(seurat.object, resolution = 0.50)
# table(Idents(seurat.object))
# To visualize

# seurat.object <- RunUMAP(object = seurat.object, reduction = "pca", dims = 1:n.pcs)
seurat.object.immunePos <- RunTSNE(object = seurat.object.immunePos, reduction = "pca", dims = 1:n.pcs, reduction.name = "tsne_cmo", reduction.key = "tsne_cmo", check_duplicates = FALSE)

# DimPlot(seurat.object, reduction = "umap", group.by = "Assignment")
# DimPlot(seurat.object, reduction = "umap", group.by = "Max_CMO", label = TRUE)
# DimPlot(seurat.object, reduction = "umap", group.by = "Max_Assignment", label = TRUE)

DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = "Assignment")
DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = "Max_CMO", label = TRUE)
DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = "Max_Assignment", label = TRUE)

seurat.object.immunePos <- HTODemux(seurat.object.immunePos, assay = "CMO", positive.quantile = 0.99)
DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = "CMO_classification.global")
DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = "CMO_classification")

VlnPlot(subset(seurat.object.immunePos, Assignment_Probability > 0.9) , features = c("eGFP"), pt.size = 0.001, ncol = 1, group.by = "Max_Assignment")
dim(subset(seurat.object.immunePos, Assignment_Probability > 0.9))

colnames(seurat.object.immunePos@meta.data)
DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = c("Assignment", "Max_Assignment")) 
FeaturePlot(seurat.object.immunePos, reduction = "tsne_cmo", features = c("Assignment_Probability", "nCount_CMO", "nCount_RNA", "nFeature_RNA"), ncol = 2)

DimPlot(seurat.object.immunePos, reduction = "tsne_cmo", group.by = c("CMO_classification", "CMO_classification.global")) 
FeaturePlot(seurat.object.immunePos, reduction = "tsne_cmo", features = c("Assignment_Probability", "nCount_CMO", "nCount_RNA", "nFeature_RNA"), ncol = 2) & scale_color_continuous(type = "viridis")

VlnPlot(seurat.object.immunePos, features = c("Assignment_Probability", "nFeature_RNA"),
        group.by = "Assignment", ncol = 2, pt.size = 0.000) / 
(VlnPlot(seurat.object.immunePos, features = c("nCount_CMO", "nCount_RNA"),
        group.by = "Assignment", ncol = 2, pt.size = 0.000) & scale_y_log10())

# save(seurat.object.immunePos, file="Robjs_ImmunePos/seurat.object.immunePos.unfiltered.Robj")
load("Robjs_ImmunePos/seurat.object.immunePos.unfiltered.Robj")

dim(seurat.object.immunePos)
DefaultAssay(seurat.object.immunePos) <- "RNA"
sum(seurat.object.immunePos$Assignment_Probability >= 0.9)
seurat.object.immunePos = subset(seurat.object.immunePos, Assignment_Probability >= 0.9)
dim(seurat.object.immunePos)

seurat.object.immunePos@meta.data
seurat.object = subset(seurat.object.immunePos, Max_Assignment %in% c("CMO301", "CMO302", "CMO303", "CMO304"))
dim(seurat.object)

DefaultAssay(seurat.object) <- "RNA"
seurat.object <- subset(seurat.object, subset = nCount_RNA >= 1000 & nFeature_RNA >= 200 & percent.mito <= 30)
VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA", "nCount_CMO", "nFeature_CMO", "eGFP"), group.by = "Max_Assignment", pt.size = 0.001)
dim(seurat.object)

# Normalize the data, then center and scale.
seurat.object <- NormalizeData(object = seurat.object)
seurat.object <- CellCycleScoring(seurat.object, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
seurat.object <- FindVariableFeatures(seurat.object) 

# Run Principal Component Analysis.
seurat.object <- ScaleData(object = seurat.object)
dim(GetAssayData(seurat.object, slot = "scale.data"))
seurat.object <- RunPCA(object = seurat.object, assay = "RNA")
ElbowPlot(object = seurat.object, ndims = 50)
n.pcs = 20
seurat.object <- FindNeighbors(seurat.object, dims = 1:n.pcs, force.recalc = TRUE)
seurat.object <- FindClusters(seurat.object, resolution = 0.10)
seurat.object <- FindClusters(seurat.object, resolution = 0.20)
seurat.object <- FindClusters(seurat.object, resolution = 0.30)
seurat.object <- FindClusters(seurat.object, resolution = 0.50)
table(Idents(seurat.object))

# To visualize
DimPlot(seurat.object, reduction = "tsne_cmo", group.by = "Max_Assignment")  | DimPlot(seurat.object, reduction = "tsne_cmo", group.by = "RNA_snn_res.0.1")
DimPlot(seurat.object, reduction = "tsne_cmo", group.by = "RNA_snn_res.0.3")

# To visualize
seurat.object <- RunTSNE(object = seurat.object, reduction = "pca", dims = 1:n.pcs)
DimPlot(seurat.object, reduction = "tsne", group.by = "Assignment")

table(seurat.object$Assignment)
Idents(seurat.object) <- seurat.object$Max_Assignment
seurat.object = RenameIdents(object = seurat.object, "CMO304" = "Tumor_DMB", 
                             "CMO303" = "Tumor_Bone", "CMO302" = "Control_DMB",
                             "CMO301" = "Control_Bone")
seurat.object$Sample <- Idents(seurat.object)
seurat.object$Sample = factor(seurat.object$Sample, levels = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB"))

DimPlot(seurat.object, reduction = "tsne", group.by = "Sample") | DimPlot(seurat.object, reduction = "tsne", group.by = "RNA_snn_res.0.1") | DimPlot(seurat.object, reduction = "tsne", group.by = "Phase")

# save(seurat.object, file="Robjs_ImmunePos/seurat.object.immunePos.0.9.Robj")
load("Robjs_ImmunePos/seurat.object.immunePos.0.9.Robj")
dim(seurat.object)
table(seurat.object$Sample)
min(seurat.object$Assignment_Probability)


library(ggrepel)
Idents(seurat.object) <- seurat.object$RNA_snn_res.0.1
markers.all.RNA = FindAllMarkers(object = seurat.object, assay = "RNA", logfc.threshold = 2.0, min.pct = 0.5, return.thresh = 0.01)
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DotPlot(seurat.object, features = unique(markers.top5.RNA$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0))

seurat.object <- subset(seurat.object, RNA_snn_res.0.1 %in% c("7"), invert = T)

Idents(seurat.object) <- seurat.object$RNA_snn_res.0.1
seurat.object <- RenameIdents(object = seurat.object, `0` = "Monocytes", `1` = "Fibroblasts", 
                              `2` = "Neutrophils", `3` = "Cancer cells", `4` = "T cells", 
                              `5` = "Smooth muscle", `6` = "Endothelial cells")
seurat.object$RNA_Celltype <- Idents(seurat.object)

Idents(seurat.object) <- seurat.object$RNA_snn_res.0.1
seurat.object <- RenameIdents(object = seurat.object, `0` = "Monocytes", `1` = "Fibroblasts", 
                              `2` = "Neutrophils", `3` = "Cancer\ncells", `4` = "T cells", 
                              `5` = "Smooth\nmuscle", `6` = "Endothelial\ncells")
seurat.object$Celltype_groups <- Idents(seurat.object)

DimPlot(seurat.object, reduction = "tsne", group.by = "Sample") | DimPlot(seurat.object, reduction = "tsne", group.by = "RNA_Celltype") | DimPlot(seurat.object, reduction = "tsne", group.by = "Phase")

DimPlot(seurat.object, reduction = "tsne", group.by = "Celltype_groups") / DimPlot(seurat.object, reduction = "tsne", group.by = "Sample") 


library(pals)
pdf(file = "./plots/immune_pos_umap.pdf", width = 2.5, height = 1.5)
DimPlot(seurat.object, reduction = "tsne", group.by = "RNA_Celltype", cols = "Dark2", pt.size = 0.001) + labs(color = "Celltypes") + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), legend.text = element_text(size = 6), 
        legend.position = 'right', 
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) 
dev.off()

dim(seurat.object)

FeaturePlot(seurat.object, features = c("eGFP"))


# save(seurat.object, file="Robjs_ImmunePos/seurat.object.immunePos.0.9.Robj")
load("Robjs_ImmunePos/seurat.object.immunePos.0.9.Robj")


library(ggrepel)
Idents(seurat.object) <- seurat.object$RNA_Celltype
markers.all.RNA = FindAllMarkers(object = seurat.object, assay = "RNA", logfc.threshold = 2.0, min.pct = 0.5, return.thresh = 0.01)
write.csv(markers.all.RNA, file = "csvs.round3/celltype_markers_matrix_immune_pos.csv")
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)

pdf(file = "./plots/immune_pos_markers.pdf", width = 6.5, height = 2.5)
DotPlot(seurat.object, features = unique(markers.top5.RNA$gene), cols = c("lightgray", "brown"), dot.min = .10, dot.scale = 3.0) + xlab("Genes") + ylab("Cell types") + theme(plot.margin = margin(1,1,1,1), axis.text.x = element_text(size = 6, angle = 90, hjust = 1.0), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), legend.position = "bottom", legend.key.size = unit(0.2,"cm"), legend.title = element_text(size = 6), legend.text = element_text(size = 6), legend.justification = "center", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10, 0, 0, 0)) + scale_color_gradient(low = "lightgray", high = "black", trans = "exp")
dev.off()

summary <- as.data.frame.matrix(table(seurat.object@meta.data[c("Sample", "Celltype_groups")]))
summary <- summary[,!colnames(summary) %in% c("NA")]
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
summary_prop$sample = factor(summary_prop$sample, levels = levels(seurat.object$Sample)) 


dim(seurat.object)

pdf(file = "./plots/immune_pos_composition.pdf", width = 4.0, height = 1.5)
ggplot(data = summary_prop) + geom_bar(mapping = aes(fill = sample, y = value, x = variable), stat = "identity", width = 0.6, position = "dodge") + labs(y = "Composition", fill = "Samples", x = "Celltypes") + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), legend.text = element_text(size = 6), axis.line = element_line(size = 0.3), legend.position = "top", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,-5,-5)) + scale_fill_brewer(palette = "Set2")
dev.off()


pdf(file = "./plots/immune_pos_dataqc.pdf", width = 6.5, height = 2.0)
VlnPlot(seurat.object, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "Sample", pt.size = 0.00) & theme_classic() & theme(plot.margin = margin(0,0,0,10),text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.text = element_text(size = 6), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-5,-5,-5,-5)) & scale_fill_brewer(palette = "Set2") & xlab("Condition")
dev.off()
