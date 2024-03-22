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

load("Robjs_ImmunePos/seurat.object.immunePos.0.9.Robj")

####################################  Monocyte cell analysis starts here! ################################################


DefaultAssay(seurat.object) <- "RNA"
neutrophils <- subset(seurat.object, subset = RNA_Celltype == "Neutrophils")
neutrophils <- NormalizeData(neutrophils) %>% FindVariableFeatures() 
neutrophils <- CellCycleScoring(neutrophils, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
neutrophils <- ScaleData(neutrophils)
neutrophils <- RunPCA(neutrophils)
ElbowPlot(object = neutrophils, ndims = 50)
DimPlot(neutrophils, reduction = "pca", group.by = "Phase")
DimPlot(neutrophils, reduction = "pca", group.by = "Sample")
neutrophils <- FindNeighbors(neutrophils, dims = 1:20, force.recalc = TRUE)
neutrophils <- FindClusters(neutrophils, resolution = 0.1)

table(Idents(neutrophils))
neutrophils <- RunTSNE(object = neutrophils, reduction = "pca", dims = 1:30)

DimPlot(neutrophils, reduction = "tsne", group.by = "RNA_snn_res.0.1") 
DimPlot(neutrophils, reduction = "tsne", group.by = "Phase")
DimPlot(neutrophils, reduction = "tsne", group.by = "RNA_snn_res.0.3") | DimPlot(neutrophils, reduction = "tsne", group.by = "Sample") 
FeaturePlot(neutrophils, features = markers.top10.neutrophils$gene[markers.top10.neutrophils$cluster == "1"][1:8], ncol = 2)
FeaturePlot(neutrophils, features = c("S100a9", "S100a8"), ncol = 2)

pdf(file = "./plots/immune_pos_neutrophils_umap.pdf", width = 2.0, height = 1.2)
DimPlot(neutrophils, reduction = "tsne", group.by = "RNA_snn_res.0.2", cols = "Set1", pt.size = 1) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Neutrophil subtypes")
dev.off()


pdf(file = "./plots/immune_pos_neutrophils_samples_umap.pdf", width = 2.0, height = 1.2)
DimPlot(neutrophils, reduction = "tsne", group.by = "Sample", cols = "Set2", pt.size = 1) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Condition")
dev.off()

# save(neutrophils, file="Robjs_ImmunePos/neutrophils.Robj")
load("Robjs_ImmunePos/neutrophils.Robj")
dim(neutrophils)

FeaturePlot(neutrophils, c("Icam1", "Cxcr2", "Sell"))

pdf(file = "./plots/tan_umap_featureplots.pdf", width = 3.3, height = 1.1)
FeaturePlot(neutrophils, c("Icam1", "Cxcr2", "Sell"), cols = c("grey", "black"), ncol = 3, pt.size = 0.5, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6), title = element_text(size = 7)) 
dev.off()


pdf(file = "./plots/tan_umap_featureplots_n1.pdf", width = 3.3, height = 1.1)
FeaturePlot(neutrophils, c(c("Il1a", "Tnf", "Ccl3" )), cols = c("grey", "black"), ncol = 3, pt.size = 0.5, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6), title = element_text(size = 7)) 
dev.off()

FeaturePlot(neutrophils, c("Ccl2", "Arg1", "Ccl5"))
pdf(file = "./plots/tan_umap_featureplots_n2.pdf", width = 3.3, height = 1.1)
FeaturePlot(neutrophils, c("Ccl2", "Arg1", "Ccl5"), cols = c("grey", "black"), ncol = 3, pt.size = 0.5, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6), title = element_text(size = 7)) 
dev.off()

Idents(neutrophils) <- neutrophils$RNA_snn_res.0.2
neutrophils.markers <- FindAllMarkers(neutrophils, logfc.threshold = 0.0, return.thresh = 0.01)
write.csv(neutrophils.markers, "csvs/ImmunePos.markers.TANs.csv", row.names = T)
markers.top20.neutrophils = neutrophils.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
markers.top10.neutrophils = neutrophils.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.neutrophils = neutrophils.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

brewer.set1(2)

pdf(file = "./plots/immune_pos_neutrophil_marker_heatmap.pdf", width = 2.5, height = 1.6)
Idents(neutrophils) <- neutrophils$RNA_snn_res.0.2
DoHeatmap(neutrophils, features = unique(markers.top10.neutrophils$gene), group.colors = brewer.set1(3), draw.lines = F, size = 0) +  
  scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6, color = "black"), axis.text.y = element_text(size = 6, color = "black"), legend.text = element_text(size = 6, colour = "black"), legend.position = "bottom", legend.margin = margin(-20,0,0,0), legend.box.margin = margin(0,0,0,0), legend.key.size = unit(6,"pt"))
dev.off()


summary <- as.data.frame.matrix(table(neutrophils@meta.data[c("Sample", "RNA_snn_res.0.2")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)

pdf(file = "./plots/immune_pos_neutrophils_composition.pdf", width = 2.0, height = 1.4)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity", position = "dodge") + 
  labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 7), axis.text = element_text(size = 6, color = "black"), axis.text.x = element_text(size = 6, angle = 0, hjust = 0.5, vjust = 0.5), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),             legend.box.margin=margin(-10,-5,-10,-10)) + scale_fill_brewer(palette = "Set1")
dev.off()

# ggplot(data = summary_prop[summary_prop$variable == "5",]) + geom_bar(mapping = aes(x = sample, y = value, fill = sample), stat = "identity") + 
#   labs(x = "Sample", y = "Composition", fill = "Sample") + theme_classic() + theme(plot.margin = margin(0,0,0,0), legend.position = "right", axis.text = element_blank(), axis.text.y = element_text(size = 7, color = "black"), axis.ticks = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.key.height = unit(6, "pt"), legend.key.width = unit(6, "pt"))
# 
# ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity") + 
#   labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic() + theme(plot.margin = margin(0,0,0,0), legend.position = "right", axis.text = element_blank(), axis.text.y = element_text(size = 7, color = "black"), axis.ticks = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), legend.key.height = unit(6, "pt"), legend.key.width = unit(6, "pt"))
# dev.off()
# 
# pdf("monocyte_cluster5.pdf", width = 1.5, height = 1.4)
# ggplot(data = summary_prop[summary_prop$sample  %in% c("Tumor Bone", "Tumor DMB") & summary_prop$variable == "5",]) + geom_bar(mapping = aes(x = sample, y = value, fill = sample), stat = "identity", width = 0.6) + labs(y = "Composition", fill = "Sample") + theme_classic() + theme(plot.margin = margin(0,0,0,0), legend.position = "none", axis.text = element_text(size = 7, colour = "black"), axis.text.y = element_text(size = 7, color = "black"), axis.title.x = element_blank(), axis.title.y = element_text(size = 10), legend.key.height = unit(6, "pt"), legend.key.width = unit(6, "pt"))
# dev.off()

pdf(file = "./plots/immune_pos_volcano_tans.pdf", width = 3.2, height = 1.6)
library(EnhancedVolcano)
tan_markers <- neutrophils.markers
EnhancedVolcano(tan_markers,
                lab = rownames(tan_markers),
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2,
                pointSize = 0.3, FCcutoff = 1.0, colAlpha = 0.8) + xlim(0.1,3.5) + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
dev.off()

DimPlot(monocytes, reduction = "umap", group.by = "RNA_snn_res.0.2")
FeaturePlot(neutrophils, c( "Tnf", "Icam1", "Ccl4"))
