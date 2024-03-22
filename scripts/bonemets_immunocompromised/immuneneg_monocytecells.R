library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(here); library(fitdistrplus); library(ggplot2)
library(EnhancedVolcano)
library(devtools); library(nichenetr)

setwd("/fs/cbsuvlaminck3/workdir/mm2937/FishbachLab/")

load("Robjs_ImmuneNeg/seurat.object.immuneNeg.exp1.Robj")
load("Robjs_ImmuneNeg/gene_names.Robj")
cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")

dim(seurat.object)

# Mouse analysis starts here
counts = GetAssayData(seurat.object, assay = "RNA", slot = "counts")[mouse_genes, seurat.object$organism == "MOUSE"]
rownames(counts) <- str_split_fixed(rownames(counts), pattern = "---", n = 2)[,2]
meta.data = seurat.object@meta.data[seurat.object$organism == "MOUSE", ]
mouse.object <- CreateSeuratObject(counts = counts, meta.data = meta.data)
table(mouse.object$RNA_Celltype)

# Getting monocyte subset
DefaultAssay(mouse.object) <- "RNA"
monocytes <- subset(mouse.object, subset = RNA_Celltype == "Monocytes")
monocytes <- NormalizeData(monocytes) %>% FindVariableFeatures() 
monocytes <- CellCycleScoring(monocytes, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
# monocytes <- ScaleData(monocytes)
monocytes <- ScaleData(monocytes, vars.to.regress = c("S.Score", "G2M.Score"))
monocytes <- RunPCA(monocytes)
ElbowPlot(object = monocytes, ndims = 50)
DimPlot(monocytes, reduction = "pca", group.by = "Phase")
DimPlot(monocytes, reduction = "pca", group.by = "Sample")
monocytes <- FindNeighbors(monocytes, dims = 1:20, force.recalc = TRUE)
monocytes <- FindClusters(monocytes, resolution = 0.3)
monocytes <- FindClusters(monocytes, resolution = 0.2)

table(Idents(monocytes))
monocytes <- RunTSNE(object = monocytes, reduction = "pca", dims = 1:20)
DimPlot(monocytes, reduction = "tsne", group.by = "Sample") 
DimPlot(monocytes, reduction = "tsne", group.by = "RNA_snn_res.0.2")
DimPlot(monocytes, reduction = "tsne", group.by = "Phase") 

monocytes <- subset(monocytes, subset = RNA_snn_res.0.2 %in% c("5", "6", "7"), invert = T)

FeaturePlot(monocytes, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"))
FeaturePlot(monocytes, features = c("Basp1", "Cd72", "Aif1", "Fcgr1", "Ms4a7", "Ly6a"), ncol = 3)
FeaturePlot(monocytes, features = c("Cd72", "Cd74", "Retnla", "Ccl8", "Ccl2", "Ccl12"), ncol = 3)
FeaturePlot(monocytes, features = c("Ccr2", "Vegfa", "Tie1", "Arg1"), ncol = 3)
grep("^Cxcr", x = rownames(monocytes), value = T)

Idents(monocytes) <- monocytes$RNA_snn_res.0.2
monocytes <- RenameIdents(object = monocytes, `0` = "Cd163+ macrophages", `1` = "Tnf+ macrophages", `2` = "Cd72+ macrophages",
                          `4` = "Cd9+ macrophages",
                          `3` = "Dendritic cells")
monocytes$subtypes <- Idents(monocytes)
table(monocytes$subtypes)

dim(monocytes)
FeaturePlot(monocytes, features = c("C1qa", "Cd72", "Mrc1", "Arg1", "Cd74", "Flt3", "Itgam", "Csf1r", "Cd14", "Cd19"), ncol = 4)

library(pals)
pdf(file = "./plots/immune_neg_monocyte_umap.pdf", width = 3.0, height = 2.0)
DimPlot(monocytes, reduction = "tsne", group.by = "subtypes", cols = "Set1", pt.size = 0.0001) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Monocyte subtypes")
dev.off()


library(pals)
pdf(file = "./plots/immune_neg_monocyte_sample_umap.pdf", width = 3.0, height = 2.2)
DimPlot(monocytes, reduction = "tsne", group.by = "Sample", cols = "Dark2", pt.size = 0.001)+ theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Monocyte subtypes")
dev.off()

# save(monocytes, file="./Robjs_ImmuneNeg/monocytes.Robj")
load("./Robjs_ImmuneNeg/monocytes.Robj")

VlnPlot(monocytes, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))

library(ggplot2)
library(viridis)
Idents(monocytes) <- monocytes$subtypes
monocyte.markers <- FindAllMarkers(monocytes, logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T)
markers.top20.monocytes = monocyte.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
markers.top10.monocytes = monocyte.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

pdf(file = "./plots/immune_neg_monocyte_marker_heatmap.pdf", width = 2.5, height = 3.8)
DoHeatmap(monocytes, features = unique(markers.top10.monocytes$gene), group.colors = brewer.set1(6), draw.lines = F, size = 0)  + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6, color = "black"), axis.text.y = element_text(size = 6, color = "black"), legend.text = element_text(size = 6, colour = "black"), legend.position = "bottom", legend.margin = margin(-20,0,0,0), legend.box.margin = margin(0,0,0,0), legend.key.size = unit(6,"pt"))
dev.off()


monocytes <- BuildClusterTree(object = monocytes)
PlotClusterTree(object = monocytes)

FeaturePlot(monocytes, c("Arg1", "Fcgr1"))

FeaturePlot(monocytes, c("Cd14", "Cd16", "C1qa", "Arg1", "Cd72", "Il1b", "Mmp12", "Cd9"))

levels(summary_prop$variable)
summary <- as.data.frame.matrix(table(monocytes@meta.data[c("Sample", "subtypes")]))
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
summary_prop$sample <- factor(summary_prop$sample, levels = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB"))
# summary_prop$variable <- factor(summary_prop$variable, levels = c("0", "1", "3", "2"))
pdf(file = "./plots/immune_neg_monocyte_composition.pdf", width = 2.7, height = 1.7)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity", position = "dodge") + 
  labs(x = "Sample", y = "Composition", fill = "subtypes") + theme_classic() + theme(plot.margin = margin(1,1,1,1), text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.text.x = element_text(size = 6, angle = 0), axis.line = element_line(size = 0.4), legend.position = "top", legend.text = element_text(size = 6), legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin = margin(-0,0,0,0), legend.box.margin = margin(0,0,0,0)) + scale_fill_manual("Monocyte clusters", values = as.vector(brewer.set1(6))) 
dev.off()

pdf(file = "./plots/immune_neg_tams_featureplots.pdf", width = 5.0, height = 1.3)
FeaturePlot(monocytes, c("Nos2", "Cd274", "Mmp13", "Cxcl3"), cols = c("grey", "black"), ncol = 4, pt.size = 0.00001, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()


FeaturePlot(monocytes, reduction = "tsne", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))

FeaturePlot(monocytes, c("Itgax", "Napsa", "H2-DMb2", "Il1b"))

FeaturePlot(monocytes, c("Cd14", "Cd68", "Cd33", "Itgam"), ncol = 4)
#TAM
FeaturePlot(monocytes, c("Nos2", "Cxcl3", "Tgm2", "Srgn", "Cd274", "Mmp13"), ncol = 3)
#M2
FeaturePlot(monocytes, c("Cd163", "Mrc1"), ncol = 1)
#M1
FeaturePlot(monocytes, c("Cd80", "Cd86", "Fcgr1", "Fcgr3", "Cd38", "Fpr2"), ncol = 3)

# save(monocytes, file="./Robjs_ImmuneNeg/monocytes.Robj")
# load("./Robjs_ImmuneNeg/monocytes.Robj")


pdf(file = "./plots/immune_neg_pos_tams_featureplots.pdf", width = 4.0, height = 1.6)
FeaturePlot(monocytes, c("Tgm2", "Nos2", "Il1b"), cols = c("grey", "firebrick"), ncol = 3, pt.size = 0.001) & theme_void(base_size = 6) & theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 10)) 
dev.off()

tam_markers <- FindMarkers(monocytes, assay = "RNA", group.by = "RNA_snn_res.0.2", logfc.threshold = 0.0, ident.1 = "2", only.pos = T)
write.csv(tam_markers, file = "./csvs/ImmuneNeg.markers.TAMs.csv", row.names = T)

pdf(file = "./plots/immune_neg_volcano_tams.pdf", width = 3.0, height = 1.6)
library(EnhancedVolcano)
EnhancedVolcano(tam_markers,
                lab = rownames(tam_markers), 
                # selectLab = tam_pos[!tam_pos %in% DMB_pos],
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2.0, 
                pointSize = 0.3, FCcutoff = 1.0, colAlpha = 0.8) + xlim(0.10,2.2)  + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
dev.off()


monocytes$Condition = str_split_fixed(monocytes$Sample, pattern = "_", n = 2)[,2]
monocytes$Tumor = str_split_fixed(monocytes$Sample, pattern = "_", n = 2)[,1]
Idents(monocytes) <- monocytes$Condition
DMB_conserved_markers <- FindConservedMarkers(monocytes, ident.1 = "DMB", ident.2 = "Bone", grouping.var = "Tumor", logfc.threshold = 0.0)
DMB_conserved_markers$product_log_FC <- DMB_conserved_markers$Tumor_avg_log2FC * DMB_conserved_markers$Control_avg_log2FC
DMB_conserved_markers$product_adj_pvals <- DMB_conserved_markers$Tumor_p_val_adj * DMB_conserved_markers$Control_p_val_adj
write.csv(DMB_conserved_markers, file = "./csvs/ImmuneNeg.markers.monocytes_DMBvsBone_conserved.csv", row.names = T)

DMB_conserved_markers = read.csv("./csvs/ImmuneNeg.markers.monocytes_DMBvsBone_conserved.csv", row.names = 1)
DMB_conserved_markers = DMB_conserved_markers[DMB_conserved_markers$product_log_FC > 0.0,]
DoHeatmap(monocytes, features = rownames(DMB_conserved_markers)[1:30], group.by = "Sample")

pdf(file = "./plots/immune_neg_volcano_monocytes_DMB_bone.pdf", width = 3.0, height = 1.6)
EnhancedVolcano(DMB_conserved_markers,
                lab = rownames(DMB_conserved_markers),
                x = 'product_log_FC', y = 'product_adj_pvals', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, labSize = 2.0,
                pointSize = 0.3, FCcutoff = 0.6, colAlpha = 0.8) + xlim(-1.0,2.7) + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 7), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
dev.off()

DMB_markers <- FindMarkers(monocytes, group.by = "Sample", ident.1 = "Control_DMB", ident.2 = "Control_Bone", logfc.threshold = 0.0)
write.csv(DMB_markers, file = "./csvs/ImmuneNeg.markers.monocytes_DMBvsBone_control.csv", row.names = T)

pdf(file = "./plots/immune_neg_volcano_monocytes_DMB_bone_control.pdf", width = 3.0, height = 1.6)
EnhancedVolcano(DMB_markers,
                lab = rownames(DMB_markers),
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, labSize = 2.0,
                pointSize = 0.3, FCcutoff = 0.6, colAlpha = 0.8) + xlim(-1.7,2.0) + ylim(1,132) + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 7), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
dev.off()


pdf(file = "./plots/immune_neg_tams_featureplots.pdf", width = 4.0, height = 1.0)
FeaturePlot(monocytes, rownames(tam_markers)[1:5], cols = c("grey", "firebrick"), ncol = 5, pt.size = 0.001) & theme_void(base_size = 6) & theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) 
dev.off()




pdf(file = "./plots/immune_neg_dmb_vlnplots.pdf", width = 2.2, height = 1.3)
FeaturePlot(monocytes, features = c("Cd72", "Isg15", "Irf7", "Ly6a", "Aif1", "Fcgr1"), ncol = 3, pt.size = 3.0, raster = T) & theme_void(base_size = 6) & theme(legend.key.size = unit(0.1,"cm")) & scale_color_viridis_c(direction = -1)
dev.off()

Hi
pdf(file = "./plots/immune_neg_tam_featureplots.pdf", width = 1.3, height = 1.3)
VlnPlot(monocytes, features = c("Ccl2", "Ccl12", "Fcrls", "Basp1"), ncol = 2, group.by = "Sample", pt.size = 0) & theme_void(base_size = 6) & theme(legend.key.size = unit(0.2,"cm"), legend.position = "none") & stat_summary(fun = "mean",
               geom = "point", size = 0.1, color = "black") & scale_fill_brewer(palette = "Dark2")
dev.off()


load(file = "./../../ess247/badframe.RData")
temp  <- as.data.frame.numeric(myframe[,1:22])
lapply(myframe, FUN = function(x) {return(as.numeric(x))})

              
ggplot(myframe) + geom_boxplot(mapping = aes(x = Sample, y = m_Tinagl1))

