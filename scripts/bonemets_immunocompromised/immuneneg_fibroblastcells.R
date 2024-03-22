library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(here); library(fitdistrplus); library(ggplot2)
library(EnhancedVolcano)
library(pals); library(RColorBrewer)
library(harmony)
library(clusterProfiler)
library(AnnotationHub)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

setwd("/fs/cbsuvlaminck3/workdir/mm2937/FishbachLab/")
load("Robjs_ImmuneNeg/seurat.object.immuneNeg.exp1.Robj")
load("Robjs_ImmuneNeg/gene_names.Robj")

dim(seurat.object)


# Mouse analysis starts here
counts = GetAssayData(seurat.object, assay = "RNA", slot = "counts")[mouse_genes, seurat.object$organism == "MOUSE"]
rownames(counts) <- str_split_fixed(rownames(counts), pattern = "---", n = 2)[,2]
meta.data = seurat.object@meta.data[seurat.object$organism == "MOUSE", ]
mouse.object <- CreateSeuratObject(counts = counts, meta.data = meta.data)
table(mouse.object$Celltype_groups)
rownames(mouse.object)


cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")


# Getting fibroblast subset
fibroblast <- subset(mouse.object, subset = Celltype_groups %in% c("Fibroblasts"))
fibroblast <- subset(fibroblast, subset = Sample %in% c("Control_Bone", "Control_DMB"))

fibroblast <- NormalizeData(fibroblast) 
fibroblast <- CellCycleScoring(fibroblast, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
fibroblast <- FindVariableFeatures(fibroblast)
fibroblast <- ScaleData(fibroblast, vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"))
fibroblast <- RunPCA(fibroblast)
ElbowPlot(fibroblast, ndims = 50)
n.pcs = 20
fibroblast <- FindNeighbors(fibroblast, reduction = "pca", dims = 1:10, force.recalc = TRUE)
fibroblast <- FindClusters(fibroblast, resolution = 0.5)
fibroblast <- FindClusters(fibroblast, resolution = 0.3)
fibroblast <- FindClusters(fibroblast, resolution = 0.2)
table(Idents(fibroblast))
fibroblast <- RunTSNE(object = fibroblast, reduction = "pca", dims = 1:10, per)

DimPlot(fibroblast, reduction = "tsne", group.by = "Sample") | DimPlot(fibroblast, reduction = "tsne", group.by = "RNA_snn_res.0.2") | DimPlot(fibroblast, reduction = "tsne", group.by = "Phase")

DimPlot(fibroblast, reduction = "pca", group.by = "Sample") | DimPlot(fibroblast, reduction = "pca", group.by = "RNA_snn_res.0.2") | DimPlot(fibroblast, reduction = "pca", group.by = "Phase")

DimPlot(fibroblast, reduction = "tsne", group.by = "RNA_snn_res.0.3")
FeaturePlot(fibroblast, reduction = "tsne", features = c("Thbs1", "Mgp", "Saa3", "Acta2", "C3", "S100a4", "Cd55"))

Idents(fibroblast) <- fibroblast$RNA_snn_res.0.2
fibroblast.markers <- FindAllMarkers(fibroblast, logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T)
write.csv(fibroblast.markers, "csvs/ImmuneNeg.markers.fibroblast_subclusters.csv", row.names = T)
markers.top10.fibroblast = fibroblast.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.fibroblast = fibroblast.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)



DotPlot(fibroblast, features = unique(markers.top10.fibroblast$gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0)) + scale_color_viridis_c(direction = -1)

VlnPlot(fibroblast, c("percent.mito"))

fibroblast <- subset(fibroblast, subset = RNA_snn_res.0.2 %in% c("0", "1", "2"))

pdf(file = "./plots/immune_neg_fibroblast_controls_umap.pdf", width = 1.6, height = 1.6)
DimPlot(fibroblast, reduction = "tsne", group.by = "RNA_snn_res.0.2", pt.size = 0.0001) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 8), 
        legend.position = 'none', 
        legend.key.size = unit(6,"pt")) + ggtitle(NULL)
dev.off()

library("pals")
library(scales)
hex = hue_pal()(4)
pdf(file = "./plots/immune_neg_fibroblast_controls_umap_sample.pdf", width = 1.6, height = 1.6)
DimPlot(fibroblast, reduction = "tsne", group.by = "Sample", cols = hex, pt.size = 0.0001) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 8), 
        legend.position = 'none', 
        legend.key.size = unit(6,"pt")) + ggtitle(NULL)
dev.off()

unique(fibroblast$Sample)
fibroblast$Sample <- factor(fibroblast$Sample, levels = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB"))
levels(fibroblast$Sample)


pdf(file = "./plots/immune_neg_fibroblast_controls_umap_features.pdf", width = 3.5, height = 1.25)
FeaturePlot(fibroblast, features = c("Thbs1", "Mgp", "Penk", "Col8a1", "C3", "C4b", "S100a4", "Aqp1"), ncol = 4, pt.size = 0.00001) & theme_void(base_size = 6) & theme(text = element_text(size = 6), plot.title = element_text(color = "black", hjust = 0.5), legend.position = "right", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0)) 
dev.off()



# save(fibroblast, file="./Robjs_ImmuneNeg/fibroblast_controls.Robj")
load("./Robjs_ImmuneNeg/fibroblast_controls.Robj")


# Getting fibroblast subset
fibroblast <- subset(mouse.object, subset = Celltype_groups %in% c("Fibroblasts"))
fibroblast <- NormalizeData(fibroblast) %>% FindVariableFeatures()
fibroblast <- CellCycleScoring(fibroblast, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
fibroblast <- ScaleData(fibroblast, vars.to.regress = c("S.Score", "G2M.Score"))
fibroblast <- RunPCA(fibroblast)
ElbowPlot(fibroblast)
n.pcs = 10
fibroblast <- FindNeighbors(fibroblast, reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
fibroblast <- FindClusters(fibroblast, resolution = 0.5)
fibroblast <- FindClusters(fibroblast, resolution = 0.3)
fibroblast <- FindClusters(fibroblast, resolution = 0.2)
table(Idents(fibroblast))
fibroblast <- RunTSNE(object = fibroblast, reduction = "pca", dims = 1:n.pcs)
DimPlot(fibroblast, reduction = "tsne", group.by = "Sample") | DimPlot(fibroblast, reduction = "tsne", group.by = "RNA_snn_res.0.2") | DimPlot(fibroblast, reduction = "tsne", group.by = "Phase")

DimPlot(fibroblast, reduction = "tsne", group.by = "RNA_snn_res.0.3")

FeaturePlot(fibroblast, features = c("Pdpn", "S100a4",  "Dcn", "Gsn", "Acta2", "Cxcl12", "Cxcl1", "Il6"))
FeaturePlot(fibroblast, features = c("Ppap2b", "Ogn",  "Timp2", "Gsn", "Acta2", "Cxcl12", "Cxcl1", "Il6", "Timp1", "Tpm1"))
FeaturePlot(fibroblast, features = c("Pdpn", "S100a4", "H2-Aa"))
FeaturePlot(fibroblast, features = c("Thbs1", "Mgp", "Penk"))

fibroblast <- subset(fibroblast, subset = RNA_snn_res.0.2 %in% c("0", "1", "2", "3"))

fibroblast <- FindNeighbors(fibroblast, reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
fibroblast <- FindClusters(fibroblast, resolution = 0.3)


# save(fibroblast, file="./Robjs_ImmuneNeg/fibroblast.Robj")
load("./Robjs_ImmuneNeg/fibroblast.Robj")

Idents(fibroblast) <- fibroblast$RNA_snn_res.0.2
fibroblast <- RenameIdents(object = fibroblast, `0` = "S100a4+ fibroblasts", `1` = "Saa3+ fibroblasts", 
                              `2` = "Col8a1+ fibroblasts", `3` = "Acta2+ fibroblasts")
fibroblast$subtypes <- Idents(fibroblast)
table(fibroblast$subtypes)


pdf(file = "./plots/immune_neg_fibroblast_umap.pdf", width = 2.5, height = 1.6)
library(pals)
DimPlot(fibroblast, reduction = "tsne", group.by = "subtypes", pt.size = 0.00001, cols = c(as.vector(alphabet()[17:18]), as.vector(alphabet()[22:23]))) + ylim(-40,35) + labs(color = "Fibroblast subtypes") +
  theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.key.size = unit(6,"pt"), legend.text = element_text(size = 6)) + ggtitle(NULL)
dev.off()

pdf(file = "./plots/immune_neg_fibroblast_sample_umap.pdf", width = 2.7, height = 2.0)
DimPlot(fibroblast, reduction = "tsne", group.by = "Sample", pt.size = 0.0001, cols = "Set2") + theme_void(base_size = 6) + labs(color = "Condition") + ylim(-40,35) + 
  theme(text = element_text(size = 6), 
        legend.position = 'right', 
        legend.key.size = unit(6,"pt"), legend.text = element_text(size = 6),legend.margin = margin(0,0,0,0)) + ggtitle(NULL)
dev.off()

fibroblast
library(ggbreak)
Idents(fibroblast) <- fibroblast$subtypes
fibroblast.markers <- FindAllMarkers(fibroblast, logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T)
write.csv(fibroblast.markers, "csvs/ImmuneNeg.markers.fibroblast_subclusters.csv", row.names = T)

markers.top10.fibroblast = fibroblast.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.fibroblast = fibroblast.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

pdf(file = "./plots/fibroblast_control_markers_dotplot.pdf", width = 2.5, height = 3.5)
DoHeatmap(fibroblast, features = unique(markers.top10.fibroblast$gene), group.colors = c(as.vector(alphabet()[17:18]), as.vector(alphabet()[22:23])), draw.lines = F, size = 0) + theme(plot.margin = margin(0,0,0,0), axis.text.y = element_text(size = 6, colour = "black"), legend.text = element_text(size = 6), legend.position = "bottom", legend.margin = margin(-20,0,0,0), legend.box.margin = margin(0,0,0,0), legend.key.size = unit(6,"pt"), text = element_text(size = 6, colour = "black")) + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

FeaturePlot(fibroblast, features = c("S100a4", "C4b"))

unique(fibroblast$Sample)
fibroblast$Sample <- factor(fibroblast$Sample, levels = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB"))


dim(fibroblast)
summary <- as.data.frame.matrix(table(fibroblast@meta.data[c("Sample", "subtypes")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$sample)
summary_prop$sample <- factor(summary_prop$sample, levels = rev(c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB")))
pdf(file = "./plots/immune_neg_fibroblast_composition.pdf", width = 1.9, height = 1.6)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity", position = "dodge") + labs(x = "Condition", y = "Composition", fill = "subtypes") + theme_classic() + theme(plot.margin = margin(1,1,1,1), text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.text.x = element_text(size = 6, angle = 0), axis.line = element_line(size = 0.4), legend.position = "none", legend.text = element_text(size = 6), legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0)) + scale_fill_manual("Fibroblast clusters", values = c(as.vector(alphabet()[17:18]), as.vector(alphabet()[22:23]))) + coord_flip()
dev.off()


library(EnhancedVolcano)
CAF_markers = FindMarkers(fibroblast, group.by = "RNA_snn_res.0.2",
                          ident.1 = c("3"), logfc.threshold = 0.0, only.pos = F, verbose = T)
write.csv(CAF_markers, "csvs/ImmuneNeg.markers.CAFs.csv", row.names = T)
CAF_markers <- read.csv("csvs/ImmuneNeg.markers.CAFs.csv", row.names = 1)
pdf(file = "./plots/immune_neg_volcano_cafs.pdf", width = 3.0, height = 1.6)
EnhancedVolcano(CAF_markers,
                lab = rownames(CAF_markers), 
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, labSize = 2,
                pointSize = 0.3, FCcutoff = 1.0, colAlpha = 0.8) + xlim(0.15,3.9) + theme_classic(base_size = 6) + theme(text = element_text(size = 6, colour = "black"), plot.margin = margin(0,0,0,0), axis.title = element_text(size = 6), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-5,-10,-10))
dev.off()

# fibroblast$Condition = str_split_fixed(fibroblast$Sample, pattern = "_", n = 2)[,2]
# fibroblast$Tumor = str_split_fixed(fibroblast$Sample, pattern = "_", n = 2)[,1]
# Idents(fibroblast) <- fibroblast$Condition
# DMB_conserved_markers <- FindConservedMarkers(fibroblast, ident.1 = "DMB", ident.2 = "Bone", grouping.var = "Tumor", logfc.threshold = 0.0)
# DMB_conserved_markers$product_log_FC <- DMB_conserved_markers$Tumor_avg_log2FC * DMB_conserved_markers$Control_avg_log2FC
# DMB_conserved_markers$product_adj_pvals <- DMB_conserved_markers$Tumor_p_val_adj * DMB_conserved_markers$Control_p_val_adj
# write.csv(DMB_conserved_markers, file = "./csvs/ImmuneNeg.markers.fibroblast_DMBvsBone_conserved.csv", row.names = T)

# pdf(file = "./plots/immune_neg_volcano_fibroblast_DMB_bone.pdf", width = 3.0, height = 1.6)
# EnhancedVolcano(DMB_conserved_markers,
#                            lab = rownames(DMB_conserved_markers),
#                            x = 'product_log_FC', y = 'product_adj_pvals', title = NULL, subtitle = NULL, caption = NULL,
#                 legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2,
#                 pointSize = 0.5, FCcutoff = 0.6, colAlpha = 0.8) + xlim(-3.2,2.3) + ylim(0, 350) + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
# dev.off()

table(fibroblast$Sample)
DMB_markers <- FindMarkers(fibroblast, group.by = "Sample", ident.1 = "Control_DMB", ident.2 = "Control_Bone", logfc.threshold = 0.0)
write.csv(DMB_markers, file = "./csvs/ImmuneNeg.markers.fibroblast_DMBvsBone_control.csv", row.names = T)

DMB_markers = read.csv("./csvs/ImmuneNeg.markers.fibroblast_DMBvsBone_control.csv", row.names = 1)
pdf(file = "./plots/immune_neg_volcano_fibroblast_DMB_bone_control.pdf", width = 4.0, height = 2.0)
EnhancedVolcano(DMB_markers,
                lab = rownames(DMB_markers),
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2,
                pointSize = 0.5, FCcutoff = 0.6, colAlpha = 0.8) + xlim(-4.0,2.5) + ylim(1, 300) + theme_classic() + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
dev.off()

marker_data = DMB_markers[DMB_markers$avg_log2FC > 1.0 & DMB_markers$p_val_adj < 0.01,]
dim(marker_data)

organism_dbs <- c("org.Mm.eg.db", "org.Hs.eg.db")
orgamisms = c("Mus musculus", "Homo sapiens")

gene.df <- bitr(rownames(marker_data), fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = organism_dbs[1])

egoRes <- enrichGO(gene        =  gene.df$ENTREZID,
                   OrgDb         = organism_dbs[1],
                   keyType       = 'ENTREZID',
                   ont           = c("BP"),
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.01, 
                   readable = T )

egoRes <- clusterProfiler::simplify(egoRes)
summaryEgoRes <- as.data.frame(egoRes)

dotplot(egoRes, showCategory = 20)

egoRes@result$minuslogpval  = -1 * log10(egoRes@result$p.adjust)
egoRes@result$minuslogpval = as.numeric(egoRes@result$minuslogpval)

pdf("./plots/ImmuneNeg_fibroblasts_DMB_go.pdf", width = 3.2, height = 2.0)
barplot(egoRes, showCategory = 15, x = "minuslogpval", label_format = 100) & theme_classic(base_size = 6) & theme(plot.margin = margin(0,0,0,0), axis.text = element_text(color= "black", size = 6), plot.title = element_text(size = 6), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt")) & scale_fill_gradient(low = "#5b5b5b", high = "#5b5b5b") & xlab("-log adjusted p-value") & ylab("Term") & ggtitle("Gene ontology terms enriched for genes upregulated on DMB")
dev.off()

#
# table(subset(fibroblast, RNA_snn_res.0.2 %in% c("3"))$Sample)
# CAF_DMB_markers <- FindMarkers(subset(fibroblast, RNA_snn_res.0.2 %in% c("3")), group.by = "Sample", ident.1 = "Tumor_DMB", ident.2 = "Tumor_Bone", logfc.threshold = 0.0)
# CAF_DMB_markers

# EnhancedVolcano(CAF_DMB_markers,
#                 lab = rownames(CAF_DMB_markers),
#                 x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
#                 legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2,
#                 pointSize = 0.5, FCcutoff = 0.2, colAlpha = 0.8)  + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10)) #+ xlim(-4.0,2.6) + ylim(1, 300)
# 



pdf(file = "./plots/immune_neg_fibroblast_dmb_vlnplots.pdf", width = 4.0, height = 1.6) 
temp = c(rownames(DMB_markers)[DMB_markers$avg_log2FC > 0.0][1:5], rownames(DMB_markers)[DMB_markers$avg_log2FC < 0.0][1:5])
VlnPlot(subset(fibroblast, Sample %in% c("Control_Bone", "Control_DMB")), features = temp, ncol = 5, group.by = "Condition", pt.size = 0, cols = brewer.set2(10)) & theme_classic(base_size = 6) & theme(legend.key.size = unit(0.2,"cm"), legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5), axis.text.x = element_text(size = 6, colour = "black"), axis.text = element_text(colour = "black"), plot.margin = margin(0,0,0,0)) & stat_summary(fun = "mean", geom = "point", size = 0.2, color = "black") & xlab("") & ylab("")
dev.off()

pdf(file = "./plots/immune_neg_fibroblast_dmb_vlnplots2.pdf", width = 4.2, height = 1.5) 
temp = c("Col8a1", "Col11a1", "Thbs1", "Mgp")
fibroblast$Sample = factor(fibroblast$Sample, levels = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB"))
VlnPlot(fibroblast, features = temp, ncol = 4, group.by = "Sample", pt.size = 0, cols = brewer.set2(4)) & theme_classic(base_size = 6) & theme(legend.key.size = unit(0.2,"cm"), legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5), axis.text.y = element_text(colour = "black", size = 5), axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1.0), plot.margin = margin(0,0,0,0)) & stat_summary(fun = "mean", geom = "point", size = 0.2, color = "black") & xlab("") & ylab("") 
dev.off() 

pdf(file = "./plots/immune_neg_fibroblast_featureplots.pdf", width = 2.0, height = 2.0)
FeaturePlot(fibroblast, features = c("S100a4", "Saa3", "Col8a1", "Mfap4"), cols = c("grey", "black"), ncol = 2, pt.size = 0.00001, max.cutoff = 'q99') & ylim(-40, 35) & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

pdf(file = "./plots/immune_neg_fibroblast_featureplots2.pdf", width = 3.0, height = 1.1)
FeaturePlot(fibroblast, features = c("Mgp", "Aspn", "Gas6"), cols = c("grey", "black"), ncol = 3, pt.size = 0.00001, max.cutoff = 'q99') & ylim(-40, 35) & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

pdf(file = "./plots/immune_neg_fibroblast_featureplots3.pdf", width = 3.0, height = 1.1)
FeaturePlot(fibroblast, features = c("Thbs1", "Mdk", "Ptn"), cols = c("grey", "black"), ncol = 3, pt.size = 0.00001, max.cutoff = 'q99') & ylim(-40, 35) & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

