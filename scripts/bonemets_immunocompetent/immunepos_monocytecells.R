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
monocytes <- subset(seurat.object, subset = RNA_Celltype == "Monocytes")
monocytes <- NormalizeData(monocytes) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
monocytes <- CellCycleScoring(monocytes, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
ElbowPlot(object = monocytes, ndims = 50)
DimPlot(monocytes, reduction = "pca", group.by = "Phase")
DimPlot(monocytes, reduction = "pca", group.by = "Sample")
monocytes <- FindNeighbors(monocytes, dims = 1:20, force.recalc = TRUE)
monocytes <- FindClusters(monocytes, resolution = 0.5)
monocytes <- FindClusters(monocytes, resolution = 0.3)
monocytes <- FindClusters(monocytes, resolution = 0.2)

table(Idents(monocytes))
monocytes <- RunTSNE(object = monocytes, reduction = "pca", dims = 1:20)
DimPlot(monocytes, reduction = "tsne", group.by = "Sample") | DimPlot(monocytes, reduction = "tsne", group.by = "RNA_snn_res.0.2") | DimPlot(monocytes, reduction = "tsne", group.by = "Phase")

monocytes$RNA_snn_res.0.2 <- factor(monocytes$RNA_snn_res.0.2, levels = c("0", "3", "2", "1"))

Idents(monocytes) <- monocytes$RNA_snn_res.0.2
monocytes <- RenameIdents(object = monocytes, `0` = "Apoe+ macrophages", `3` = "Arg1+ macrophages", 
                           `2` = "Nos2+ Macrophages", `1` = "Dendritic cells")
monocytes$subtypes <- Idents(monocytes)
table(monocytes$subtypes)


library(pals)
pdf(file = "./plots/immune_pos_monocyte_umap.pdf", width = 2.5, height = 1.6)
DimPlot(monocytes, reduction = "tsne", group.by = "subtypes", cols = as.vector(alphabet()[12:16]), pt.size = 0.001) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Monocyte subtypes")
dev.off()

library(pals)
pdf(file = "./plots/immune_pos_monocyte_sample_umap.pdf", width = 2.3, height = 1.5)
DimPlot(monocytes, reduction = "tsne", group.by = "Sample", cols = "Set2", pt.size = 0.001) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Condition")
dev.off()


(DimPlot(monocytes, reduction = "tsne", group.by = "Sample", cols = "Dark2", pt.size = 0.01) | FeaturePlot(monocytes, c("C1qa", "Arg1"), cols = c("grey", "firebrick"), pt.size = 0.01)) & theme_void(base_size = 6) & theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) 

FeaturePlot(monocytes, c("Nos2", "Arg1"), cols = c("grey", "firebrick"), pt.size = 0.01, split.by = "Sample") & theme_void(base_size = 6)

VlnPlot(monocytes, features = c("Arg1", "Nos2"), ncol = 5, group.by = "Sample", pt.size = 0.001) & theme_void(base_size = 6) & theme(legend.key.size = unit(0.2,"cm"), legend.position = "none", plot.title = element_text(hjust = 0.5))  & scale_fill_brewer(palette = "Dark2") 

FeatureScatter(monocytes, feature1 = "Il1b", feature2 = "Nos2", group.by = "Sample")
monocytes

FeaturePlot(monocytes, c("Nos2", "Il10", "Tnf", "Ifng"), cols = c("grey", "firebrick"), pt.size = 0.01, split.by = "Sample")

# save(monocytes, file="Robjs_ImmunePos/monocytes.Robj")
load("Robjs_ImmunePos/monocytes.Robj")

Idents(monocytes) <- monocytes$RNA_snn_res.0.2
monocyte.markers <- FindAllMarkers(monocytes, logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T)
write.csv(monocyte.markers, "./csvs/ImmunePos.markers.monocyte_subclusters.csv", row.names = T)
monocyte.markers = read.csv("./csvs/ImmunePos.markers.monocyte_subclusters.csv", row.names = 1)
markers.top10.monocytes = monocyte.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
markers.top5.monocytes = monocyte.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

pdf(file = "./plots/immune_pos_monocyte_marker_heatmap.pdf", width = 2.5, height = 3.5)
Idents(monocytes) <- monocytes$RNA_snn_res.0.2
DoHeatmap(monocytes, features = unique(markers.top10.monocytes$gene), group.colors = as.vector(alphabet()[12:16]), draw.lines = F, size = 0) + scale_color_manual(
  values = as.vector(alphabet()[12:16])) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6, color = "black"), axis.text.y = element_text(size = 6, color = "black"), legend.text = element_text(size = 6, colour = "black"), legend.position = "bottom", legend.margin = margin(-20,0,0,0), legend.box.margin = margin(0,0,0,0), legend.key.size = unit(6,"pt"))
dev.off()

summary <- as.data.frame.matrix(table(monocytes@meta.data[c("Sample", "subtypes")]))
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
summary_prop$sample <- factor(summary_prop$sample, levels = rev(c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB")))
# summary_prop$variable <- factor(summary_prop$variable, levels = c("0", "1", "3", "2"))
pdf(file = "./plots/immune_pos_monocyte_composition.pdf", width = 1.9, height = 2.0)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity", position = "dodge") + 
  labs(x = "Sample", y = "Composition", fill = "subtypes") + theme_classic() + theme(plot.margin = margin(1,1,1,1), text = element_text(size = 6), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.text.x = element_text(size = 6, angle = 0), axis.line = element_line(size = 0.4), legend.position = "none", legend.text = element_text(size = 6), legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin = margin(-0,0,0,0), legend.box.margin = margin(0,0,0,0)) + scale_fill_manual("Fibroblast clusters", values = as.vector(alphabet()[12:16])) + coord_flip()
dev.off()

tam_markers <- FindMarkers(subset(monocytes, RNA_snn_res.0.2 %in% c("0", "2", "3")), assay = "RNA", group.by = "RNA_snn_res.0.2", logfc.threshold = 0.0, ident.1 = "2", only.pos = T)
write.csv(tam_markers, file = "./csvs/ImmunePos.markers.TAMs.csv", row.names = T)

tam_markers = read.csv("./csvs/ImmunePos.markers.TAMs.csv", row.names = 1)
pdf(file = "./plots/immune_pos_volcano_tams.pdf", width = 3.5, height = 1.8)
library(EnhancedVolcano)
EnhancedVolcano(tam_markers,
                lab = rownames(tam_markers), 
                # selectLab = tam_pos[!tam_pos %in% DMB_pos],
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2.0, 
                pointSize = 0.3, FCcutoff = 1.2, colAlpha = 0.8) + xlim(0.15,3.65)  + theme_classic() + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6), axis.title = element_text(size = 6, colour = "black"), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),   legend.box.margin=margin(-10,-5,-10,-10))
dev.off()


# monocytes$Condition = str_split_fixed(monocytes$Sample, pattern = "_", n = 2)[,2]
# monocytes$Tumor = str_split_fixed(monocytes$Sample, pattern = "_", n = 2)[,1]
# Idents(monocytes) <- monocytes$Condition
# DMB_conserved_markers <- FindConservedMarkers(monocytes, ident.1 = "DMB", ident.2 = "Bone", grouping.var = "Tumor", logfc.threshold = 0.0)
# DMB_conserved_markers$product_log_FC <- DMB_conserved_markers$Tumor_avg_log2FC * DMB_conserved_markers$Control_avg_log2FC
# DMB_conserved_markers$product_adj_pvals <- DMB_conserved_markers$Tumor_p_val_adj * DMB_conserved_markers$Control_p_val_adj
# write.csv(DMB_conserved_markers, file = "./csvs/ImmunePos.markers.monocytes_DMBvsBone_conserved.csv", row.names = T)
# DMB_conserved_markers = DMB_conserved_markers[DMB_conserved_markers$product_log_FC > 0.0,]
# DMB_conserved_markers = DMB_conserved_markers[DMB_conserved_markers$Control_avg_log2FC > 0.0,]
# 
# pdf(file = "./plots/immune_pos_volcano_monocytes_DMB_bone.pdf", width = 3.0, height = 1.6)
# EnhancedVolcano(DMB_conserved_markers,
#                 lab = rownames(DMB_conserved_markers),
#                 x = 'product_log_FC', y = 'product_adj_pvals', title = NULL, subtitle = NULL, caption = NULL,
#                 legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, labSize = 2,
#                 pointSize = 0.5, FCcutoff = 0.5, colAlpha = 0.8) + xlim(0, 1.4)  + theme_classic() + theme(text = element_text(size = 6), axis.title = element_text(size = 7), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
# dev.off()


DMB_markers <- FindMarkers(monocytes, group.by = "Sample", ident.1 = "Control_DMB", ident.2 = "Control_Bone", logfc.threshold = 0.0)
write.csv(DMB_markers, file = "./csvs/ImmunePos.markers.monocytes_DMBvsBone_control.csv", row.names = T)

DMB_markers = read.csv("./csvs/ImmunePos.markers.monocytes_DMBvsBone_control.csv", row.names = 1)
pdf(file = "./plots/immune_pos_volcano_monocytes_DMB_bone_control.pdf", width = 4.0, height = 2.0)
EnhancedVolcano(DMB_markers,
                lab = rownames(DMB_markers),
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2,
                pointSize = 0.5, FCcutoff = 0.5, colAlpha = 0.8) + xlim(-1.85, 1.85) + ylim(1.0,100) + theme_classic() + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 5), axis.title = element_text(size = 6), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10))
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

pdf("./plots/ImmunePos_monocytes_DMB_go.pdf", width = 4.0, height = 2.0)
barplot(egoRes, showCategory = 15, x = "minuslogpval", label_format = 100) & theme_classic(base_size = 6) & theme(plot.margin = margin(0,0,0,0), axis.text = element_text(color= "black", size = 6), plot.title = element_text(size = 6), legend.position = "none", legend.text = element_text(margin = margin(0,0,0,0)), legend.spacing = unit(0, "pt")) & scale_fill_gradient(low = "#5b5b5b", high = "#5b5b5b") & xlab("-log adjusted p-value") & ylab("Term") & ggtitle("Gene ontology terms enriched for genes upregulated on DMB")
dev.off()


pdf(file = "./plots/immune_pos_dmb_vlnplots1.pdf", width = 4.2, height = 1.5) 
temp = c("Arg1", "Cd163", "Mmp9", "Fn1")
monocytes$Sample = factor(monocytes$Sample, levels = c("Control_Bone", "Control_DMB", "Tumor_Bone", "Tumor_DMB"))
VlnPlot(monocytes, features = temp, ncol = 4, group.by = "Sample", pt.size = 0, cols = brewer.set2(4)) & theme_classic(base_size = 6) & theme(legend.key.size = unit(0.2,"cm"), legend.position = "none", plot.title = element_text(size = 6, hjust = 0.5), axis.text.y = element_text(colour = "black", size = 5), axis.text.x = element_text(size = 6, colour = "black", angle = 45, hjust = 1.0), plot.margin = margin(0,0,0,0)) & stat_summary(fun = "mean", geom = "point", size = 0.2, color = "black") & xlab("") & ylab("") 
dev.off() 

tam_genes = markers.top10.monocytes$gene[markers.top10.monocytes$cluster == "2"]
tam_genes
tam_genes_temp = c("Nos2", "Arg1", "Ccl5", "Il1b", "Mmp12",  "Mmp13", "Cxcl3", "Cxcl2")


pdf(file = "./plots/immune_pos_tams_featureplots.pdf", width = 3.0, height = 1.6)
FeaturePlot(monocytes, tam_genes_temp, cols = c("grey", "black"), ncol = 4, pt.size = 0.0000001, max.cutoff = 'q99') & theme_void(base_size = 6) & theme(plot.margin = margin(0,0,0,0), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 6)) 
dev.off()

# DMB_conserved_markers = DMB_conserved_markers[order(DMB_conserved_markers$product_adj_pvals),]
# rownames(DMB_conserved_markers)[1:10]
# pdf(file = "./plots/immune_pos_DMB_featureplots.pdf", width = 4.0, height = 1.8)
# FeaturePlot(monocytes, rownames(DMB_conserved_markers)[1:10], cols = c("grey", "firebrick"), ncol = 5, pt.size = 0.001) & theme_void(base_size = 10) & theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) 
# dev.off()

# 
# pdf(file = "./plots/immune_pos_tam_vlmplots.pdf", width = 3.0, height = 1.5)
# VlnPlot(monocytes, features = rownames(DMB_conserved_markers)[1:10], ncol = 5, group.by = "Condition", pt.size = 0) & theme_void(base_size = 6) & theme(legend.key.size = unit(0.2,"cm"), legend.position = "none", plot.title = element_text(hjust = 0.5))  & scale_fill_brewer(palette = "Dark2") 
# dev.off()




gene_list1 = read.csv("csvs/ImmuneNeg.tam_markers.csv", row.names = 1)
EnhancedVolcano(gene_list1,
                lab = rownames(gene_list1), selectLab = temp,
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, labSize = 4,
                pointSize = 2.0, FCcutoff = 1.0, colAlpha = 0.8, drawConnectors = T) + xlim(0.9,2.2)

gene_list2 = read.csv("csvs/ImmunePos.tam_markers.csv", row.names = 1)

gene_list1_sig <- gene_list1[(gene_list1$avg_log2FC > 1.0 & gene_list1$p_val_adj<0.01),]
gene_list2_sig <- gene_list2[(gene_list2$avg_log2FC > 1.0 & gene_list2$p_val_adj<0.01),]

temp = intersect(rownames(gene_list1_sig), rownames(gene_list2_sig))


FeaturePlot(tams, c("Dcn", "Postn", "Cd8a", "S100a8"))

FeaturePlot(monocytes, c("Cd14", "Cd68", "Cd33", "Itgam"), ncol = 4)
#TAM
FeaturePlot(monocytes, c("Nos2", "Cxcl3", "Tgm2", "Srgn", "Cd274", "Mmp13"), ncol = 3)
#M2
FeaturePlot(monocytes, c("Cd163", "Mrc1"), ncol = 1)
#M1
FeaturePlot(monocytes, c("Cd80", "Cd86", "Fcgr1", "Fcgr3", "Cd38", "Fpr2"), ncol = 3)


tam_pos <- subset(monocytes, RNA_snn_res.0.2 == "1")
dim(tam_pos)

save(tam_pos, file = "Robjs_ImmunePos/tam_pos.Robj")

load("Robjs_ImmunePos/monocytes.Robj")
monocytes_pos <- monocytes
load("./Robjs_ImmuneNeg/monocytes.Robj")
monocytes_neg <- monocytes

DimPlot(monocytes_pos, reduction = "tsne", group.by = "RNA_snn_res.0.2")
monocytes_pos$TAM <- "Non-TAM"
monocytes_pos$TAM[monocytes_pos$RNA_snn_res.0.2 == "1"] <- "TAM"
DimPlot(monocytes_neg, reduction = "tsne", group.by = "RNA_snn_res.0.2") | DimPlot(monocytes_neg, reduction = "tsne", group.by = "Sample")
monocytes_neg$TAM <- "Non-TAM"
monocytes_neg$TAM[monocytes_neg$RNA_snn_res.0.2 == "2"] <- "TAM"

monocytes_merged = merge(monocytes_neg, c(monocytes_pos))
dim(monocytes_merged)

# Normalize the data, then center and scale.
monocytes_merged <- NormalizeData(object = monocytes_merged)
monocytes_merged <- CellCycleScoring(monocytes_merged, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)
monocytes_merged <- FindVariableFeatures(object = monocytes_merged, selection.method = "vst", nfeatures = 2000)

# Run Principal Component Analysis.
monocytes_merged <- ScaleData(object = monocytes_merged)
dim(GetAssayData(monocytes_merged, slot = "scale.data"))
monocytes_merged <- RunPCA(object = monocytes_merged, npcs = 50)

monocytes_merged <- RunHarmony(monocytes_merged, group.by.vars = "orig.ident")

ElbowPlot(object = monocytes_merged, reduction = "harmony", ndims = 50)
n.pcs = 50
monocytes_merged <- FindNeighbors(monocytes_merged, dims = 1:n.pcs, reduction = "harmony", force.recalc = TRUE)
monocytes_merged <- FindClusters(monocytes_merged, resolution = 0.20)
monocytes_merged <- FindClusters(monocytes_merged, resolution = 0.30)

table(Idents(monocytes_merged))
# To visualize

monocytes_merged <- RunTSNE(object = monocytes_merged, reduction = "harmony", dims = 1:n.pcs)

((DimPlot(monocytes_merged, reduction = "tsne", group.by = "orig.ident") |
DimPlot(monocytes_merged, reduction = "tsne", group.by = "RNA_snn_res.0.3")) / 
(DimPlot(monocytes_merged, reduction = "tsne", group.by = "Sample") |
   DimPlot(monocytes_merged, reduction = "tsne", group.by = "TAM")) ) | FeaturePlot(monocytes_merged, features = c("Ccr2"))

Idents(monocytes_merged) <- monocytes_merged$TAM
library(metap)
TAM_conserved_markers <- FindConservedMarkers(monocytes_merged, ident.1 = "TAM", grouping.var = "orig.ident",
                                              logfc.threshold = 0.5, only.pos = TRUE, meta.method = minimump)
write.csv(TAM_conserved_markers, "./csvs/tam_conserved_markers.csv")


library(pals)

monocytes_neg$TAM[monocytes_neg$TAM == "Non-TAM"] <- "Non\nTAM"

pdf(file = "./plots/vln_tam_conserved_neg.pdf", width = 5.0, height = 1.2)
Idents(monocytes_neg) <- monocytes_neg$TAM
VlnPlot(monocytes_neg, features = rownames(TAM_conserved_markers[1:5,]), pt.size = 0, ncol = 5) & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black") &  theme_classic() & theme(text = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-5,-10,-10), plot.title = element_text(hjust = 0.5, size = 8)) & scale_fill_brewer(palette = "Dark2") & xlab(NULL) & ylab(NULL)
dev.off()

levels(monocytes_neg)


monocytes_pos$TAM <- factor(monocytes_pos$TAM, levels = c("TAM", "Non-TAM"))
pdf(file = "./plots/vln_tam_conserved_pos.pdf", width = 5.0, height = 1.2)
Idents(monocytes_pos) <- monocytes_pos$TAM
VlnPlot(monocytes_pos, features = rownames(TAM_conserved_markers[1:5,]), pt.size = 0, ncol = 5) & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black") &  theme_classic() & theme(text = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-5,-10,-10), plot.title = element_text(hjust = 0.5, size = 8)) & scale_fill_brewer(palette = "Dark2") & xlab(NULL) & ylab(NULL)
dev.off()


FeaturePlot(monocytes_neg, reduction = "tsne", features = "Cd274")


