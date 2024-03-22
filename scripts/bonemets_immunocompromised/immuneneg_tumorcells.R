# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(here); library(fitdistrplus); library(ggplot2)
library(future)
library(harmony)
library(EnhancedVolcano)
library(pals); library(RColorBrewer)
setwd("/workdir/mm2937/FishbachLab/")
library(NMF)

load(file = "cc.genes.rda") 
cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
samples_round1= c("Tumor_Bone_111425_HNYY2BGXC_1M", "Tumor_DMB_111426_HNYY2BGXC_2C")
samples_round2 = c("Tumor_DMB_116340_H7HWFBGXF_1A", "Tumor_Bone_116341_H7HWFBGXF_2B", "Control_DMB_116342_H7HWFBGXF_3C", "Control_Bone_116343_H7HWFBGXF_4D")
cc.genes.humans <- list("s.genes" = c(paste("GRCh38-", cc.genes.updated.2019$s.genes, sep = "")), "g2m.genes" = c(paste("GRCh38-", cc.genes.updated.2019$g2m.genes, sep = ""))) 
cc.genes.mouse <- list("s.genes" = c(paste("mm10---", cc.genes.mouse$s.genes, sep = "")), "g2m.genes" = c(paste("mm10---", cc.genes.mouse$g2m.genes, sep = ""))) 


load("Robjs_ImmuneNeg/seurat.object.immuneNeg.exp1.Robj")
dim(seurat.object)

load("Robjs_ImmuneNeg/gene_names.Robj")

# Human cancer cell analysis starts here!
counts = GetAssayData(seurat.object, assay = "RNA", slot = "counts")[human_genes, seurat.object$organism == "HUMAN"]
rownames(counts) <- str_split_fixed(rownames(counts), pattern = "-", n = 2)[,2]
meta.data = seurat.object@meta.data[seurat.object$organism == "HUMAN", c("orig.ident", "percent.mito", "Sample")]


human.object <- CreateSeuratObject(counts = counts, meta.data = meta.data)
dim(human.object)

human.object <- NormalizeData(human.object) 
human.object <- CellCycleScoring(human.object, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
human.object <- FindVariableFeatures(human.object)

VariableFeaturePlot(human.object)

# human.object <- ScaleData(human.object)
human.object <- ScaleData(human.object, vars.to.regress = c("S.Score", "G2M.Score", "percent.mito"))
human.object <- RunPCA(human.object)
DimPlot(human.object, reduction = "pca", group.by = "Experiment")
DimPlot(human.object, reduction = "pca", group.by = "orig.ident")
ElbowPlot(object = human.object, ndims = 50)
n.pcs = 20

human.object <- FindNeighbors(human.object, reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
human.object <- FindClusters(human.object, resolution = 0.3)
human.object <- FindClusters(human.object, resolution = 0.2)
human.object <- FindClusters(human.object, resolution = 0.5)

table(Idents(human.object))
human.object <- RunTSNE(object = human.object, reduction = "pca", dims = 1:n.pcs)
DimPlot(human.object, reduction = "tsne", group.by = "Phase")
DimPlot(human.object, reduction = "tsne", group.by = "RNA_snn_res.0.3") | DimPlot(human.object, reduction = "tsne", group.by = "Sample")| DimPlot(human.object, reduction = "tsne", group.by = "Phase")
VlnPlot(human.object, features = c("nFeature_RNA", "nCount_RNA"))


# save(human.object, file="Robjs_ImmuneNeg/tumor.object.ImmuneNeg.exp1.Robj")
load("Robjs_ImmuneNeg/tumor.object.ImmuneNeg.exp1.Robj")
human.object

dim(human.object)
var_data <- human.object@assays$RNA@meta.features
temp = var_data[var_data$vst.variance != 0,]
ggplot(temp, mapping = aes(x = vst.variance.standardized)) + geom_histogram(position = "dodge", bins = 200) + theme_classic()
sd(temp$vst.variance)

var_data <- mouse_tumor@assays$RNA@meta.features
temp = var_data[var_data$vst.variance != 0,]
ggplot(temp, mapping = aes(x = vst.variance.standardized)) + geom_histogram(position = "dodge", bins = 200) + theme_classic()
sd(temp$vst.variance)


load("Robjs_ImmuneNeg/tumor.object.ImmuneNeg.exp1.Robj")



## NMF
library(NMF)
source('./PanCancer-main/seurat_functions_public.R')
# Run NMF
# 
data = as.matrix(GetAssayData(human.object, assay = 'RNA', slot = 'scale.data'))

dim(data)
#if (unique(srt$author) == ''){
data = data[VariableFeatures(human.object),]
#}
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]
print(dim(data))

range = 2:25
res.list = mclapply(range, function(r){
  nmf(data, nrun = 1, rank = r, seed = 'ica', method = 'nsNMF')
}, mc.cores = 20)
names(res.list) = range
# Select rank
gmin = 5
modules.list = lapply(res.list, NMFToModules, gmin = gmin)
print(sapply(modules.list,length))
comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
mi = min(comp)
r = names(which(comp == mi))
r = r[length(r)]
print(r)
res = res.list[[r]]

# Process output
modules = NMFToModules(res, gmin = gmin)
scores = basis(res)
dim(scores)
colnames(scores) = names(modules)
coefs = coefficients(res)
dim(coefs)

rownames(coefs) = names(modules)
# Order modules
h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
o = row_order(h)
scores = scores[, o]
coefs = coefs[o, ]
modules = modules[o]
print(modules)
length(modules)s

human.object = AddMetaData(human.object, t(coefs), col.name = rownames(coefs))
# Cluster NMF output
# h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
hcl = as.hclust(column_dend(h))
sig = cutree(hcl, k = length(modules))
nmf = c(by(t(coefs), INDICES = sig, FUN = function(x){names(modules)[which.max(colMeans(x))]}))[sig]
human.object$nmf = factor(nmf, levels = names(modules))
col_nmf = c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Set1'))[1:nlevels(human.object$nmf)]
names(col_nmf) = levels(human.object$nmf)
## Saving
human.object
save(human.object, file = 'Robjs_ImmuneNeg/human.object_nmf.RData')
save(res.list, file = 'Robjs_ImmuneNeg/res.list_nmf.RData')

load('Robjs_ImmuneNeg/human.object_nmf.RData')
load('Robjs_ImmuneNeg/res.list_nmf.RData')

for (reduction in c('pca','umap')){
  h = DimPlot(human.object, pt.size = 2, reduction = reduction, group.by = 'nmf', cols = col_nmf)
  print(h)
  h = DimPlot(human.object, pt.size = 2, reduction = reduction, group.by = 'Sample', cols = col_cluster)
  print(h)
  h = FeaturePlot(human.object, reduction = reduction, features = names(modules))
  print(h)
}

colnames(human.object@meta.data)


DimPlot(human.object, pt.size = 2, reduction = "tsne", group.by = 'nmf') | DimPlot(human.object, pt.size = 2, reduction = "tsne", group.by = 'Sample')



library(pals)
pdf(file = "./plots/immune_neg_umap_tumor_cells.pdf", width = 1.5, height = 1.4)
DimPlot(human.object, reduction = "tsne", group.by = "RNA_snn_res.0.3", cols = "Dark2", pt.size = 0.5)  + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Clusters")
dev.off()

library(pals)
pdf(file = "./plots/immune_neg_umap_tumor_cells_by_samples.pdf", width = 1.75, height = 1.4)
DimPlot(human.object, reduction = "tsne", group.by = "Sample", cols = "Set2", pt.size = 0.5) + scale_color_manual(values = brewer.set2(6)[3:4]) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 6), plot.margin = margin(0,0,0,0),
        legend.position = 'right', legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL) + labs(color = "Condition")
dev.off()

dim(human.object)

library(ggrepel)
Idents(human.object) <- human.object$RNA_snn_res.0.3
markers.all.RNA = FindAllMarkers(object = human.object, assay = "RNA", return.thresh = 0.01, only.pos = T)
# write.csv(markers.all.RNA, file = "./csvs/ImmuneNeg.markers.tumor_subclusters.csv")
markers.top30.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(30, avg_log2FC)
markers.top5.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(5, avg_log2FC)go_results

library(RColorBrewer)
DoHeatmap(human.object, features = unique(markers.top30.RNA$gene), group.colors = brewer.dark2(3)) + scale_color_manual(
  values = brewer.dark2(3), limits = c('0', '1', '2')) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) 

bone_dmb <- FindMarkers(human.object, assay = "RNA", group.by = "Sample", ident.1 = "Tumor_DMB", ident.2 = "Tumor_Bone", test.use = "wilcox")
pdf(file = "./plots/immune_neg_volcano_tumor_dmb.pdf", width = 3.0, height = 1.8)
library(EnhancedVolcano)
EnhancedVolcano(bone_dmb,
                lab = rownames(bone_dmb),
                x = 'avg_log2FC', y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL,
                legendPosition = "none", pCutoff = 10^-2, axisLabSize = 6, labSize = 2.0, 
                pointSize = 0.3, FCcutoff = 1.2, colAlpha = 0.8) + xlim(-5,5)  + theme_classic() + theme(plot.margin = margin(0,0,0,0), text = element_text(size = 6), axis.title = element_text(size = 6, colour = "black"), axis.text = element_text(size = 5, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm") )
dev.off()

write.csv(bone_dmb, file = "csvs/human_tumorbone_tumor_dmb.csv", row.names = T)

tumor_2 <- FindMarkers(human.object, assay = "RNA", group.by = "RNA_snn_res.0.5", ident.1 = "2")
EnhancedVolcano(tumor_2,
                lab = rownames(tumor_2),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 4, pointSize = 2.0, FCcutoff = 0.5, colAlpha = 0.5) + xlim(-1.5,1.5) + ylim(0,30)
write.csv(tumor_1, file = "csvs/human_tumor_2.csv", row.names = T)

summary <- as.data.frame.matrix(table(human.object@meta.data[c("Sample", "RNA_snn_res.0.3")]))
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)

pdf(file = "./plots/immune_neg_tumor_composition.pdf", width = 1.7, height = 1.5)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity") + 
  labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic() + theme(text = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "right", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10)) + scale_fill_brewer(palette = "Dark2")
dev.off()


tumor_2_1 <- FindMarkers(human.object, assay = "RNA", group.by = "RNA_snn_res.0.5", logfc.threshold = 0.0, ident.1 = "2", ident.2 = "1")
pdf(file = "./plots/immune_neg_volcano_tumor_2_1.pdf", width = 3.5, height = 1.7)
EnhancedVolcano(tumor_2_1,
                lab = rownames(tumor_2_1),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 2, pointSize = 0.5, FCcutoff = 0.5, colAlpha = 0.5) + xlim(-1.6,1.5) + ylim(0,25) + theme_classic() + theme(text = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10)) + scale_fill_brewer(palette = "Dark2")
dev.off()
write.csv(tumor_2_1, file = "csvs/ImmuneNeg.markers.tumor_2vs1", row.names = T)



# Using GO terms to visualize genes of interest

Genes = intersect(str_split(egoResSummary["GO:0052547","geneID"], "/")[[1]], str_split(egoResSummary["GO:0061041","geneID"], "/")[[1]])
Genes = c("FN1", "NR2F1", "HAS2", "SOX4")  
FeaturePlot(human.object, features = Genes)

VlnPlot(human.object, features = str_split(egoResSummary["GO:0048762","geneID"], "/")[[1]], pt.size = 0.0)  & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")
VlnPlot(human.object, features = str_split(egoResSummary["GO:0048762","geneID"], "/")[[1]], pt.size = 0.0, group.by = "Sample") & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")

VlnPlot(human.object, features = str_split(egoResSummary["GO:0031589","geneID"], "/")[[1]], pt.size = 0.0)  & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")
VlnPlot(human.object, features = str_split(egoResSummary["GO:0031589","geneID"], "/")[[1]], log = T, pt.size = 0.0, group.by = "Sample") & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")

VlnPlot(human.object, features = Genes, pt.size = 0.0, ncol = 4)  & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black") 
VlnPlot(human.object, features = Genes, pt.size = 0.0, group.by = "Sample", ncol = 4) & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black") & theme(axis.text.x = element_blank())



bulk_markers <- read.csv2("20190619Top50Genes.txt", sep = "\t", header = TRUE, row.names = 1)
markers.top10.bulk = bulk_markers %>% top_n(10, log2FoldChange)
markers.bottom10.bulk = bulk_markers %>% top_n(-10, log2FoldChange)

Idents(human.object) <- human.object$orig.ident
genesTest = paste("HUMAN-", bulk_markers$geneName, sep = "")
genesTest = genesTest[genesTest %in% rownames(human.object)]
cm.markers.bulk = FindAllMarkers(object = human.object, features = genesTest, min.pct = 0.0, logfc.threshold = 0.25, return.thresh = 0.1, only.pos = T)
EnhancedVolcano::EnhancedVolcano(toptable = cm.markers.bulk, x = "avg_logFC", y = "p_val_adj", lab = cm.markers.bulk$gene, pCutoff = 0.01, FCcutoff = 0.25)
cm.markers = FindAllMarkers(object = human.object, logfc.threshold = 1.0, return.thresh = 0.01)



gene_list1 <- read.csv(file = "./csvs/ImmuneNeg.conservedmarkers_tumor_c1.csv", row.names = 1)
EnhancedVolcano(gene_list1,
                lab = rownames(gene_list1),selectLab = temp,
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 4, pointSize = 1.5, FCcutoff = 0.5, colAlpha = 0.5, drawConnectors = T) + xlim(-1.5,1.6) + ylim(0,80)

gene_list2 <- read.csv(file = "./csvs/ImmunePos.markers_tumor_c2c0.csv", row.names = 1)
EnhancedVolcano(gene_list2,
                lab = rownames(gene_list2), selectLab = temp,
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 4, pointSize = 1.5, FCcutoff = 0.5, colAlpha = 0.5, drawConnectors = T) + xlim(-1.65,1.65) + ylim(0,37)

gene_list1_sig <- gene_list1[(gene_list1$avg_log2FC > 0.5 & gene_list1$p_val_adj<0.01) | (gene_list1$avg_log2FC < -0.5 & gene_list1$p_val_adj<0.01),]
gene_list2_sig <- gene_list2[(gene_list2$avg_log2FC > 0.5 & gene_list2$p_val_adj<0.01) | (gene_list2$avg_log2FC < -0.5 & gene_list2$p_val_adj<0.01),]

temp = intersect(convert_human_to_mouse_symbols(rownames(gene_list1_sig)), rownames(gene_list2_sig))
temp = intersect(convert_mouse_to_human_symbols(rownames(gene_list2_sig)), rownames(gene_list1_sig))
length(temp)



############################################## 

library(nichenetr)
temp <- convert_human_to_mouse_symbols(rownames(human.object))
human.object.mouse <- human.object[names(temp[!is.na(temp)]),]
dim(human.object.mouse)

new_gene_names <- convert_human_to_mouse_symbols(rownames(human.object.mouse))
counts = GetAssayData(human.object.mouse, assay = "RNA", slot = "counts")
meta.data <- human.object.mouse@meta.data
rownames(counts) <- as.vector(new_gene_names)

human.object.mouse <- CreateSeuratObject(counts = counts, meta.data = meta.data)
human.object.mouse$RNA_Celltype = "Cancer cells"
human.object.mouse$Celltype_groups = "Cancer cells"
seurat.object$Celltype_groups

# Mouse analysis starts here
counts = GetAssayData(seurat.object, assay = "RNA", slot = "counts")[mouse_genes, seurat.object$organism == "MOUSE"]
rownames(counts) <- str_split_fixed(rownames(counts), pattern = "---", n = 2)[,2]
meta.data = seurat.object@meta.data[seurat.object$organism == "MOUSE", ]
mouse.object <- CreateSeuratObject(counts = counts, meta.data = meta.data)
table(mouse.object$RNA_Celltype)
dim(mouse.object)

seurat.object.hgtomm <- merge(human.object.mouse, mouse.object)
dim(seurat.object.hgtomm)

# save(seurat.object.hgtomm, file = "./Robjs_ImmuneNeg/seurat.object.immuneNeg.exp1.hgtomm.Robj")

