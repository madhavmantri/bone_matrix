# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat);
packageVersion("Seurat")
library(Matrix); library(stringr); library(dplyr); library(reticulate)
library(readr); library(fitdistrplus); library(ggplot2)
library(nichenetr)
library(pals); library(RColorBrewer)
library(EnhancedVolcano)
library(clusterProfiler) 

setwd("/fs/cbsuvlaminck3/workdir/mm2937/FishbachLab/")

cc.genes.mouse <- readRDS("mouse_cell_cycle_genes/mouse_cell_cycle_genes.rds")
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")

# save(seurat.object.immunePos, file="Robjs_ImmunePos/seurat.object.immunePos.Robj")
load("Robjs_ImmunePos/seurat.object.immunePos.0.9.Robj")
dim(seurat.object)


####################################  Cancer cell analysis starts here! ################################################
DefaultAssay(seurat.object) <- "RNA"
head(seurat.object@meta.data)
human.object <- subset(seurat.object, subset = RNA_Celltype == "Cancer cells")
table(human.object$Sample)
human.object <- subset(human.object, subset = Sample %in% c("Tumor_Bone", "Tumor_DMB"))

VariableFeaturePlot(human.object)
dim(human.object)


dim(human.object)
human.object <- NormalizeData(human.object) 
human.object <- CellCycleScoring(human.object, s.features = cc.genes.mouse$s.genes, g2m.features = cc.genes.mouse$g2m.genes)

# RUn again from here after subsetting the clusters
human.object <- FindVariableFeatures(human.object)

# human.object <- ScaleData(human.object)
human.object <- ScaleData(human.object, vars.to.regress = c("S.Score", "G2M.Score", "percent.mito"))
human.object <- RunPCA(human.object)
DimPlot(human.object, reduction = "pca", group.by = "Phase")
DimPlot(human.object, reduction = "pca", group.by = "Sample")
FeaturePlot(human.object, reduction = "pca", features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "eGFP"))

ElbowPlot(object = human.object, ndims = 50)
human.object <- FindNeighbors(human.object, dims = 1:20, force.recalc = TRUE)
human.object <- FindClusters(human.object, resolution = 0.1)
human.object <- FindClusters(human.object, resolution = 0.2)
human.object <- FindClusters(human.object, resolution = 0.3)
human.object <- FindClusters(human.object, resolution = 0.5)
table(human.object$RNA_snn_res.0.3)

human.object <- RunTSNE(object = human.object, reduction = "pca", dims = 1:20)
DimPlot(human.object, reduction = "tsne", group.by = "RNA_snn_res.0.3") | DimPlot(human.object, reduction = "tsne", group.by = "Sample") | DimPlot(human.object, reduction = "tsne", group.by = "Phase")
VlnPlot(human.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))

DimPlot(human.object, reduction = "tsne", group.by = "RNA_snn_res.0.3") | DimPlot(human.object, reduction = "tsne", group.by = "Sample") | DimPlot(human.object, reduction = "tsne", group.by = "Phase") | FeaturePlot(human.object, features = c("eGFP"))

human.object <- subset(human.object, subset = RNA_snn_res.0.3 %in% c("3"), invert = T)
# Run preprocessing again after this step

library(pals)
pdf(file = "./plots/immune_pos_umap_tumor_cells.pdf", width = 1.5, height = 1.5)
DimPlot(human.object, reduction = "tsne", group.by = "RNA_snn_res.0.3", cols = "Dark2", pt.size = 0.5) + theme_void(base_size = 6) + 
  theme(text = element_text(size = 8), 
        legend.position = 'right', 
        legend.key.size = unit(0.3,"cm")) + ggtitle(NULL)
dev.off()

# save(human.object, file="Robjs_ImmunePos/tumor.object.ImmunePos.Robj")
load("Robjs_ImmunePos/tumor.object.ImmunePos.Robj")

head(human.object@meta.data)
summary <- as.data.frame.matrix(table(human.object@meta.data[c("Sample", "RNA_snn_res.0.3")]))
summary <- summary[rowSums(summary[])>0,]
summary <- summary[,colSums(summary[])>0]
summary
library(reshape)
summary_prop <- summary/rowSums(summary)
summary_prop$sample <- rownames(summary_prop)
summary_prop <- melt(as.data.frame(summary_prop), id.vars = c("sample"))
levels(summary_prop$variable)
ggplot(data = summary_prop) + geom_bar(mapping = aes(x = sample, y = value, fill = variable), stat = "identity") + 
  labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic() + theme(text = element_text(size = 16))

pdf(file = "./plots/immune_pos_tumor_composition.pdf", width = 1.7, height = 1.5)
ggplot(data = summary_prop, mapping = aes(x = sample, y = value, fill = variable)) + geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Composition", fill = "Cell types") + theme_classic() + theme(text = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "right", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10)) + scale_fill_brewer(palette = "Dark2")
dev.off()

library(ggrepel)
Idents(human.object) <- human.object$RNA_snn_res.0.3
markers.all.RNA = FindAllMarkers(object = human.object, assay = "RNA", return.thresh = 0.01, only.pos = T)
write.csv(markers.all.RNA, file = "./csvs/ImmunePos.markers.tumor_subclusters.csv")
markers.top20.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(20, avg_log2FC)
markers.top10.RNA = markers.all.RNA %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DotPlot(human.object, features = markers.top20.RNA$gene) + theme(axis.text.x = element_text(angle = 90, hjust = 1.0))

DoHeatmap(human.object, features = unique(markers.top20.RNA$gene), group.colors = brewer.dark2(3)) + scale_color_manual(
  values = brewer.dark2(3), limits = c('0', '1', '2')) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) 


FeaturePlot(human.object, features = c("Plau", "Mmp3"))

tumor_1 <- FindMarkers(human.object, assay = "RNA", group.by = "RNA_snn_res.0.3", ident.1 = "1")
EnhancedVolcano(tumor_1,
                lab = rownames(tumor_1),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 4, pointSize = 2.0, FCcutoff = 0.5, colAlpha = 0.5) + xlim(-1.6,2.2) + ylim(0,45)

tumor_2 <- FindMarkers(human.object, assay = "RNA", group.by = "RNA_snn_res.0.3", ident.1 = "2")
EnhancedVolcano(tumor_2,
                lab = rownames(tumor_2),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 4, pointSize = 2.0, FCcutoff = 0.8, colAlpha = 0.5) + xlim(-2.0,2.0)+ ylim(0,55)

library(ggbreak)
tumor_1_2 <- FindMarkers(human.object, assay = "RNA", group.by = "RNA_snn_res.0.3", logfc.threshold = 0.0, ident.1 = "1", ident.2 = "2")
write.csv(tumor_1_2, file = "csvs/ImmunePos.markers_tumor_1vs2.csv")
pdf(file = "./plots/immune_pos_volcano_tumor_1_2.pdf", width = 3.5, height = 1.7)
EnhancedVolcano(tumor_1_2,
                lab = rownames(tumor_1_2),
                x = 'avg_log2FC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-2, axisLabSize = 10, 
                labSize = 2, pointSize = 0.5, FCcutoff = 0.5, colAlpha = 0.5) + xlim(-1.65,3.3) + ylim(0,37) + theme_classic() + theme(text = element_text(size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 6, color = "black"), axis.line = element_line(size = 0.4), legend.position = "none", legend.key.size = unit(0.2,"cm"), legend.justification = "center", legend.margin=margin(0,0,0,0),                                                                                     legend.box.margin=margin(-10,-5,-10,-10)) + scale_fill_brewer(palette = "Dark2")
dev.off()


load("Robjs_ImmunePos/tumor.object.ImmunePos.Robj")
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

#Rename modules

#top gene in each module
module.genes <- c("Hist1h3c", "Plau", "Ifih1", "Ccnb2", "H2-Ea-ps", "Cyp51", "Ccn5", "Nsd2", "Hilpda", "Ankrd33b", "Serpinc1", "Senp2", "Ano4", "Arg1", "Rac2", "Epcam", "Fbxo32", "Scarb1", "Cst6", "Mad2l1", "Col6a3", "Tinagl1")
#corresponding module name for genes in module.genes (make sure indexes match)
module.names <- c("Cycle1", "Mechanoactivated.pEMT", "Interferon", "Cycle2", "Immunoregulatory", "Fat.metabolism", "Glandular", "Stress1", "Hypoxia", "Interferon2", "Non.canonical.wnt", "EMT.Bone", "Interferon1", "Macrophage.like", "m_Rac2", "Epithelial", "m_Fbxo32", "Wound.healing", "m_Cst6", "Cycle3", "Mesenchymal", "m_Tinagl1")

temp.mat <- as.matrix(modules)
for (i in 1:length(module.genes)) {
  for (j in 1:length(names(modules))) {
    if (module.genes[i] %in% temp.mat[[j]]) {
      names(modules)[j] <- module.names[i]
    }
  }
}

rownames(coefs) <- names(modules)

# Order modules
h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
o = row_order(h)
scores = scores[, o]
coefs = coefs[o, ]
modules = modules[o]
print(modules)
length(modules)

lapply(modules, write, "./csvs/ImmunePos.NMFmodules.txt", append=TRUE, ncolumns=1000)

human.object = AddMetaData(human.object, t(coefs), col.name = rownames(coefs))
# Cluster NMF output
# h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
hcl = as.hclust(column_dend(h))
sig = cutree(hcl, k = length(modules))
nmf = c(by(t(coefs), INDICES = sig, FUN = function(x){names(modules)[which.max(colMeans(x))]}))[sig]
length(nmf)
dim(human.object)

human.object$nmf = factor(nmf, levels = names(modules))
library("RColorBrewer")
col_nmf = c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Set1'))[1:nlevels(human.object$nmf)]
names(col_nmf) = levels(human.object$nmf)
## Saving
human.object
save(human.object, file = 'Robjs_ImmunePos/tumor.object_nmf.RData')
save(res.list, file = 'Robjs_ImmunePos/res.list_nmf.RData')

load('Robjs_ImmunePos/human.object_nmf.RData')
load('Robjs_ImmunePos/res.list_nmf.RData')

DimPlot(human.object, pt.size = 2, reduction = reduction, group.by = 'nmf') | DimPlot(human.object, pt.size = 2, reduction = reduction, group.by = 'Sample')

for (reduction in c('tsne')){
  h = DimPlot(human.object, pt.size = 2, reduction = reduction, group.by = 'nmf')
  print(h)
  h = DimPlot(human.object, pt.size = 2, reduction = reduction, group.by = 'Sample')
  print(h)
  #h = FeaturePlot(human.object, reduction = reduction, features = names(modules))
  #print(h)
}

VlnPlot(human.object, features = "Serpinc1", group.by = "Sample")


Genes = str_split(summaryEgoRes["GO:0006119","geneID"], "/")[[1]][1:4]
Genes = c("Pgk1", "Higd1a", "Gas5", "Ndrg1")  
FeaturePlot(human.object, features = Genes, ncol = 4)

VlnPlot(human.object, features = str_split(egoResSummary["GO:0048762","geneID"], "/")[[1]], pt.size = 0.0)  & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")
VlnPlot(human.object, features = str_split(egoResSummary["GO:0048762","geneID"], "/")[[1]], pt.size = 0.0, group.by = "Sample") & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")

VlnPlot(human.object, features = str_split(egoResSummary["GO:0031589","geneID"], "/")[[1]], pt.size = 0.0)  & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black")
VlnPlot(human.object, features = str_split(egoResSummary["GO:0031589","geneID"], "/")[[1]], log = T, pt.size = 0.0, group.by = "Sample") & stat_summary(fun = "mean",  geom = "crossbar", width = 0.5, colour = "black")

VlnPlot(human.object, features = Genes, pt.size = 0.0, ncol = 4)  & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black") 
VlnPlot(human.object, features = Genes, pt.size = 0.0, group.by = "Sample", ncol = 4) & stat_summary(fun = "mean",  geom = "crossbar",  width = 0.5, colour = "black") & theme(axis.text.x = element_blank())




# save(human.object, file=here("robjs.round2", "human.object.round2.Robj"))
load(here("robjs.round2", "human.object.round2.Robj"))

bulk_markers <- read.csv2("20190619Top50Genes.txt", sep = "\t", header = TRUE, row.names = 1)
markers.top10.bulk = bulk_markers %>% top_n(10, log2FoldChange)
markers.bottom10.bulk = bulk_markers %>% top_n(-10, log2FoldChange)

Idents(human.object) <- human.object$Sample
genesTest = bulk_markers$geneName
genesTest = genesTest[genesTest %in% rownames(human.object)]
cm.markers.bulk = FindAllMarkers(object = human.object, features = genesTest)
EnhancedVolcano::EnhancedVolcano(toptable = cm.markers.bulk, x = "avg_log2FC", y = "p_val_adj", lab = cm.markers.bulk$gene, pCutoff = 0.01, FCcutoff = 0.25)
cm.markers = FindAllMarkers(object = human.object, logfc.threshold = 1.0, return.thresh = 0.01)

