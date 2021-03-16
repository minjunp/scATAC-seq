#!/usr/bin/env Rscript
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(hdf5r)
#Signac is an extension of Seurat for the analysis, interpretation, and exploration of single-cell chromatin datasets.
library(Signac)
library(patchwork)
library(harmony)
# for integration
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

load('./obj/wh_harmony.Robj')
# Plot a legend to map colors to expression levels
png('umap_peak_region_fragments.png', width=800, height=600)
FeaturePlot(obj, features = "peak_region_fragments")
dev.off()

png('blacklist_ratio_peak_region_fragments.png', width=800, height=600)
FeaturePlot(obj, features = "blacklist_ratio")
dev.off()

png('nucleosome_signal_peak_region_fragments.png', width=800, height=600)
FeaturePlot(obj, features = "nucleosome_signal")
dev.off()

baseplot <- DimPlot(obj, reduction = "umap")

# Add custom labels and titles
png('umap_baseplot.png', width=800, height=600)
baseplot + labs(title = "P1 & P8 Whole Heart")
dev.off()

# Use community-created themes, overwriting the default Seurat-applied theme Install ggmin with
remotes::install_github('sjessa/ggmin')
png('umap_ggmin_baseplot.png', width=800, height=600)
baseplot + ggmin::theme_powerpoint()
dev.off()

# Seurat also provides several built-in themes, such as DarkTheme; for more details see
# ?SeuratTheme
png('umap_dark_baseplot.png', width=800, height=600)
baseplot + DarkTheme()
dev.off()

#dtype
#P1_r1
#P1_r2
#P8_w1_r1
#P8_w1_r2
#P8_w2_r2
obj@meta.data['data'] <- 'na'
obj@meta.data$data[obj@meta.data$dtype == "P1_r1"] <- "P1"
obj@meta.data$data[obj@meta.data$dtype == "P1_r2"] <- "P1"
obj@meta.data$data[obj@meta.data$dtype == "P8_w1_r1"] <- "P8"
obj@meta.data$data[obj@meta.data$dtype == "P8_w1_r2"] <- "P8"
obj@meta.data$data[obj@meta.data$dtype == "P8_w2_r2"] <- "P8"

Idents(obj) <- 'data'
png('wh_umap_by_data.png', width=800, height=600)
DimPlot(object = obj, group.by = 'data', repel = TRUE, label = TRUE) + NoLegend()
dev.off()
