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

load('wh_harmony.Robj')
DefaultAssay(obj) <- 'peaks'
Idents(obj) <- 'seurat_clusters'

# Remove cluster 3 -- which shows low sequencing depth
obj <- subset(x = obj, idents = 3, invert = TRUE)

pdf('wh_after_filter.pdf')
DimPlot(object = obj, group.by = 'predicted.id', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

umap_coords <- obj[["umap"]]@cell.embeddings
write.csv(as.data.frame(umap_coords), 'umap_coordinates.csv')

save(obj,file='wh_harmony_filtered.Robj')
