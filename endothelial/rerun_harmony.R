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

load('ec_harmony.Robj')
snRNA <- readRDS("/mount/samee/hali_data/ec_cd31_scrnaseq/ec_integrated_v3.rds")

endo_cells <- subset(obj, subset = (predicted.id == "Endothelial" | predicted.id == "Endocardium"))

endo_cells <- RunUMAP(object = endo_cells, dims= 2:30, reduction = 'harmony')
endo_cells <- FindNeighbors(object = endo_cells, reduction = 'harmony', dims = 2:50, graph.name = 'harmony_snn')
endo_cells <- FindClusters(object = endo_cells, graph.name = 'harmony_snn', verbose = FALSE, algorithm = 3, resolution = 1.0)
umapCoord <- as.data.frame(Embeddings(object = endo_cells[["umap"]]))
endo_cells@meta.data <- cbind(endo_cells@meta.data, umapCoord)
DefaultAssay(endo_cells) <- 'RNA'

transfer.anchors <- FindTransferAnchors(reference = snRNA, query = endo_cells, reduction = 'cca', dims = 2:50)
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = snRNA$celltypes_v3, weight.reduction = endo_cells[['harmony']], dims = 2:50)
endo_cells <- AddMetaData(object = endo_cells, metadata = predicted.labels)

png('ec_rerun__harmony_predicted_id.png', width=800, height=600)
DimPlot(object = endo_cells, group.by = 'predicted.id', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

pdf('ec_rerun__harmony_predicted_id.pdf')
DimPlot(object = endo_cells, group.by = 'predicted.id', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

# save the object
save(endo_cells,file='endo_cells.Robj')
