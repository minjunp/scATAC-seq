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

#predicted.ids
#P1_r1
#P1_r2
#P8_w1_r1
#P8_w1_r2
#P8_w2_r2
cm <- subset(obj, subset = predicted.id == "Cardiomyocytes")
Idents(cm) <- 'dtype'
da_cm_genes <- FindMarkers(
  object = cm,
  ident.1 = c('P1_r1', 'P1_r2'),
  ident.2 = c('P8_w1_r1', 'P8_w1_r2', 'P8_w2_r2'),
  min.pct = 0.2,
  test.use = 'LR'
)
write.csv(da_cm_genes, 'da_cm_genes.csv')

DefaultAssay(cm) <- 'peaks'
da_cm_genes <- FindMarkers(
  object = cm,
  ident.1 = c('P1_r1', 'P1_r2'),
  ident.2 = c('P8_w1_r1', 'P8_w1_r2', 'P8_w2_r2'),
  min.pct = 0.2,
  test.use = 'LR'
)
write.csv(da_cm_genes, 'da_cm_peaks.csv')
