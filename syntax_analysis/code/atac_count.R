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

load('wh_harmony_filtered.Robj')

obj@meta.data['data'] <- 'na'
obj@meta.data$data[obj@meta.data$dtype == "P1_r1"] <- "P1"
obj@meta.data$data[obj@meta.data$dtype == "P1_r2"] <- "P1"
obj@meta.data$data[obj@meta.data$dtype == "P8_w1_r1"] <- "P8"
obj@meta.data$data[obj@meta.data$dtype == "P8_w1_r2"] <- "P8"
obj@meta.data$data[obj@meta.data$dtype == "P8_w2_r2"] <- "P8"

P1 <- subset(obj, subset = data == "P1")
P8 <- subset(obj, subset = data == "P8")
Idents(P1) <- "data"
Idents(P8) <- "data"

P1_count <- Seurat::AverageExpression(P1, slot="counts")
P1_atac_count <- P1_count[1]
write.csv(P1_atac_count, "P1_atac_count.csv")

P8_count <- Seurat::AverageExpression(P8, slot="counts")
P8_atac_count <- P8_count[1]
write.csv(P8_atac_count, "P8_atac_count.csv")

cardiomyocytes <- subset(obj, subset = predicted.id == "Cardiomyocytes")
endothelial <- subset(obj, subset = predicted.id == "Endothelial")
fibroblasts <- subset(obj, subset = predicted.id == "Fibroblasts")
Idents(cardiomyocytes) <- "predicted.id"
Idents(endothelial) <- "predicted.id"
Idents(fibroblasts) <- "predicted.id"

cardiomyocytes_count <- Seurat::AverageExpression(cardiomyocytes, slot="counts")
cardiomyocytes_atac_count <- cardiomyocytes_count[1]
write.csv(cardiomyocytes_atac_count, "cardiomyocytes_atac_count.csv")

endothelial_count <- Seurat::AverageExpression(endothelial, slot="counts")
endothelial_atac_count <- endothelial_count[1]
write.csv(endothelial_atac_count, "endothelial_atac_count.csv")

fibroblasts_count <- Seurat::AverageExpression(fibroblasts, slot="counts")
fibroblasts_atac_count <- fibroblasts_count[1]
write.csv(fibroblasts_atac_count, "fibroblasts_atac_count.csv")
