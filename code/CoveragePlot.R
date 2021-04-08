#!/usr/bin/env Rscript

library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(hdf5r)
#Signac is an extension of Seurat for the analysis, interpretation, and exploration of single-cell chromatin datasets.
library(Signac)
library(patchwork)
# for integration
library(GenomicRanges)
library(harmony)
set.seed(1234)

load("/mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/obj/wh_harmony_filtered.Robj")

obj
DefaultAssay(obj) <- "peaks"
#Idents(obj) <- "predicted.id"
Idents(obj) <- "dtype"

obj@meta.data['data'] <- 'na'
obj@meta.data$data[obj@meta.data$dtype == "P1_r1"] <- "P1"
obj@meta.data$data[obj@meta.data$dtype == "P1_r2"] <- "P1"
obj@meta.data$data[obj@meta.data$dtype == "P8_w1_r1"] <- "P8"
obj@meta.data$data[obj@meta.data$dtype == "P8_w1_r2"] <- "P8"
obj@meta.data$data[obj@meta.data$dtype == "P8_w2_r2"] <- "P8"

Idents(obj) <- "data"
## DE genes for CM
Actb <- "chr5-142803654-143003654"
Actc1 <- "chr2-113952887-114152887"
Myl2 <- "chr5-122000951-122200951"
Myl3 <- "chr9-110641861-110841861"
Slc25a4 <- "chr8-46111030-46311030"
Sparc <- "chr11-55299574-55499574"
Tnnt2 <- "chr1-135736354-135936354"
Tpm1 <- "chr9-66932825-67132825"
Tnni1 <- "chr1-135679434-135879434"
Fabp3 <- "chr4-130208595-130408595"

pdf('Actb_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Actb
)
cov_plot
dev.off()

pdf('Actc1_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Actc1
)
cov_plot
dev.off()

pdf('Myl2_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Myl2
)
cov_plot
dev.off()

pdf('Myl3_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Myl3
)
cov_plot
dev.off()

pdf('Slc25a4_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Slc25a4
)
cov_plot
dev.off()

pdf('Sparc_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Sparc
)
cov_plot
dev.off()

pdf('Tnnt2_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Tnnt2
)
cov_plot
dev.off()

pdf('Tpm1_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Tpm1
)
cov_plot
dev.off()

pdf('Tnni1_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Tnni1
)
cov_plot
dev.off()

pdf('Fabp3_covplot.pdf', width=20, height=10)
cov_plot <- CoveragePlot(
  object = obj,
  region = Fabp3
)
cov_plot
dev.off()

"""
Code not working as expected...

for (gene in c(Actb, Actc1, Myl2, Myl3, Slc25a4, Sparc, Tnnt2, Tpm1, Tnni1, Fabp3)){
  cov_plot <- CoveragePlot(
    object = obj,
    region = gene,
    annotation = TRUE,
    peaks = TRUE,
    tile = TRUE,
    links = TRUE
  )

  gene.name <- deparse(substitute(gene))
  pdf_name <- sprintf("%s_cov_plot.pdf", gene.name)
  pdf('cov_plot.pdf', width=20, height=10)
  cov_plot
  dev.off()

}
"""
