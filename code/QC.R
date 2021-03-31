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

load("LAPV_detailed.Robj")
Idents(obj) <- 'dtype'

library(gridExtra)
pdf("scATAC_QC.pdf", width=10, height=10, useDingbats=FALSE)
options(repr.plot.width=10, repr.plot.height=10)
temp_combo <- obj@meta.data
g <- ggplot(temp_combo, aes(x = dtype))


p1 <- g + geom_violin(aes(y = pct_reads_in_peaks, fill = dtype)) +
        theme(legend.position = "none") +
        geom_boxplot(aes(y = pct_reads_in_peaks), width = 0.1) + theme_minimal()

p2 <- g + geom_violin(aes(y = peak_region_fragments, fill = dtype)) +
        theme(legend.position = "none") +
        geom_boxplot(aes(y = peak_region_fragments), width = 0.1) + theme_minimal()

p3 <- g + geom_violin(aes(y = nucleosome_signal, fill = dtype)) +
        theme(legend.position = "none") +
        geom_boxplot(aes(y = nucleosome_signal), width = 0.1) + theme_minimal()

p4 <- g + geom_violin(aes(y = blacklist_ratio, fill = dtype)) +
        theme(legend.position = "none") +
        geom_boxplot(aes(y = blacklist_ratio), width = 0.1) + theme_minimal()

p5 <- g + geom_violin(aes(y = TSS.enrichment, fill = dtype)) +
        theme(legend.position = "none") +
        geom_boxplot(aes(y = TSS.enrichment), width = 0.1) + theme_minimal()

grid.arrange(p1, p2, p3, p4, p5, ncol = 1)

dev.off()
