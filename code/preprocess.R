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

count <- '/mount/samee/hali_data/aggr_results_whole_heart/agg_whole_heart/outs/filtered_peak_bc_matrix.h5'
meta <- '/mount/samee/hali_data/aggr_results_whole_heart/minjun_analysis/preprocess/singlecell_v2.csv'
frag <- '/mount/samee/hali_data/aggr_results_whole_heart/agg_whole_heart/outs/fragments.tsv.gz'
# Load the pre-processed scRNA-seq data
snRNA <- readRDS("/mount/samee/hali_data/seurat_whole_heart/RNAseq/integrated_v2.rds")

counts <- Read10X_h5(count)
metadata <- read.csv(
    file = meta,
    header = TRUE,
    row.names = 1
)

chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = frag,
    min.cells = 1
)

obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = 'peaks',
    meta.data = metadata
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(obj) <- annotations

# compute nucleosome signal score per cell
obj <- NucleosomeSignal(object = obj)
# compute TSS enrichment score per cell
obj <- TSSEnrichment(object = obj, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$passed_filters * 100
obj$blacklist_ratio <- obj$blacklist_region_fragments / obj$peak_region_fragments

# I will use the threshold:
obj <- subset(
	         x = obj,
		   subset = peak_region_fragments > 1500 &
			       peak_region_fragments < 30000 &
			           pct_reads_in_peaks > 15
			   )

png('vln_after_filter.png', width=800, height=600)
VlnPlot(
 object = obj,
 features = c('pct_reads_in_peaks', 'peak_region_fragments',
              'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
 pt.size = 0.1,
 ncol = 5
)
dev.off()

Idents(obj) <- 'dtype'

# make sure to change to the assay containing common peaks
DefaultAssay(obj) <- "peaks"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 20)
obj <- RunSVD(
	         obj,
		   reduction.key = 'LSI_',
		   reduction.name = 'lsi',
		     irlba.work = 400
		   )
obj <- RunUMAP(obj, dims = 2:30, reduction = 'lsi')

pdf('depthcor.pdf')
DepthCor(obj)
dev.off()

png('wh_before_harmony_label.png', width=800, height=600)
DimPlot(object = obj, group.by = 'dtype', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

# Get ready for snRNA-seq integration
gene.activities <- GeneActivity(obj)

# convert rownames from chromsomal coordinates into gene names
#gene.key <- genebodyandpromoter.coords$gene_name
#names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
#rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
#gene.activities <- gene.activities[rownames(gene.activities)!="",]

#Add the gene activity matrix to the Seurat object as a new assay, and normalize it
obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
obj <- NormalizeData(
		        object = obj,
			  assay = 'RNA',
			  normalization.method = 'LogNormalize',
			    scale.factor = median(obj$nCount_RNA)
			  )

# save the object
save(obj,file="wh_before_harmony.Robj")

# run harmony
obj <- RunHarmony(
  object = obj,
  group.by.vars = "types",
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

obj <- RunUMAP(object = obj, dims= 2:30, reduction = 'harmony')
obj <- FindNeighbors(object = obj, reduction = 'harmony', dims = 2:50, graph.name = 'harmony_snn')
obj <- FindClusters(object = obj, graph.name = 'harmony_snn', verbose = FALSE, algorithm = 3, resolution = 1.0)
DefaultAssay(obj) <- 'RNA'

transfer.anchors <- FindTransferAnchors(reference = snRNA, query = obj, reduction = 'cca', dims = 2:50)

predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = snRNA$celltypes_v1, weight.reduction = obj[['harmony']], dims = 2:50)

obj <- AddMetaData(object = obj, metadata = predicted.labels)

png('wh_harmony_predicted_id.png', width=800, height=600)
DimPlot(object = obj, group.by = 'predicted.id', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

pdf('wh_harmony_predicted_id.pdf')
DimPlot(object = obj, group.by = 'predicted.id', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

png('wh_harmony_dtype.png', width=800, height=600)
DimPlot(object = obj, group.by = 'dtype', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

png('wh_harmony_seurat_clusters.png', width=800, height=600)
DimPlot(object = obj, group.by = 'seurat_clusters', repel = TRUE, label = TRUE) + NoLegend()
dev.off()

write.csv(obj@meta.data, 'mdata.csv')

# save the object
save(obj,file='wh_harmony.Robj')
