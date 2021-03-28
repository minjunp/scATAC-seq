#!/usr/bin/env Rscript
library(ArchR)
set.seed(1)

# Load the saved model (Up to doublets filtered)
proj <- loadArchRProject(path = "save-projWH")

# Dimensionality Reduction and Clustering
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# Add Harmony for batch correction
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# To call clusters in this reduced dimension sub-space, we use the addClusters() function which uses Seurat’s graph clustering as the default clustering method.
proj <- addClusters(input = proj, reducedDims = "Harmony")

# Visualizing in a 2D UMAP Embedding
proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony")

### SAVE Project ###
saveArchRProject(ArchRProj = proj, outputDirectory = "save-projWH-harmony", load = FALSE)

# We can color by “Sample” or "Clusters":
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

png('wh_archR_harmony_umap.png', width=800, height=600)
ggAlignPlots(p1, p2, type = "h")
dev.off()

# To save an editable vectorized version of this plot, we use the plotPDF() function.
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
