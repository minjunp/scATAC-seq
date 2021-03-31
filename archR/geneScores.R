#!/usr/bin/env Rscript

library(ArchR)
set.seed(1)

# Load the saved model (Up to doublets filtered)
proj <- loadArchRProject(path = "save-projWH-harmony")
#proj <- loadArchRProject(path = "save-projWH-v2")

markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8",
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

# 7.4 - 7.5 Visualizing Marker Genes on an Embedding - Marker Genes Imputation with MAGIC
proj <- addImputeWeights(proj)

markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

png('CD14.png', width=800, height=600)
p$CD14
dev.off()

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p,
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)


p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "Clusters",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
)

#grid::grid.newpage()
#grid::grid.draw(p$CD14)

plotPDF(plotList = p,
    name = "Plot-Tracks-Marker-Genes.pdf",
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)

### SAVE Project ###
saveArchRProject(ArchRProj = proj, outputDirectory = "save-projWH-harmony-v2", load = FALSE)
