#!/usr/bin/env Rscript

library(ArchR)
set.seed(1)

proj <- loadArchRProject(path = "save-projWH-harmony-v2")

scRNA <- readRDS("/mount/samee/hali_data/seurat_whole_heart/RNAseq/integrated_v2.rds")

# The metadata contains a column called BioClassification which contains the cell type classifications for each cell in the scRNA-seq dataset.
#colnames(colData(scRNA))
# Using table() we can see how many cells are in each of the scRNA-seq cell type classifications.
#table(colData(scRNA)$celltypes_v1)

# 8.1.2 Constrained Integration
#cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
#preClust <- colnames(cM)[apply(cM, 1 , which.max)]
#cbind(preClust, rownames(cM)) #Assignments

# First, lets look at the cell type labels from our scRNA-seq data that were used in our unconstrained integration:
#unique(unique(proj$predictedGroup_Un))

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = scRNA,
    addToArrow = FALSE,
    groupRNA = "celltypes_v1",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

p1 <- plotEmbedding(
    proj,
    colorBy = "cellColData",
    name = "predictedGroup_Un",
    pal = pal
)

plotPDF(plotList = p1,
    name = "Plot-UMAP-RNA-Integration.pdf.pdf",
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = "save-projWH-harmony-v3", load = FALSE)
