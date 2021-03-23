#!/usr/bin/env Rscript

library(ArchR)
set.seed(1)
#addArchRThreads(threads = 2)

#inputFiles <- getTutorialData("Hematopoiesis")
inputFiles <- read.csv('configs.txt', header = TRUE, sep = ",")
df <- as.data.frame(inputFiles)
df_chr <- unlist(df)

addArchRGenome("mm10")

ArrowFiles <- createArrowFiles(
  inputFiles = df_chr,
  sampleNames = names(df_chr),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

# Inferring Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

# Create ArchR Project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "outs",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
