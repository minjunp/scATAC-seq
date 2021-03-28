#!/bin/bash

# -------
# We need only one master set -- they are all same!
# -------
#cat P1_atac_count.csv | sed '1d' | cut -f1 -d "," | sed 's|"||g' | sed 's|-|\t|g' > ATAC_peaks.bed
cat P8_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|g' | sed 's|,|\t|g' > P8_peaks.bed
#cat cardiomyocytes_atac_count.csv | sed '1d' | cut -f1 -d "," | sed 's|"||g' | sed 's|-|\t|g' > cardiomyocytes_peaks.bed
#cat endothelial_atac_count.csv | sed '1d' | cut -f1 -d "," | sed 's|"||g' | sed 's|-|\t|g' > endothelial_peaks.bed
#cat fibroblasts_atac_count.csv | sed '1d' | cut -f1 -d "," | sed 's|"||g' | sed 's|-|\t|g' > fibroblasts_peaks.bed
