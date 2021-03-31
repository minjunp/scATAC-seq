#!/bin/bash

# -------
# We need only one master set -- they are all same!
# To include the count -- they are all different!
# -------

#cat P1_atac_count.csv | sed '1d' | cut -f1 -d "," | sed 's|"||g' | sed 's|-|\t|g' > ATAC_peaks.bed

#cat P1_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P1_peaks.bed
#cat P8_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P8_peaks.bed
#cat cardiomyocytes_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > cardiomyocytes_peaks.bed
#cat endothelial_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > endothelial_peaks.bed
#cat fibroblasts_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > fibroblasts_peaks.bed

cat P1_cardiomyocytes_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P1_cardiomyocytes_peaks.bed
cat P1_endothelial_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P1_endothelial_peaks.bed
cat P1_fibroblasts_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P1_fibroblasts_peaks.bed

cat P8_cardiomyocytes_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P8_cardiomyocytes_peaks.bed
cat P8_endothelial_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P8_endothelial_peaks.bed
cat P8_fibroblasts_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|' | sed 's|-|\t|' | sed 's|,|\t|g' > P8_fibroblasts_peaks.bed
#cat P1_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|g' | sed 's|,|\t|g' > P1_peaks.bed
#cat P8_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|g' | sed 's|,|\t|g' > P8_peaks.bed
#cat cardiomyocytes_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|g' | sed 's|,|\t|g' > cardiomyocytes_peaks.bed
#cat endothelial_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|g' | sed 's|,|\t|g' > endothelial_peaks.bed
#cat fibroblasts_atac_count.csv | sed '1d' | cut -f1,2,3 -d "," | sed 's|"||g' | sed 's|-|\t|g' | sed 's|,|\t|g' > fibroblasts_peaks.bed

#cat AllGenes_100kb_flanking_regions_dropped_duplicates.bed | awk -F"\t" '{print $2,$3,$4}' | sed 's/ /\t/g'  > AllGenes_100kb_flanking_regions_dropped_duplicates_v2.bed
