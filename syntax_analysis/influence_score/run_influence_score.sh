#!/bin/bash

# Fimo processed files
fimo_files=$(ls /mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/fimo_counts/processed*)
cell_type='cardiomyocytes_peaks.bed'
allGenes='AllGenes_100kb_flanking_regions_dropped_duplicates_v2.bed'

for fimo_file in ${fimo_files}
  do
    motif=$(echo ${fimo_file} | cut -f9 -d"_")
    python calculate_influence_score.py ${fimo_file} ${cell_type} ${allGenes} ${motif}_influence_score.csv
  done

#python calculate_influence_score.py processed_ETS1_fimo.tsv cardiomyocytes_peaks.bed AllGenes_100kb_flanking_regions_dropped_duplicates.bed ETS1_influence_score.csv
