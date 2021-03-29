#!/bin/bash

## Find the intersection with up- & down-regulated genes for LAcm1, PVcm1, PVcm2 data
module load bedtools

TF_dir="/mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/tiled_heatmap/TFs"
cd ${TF_dir}
WORKDIR="/mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/tiled_heatmap"

all_genes="${WORKDIR}/AllGenes_100kb_flanking_regions_dropped_duplicates_v2.bed"
peaks="${WORKDIR}/ATAC_peaks.bed"
bedtools intersect -u -a ${peaks} -b ${all_genes} > ${WORKDIR}/genes_peaks_subset.bed

for motif in TEAD* FOX* JUN* ETS* FOS*
do
	motif_name=$(echo ${motif} | cut -f1 -d'.')
	mkdir -p ${WORKDIR}/outs/${motif_name}_outs

	# Get the subset of peaks that 100 kb overlapping with all genes --> later filter by up- and down-regulated genes

	bedtools intersect -u -a ${motif} -b ${WORKDIR}/genes_peaks_subset.bed > ${WORKDIR}/motif_gene.bed

	for tf in $(ls *.bed)
	do
		bedtools intersect -u -a $tf -b ${WORKDIR}/genes_peaks_subset.bed > ${WORKDIR}/${motif_name}_gene.bed
		tf_name=$(echo $tf | cut -f10 -d"/" | cut -f1 -d".")
		echo $tf_name
		cat ${WORKDIR}/motif_gene.bed | wc -l
		cat ${WORKDIR}/${motif_name}_gene.bed | wc -l
		bedtools intersect -c -a ${WORKDIR}/motif_gene.bed -b ${WORKDIR}/${motif_name}_gene.bed > ${WORKDIR}/outs/${motif_name}_outs/${tf_name}.bed
	done

	rm ${WORKDIR}/motif_gene.bed
	rm ${WORKDIR}/genes_peaks_subset.bed
	rm ${WORKDIR}/${motif_name}_gene.bed
done
