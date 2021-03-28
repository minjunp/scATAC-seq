#!/bin/bash

module load meme
module load bedtools

fasta_file="ATAC_peaks.fa"
# $1: subgroup, e.g. g1 g2 g3...
# ./run_fimo.sh g1 &

#files=$(ls ./indiv_motifs)
files=$(ls ./cisBP_separate_subfolders/$1)
mkdir -p cisBP_motif_outputs

for i in $files
do
	motif_name=$(echo $i | cut -f1 -d".")
	mkdir -p ./cisBP_motif_outputs/$motif_name
	echo $motif_name

	thr="1e-4"
	fimo --max-stored-scores 1600000 --verbosity 1 --thresh ${thr} --oc ./cisBP_motif_outputs/${motif_name} ./cisBP_separate_subfolders/$1/${i} ${fasta_file}
done
