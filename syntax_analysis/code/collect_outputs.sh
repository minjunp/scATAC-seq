#!/bin/bash

mkdir fimo_outs
files=$(ls ./cisBP_motif_outputs)
echo $files

for i in $files
do
	cp ./cisBP_motif_outputs/${i}/fimo.tsv fimo_outs/${i}_fimo.tsv
done
