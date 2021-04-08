#!/bin/bash

## Get the peak information for all TFs
# example: ./data_preprocess.sh
workdir='/mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/tiled_heatmap'
mkdir -p ${workdir}/TFs

cd /mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/fimo_outs

for tf in $(ls *.tsv)
do
	tf_name=$(echo $tf | cut -f1 -d"_")
	cat $tf | awk '{print $3}' | sed 's/ /\t/g' | sed 's/:/\t/' | sed 's/-/\t/' | sed '1d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' > ${workdir}/TFs/${tf_name}.bed
	echo $tf_name
done
