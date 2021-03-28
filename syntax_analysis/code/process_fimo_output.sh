#!/bin/bash
cd /mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/fimo_outs

WORKDIR='/mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo/fimo_counts'
BASEDIR='/mount/samee/hali_data/aggr_results_whole_heart/seurat_analysis/syntax_analysis/fimo'
mkdir -p $WORKDIR
for f in TEAD* FOX* JUN* ETS* FOS*
do
	head -n -4 ${f} | sed '1d' > ${WORKDIR}/${f}
	python ${BASEDIR}/process_fimo_output.py ${WORKDIR}/${f} ${WORKDIR}/processed_${f}
done
