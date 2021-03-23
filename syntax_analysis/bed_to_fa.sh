#!/bin/bash

cat P1_atac_count.csv | sed '1d' | cut -f1 -d "," | sed 's|"||g' | sed 's|-|\t|g' > ATAC_peaks.bed
