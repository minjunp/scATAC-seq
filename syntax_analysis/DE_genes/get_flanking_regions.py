import pandas as pd
import numpy as np
import sys

# STEP 1: get the flanking region of the DE genes
def flanking_region(DEfile, bedfile, peakfile, threshold, pval_thres, direction):
    """
    Set flanking threshold: 100000 is 100Kb
    """

    ## Load DE genes output
    #DE_file = open(DEfile, 'r').read().split('\n') # make into vector
    #DE_file = pd.read_csv(DEfile, sep=',')

    ## Load for ALL genes not just DE genes
    gene_names = pd.read_csv(bedfile, sep='\t', header=None)
    gene_names = gene_names.iloc[:, 3].values

    ## Load TSS sites bed file
    TSS_sites = pd.read_csv(bedfile, sep='\t', header=None)

    """
    ## Filter DE genes based on direction & adjusted p-value
    if direction == 'up-regulated':
        DE_genes = DE_file[DE_file.p_val_adj < pval_thres]
        DE_genes = DE_file[DE_file.avg_log2FC > 0]
        gene_names = DE_genes.iloc[:, 0].values
    elif direction == 'down-regulated':
        DE_genes = DE_file[DE_file.p_val_adj < pval_thres]
        DE_genes = DE_file[DE_file.avg_log2FC < 0]
        gene_names = DE_genes.iloc[:, 0].values
    else:
        raise Error('Wrong direction name')
    """
    TSS_genes = TSS_sites.loc[:, 3].values

    count = 0
    with open(peakfile, 'w') as f:
        for gene in gene_names:
            if gene in TSS_genes:
                count += 1 # to see how many genes are present
                gene_index = np.where(TSS_genes == gene) # index of a gene

                gene_info = TSS_sites.iloc[gene_index]
                chr = gene_info[0].values.tolist()[0]
                tss_start = gene_info[1].values.tolist()[0]
                tss_end = gene_info[2].values.tolist()[0] # type: int64
                g = gene_info[3].values.tolist()[0] # gene name

                expand_tss_start = tss_start - threshold
                expand_tss_end = tss_end + threshold

                peak = g + '\t' + chr + '\t' + str(expand_tss_start) + '\t' + str(expand_tss_end) + '\n'
                f.write(peak)


"""
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Cardiomyocytes.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'CM_up_DEgenes_100kb_flanking_regions.bed', threshold=100000, pval_thres=0.01, direction='up-regulated')
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Cardiomyocytes.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'CM_down_DEgenes_100kb_flanking_regions.bed', threshold=100000, pval_thres=0.01, direction='down-regulated')
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Endothelial_1.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'CM_up_DEgenes_100kb_flanking_regions.bed', threshold=100000, pval_thres=0.01, direction='up-regulated')
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Endothelial_1.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'CM_down_DEgenes_100kb_flanking_regions.bed', threshold=100000, pval_thres=0.01, direction='down-regulated')
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Fibroblasts.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'CM_up_DEgenes_100kb_flanking_regions.bed', threshold=100000, pval_thres=0.01, direction='up-regulated')
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Fibroblasts.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'CM_down_DEgenes_100kb_flanking_regions.bed', threshold=100000, pval_thres=0.01, direction='down-regulated')
"""

## Cell-type does not matter anymore...
flanking_region('DEG_whole_heart_P1_vs_P8_MI_RNA_Cardiomyocytes.csv', '../data/mm10-3.0.0.premrna.tss.TAD.txt',
                'AllGenes_200kb_flanking_regions.bed', threshold=200000, pval_thres=0.01, direction='up-regulated')

df = pd.read_csv('AllGenes_200kb_flanking_regions.bed', sep='\t', header=None)
df = df.drop_duplicates()
df.to_csv('AllGenes_200kb_flanking_regions_dropped_duplicates.bed', sep='\t', header=None, index=False)
