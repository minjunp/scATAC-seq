import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation

DEfile = '/Users/minjunpark/Documents/hali/syntax_analysis/DE_genes/DEG_whole_heart_P1_vs_P8_MI_RNA_Cardiomyocytes.csv'
## Load DE genes output
DE_file = open(DEfile, 'r').read().split('\n') # make into vector
DE_file = pd.read_csv(DEfile, sep=',')

## Split into up-regulated, down-regulated, and non-DE genes (adj p-val < 0.01 & abs(log2FC) > 1.25)
p_val_thres = 0.01
log2FC_thres = 1.25

up_genes = DE_file[DE_file.p_val_adj < p_val_thres]
up_genes = up_genes[up_genes.avg_log2FC > log2FC_thres]
up_gene_names = up_genes.iloc[:, 0].values

down_genes = DE_file[DE_file.p_val_adj < p_val_thres]
down_genes = down_genes[down_genes.avg_log2FC < -log2FC_thres]
down_gene_names = down_genes.iloc[:, 0].values

subset = 'P8_cm'

motifs = ['ETS1', 'ETS2',
          'FOS', 'FOSB', 'FOSB', 'FOSL1', 'FOSL2', 'FOXA1', 'FOXA2', 'FOXA3', 'FOXB1', 'FOXC1', 'FOXC2', 'FOXD1', 'FOXD3',
          'FOXG1', 'FOXI1', 'FOXJ1', 'FOXJ3', 'FOXK1', 'FOXL1', 'FOXL2', 'FOXN1', 'FOXN4', 'FOXO1', 'FOXO3', 'FOXO4', 'FOXO6', 'FOXP1', 'FOXP2', 'FOXQ1',
          'JUN', 'JUNB', 'JUND',
          'TEAD1', 'TEAD2', 'TEAD4']
for motif in motifs:
    plt.cla()   # Clear axis
    plt.clf()   # Clear figure

    #motif = 'ETS1'
    influence_score = '/Users/minjunpark/Documents/hali/syntax_analysis/influence_score/{}/{}_influence_score.csv'.format(subset, motif)
    df = pd.read_csv(influence_score, index_col=0)
    df['direction'] = 'void'

    for i in range(len(df)):
        if df.iloc[i,0] in up_gene_names:
            df.iloc[i,2] = 'up-regulated'
        elif df.iloc[i,0] in down_gene_names:
            df.iloc[i,2] = 'down-regulated'
        else:
            df.iloc[i,2] = 'non-DE'

    sns.set(style="whitegrid")
    x = "direction"
    y = "influence score"
    order = ['up-regulated', 'down-regulated', 'non-DE']
    ax = sns.boxplot(data=df, x=x, y=y, order=order)
    ax.set_ylabel('Normalized Motif Score')
    ax.set_title('{} box plot'.format(motif))
    add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                        box_pairs=[("up-regulated", "down-regulated"), ("up-regulated", "non-DE"), ("down-regulated", "non-DE")],
                        test='Mann-Whitney', text_format='full', loc='inside', verbose=2)
    plt.savefig('../figs/{}/{}_NMS.pdf'.format(subset, motif))

"""
df_up_gene = df.loc[df['gene name'].isin(up_gene_names)]
df_down_gene = df.loc[df['gene name'].isin(down_gene_names)]
df_nonDE_gene = df.loc[~df['gene name'].isin(np.concatenate((up_gene_names, down_gene_names)))]

df_up_score = df_up_gene['influence score']
df_down_score = df_down_gene['influence score']
df_nonDE_score = df_nonDE_gene['influence score']


### Prettier way
combined = np.array([df_up_score, df_down_score, df_nonDE_score])
labels = ['up-regulated genes', 'down-regulated genes', 'non-DE genes'
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
# rectangular box plot
bplot = ax.boxplot(combined,
                     vert=True,  # vertical box alignment
                     patch_artist=True,  # fill with color
                     labels=labels)  # will be used to label x-ticks

ax.set_title('{} box plot'.format(motif))
# adding horizontal grid lines
ax.yaxis.grid(True)
ax.set_xlabel('Direction')
ax.set_ylabel('Normalized Motif Score')
# fill with colors
colors = ['pink', 'lightblue', 'lightgreen']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
plt.show()
"""
