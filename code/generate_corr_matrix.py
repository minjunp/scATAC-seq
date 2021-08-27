from os import listdir
from os.path import isdir,isfile, join
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy

def corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName):
    #example p1 cm peaks
    all_peaks = pd.read_csv('/Users/minjunpark/Documents/Pitx2/software/co-occupancy/peaks.bed', sep='\t', header=None)
    all_peaks.head()

    # get the subset of interest
    vec = []
    for i in range(len(rawData)):
        dff = all_peaks[(all_peaks[0] == rawData.iloc[i, 0]) & (all_peaks[1] == rawData.iloc[i, 1]) & (all_peaks[2] == rawData.iloc[i, 2])]
        if len(dff) == 1:
            vec.append('exist')
        else:
            vec.append('null')

    rawData['peak_info'] = vec
    rawData = rawData[rawData.peak_info == 'exist']
    rawData = rawData.drop(['peak_info'], axis=1)

    end = rawData[2]
    start = rawData[1]
    peak_length = np.subtract(end, start).values
    peak_length = pd.Series(peak_length)

    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    df = pd.DataFrame()
    TFs = []
    for i in range(len(onlyfiles)):
        f = pd.read_csv(f'{mypath}/{onlyfiles[i]}', sep='\t', header=None)
        TF = onlyfiles[i].split('.')[0]
        df[i] = f[3]
        TFs.append(TF)
    df.columns = TFs

    df['peak_info'] = vec
    df = df[df.peak_info == 'exist']
    df = df.drop(['peak_info'], axis=1)

    # CM1 or CM2
    s1 = merged[merged[celltype] == True]
    s1 = s1[s1["Family"] != 'Homeodomain']
    s1 = s1[s1["Family"] != 'Homeodomain,Paired']
    s1 = s1[s1["Family"] != 'Homeodomain,POU']
    s1 = s1[s1["Family"] != 'CUT,Homeodomain']
    subset_genes = s1['genes'].values

    subs = df[subset_genes].copy()
    subs[subs > 1] = 1
    thres = len(subs.iloc[:,1]) * threshold
    temp = pd.DataFrame(subs.sum(axis=0))
    df_core = temp.loc[temp[0] > thres].index

    df_subs = df[df_core].copy()
    # normalize by peak lengths
    subs_normalized = df_subs.apply(lambda x: np.asarray(x) / np.asarray(peak_length))
    subs_normalized = subs_normalized * 1000

    corr = subs_normalized.corr()
    corr = corr.fillna(0)

    corr.to_csv(outputName)
    print('Correlation matrix created')

mypath = '/Users/minjunpark/Documents/Pitx2/manuscript/fig3/data/LAdown_bicluster_outs'
filename = f'{mypath}/A1JVI6.bed'
rawData = pd.read_csv(filename, sep='\t', header=None)
merged = pd.read_csv('df_merged.csv', index_col=0)
celltype = 'RNA.CM_Atrial'
#celltype = 'RNA.CM_Myh7b...'
threshold = 0.25
outputName = 'LAdown_CM1_corr_25.csv'
corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName)

mypath = '/Users/minjunpark/Documents/Pitx2/manuscript/fig3/data/LAup_bicluster_outs'
filename = f'{mypath}/A1JVI6.bed'
rawData = pd.read_csv(filename, sep='\t', header=None)
merged = pd.read_csv('df_merged.csv', index_col=0)
celltype = 'RNA.CM_Atrial'
#celltype = 'RNA.CM_Myh7b...'
threshold = 0.25
outputName = 'LAup_CM1_corr_25.csv'
corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName)

mypath = '/Users/minjunpark/Documents/Pitx2/manuscript/fig3/data/PVdown_bicluster_outs'
filename = f'{mypath}/A1JVI6.bed'
rawData = pd.read_csv(filename, sep='\t', header=None)
merged = pd.read_csv('df_merged.csv', index_col=0)
celltype = 'RNA.CM_Atrial'
#celltype = 'RNA.CM_Myh7b...'
threshold = 0.25
outputName = 'PVdown_CM1_corr_25.csv'
corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName)

mypath = '/Users/minjunpark/Documents/Pitx2/manuscript/fig3/data/PVup_bicluster_outs'
filename = f'{mypath}/A1JVI6.bed'
rawData = pd.read_csv(filename, sep='\t', header=None)
merged = pd.read_csv('df_merged.csv', index_col=0)
celltype = 'RNA.CM_Atrial'
#celltype = 'RNA.CM_Myh7b...'
threshold = 0.25
outputName = 'PVup_CM1_corr_25.csv'
corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName)

mypath = '/Users/minjunpark/Documents/Pitx2/manuscript/fig3/data/PVdown_bicluster_outs'
filename = f'{mypath}/A1JVI6.bed'
rawData = pd.read_csv(filename, sep='\t', header=None)
merged = pd.read_csv('df_merged.csv', index_col=0)
#celltype = 'RNA.CM_Atrial'
celltype = 'RNA.CM_Myh7b...'
threshold = 0.25
outputName = 'PVdown_CM2_corr_25.csv'
corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName)

mypath = '/Users/minjunpark/Documents/Pitx2/manuscript/fig3/data/PVup_bicluster_outs'
filename = f'{mypath}/A1JVI6.bed'
rawData = pd.read_csv(filename, sep='\t', header=None)
merged = pd.read_csv('df_merged.csv', index_col=0)
#celltype = 'RNA.CM_Atrial'
celltype = 'RNA.CM_Myh7b...'
threshold = 0.25
outputName = 'PVup_CM2_corr_25.csv'
corr_matrix(mypath, filename, rawData, merged, celltype, threshold, outputName)
