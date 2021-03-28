import sys

# STEP 2: Join the two data (fimo output + count matrix)
## first generate LAcm_atac_count_ver2.csv and PVcm_atac_count_ver2.csv
def merge_fimo_count(fimo_input, count_file):
    fimo = pd.read_csv(fimo_input, sep='\t', header=None)
    count = pd.read_csv(count_file, sep='\t', header=None)

    # merge on fimo output
    df = pd.merge(fimo, count, how='left', left_on=[0,1,2], right_on=[0,1,2])
    return df

# STEP 3: Calculate influence score for each gene
def influence_score(merge_df, genes_bed, output_name):
    genes = pd.read_csv(genes_bed, sep='\t', header=None)
    #print(genes)
    dff = pd.DataFrame()
    gene_scores = []
    gene_names = []
    gene_score = []
    accessibility_scores = []

    ## first get influence scores for all peaks
    for j in merge_df.index:
        df_chr = merge_df.iloc[j,0]
        df_start = merge_df.iloc[j,1]
        df_end = merge_df.iloc[j,2]

        peak_length = df_end - df_start
        pitx2_occurrence = merge_df.iloc[j,3]
        accessibility_score = merge_df.iloc[j,4]

        ## not divide by accessibility score
        influence_score = accessibility_score * pitx2_occurrence / peak_length

        ## not divide by peak_length
        #influence_score = accessibility_score * pitx2_occurrence

        #print('accessibility score is {}'.format(accessibility_score))
        #print('Pitx2 occurrence is {}'.format(pitx2_occurrence))
        #print('Peak length is {}'.format(peak_length))
        #print('influence score is {}'.format(influence_score))
        gene_score.append(influence_score)

    merge_df['influence score'] = gene_score

    for i in genes.index:
        """
        ## Take a unique values from the gene list (Use gene once)
        """
        gene_info = genes.loc[i]
        # get gene information
        gene = gene_info[0]
        chr = gene_info[1]
        start = gene_info[2]
        end = gene_info[3]

        # get the subset of interest
        df = merge_df[(merge_df[0] == chr) & (merge_df[1] > start) & (merge_df[2] < end)]
        ##try 2: maybe try without peak_length
        ##try 1: try without dividing by accessibility_scores
        sum_influence_score = np.sum(df.loc[:,'x score'].values)

        gene_names.append(gene)
        gene_scores.append(sum_influence_score)
        print(i)

    dff['gene name'] = gene_names
    dff['influence score'] = gene_scores
    dff.to_csv(output_name)

def main(argv):
	if argv is None:
		argv = sys.argv
	else:
		argv = argv.split()

	input1 = argv[1] # e.g. 'processed_ETS1_fimo.tsv'
    input2 = argv[2] # e.g. 'cardiomyocytes_peaks.bed'

    input3 = argv[3] # e.g. 'AllGenes_100kb_flanking_regions_dropped_duplicates.bed'
    input4 = argv[4] # e.g. 'ETS1_influence_score.csv'
    merged_df = merge_fimo_count(input1, input2)
    influence_score(merged_df, input3, input4)
	return 0

if __name__ == "__main__":
	sys.exit(main(None))
