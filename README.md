## Usage

Install conda environment

```sh
conda env create -f environment.yml
```

Activate conda environment

```sh
conda activate museam
```

Run model on your data

```sh
python train.py sequences.fa wt_readout.dat parameters.txt
```

## Demo




# Normalized Motif Score
### atac_count.R, csv_to_bed.sh, getfasta_script.sh
* STEP 1: Get ATAC counts, Prepare DE genes list
* STEP 2: ATAC counts (.csv) --> BED format --> FASTA format

### run_fimo.sh, collect_outputs.sh, process_fimo_output.sh
* STEP 3-a: Apply FIMO to find motif hits using cisBP motif library (remove first and last 3 lines of the output for the next step!)
* STEP 3-b: save flanking regions of DE genes in BED format

### calclulate_influence_score.py
* STEP 4-a: Calculate influence score for ALL genes. Identify the motifs of interest and get the processed FIMO output. We can use them to calculate motif scores.
* STEP 4-b: Split into up-regulated, down-regulated, and non-DE genes (adj p-val < 0.01 & abs(log2FC) > 1.25)

### viz.py
* STEP 5: Visualize using a box plots

# Syntax Analysis
### data_preprocess.sh
* STEP 1: Get all the peak information for each TFs from fimo outputs (all .tsv files together)

### preprocess.sh
* STEP 2:





################## Application Note Paper #######################
Input: ATAC counts, DE genes list
