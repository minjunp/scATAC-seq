# !/bin/bash

# Script description:
# Gets the fasta files output given an input BED file

# Load config
source script_config.txt

# Sanity check for variables
echo "Input file list: "$Blist
echo "Genome file: "$Genome

# Script body
for Bedin in $Blist;
do
Bedout=$(echo $Bedin | cut -d '.' -f1)"_v2.fa"
bedtools getfasta -fi $Genome -bed $Bedin |\
    sed 's/>c/>::c/' > $Bedout
done

echo "All input files processed."
