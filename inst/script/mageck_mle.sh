#!/bin/bash
# The MAGeCK RRA tutorial #4 was followed almost exactly:
# https://sourceforge.net/p/mageck/wiki/demo/#the-fourth-tutorial-using-mageck-mle-module

# It is roughly outlined below.

# Download the counts table file (leukemia.new.csv) from:
# http://https//sourceforge.net/projects/mageck/files/example/leukemia.new.csv/download

# Counts for 1000 randomly selected genes were selected.
# Extract unique genes
awk -F',' 'NR>1 {print $2}' leukemia.new.csv | sort -u > unique_genes.txt

# Randomly select 1000 genes (so package data is not too large)
shuf -n 1000 unique_genes.txt > random_1000_genes.txt

# Add header
touch leukemia.new_1000.csv
echo "sgRNA,Gene,HL60.initial,KBM7.initial,HL60.final,KBM7.final" > leukemia.new_1000.csv

# Extract lines for the 1000 randomly selected genes
awk -F',' 'NR==FNR {genes[$1]; next} $2 in genes' random_1000_genes.txt leukemia.new.csv >> leukemia.new_1000.csv

# Create the design matrix file.
echo "Samples        baseline        HL60        KBM7" > designmat.txt
echo "HL60.initial   1               0           0" >> designmat.txt
echo "KBM7.initial   1               0           0" >> designmat.txt
echo "HL60.final     1               1           0" >> designmat.txt
echo "KBM7.final     1               0           1" >> designmat.txt

# Run MAGeCK mle
mageck mle -k leukemia.new_1000.csv -d designmat.txt -n beta_leukemia --threads 4 --max-sgrnapergene-permutation 15

# The beta_leukemia.gene_summary.txt, beta_leukemia.sgrna_summary.txt, and leukemia.new_1000.csv files were then included in the package.