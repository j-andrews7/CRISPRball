#!/bin/bash
# The MAGeCK RRA tutorial #3 was followed almost exactly:
# https://sourceforge.net/p/mageck/wiki/demo/#the-third-tutorial-going-through-a-public-crisprcas9-screening-dataset

# It is roughly outlined below:
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR376/ERR376998/ERR376998.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR376/ERR376999/ERR376999.fastq.gz

gunzip ERR376998.fastq.gz
gunzip ERR376999.fastq.gz

# Download the library file and unzip it from:
# https://sourceforge.net/projects/mageck/files/libraries/yusa_library.csv.zip/download

# sgRNAs for 1000 randomly selected genes were selected.
# Extract unique genes
awk -F',' 'NR>1 {print $3}' yusa_library.csv | sort -u > unique_genes.txt

# Randomly select 1000 genes
shuf -n 1000 unique_genes.txt > random_1000_genes.txt

# Add header
touch yusa_library_1000.csv
echo "id,gRNA.sequence,gene" > yusa_library_1000.csv

# Extract lines for the 1000 randomly selected genes
awk -F',' 'NR==FNR {genes[$1]; next} $3 in genes' random_1000_genes.txt yusa_library.csv >> yusa_library_1000.csv


# Run MAGeCK count
mageck count -l yusa_library_1000.csv -n demo -n escneg --sample-label "plasmid,ESC1" --fastq ERR376998.fastq  ERR376999.fastq

# Run the RRA tests
mageck test -k escneg.count.txt -t ESC1 -c plasmid -n esc1
mageck test -k escneg.count.txt -t plasmid -c ESC1 -n plasmid

# After running the tests, the "neg|fdr" column for DDX27 the "esc1.gene_summary.txt" file was manually set to 1 as a sanity check.