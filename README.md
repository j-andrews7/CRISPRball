# CRISPRball - Lift the Fog from Your CRISPR Analyses

CRISPR screens are becoming more and more common, and as such, so are analysis tools. Perhaps the most popular is [MAGeCK](https://sourceforge.net/projects/mageck/), which is both quite easy to run and interpret.

**CRISPRball** is a Shiny application to explore, visualize, filter, and integrate CRISPR screens with public data and multiple datasets. In particular, it allows for publication-quality figure generation including full aesthetic customization and interactive labeling, filtering of results using DepMap genes, simple comparisons between datasets/timepoints/treatments, etc.

It is designed for end users and may be particularly useful for bioinformatics/genome editing cores that perform MAGeCK analyses before returning results to users. Pointing users to the online version of the app (or a hosted one) will allow them to quickly wade through and interpret their data.

## Installation

This package is in development - it may break at any time and contain unstable or untested features. A stable version will be submitted to Bioconductor once the initially planned features are completed.

To install the package via Github:

```
install.packages("devtools")
devtools::install_github("j-andrews7/CRISPRball")
```

## Usage

**CRISPRball** is simple to use after running a typical MAGeCK RRA analysis.

### Quick Start


```
library("CRISPRball")

# Create lists of results summaries for each dataset.
d1.genes <-  read.delim("d1.gene_summary.txt", check.names = FALSE)
d2.genes <-  read.delim("d2.gene_summary.txt", check.names = FALSE)

d3.genes <-  read.delim("d3.gene_summary.txt", check.names = FALSE)
d4.genes <-  read.delim("d4.gene_summary.txt", check.names = FALSE)

d1.grnas <-  read.delim("d1.sgrna_summary.txt", check.names = FALSE)
d2.grnas <-  read.delim("d2.sgrna_summary.txt", check.names = FALSE)

d3.grnas <-  read.delim("d3.sgrna_summary.txt", check.names = FALSE)
d4.grnas <-  read.delim("d4.sgrna_summary.txt", check.names = FALSE)

count_summ <- read.table("counts.countsummary.txt", sep = "\t", header = TRUE)
norm_counts <- read.table("counts.count_normalized.txt", sep = "\t", header = TRUE)

genes <- list(D1 = d1.genes, D2 = d2.genes, D3 = d3.genes, D4 = d4.genes)
grnas <- list(D1 = d1.grnas, D2 = d2.grnas, D3 = d3.grnas, D4 = d4.grnas)

CRISPRball(gene.data = genes, sgrna.data = grnas, count.summary = count_summ, norm.counts = norm_counts)
```



