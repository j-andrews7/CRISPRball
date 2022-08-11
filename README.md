# CRISPRball - Lift the Fog from Your CRISPR Analyses

CRISPR screens are becoming more and more common, and as such, so are analysis tools. Perhaps the most popular is [MAGeCK](https://sourceforge.net/projects/mageck/), which is both quite easy to run and interpret.

**CRISPRball** is a Shiny application to explore, visualize, filter, and integrate CRISPR screens with public data and multiple datasets. In particular, it allows for publication-quality figure generation including full aesthetic customization and interactive labeling, filtering of results using DepMap genes, simple comparisons between datasets/timepoints/treatments, etc.

It is designed for end users and may be particularly useful for bioinformatics/genome editing cores that perform MAGeCK analyses before returning results to users. Pointing users to the online version of the app (or a hosted one) will allow them to quickly wade through and interpret their data.

## Installation

This package is in development - it may break at any time and contain unstable or untested features. A stable version will be submitted to Bioconductor once the initially planned features are completed.

To install the package via Github:

```r
install.packages("devtools")
devtools::install_github("j-andrews7/CRISPRball")
```

## Usage

Starting the app is as simple as calling the `CRISPRball` function, which will return a shiny app.

```r
library("CRISPRball")
CRISPRball()
```

Users can upload their data within the app very easily.

One can also pass their input data directly as input - all that are needed are file paths to the MAGeCK RRA output files. In this case, we'll use the example data from the [third MAGeCK tutorial](https://sourceforge.net/p/mageck/wiki/demo/#the-third-tutorial-going-through-a-public-crisprcas9-screening-dataset). In this example, the two datasets are just reverse comparisons (ESC vs plasmid & plasmid vs ESC).

```r
# Create lists of results summaries for each dataset.
d1.genes <-  read.delim(system.file("extdata", "esc1.gene_summary.txt", package = "CRISPRball"), check.names = FALSE)
d2.genes <-  read.delim(system.file("extdata", "plasmid.gene_summary.txt", package = "CRISPRball"), check.names = FALSE)

d1.sgrnas <-  read.delim(system.file("extdata", "esc1.sgrna_summary.txt", package = "CRISPRball"), check.names = FALSE)
d2.sgrnas <-  read.delim(system.file("extdata", "plasmid.sgrna_summary.txt", package = "CRISPRball"), check.names = FALSE)

count.summ <- read.delim(system.file("extdata", "escneg.countsummary.txt", package = "CRISPRball"), check.names = FALSE)
norm.counts <- read.delim(system.file("extdata", "escneg.count_normalized.txt", package = "CRISPRball"), check.names = FALSE)

genes <- list(ESC = d1.genes, plasmid = d2.genes)
sgrnas <- list(ESC = d1.sgrnas, plasmid = d2.sgrnas)

CRISPRball(gene.data = genes, sgrna.data = sgrnas, count.summary = count.summ, norm.counts = norm.counts)
```

Passing data directly can be useful when hosting the app on a local Shiny server where having pre-loaded data for the user is wanted. 

### Adding Genesets

Often, it can be useful to highlight a set of genes on the plots. This can be done by passing a named list of gene identifiers to the `genesets` argument.

```r
library("msigdbr")

# Retrieve MSigDB Hallmark gene sets and convert to a named list.
gene.sets <- msigdbr(species = "Homo sapiens", category = "H")
gene.sets <- hs.gene.sets %>% split(x = .$gene_symbol, f = .$gs_name)

# Can also add genesets manually.
gene.sets["my_fav_genes"] <- c("TOP2A", "FECH", "SOX2", "DUT", "RELA")

CRISPRball(gene.data = genes, sgrna.data = sgrnas, count.summary = count.summ, norm.counts = norm.counts, genesets = gene.sets)
```

Such genesets can then be highlighted very easily on the plots in the **Gene (Overview)** tab using the "Highlight Gene(sets)" inputs in the sidebar.

### Incorporating DepMap Data

[DepMap](https://depmap.org/portal/) contains a multitude of data for hundreds of cell lines, including CRISPR/RNAi screen results, gene expression data, copy number variation, and more. This data can be extremely useful to remove common dependencies from your own results, look at expression or copy number for a given gene across various lineages or diseases, and to identify hits that are selective to a given (sub)lineage.

For fast access to this data, **CRISPRball** includes a function (`build_depmap_db()`) to build a SQLite database using the [depmap](https://bioconductor.org/packages/release/data/experiment/html/depmap.html) R package. **This database will be large, >4 GB.**

This SQLite database can then be passed to the app and the data contained therein easily explored in the **DepMap** tab.

```r
library("depmap")
library("pool")
library("RSQLite")

# This will likely take several minutes to run.
# The database will be named "depmap_db.sqlite" and place in the working directory.
build_depmap_db()

CRISPRball(gene.data = genes, sgrna.data = sgrnas, count.summary = count.summ, 
           norm.counts = norm.counts, genesets = gene.sets, depmap.db = "depmap_db.sqlite")
```

### Additional Help

Almost every input in the app will display a helpful tooltip explaining its function on hover. Plots also have an information icon that will explain the plot contents when hovered or clicked.
