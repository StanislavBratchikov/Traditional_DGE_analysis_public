# Single cell RNAseq analysis using traditional bioinformatics tools
This repo highlights parts of analysis perfomed as part of larger project. The project aims to answer whether pathogen response is predictive of pneuomonia patients outcomes.

## Methods:

To find differentially-expressed genes between pathogen groups (7) and control groups (2), we use pseudobulk approach with certain modification, detailed below. Pseudobulk approach shows better agreement with bulk approach during benchmarking (Squair et al., Nat Commun, 2021) compared to tests that do not account for cells’ belonging to different samples, and does not suffer from inflated p-values.

First, we prepared pseudobulk samples from our single-cell data. For each annotated cell type *T*, we excluded samples that didn’t have 50 cells. Each sample belongs to strictly one pathogen or control group *G*. Next, we excluded genes that are not consistently expressed in the data from differential gene expression analysis: if a gene was detected (count > 0) in less than 80% of pseudobulk samples after filtering from sample group *G*, we excluded that gene. We also excluded

- mitochondrial genes
- ribosomal genes
- genes with certain name patterns (see [here](https://www.biostars.org/p/9553891/))

## NB:

is it right to exclude these genes before running DESeq2? It increases power but affects counts of the pseudobulks. Maybe a more correct procedure would be to:

- Sum all detected genes for pseudobulks
- Construct DESeq2 object with normalization
- Run DE tests
- Remove genes from test results
- Redo FDR correction on the set of selected genes


Then we summed up raw gene counts for all cells in that sample that belong to the annotated cell type  (see here). Next, we constructed DESeq2 object will all pathogen and control groups (9 total) together with model expression *'~ group + sex'* and *'fitType'* set to *'local'*. Finally, we performed pairwise differential gene expression testing between each pathogen group and each control group (14 tests maximum) where we found at least 3 samples on both sides. We discarded genes that had baseMean expression of less than 50. We used adjusted p-values from DESeq2 to determine significant up- or down-regulated genes at 0.05 alpha threshold. Code is in this repo.
