# Yeast RNA-seq Replicate Sensitivity Analysis

A reanalysis of the highly replicated yeast RNA-seq experiment published by *Schurch et al. (2016)*, using DESeq2

## Project aim

This project investigates how the number of biological replicates affects the recovery of differentially expressed genes.

The initial analysis useS the processed RNA-seq count matrix from the *Schurch et al*. dataset.

## Research questions

- What visualisations can be made from the orignal *Schurch et al. (2016)* dataset?
- How distinct are the overall gene-expression profiles of wild-type and knockout samples?
- How does differential-expression sensitivity change as the number of biological replicates increase?
- How variable are the results obtained during different random selections of samples?

## Dataset

- Organism: *Saccharomyces* *cerevisiae (Yeast)*
- Conditions: wild type and `snf2` knockout
- Study accession: `PRJEB5348`
- Original Publication: Schurch et al. (2016), *How many biological replicates are needed in an RNA-seq experiment?*

## The workflow

1. Inspected and reproduced the count matrix and sample metadata using the author's HighlyReplicatedRNASeq BioConductor package.
2. Performed exploratory analysis, visualisations (PCA, sample distance heatmap etc).
3. Ran a full-data DESeq2 analysis.
4. Randomly subsampled different numbers of replicates. Tested out n=3 before n=3-10.
5. Calculated DEG recovery relative to the full-data analysis.
6. Plotted sensitivity against replicate number.

## Repository structure

```text
yeast-rnaseq-replicate-sensitivity/
├── data/
│   ├── metadata/
│   └── processed/
├── results/
│   ├── figures/
│   │   ├── full_data_qc/
│   │   └── subsampling/
│   ├── full_data/
│   └── subsampling/
│       ├── n3_single_run/
│       └── repeated_runs/
└── scripts/
    └── R/
```
