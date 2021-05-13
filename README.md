# R Single Cell RNA-seq Analysis
R Mardown using the Seurat scRNA-seq package to analyse scRNA-seq cellranger
data consisting of six stages.

1. Data Import
2. Quality Control
3. Dimensionality Reduction
4. Miscellaneous plots
5. Differential Expression (DE) analysis
6. Gene Enrichment Analysis

<b>IMPORTANT: </b> Program can run in single-sample or multi-sample mode and
can compare samples between eachother (multi-sample only) or compare clusters identified post dimensionality reduction.

Special attention should be payed to QC chunks (stage 2). It expects the user
to pass qc cutoffs based on the qc plots generated.

Output from differential gene expression testing and GO enrichment are written
to .xlsx (excel) workbooks. These tests are pairwise between two different
samples and are added to the workbook as a new sheet.

### Required
`sample_names`: List of filtered feature-barcode matrix names

### Settings
| Setting     | Description | Default    |
| ----------- | ----------- | -----------|
| remove_ribo     | Remove all ribosomal genes | False
| remove_mito   | Remove all mitochondiral genes       | False

<br>

# Python Single Cell RNA-seq Analysis

Interactive python script using the ScanPy scRNA-seq package to analyse scRNA-seq data in cellranger feature-barcode format (https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

The workflow includes a QC stage where QC plots are generated and the user is promted for qc cutoffs to be entered based on the plots. 

All intermediate plots will be written to a ```/figures``` folder.

Visualisation follows, with PCA dimensionality reduction and UMAP plots. Finally cells are clustered and differential expression testing carried out listing the marker genes for each cluster created.

## Getting Started
This is an example of running the analysis

### Install ScanPy

```python
pip install scanpy
```
As a minimum the script requires one or more samples to be passed with the ```--samples``` option. For more than one sample, seperate with commas. For example ```--samples A,B,C```

```unix
python scRNA-seq.py --samples <samples>
```

### Optional Arguments
| Option      | Description | Default    |
| ----------- | ----------- | -----------|
| --removeRibo      | Remove all ribosomal genes | False
| --show   | Show each plot upon generation        | False
