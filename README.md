## Single Cell RNA-seq Analysis

Interactive python script using the ScanPy scRNA-seq package to analyse scRNA-seq data in cellranger feature-barcode format (https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

The workflow includes a QC stage where QC plots are generated and the user is promted for qc cutoffs to be entered based on the plots. 

All intermediate plots will be written to a ```/figures``` folder.

Visualisation follows, with PCA dimensionality reduction and UMAP plots. Finally cells are clustered and differential expression testing carried out listing the marker genes for each cluster created.

## Getting Started
This is an example of running the analysis

## Install ScanPy

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
| ```--removeRibo```      | Remove all ribosomal genes | False
| ```--show```   | Show each plot upon generation        | False
