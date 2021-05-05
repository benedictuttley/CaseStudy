import argparse
import pandas as pd
import scanpy as sc

'''
scRNA-seq analysis of single sample
'''

# Config variables
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = None
show_plots = False
root_dir = "C:/Users/bened/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/benedict/CaseStudy/data/united/RNA/"

# Entry point
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', help='Root dir name of sample')
    parser.add_argument('--show', help='Root dir name of sample', default=False)
    args = parser.parse_args()
    root_dir += args.sample
    results_file = (args.sample + ".h5ad")
    show_plots = args.show


# Read in count matrix
data = sc.read_10x_mtx(
    (root_dir + "/filtered_feature_bc_matrix"), # Directory containing feature-barcode .mtx file
    var_names='gene_symbols', # Use gene symbols for the variable names
    cache=True) # Create cache file for faster loading of file in the future
data.var_names_make_unique()
print(data)


# Preprocessing

# Show genes with highest fraction of counts across all cells
sc.pl.highest_expr_genes(data, n_top=20, save="highest_expressed_genes.png", show=show_plots)

# Broad filtration step
sc.pp.filter_cells(data, min_genes=200)
sc.pp.filter_genes(data, min_cells=3)
print("Cells with < 200 genes removed")
print("Genes in < 3 cells are removed")


# Compute QC metrics
data.var['mt'] = data.var_names.str.startswith('mt-') # Annotate the group of mitochondiral genes as 'mt'
sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# A violin plot of computed quality measures:
# Number of genes expressed in the count matrix
# Total counts per cell
# Percentage of counts in mitochondrial genes
sc.pl.violin(data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save='_qc.png', show=show_plots)
sc.pl.scatter(data, x='total_counts', y='pct_counts_mt', save='_total_counts_vs_pct_counts_mt.png', show=show_plots)
sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts', save='_total_counts_vs_n_genes_by_counts.png', show=show_plots)

print("Please review QC plots in QC folder \n")
min_n_genes_by_counts = int(input("Enter Min number of counts per gene: ") or 0)
max_n_genes_by_counts = int(input("Enter Max number of counts per gene: ") or float('Inf'))
max_pct_counts_mt = int(input("Enter maximum percentage of reads permitted to be mitochondrial: ") or 100)

data = data[(data.obs.n_genes_by_counts > min_n_genes_by_counts) & (data.obs.n_genes_by_counts < max_n_genes_by_counts), :]
data = data[data.obs.pct_counts_mt < max_pct_counts_mt, :]
print(data) #TODO: Print a nicer formatted summary


# Total-count normalise the data matrix to 10,000 reads per cell so that counts become comparable between cells
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data) # Logarithmise the data

# Identify highly variable genes
sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(data, save='_highly_variable_genes.png', show=show_plots)
data.raw = data

# Filter down to highly variable genes only
data = data[:, data.var.highly_variable]

# Regress out effects of total counts per cell and pct mitochondrial genes expressed
sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'])

# Scale each gene to unit variance
sc.pp.scale(data, max_value=10)

# Principal Components Analysis
# Reduce dimensionality of data to reveal main axes of variation and denoise the data
sc.tl.pca(data, svd_solver='arpack')

# Scatter plot in PCA coordiantes
sc.pl.pca(data, save='_pca.png', show=show_plots)
sc.pl.pca_variance_ratio(data, log=True, save='_pca_variance.png', show=show_plots) # See contribution of each PC to toal variance

# Save result
data.write(results_file)

# Computing Neighbourhood Graph
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
# Embed neighbourhood graph in 2 dimensions with UMAP
sc.tl.umap(data)
sc.pl.umap(data, color=["Sell", "Adam17"], save='_sell_adam17_umap.png', show=show_plots)


# Previous plot shows raw gene expression (normalised, logarithmised and uncorrected)
# Instead the scaled and corrected gene expression can be plotted
# sc.pl.umap(data, color=["Sell", "Adam17"], use_raw=False)

# Cluster the neighborhood graph
sc.tl.leiden(data)
# Plot clusters
sc.pl.umap(data, color=['leiden'], save="_umap_clusters.png", show=show_plots)
exit(0)

# Finding marker genes
sc.tl.rank_genes_groups(data, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(data, n_genes=25, sharey=False)

# Top 5 genes per cluster
pd.DataFrame(data.uns['rank_genes_groups']['names']).head(5)

# Comparison plots
sc.pl.violin(data, ['Sell'], groupby='leiden')
sc.pl.dotplot(data, ['Sell', 'Ccr7', 'Cd44', 'Fas', 'Ccl5'], groupby='leiden')