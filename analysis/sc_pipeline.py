# this pipeline assumes you already have scRNA-seq data stored in an Anndata object with cell IDs as your observations (rows) and genes are your variables (columns). 

import scanpy as sc

# adjust filtering criteria as needed
def filter_cells_and_genes(adata, min_genes=200, min_cells=20):
    # Filter out cells that have fewer than "min_genes" detected genes (ex, 200)
    sc.pp.filter_cells(adata, min_genes=min_genes) 
    # Filter out genes that appear in fewer than "min_cells" cells (ex, 20)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def normalize_and_log_transform(adata, target_sum=1e4):
    # Apply global-scaling normalization
    sc.pp.normalize_total(adata, target_sum=target_sum)
    # Log-transform the data
    sc.pp.log1p(adata)
    return adata

# adjust # of highly variable genes to retain
def find_highly_variable_genes(adata, n_top_genes=2000):
    # Find the most highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)
    return adata

def scale_data(adata, zero_center=True):
    # Apply z-transformation
    sc.pp.scale(adata, zero_center=zero_center)
    return adata

def perform_pca(adata, svd_solver='arpack'):
    # Perform dimensionality reduction via PCA
    sc.tl.pca(adata, svd_solver=svd_solver)
    return adata

# adjust number of nearest neighbors and principal components to retain
def construct_nearest_neighbors_graph(adata, n_neighbors=10, n_pcs=10):
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata

# adjust cluster resolution and number of iterations as needed
def apply_leiden_clustering(adata, resolution=0.5, n_iterations=3, flavor='igraph', directed=False):
    sc.tl.leiden(adata, key_added='clusters', resolution=resolution, n_iterations=n_iterations, flavor=flavor, directed=directed)
    return adata

def create_and_visualize_umap(adata, color='clusters'):
    sc.tl.umap(adata)
        sc.pl.umap(adata, color=color, add_outline=True, legend_loc='on data', legend_fontsize=12, legend_fontoutline=2, frameon=True)

def process_single_cell_data(adata):
    # Step 1: Filter cells and genes
    adata = filter_cells_and_genes(adata)
    print(f"Number of Genes: {adata.n_vars}")
    print(f"Number of cells: {adata.n_obs}")

    # Step 2: Normalize and log-transform
    adata = normalize_and_log_transform(adata)
    
    # Step 3: Find highly variable genes
    adata = find_highly_variable_genes(adata)
    print(adata)
    
    # Step 4: Scale the data
    adata = scale_data(adata)
    
    # Step 5: Perform PCA
    adata = perform_pca(adata)
    
    # Step 6: Construct nearest neighbors graph
    adata = construct_nearest_neighbors_graph(adata)
    
    # Step 7: Apply Leiden clustering
    adata = apply_leiden_clustering(adata)
    
    # Step 8: Create and visualize UMAP
    create_and_visualize_umap(adata)

# Example usage with adata_transposed (AnnData object name in analysis notebook)
process_single_cell_data(adata_transposed)
