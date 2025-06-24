import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

def plot_violin(adata, genes, groupby, output_path=None):

    plt.close('all')
    # 自动将 groupby 列转为类别型（如不是的话）
    if groupby in adata.obs and not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        # 可以根据实际需求分箱，比如分10组
        adata.obs[groupby] = pd.qcut(adata.obs[groupby], q=10, duplicates='drop').astype("category")
    sc.pl.violin(adata, keys=genes, groupby=groupby, jitter=0.4, rotation=45, show=False)
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

def plot_qc(adata, output_path=None):
    plt.close('all')
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

def plot_umap(adata, color=None, output_path=None):
    import umap
    import inspect
    plt.close('all')
    sc.pp.neighbors(adata)
    umap_params = {}
    if "backend" in inspect.signature(umap.UMAP).parameters:
        try:
            import cuml
            umap_params["backend"] = "cuml"
        except ImportError:
            umap_params["backend"] = "sklearn"
    reducer = umap.UMAP(**umap_params)
    adata.obsm['X_umap'] = reducer.fit_transform(adata.X)
    sc.pl.umap(adata, color=color, show=False)
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

def check_adata_structure(adata, output_path=None):
    info = str(adata)
    if output_path:
        with open(output_path.replace('.png', '.txt'), 'w') as f:
            f.write(info)
    return info

def plot_dotplot(adata, genes, groupby, output_path=None):
    plt.close('all')
    sc.pl.dotplot(adata, var_names=genes, groupby=groupby, show=False)
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

def plot_pseudotime(adata, color=None, output_path=None):
    plt.close('all')
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.dpt(adata)
    sc.pl.pca(adata, color='dpt_pseudotime', show=False)
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()
    