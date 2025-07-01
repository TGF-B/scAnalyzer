import scanpy as sc
import matplotlib.pyplot as plt

def plot_violin(adata, genes, groupby, output_path=None):
    import pandas as pd
    plt.close('all')
    
    # 确保基因存在于数据集中
    existing_genes = [gene for gene in genes if gene in adata.var_names]
    if not existing_genes:
        available_genes = list(adata.var_names[:10])  # 列出前10个基因作为示例
        return {"error": f"请求的基因不在数据集中。可用基因示例: {available_genes}"}
    
    # 检查分组变量是否存在
    if groupby not in adata.obs.columns:
        available_groups = list(adata.obs.columns)
        return {"error": f"分组变量 '{groupby}' 不在数据集中。可用分组: {available_groups}"}
    
    # 自动将连续型分组变量转为分类型
    if not pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        # 如果是二元对照组（如Ctrl vs CP），直接转换为分类
        if len(adata.obs[groupby].unique()) <= 5:
            adata.obs[groupby] = adata.obs[groupby].astype("category")
        else:
            # 否则进行分箱
            try:
                adata.obs[f"{groupby}_binned"] = pd.qcut(
                    adata.obs[groupby], q=5, duplicates='drop').astype("category")
                groupby = f"{groupby}_binned"
            except:
                # 如果无法分箱，直接转为分类
                adata.obs[groupby] = adata.obs[groupby].astype("category")
    
    # 使用更好的图形设置
    sc.pl.violin(
        adata, 
        keys=existing_genes, 
        groupby=groupby, 
        jitter=0.4, 
        rotation=45, 
        show=False,
        multi_panel=True if len(existing_genes) > 3 else False  # 多个基因时分面展示
    )
    
    if output_path:
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
    
    return {"success": f"已绘制 {len(existing_genes)} 个基因的小提琴图，按 {groupby} 分组"}
def plot_qc(adata, output_path=None):
    plt.close('all')
    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

def plot_umap(adata, color=None, output_path=None):
    plt.close('all')
    
    # 检查是否已计算 UMAP
    if 'X_umap' not in adata.obsm:
        # 检查是否已有 PCA 结果
        if 'X_pca' not in adata.obsm:
            sc.pp.pca(adata)
        # 计算邻居图
        sc.pp.neighbors(adata)
        # 使用标准 scanpy UMAP 实现
        sc.tl.umap(adata)
    
    # 检查 cell_type 列，优先使用带下划线的格式
    if color is None:
        if 'cell_type' in adata.obs.columns:
            color = 'cell_type'
        elif 'celltype' in adata.obs.columns:
            color = 'celltype'
    
    # 用合适的参数绘图
    if color in adata.obs.columns and adata.obs[color].dtype.name == 'category':
        # 如果是分类变量，使用 'on data' 图例
        sc.pl.umap(adata, color=color, legend_loc='on data', 
                  frameon=False, title=f'UMAP by {color}', show=False)
    else:
        # 其他情况使用默认设置
        sc.pl.umap(adata, color=color, show=False)
    
    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

def check_adata_structure(adata, output_path=None):
    """
    返回 AnnData 对象的详细结构信息，以 JSON 格式展示
    """
    import json
    
    # 构建包含 AnnData 结构的字典
    structure = {
        "shape": {"n_obs": adata.n_obs, "n_vars": adata.n_vars},
        "obs_keys": list(adata.obs.columns),
        "var_keys": list(adata.var.columns),
        "obsm_keys": list(adata.obsm.keys()) if hasattr(adata, 'obsm') else [],
        "layers_keys": list(adata.layers.keys()) if hasattr(adata, 'layers') else [],
        "uns_keys": list(adata.uns.keys()) if hasattr(adata, 'uns') else [],
        "obsp_keys": list(adata.obsp.keys()) if hasattr(adata, 'obsp') else [],
        "varm_keys": list(adata.varm.keys()) if hasattr(adata, 'varm') else [],
        "varp_keys": list(adata.varp.keys()) if hasattr(adata, 'varp') else []
    }
    
    # 添加 obs 中重要列的类型和值的统计信息
    structure["obs_sample"] = {}
    for col in adata.obs.columns:
        col_type = str(adata.obs[col].dtype)
        if pd.api.types.is_categorical_dtype(adata.obs[col]):
            # 对于分类变量，列出类别
            structure["obs_sample"][col] = {
                "type": col_type,
                "categories": list(adata.obs[col].cat.categories),
                "n_categories": len(adata.obs[col].cat.categories)
            }
        else:
            # 对于数值变量，给出统计摘要
            try:
                if pd.api.types.is_numeric_dtype(adata.obs[col]):
                    structure["obs_sample"][col] = {
                        "type": col_type,
                        "min": float(adata.obs[col].min()),
                        "max": float(adata.obs[col].max()),
                        "mean": float(adata.obs[col].mean())
                    }
                else:
                    # 非数值非分类，列出唯一值（最多10个）
                    unique_values = list(adata.obs[col].unique())
                    structure["obs_sample"][col] = {
                        "type": col_type,
                        "unique_values": unique_values[:10],
                        "n_unique": len(unique_values)
                    }
            except:
                structure["obs_sample"][col] = {"type": col_type}
    
    # 转换为 JSON
    json_result = json.dumps(structure, indent=2)
    
    # 如果需要输出到文件
    if output_path:
        with open(output_path.replace('.png', '.json'), 'w') as f:
            f.write(json_result)
    
    return json_result

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
