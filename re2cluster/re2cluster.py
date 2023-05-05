#!/usr/bin/env python 

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

from typing import List, Tuple, Iterable, Literal


def _split_adata_on_identity(adata, identity_column: str) -> List[anndata.AnnData]:
    """ Subset an anndata.AnnData object based on identity """

    unique_identities = pd.unique(adata.obs[identity_column])
    adatas = [adata[adata.obs[identity_column] == identity] for identity in unique_identities]
    return adatas


def quality_control(adata):
    """ Calculate qc parameters """

    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])[0]
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=['mt'],
                               percent_top=None,
                               log1p=False,
                               inplace=True
                               )
    
    return adata


def plot_quality_control(adata: anndata.AnnData, save: str = None):
    """ Plot quality control parameters for initial object """

    p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes',
                       color='pct_counts_mt'
                       )
    
    p2 = sc.pl.scatter(adata[adata.obs['n_counts']<5000], 'n_counts', 'n_genes',
                       color='pct_counts_mt'
                       )
    
    t1 = sc.pl.violin(adata, 
                      keys=['n_counts','n_genes','pct_counts_mt'],
                      size=2, log=True, cut=0
                      )
    
    if save is not None: 
        plt.savefig(save, bbox_inches='tight')

    return fig, axs


def normalize_hvg_pearson(adata: anndata.AnnData, 
                          min_cells_per_gene: int,
                          n_hvg: int
                          ) -> anndata.AnnData:
    """ Normalize with pearson residuals and find highly variable genes (hvg) 
    
    Parameters
    ----------
    adata : anndata.Anndata 
        AnnData object 
    min_cells_per_gene : int 
        Passed to `sc.pp.filter_genes` `min_counts`. Remove genes that measured
        less cells
    n_hvg : int 
        Number of highly variable genes that are considered. Passed to 
        `sc.experimental.pp.highly_variable_genes` `n_top_genes` parameter

    Returns 
    -------
    anndata.AnnData 
        With highly variable genes

    """

    adata.raw = adata 
    sc.experimental.pp.normalize_pearson_residuals(adata)
    sc.pp.filter_genes(adata, min_counts=min_cells_per_gene)
    sc.experimental.pp.highly_variable_genes(adata, 
                                             flavor='pearson_residuals', 
                                             n_top_genes=n_hvg
                                             )
    
    sc.experimental.pp.normalize_pearson_residuals(adata)

    return adata 


def normalize_hvg_tpm(adata: anndata.AnnData, 
                        n_hvg: int) -> anndata.AnnData: 
    """ Normalize to transcripts per million, log1p normalize and find 
    highly variable genes (hvg) 
    """

    # normalize data with TPM + log1p
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=10**4)
    sc.pp.log1p(adata)

    # get highly variable genes 
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)

    return adata


def select_ideal_pcs(adata: anndata.AnnData): 
    """ Select ideal number of principal components of sc data 
    
    Current criterium: Change in explained variance is in the lowest 15 % quantile

    Parameters
    ----------
    adata : anndata.Anndata 
        anndata object. Is expected to contain information about pca
    """
    variance_ratio = adata.uns['pca']['variance_ratio']
    variance_ratio_diff = -1*np.diff(variance_ratio)

    # Number of PCs that fullfill criterium
    n_pcs = (variance_ratio_diff > np.quantile(variance_ratio_diff, 0.85)).sum()

    return n_pcs


def silhouette_param_scan(adata: anndata.AnnData,
                          leiden_resolution_min: float,
                          leiden_resolution_max: float,
                          steps: int) -> Tuple[Iterable, Iterable, Iterable]:
    """ Optimize resolution of leiden clustering by maximizing
        silhouette score.

        Parameters
        ----------
        adata : anndata.Anndata 
            Anndata object 
        leiden_resolution_min : float 
            Minimal resolution parameter to test 
        leiden_resolution_max : float 
            Maximal resolution parameter to test 
        steps : int 
            Number of steps for which clustering should be tested

        Returns 
        -------
        leiden_resolution_range : np.array
            All resolution values that were tested
        silhouette_scores : np.array 
            Silhouette scores of clustering in the range of -1 
            (incorrect clustering) to 1 (perfect assignment). 
            Is np.nan if `sklearn.metrics.silhouette_score` failes
        n_clusters : np.array 
            Number of clusters obtained for the respective resolution.  

    """
    # iterate over these resolutions 
    leiden_resolution_range = np.linspace(leiden_resolution_min, 
                                          leiden_resolution_max, 
                                          steps)
    
    # initialize result arrays 
    silhouette_scores = np.zeros(leiden_resolution_range.shape)
    n_clusters = np.zeros(leiden_resolution_range.shape)

    for idx, resolution in enumerate(leiden_resolution_range):
        adata_tmp = sc.tl.leiden(adata,
                                resolution=resolution,
                                key_added='leiden',
                                copy=True
                                )
        
        n_clusters[idx] = adata_tmp.obs['leiden'].unique().size

        try:
            silhouette_scores[idx] = silhouette_score(X=adata_tmp.obsm['X_pca'],
                                                      labels=adata_tmp.obs['leiden']
                                                      )
        # silhouette score might fail if ther is only one cluster 
        except:
            silhouette_scores[idx] = np.nan

        del adata_tmp

    return leiden_resolution_range, silhouette_scores, n_clusters


def plot_silhouette_param_scan(leiden_resolution_range: Iterable, 
                               silhouette_scores: np.array, 
                               n_clusters: np.array, 
                               save: str = None):
    """ Plot performance of leiden clustering for different resolutions 

        Parameters
        ----------
        leiden_resolution_range : np.array 
            Range of leiden parameters 
        silhouette_scores : np.array 
            Performance measure. Same indexing as `leiden_resolution_range` 
        n_clusters : np.array 
            Number of clusters. Same indexing as `leiden_resolution_range` 
        save : None|str 
            If string, save plot at indicated path. Default is None, which 
            does not save

        Returns
        -------
        fig, axs
    """

    ideal_idx = np.argmax(silhouette_scores)

    plot_params = dict(
        color='#13688D',
        linestyle='',
        marker='o',
        markersize=5
    )

    line_params = dict(
        color='#E74C3C',
        linewidth=2
    )

    # Initialize plot 
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    axs[0].axvline(leiden_resolution_range[ideal_idx], **line_params)
    axs[0].plot(leiden_resolution_range, silhouette_scores, **plot_params)
    axs[0].set_xlabel('Leiden Resolution')
    axs[0].set_ylabel('Silhouette Scores')
    axs[0].set_title('Resolution vs. Average Silhouette')

    axs[1].plot(leiden_resolution_range, n_clusters, **plot_params)
    axs[1].axvline(leiden_resolution_range[ideal_idx], **line_params)
    axs[1].set_xlabel('Leiden Resolution')
    axs[1].set_ylabel('#Clusters')
    axs[1].set_title('Resolution vs. Number of clusters')

    axs[2].plot(n_clusters, silhouette_scores, **plot_params)
    axs[2].axvline(n_clusters[ideal_idx], **line_params)
    axs[2].set_xlabel('#Clusters')
    axs[2].set_ylabel('Silhouette Scores')
    axs[2].set_title('Number of clusters vs. Silhouette Scores')

    plt.tight_layout()

    if save is not None:
        plt.savefig(save, bbox_inches='tight')

    return fig, axs

#%%

# if __name__ == '__main__':

# Parameters 
adata = sc.datasets.pbmc3k()
leiden_resolution_min = 0.01
leiden_resolution_max = 1.4
steps = 30
normalization_method = 'tpm'
n_hvg = 2000 
min_cells_per_gene = 3
save_qc = None
save_deg = None
save_param_scan = None

# QC 
adata = quality_control(adata)

# ## Plotting 
fig, axs = plot_quality_control(adata, save=save_qc)

if normalization_method == 'pearson':
    adata = normalize_hvg_pearson(adata,
                                  min_cells_per_gene=min_cells_per_gene, 
                                  n_hvg=n_hvg
                                  )
elif normalization_method == 'tpm':
    adata = normalize_hvg_tpm(adata, n_hvg=n_hvg)
else:
    raise ValueError(f'Normalization method {normalization_method} not available. Select "pearson" or "tpm"')

# PCA
sc.pp.pca(adata, n_comps=50, use_highly_variable= True, svd_solver='arpack')
n_pcs = select_ideal_pcs(adata)

# Neighbors
# Select n_neighbors parameter for neighbor graph generation
# at most 30
n_neighbors = int(min(0.5*np.sqrt(adata.n_obs), 30))

sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

# compute optimal clustering 
leiden_resolution_range, silhouette_scores, n_clusters = silhouette_param_scan(adata, 
                                                                               leiden_resolution_min, 
                                                                               leiden_resolution_max, 
                                                                               steps
                                                                               )                                                                             )
 
optimal_resolution = leiden_resolution_range[np.argmax(silhouette_scores)]

fig, axs = plot_silhouette_param_scan(leiden_resolution_range, 
                                      silhouette_scores,
                                      n_clusters,
                                      save=save_param_scan
                                      )

# run again with optimal parameters 
sc.tl.leiden(adata,
             resolution=optimal_resolution,
             key_added='leiden',
             copy=False
            )

sc.tl.rank_genes_groups(adata, 
                        groupby='leiden', 
                        method='wilcoxon', 
                        corr_method='benjamini-hochberg'
                        )


sc.pl.rank_genes_groups_heatmap(adata, save=save_deg)


# UMAP plotting 
sc.tl.umap(adata, random_state=42)
sc.pl.umap(adata, color=['leiden'])
