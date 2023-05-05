#!/usr/bin/env python 

#%%
import numpy as np
import pandas as pd

import anndata
import scanpy as sc

from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

from collections import namedtuple

from tqdm import tqdm

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

    return p1, p2, t1


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
                          leiden_steps: int) -> Tuple[Iterable, Iterable, Iterable]:
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
        leiden_steps : int 
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
    leiden_resolution_range = np.logspace(np.log10(leiden_resolution_min), 
                                          np.log10(leiden_resolution_max), 
                                          leiden_steps)
    
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
                               silhouette_scores: Iterable, 
                               n_clusters: Iterable, 
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


def cluster(adata: anndata.AnnData, 
            leiden_resolution_min: float = 0.01,
            leiden_resolution_max: float = 1.4, 
            leiden_steps: int = 30, 
            save_deg: str = None, 
            save_param_scan: str = None, 
            save_umap: str = None):
    """ Wrapper function for workflow based on arbolpy. 

    Performs automated PCA, neighbor graph building, clustering 
    and UMAP visualization. 

    Returns information about the clustering as named tuple. Optimal parameters
    (nr of pcs, nr of neighbors, optimal resolution), cluster assignment 
    and marker genes. 
    
    Parameter
    ---------
    
    adata : anndata.Anndata 
        anndata Object. 
    leiden_resolution_min : float 
        Minimal resolution tested during leiden clustering
    leiden_resolution_max : float 
        Minimal resolution tested
    leiden_steps : int
        Number of different resolutions to test. Evenly (linearly) spaced 
        between min and max resolution
    save_deg : None|str 
        Whether to safe heatmap of differentially expressed genes. 
        Provide path if yes, default is None (no saving)
    save_param_scan : None|str 
        Whether to to safe parameter scan results. 
        Provide path if yes, default is None (no saving)
    save_umap : None|str
        Whether to safe umap plot with cluster annotation
        Provide path if yes, default is None (no saving)

    
    Returns 
    -------
    n_pcs : int 
        Number of PCs that where considered
    n_neighbors : int
        Number of neighbors that where considered 
    optimal_resolution : float 
        Optimal resolution for leiden clustering
    leiden_clusters : pandas.Series
        Assignment of clusters to cells 
    gene_groups : pandas.DataFrame
        

    """
    cluster_info = namedtuple('cluster_info', ['n_pcs', 
                                              'n_neighbors', 
                                              'optimal_resolution', 
                                              'leiden_clusters', 
                                              'markers']
                                              )

    # PCA
    # number of PCs can be at most the number of observations in dataset
    n_comps = min(adata.n_obs - 1, 50)
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable= True, svd_solver='arpack')
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
                                                                                   leiden_steps
                                                                                   )                                                                             
    
    optimal_resolution = leiden_resolution_range[np.nanargmax(silhouette_scores)]

    fig, axs = plot_silhouette_param_scan(leiden_resolution_range, 
                                          silhouette_scores,
                                          n_clusters,
                                          save=save_param_scan
                                          )
    
    # run again with optimal parameters 
    sc.tl.leiden(adata, resolution=optimal_resolution, key_added='leiden', copy=False)

    # Get marker genes
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', corr_method='benjamini-hochberg')

    # Plotting
    sc.pl.rank_genes_groups_heatmap(adata, save=save_deg)


    # UMAP plotting 
    sc.tl.umap(adata, random_state=42)
    sc.pl.umap(adata, color=['leiden'], save=save_umap)

    return cluster_info(n_pcs, n_neighbors, optimal_resolution, adata.obs['leiden'], sc.get.rank_genes_groups_df(adata, group=None))


def re2cluster(adata: anndata.AnnData, 
               leiden_resolution_min: float = 0.01,
               leiden_resolution_max: float = 1.4,
               leiden_steps: int = 30,
               normalization_method: Literal['tpm'] = 'tpm',
               n_hvg: int = 2000,
               min_cluster_size: int = 50,
               save_qc: str = None,
               save_deg: str = None,
               save_param_scan: str = None,
               save_umap: str = None,
               n_tiers: int = 3) -> anndata.AnnData: 
    """ Run automated QC and reclustering algorithm """

    # QC 
    adata = quality_control(adata)

    # ## Plotting 
    a, b, c = plot_quality_control(adata, save=save_qc)

    # Normalization 
    if normalization_method == 'tpm':
        adata = normalize_hvg_tpm(adata, n_hvg=n_hvg)
    # for pearson normalization, currently not implemented 
    # elif normalization_method == 'pearson':
    #     adata = normalize_hvg_pearson(adata, 
    #                                   min_cells_per_gene=min_cells_per_gene, 
    #                                   n_hvg=n_hvg
    #                                   )
    else:
        raise ValueError(f'Normalization method {normalization_method} not available. Select "tpm"')


    # Root node (includes all cells)
    adata.obs[f'leiden_tier_0'] = '0'

    # Parameter dict is supposed to have the format 
    # { tier: 
    #   {node : 
    #       {n_pcs: XXX, 
    #        n_neighbors: XXX, 
    #        optimal_resolution: XXX
    #        }
    #    ...
    #   }
    #   ...
    # }
    adata.uns['re2cluster_parameters'] = dict()

    for tier in tqdm(range(1, n_tiers+1)):

        parameter_dict_tier = dict()

        leiden_list = list()

        # adata.obs is just a dataframe 
        # groupby all previous cluster assignments of tiers to get unique cluster index per cell
        # nodes has the form dict{(idx0, .., idxn): [x0, ..., xn]} where idx is the node identity
        # and x is the row index
        nodes = adata.obs.groupby(by=[f'leiden_tier_{tier_idx}' 
                                        for tier_idx in range(tier)]).groups

        for node, indices in nodes.items():

            adata_tmp = adata[indices, :].copy()

            # Stop conditions (no clustering)

            # Minimal cluster size reached - do not cluster
            # Cluster assignment will be nan (implicitly handled by pd.concat)
            # In the later steps, np.nan is not a valid groupby key so that all subsequent cluster assignments will also be nan (endpoint reached)
            if adata_tmp.n_obs < min_cluster_size: 
                parameter_dict_tier[node] = dict(n_pcs=None, n_neighbors=None, optimal_resolution=None)

                continue

            cluster_data = cluster(adata_tmp, 
                                    leiden_resolution_min=leiden_resolution_min, leiden_resolution_max=leiden_resolution_max, 
                                    leiden_steps=leiden_steps, 
                                    save_deg=None, 
                                    save_param_scan=None, 
                                    save_umap=None)
            
            parameter_dict_tier[node] = dict(n_pcs=cluster_data.n_pcs, n_neighbors=cluster_data.n_neighbors, optimal_resolution=cluster_data.optimal_resolution)

            leiden_list.append(pd.Series(cluster_data.leiden_clusters, name=f'leiden_tier_{tier}'))

            del adata_tmp
        
        # concatenate all cluster assignments of current tier
        leiden = pd.concat(leiden_list)

        # Store info in adata object
        adata.obs = pd.concat([adata.obs, leiden], axis=1)
        adata.uns['re2cluster_parameters'][tier] = parameter_dict_tier

    return adata 
    

if __name__ == '__main__': 

    adata = sc.datasets.pbmc3k()
    re2cluster(adata)
