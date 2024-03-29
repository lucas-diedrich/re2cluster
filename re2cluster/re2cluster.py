#!/usr/bin/env python 
#%%
import os 

import numpy as np
import pandas as pd
from scipy import stats

import anndata
import scanpy as sc 

from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

from collections import namedtuple

from typing import List, Tuple, Iterable, Literal


def quality_control(adata):
    """ Calculate qc parameters """

    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])[0]
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    sc.pp.calculate_qc_metrics(adata,
                               qc_vars=['mt'],
                               log1p=True,
                               percent_top=[20],
                               inplace=True
                               )
    
    return adata


def flag_outliers(adata, metric: str, nmads: float = 5):
    """ Define outliers as cells that deviate by more than N MADs from 
    the overall cell population. Proposed by sc-best-practices. 
    
    The MAD is defined as 

    .. math 
        MAD = median( |X_i - median(X)| ), X_i \in X    

    Parameters
    ----------
    adata: anndata.Adata
        Dataset for which quality metrics shall be calculated
    metric: str
        Column in adata.obs that represents the desired quality metric
        (e.g. mitochondrial%, ribosomal% etc.)
    nmads: float
        Minimal number of median absolut deviations to define a cell as 
        outlier. Defaults to 5, which is a permissive value according to 
        sc-best-practices. 

    Returns 
    -------
    outlier : pd.Series
        pd.Series that flags outliers 

    Citation 
    --------
    Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities. 
    Nat Rev Genet (2023). https://doi.org/10.1038/s41576-023-00586-w

    Germain, PL., Sonrel, A. & Robinson, M.D. pipeComp, a general framework for the evaluation 
    of computational pipelines, reveals performant single cell RNA-seq preprocessing tools. 
    Genome Biol 21, 227 (2020). https://doi.org/10.1186/s13059-020-02136-7
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * stats.median_abs_deviation(M)) | (
        np.median(M) + nmads * stats.median_abs_deviation(M) < M
    )
    return outlier


def plot_quality_control(adata: anndata.AnnData, save: str = None):
    """ Plot quality control parameters for initial object """

    p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes',
                       color='pct_counts_mt'
                       )
    
    p2 = sc.pl.scatter(adata, 'n_counts', 'n_genes_by_counts',
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
    """ Normalize to transcripts per 10000, log1p normalize and find 
    highly variable genes (hvg) 
    """

    # normalize data with TPM + log1p
    adata.raw = adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # get highly variable genes 
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)

    return adata


def select_ideal_pcs(adata: anndata.AnnData, variance_ratio_quantile: float = 0.80): 
    """ Select ideal number of principal components of sc data 
    
    Current criterium: Change in explained variance is in the lowest 15 % quantile

    Parameters
    ----------
    adata : anndata.Anndata 
        anndata object. Is expected to contain information about pca
    variance_ratio_quantile : float 
        Quantile of accepted variance ratio diff. Arbol: 0.85, here we are more permissive
    """
    variance_ratio = adata.uns['pca']['variance_ratio']
    variance_ratio_diff = -1*np.diff(variance_ratio)

    # Number of PCs that fullfill criterium
    n_pcs = (variance_ratio_diff > np.quantile(variance_ratio_diff, variance_ratio_quantile)).sum()

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


def _cluster(adata: anndata.AnnData,
            variance_ratio_quantile: float,  
            leiden_resolution_min: float,
            leiden_resolution_max: float, 
            leiden_steps: int, 
            save_deg: str, 
            save_param_scan: str, 
            save_umap: str,
            show: bool = False 
            ):
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
    variance_ratio_quantile : float 
        Parameter for the selection of the ideal number of principal components
        during the PCA step. The lower, the more permissive. ARBOL defaults to 
        0.85, while the default here is 0.80. 
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
    show: 
        Default False

    
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
    n_pcs = select_ideal_pcs(adata, variance_ratio_quantile)

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
    fig = sc.pl.rank_genes_groups_heatmap(adata, return_fig=True)
    if save_deg is not None:
        plt.savefig(save_deg, bbox_inches='tight')


    # UMAP plotting 
    sc.tl.umap(adata, random_state=42)

    fig = sc.pl.umap(adata, color=['leiden'], return_fig=True)
    if save_umap is not None: 
        plt.savefig(save_umap, bbox_inches='tight')

    leiden_clusters = adata.obs['leiden']
    markers = sc.get.rank_genes_groups_df(adata, group=None)

    return cluster_info(n_pcs, 
                        n_neighbors, 
                        optimal_resolution, 
                        leiden_clusters, 
                        markers)


def re2cluster(adata: anndata.AnnData, 
              variance_ratio_quantile : float = 0.8,  
               leiden_resolution_min: float = 0.2,
               leiden_resolution_max: float = 1.4,
               leiden_steps: int = 30,
               normalization_method: Literal['tenk'] = 'tenk',
               n_hvg: int = 2000,
               min_cluster_size: int = 50,
               outdir = None,
               save_deg: str = None,
               save_param_scan: str = None,
               save_umap: str = None,
               n_tiers: int = 3, 
               show: bool = False
               ) -> anndata.AnnData: 
    """ Run automated QC and reclustering algorithm 
    Expects quality controlled anndata object and returns anndata object with assigned subclusters, 

    Parameters 
    ----------

    adata : anndata.Anndata 
        anndata Object.
    variance_ratio_quantile : float 
        Parameter for the selection of the ideal number of principal components
        during the PCA step. The lower, the more permissive. ARBOL defaults to 
        0.85, while the default here is 0.80. 
    leiden_resolution_min : float 
        Minimal resolution tested during leiden clustering
    leiden_resolution_max : float 
        Minimal resolution tested
    leiden_steps : int
        Number of different resolutions to test. Evenly (linearly) spaced 
        between min and max resolution
    normalization_method : Literal['tenk']
        Choose between normalization methods. Currently, only normalization to 10k counts is implemented
    min_cluster_size : int 
        Minimal cluster size to compute PCA for. Else, cluster identity nan is assigned to
        cluster in the current tier. 
    save_deg : None|str 
        Whether to safe heatmap of differentially expressed genes. 
        Provide path if yes, default is None (no saving)
    save_param_scan : None|str 
        Whether to to safe parameter scan results. 
        Provide path if yes, default is None (no saving)
    save_umap : None|str
        Whether to safe umap plot with cluster annotation
        Provide path if yes, default is None (no saving)
    n_tiers : int 
        Number of clustering steps


    Returns 
    -------
    adata : anndata.Anndata
        Modified anndata object that contains cluster assignment in adata.obs columns leiden_tier_{1..n_tier}
    re2cluster_paramters : pd.DataFrame 
        Optimal cluster parameters for every iteration (tier: 1...n_tier, subclusters: X.Y.xxx)
    re2cluster_markers : pd.DataFrame
        Markers (currently cluster vs. the respective subclusters, which is likely non-desirable)
    
    """

    if outdir is None: 
        outdir = os.getcwd()

    # Normalization 
    if normalization_method == 'tenk':
        adata = normalize_hvg_tpm(adata, n_hvg=n_hvg)
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
    # Stores optimal clustering parameters
    # adata.uns['re2cluster_parameters'] = list()
    re2cluster_parameters = list()

    # For every marker genes per cluster 
    # adata.uns['re2cluster_markers'] 
    re2cluster_markers = list()


    for tier in range(1, n_tiers+1):

        parameters_tier = list()
        markers_tier = list()

        leiden_list = list()

        # adata.obs is just a dataframe 
        # groupby all previous cluster assignments of tiers to get unique cluster index per cell
        # nodes has the form dict{(idx0, .., idxn): [x0, ..., xn]} where idx is the node identity
        # and x is the row index
        nodes = adata.obs.groupby(by=[f'leiden_tier_{tier_idx}' 
                                        for tier_idx in range(tier)]).groups

        for node, indices in nodes.items():
            
            node_name = '.'.join([str(nidx) for nidx in node])

            if save_deg is not None: 
                save_deg_path = os.path.join(outdir, f'{node_name}_{save_deg}')
            if save_param_scan is not None:
                save_param_scan_path = os.path.join(outdir, f'{node_name}_{save_param_scan}')
            if save_umap is not None:
                save_umap_path = os.path.join(outdir, f'{node_name}_{save_umap}')

            adata_tmp = adata[indices, :].copy()

            # Stop conditions (no clustering)

            # Minimal cluster size reached - do not cluster
            # Cluster assignment will be nan (implicitly handled by pd.concat)
            # In the later steps, np.nan is not a valid groupby key so that all subsequent cluster assignments will also be nan (endpoint reached)
            if adata_tmp.n_obs < min_cluster_size: 
                parameters_tier.append(pd.Series(dict(n_pcs=None, 
                                                       n_neighbors=None, 
                                                       optimal_resolution=None
                                                       ), 
                                             name=node_name)
                                    )
                
                
                del adata_tmp  
                continue

            cluster_data = _cluster(adata_tmp, 
                                    variance_ratio_quantile=variance_ratio_quantile,
                                    leiden_resolution_min=leiden_resolution_min, leiden_resolution_max=leiden_resolution_max, 
                                    leiden_steps=leiden_steps, 
                                    save_deg=save_deg_path, 
                                    save_param_scan=save_param_scan_path, 
                                    save_umap=save_umap_path, 
                                    show = show
                                    )

            leiden_list.append(pd.Series(cluster_data.leiden_clusters, name=f'leiden_tier_{tier}'))
            
            parameters_tier.append(pd.Series(dict(n_pcs=cluster_data.n_pcs, 
                                                       n_neighbors=cluster_data.n_neighbors, 
                                                       optimal_resolution=cluster_data.optimal_resolution
                                                       ), 
                                             name='.'.join(node)
                                             )
                                  )

            df_markers = cluster_data.markers
            df_markers['tier'] = tier
            df_markers['node'] = '.'.join(node)
            markers_tier.append(df_markers)

            del adata_tmp
        
        # concatenate all cluster assignments of current tier
        leiden = pd.concat(leiden_list)

        # Store info in adata object

        # Cluster assignments in current tier (only unique in combination with prior cluster assignments)
        adata.obs = pd.concat([adata.obs, leiden], axis=1)

        # Optimal cluster parameters
        df_parameters_tier = pd.concat(parameters_tier, axis=1).T
        df_parameters_tier['tier'] = tier
        re2cluster_parameters.append(df_parameters_tier)

        # Marker genes 
        re2cluster_markers.append(pd.concat(markers_tier, axis=0, join='inner'))

    re2cluster_parameters = pd.concat(re2cluster_parameters, axis=0, join='inner')
    re2cluster_markers = pd.concat(re2cluster_markers, axis=0, join='inner')

    return adata, re2cluster_parameters, re2cluster_markers