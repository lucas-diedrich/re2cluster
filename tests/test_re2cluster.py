import shutil
import scanpy as sc

import pytest
from pathlib import Path 

from re2cluster import re2cluster 




@pytest.fixture(scope='module')
def tmp_dataset_dir(tmpdir_factory):
    """ Store datasets in temporary directory. Method from 
    scanpy test suite (Heumos et al)
    """ 
    new_dir =Path(tmpdir_factory.mktemp('data'))
    old_dir = sc.settings.datasetdir 
    sc.settings.datasetdir = new_dir 
    yield sc.settings.datasetdir 
    sc.settings.datasetdir = old_dir


def test_re2cluster(tmp_dataset_dir, tmp_path): 

    adata = sc.datasets.pbmc3k()
    adata, re2cluster_parameters, re2cluster_markers = re2cluster(adata, 
                                                                  outdir=tmp_path, 
                                                                  save_deg='deg.png', 
                                                                  save_param_scan='param-scan.png', 
                                                                  save_umap='umap.png'
                                                                  )

    assert 'leiden_tier_0' in adata.obs.columns
