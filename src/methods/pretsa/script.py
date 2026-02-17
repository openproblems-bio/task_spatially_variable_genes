import warnings
warnings.filterwarnings('ignore')

import anndata as ad
import numpy as np
import pandas as pd

# VIASH START
par = {
    'input_data': 'resources_test/task_spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'name': 'pretsa'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run PreTSA', flush=True)
from pretsa_spatial import spatialTest

coords = adata.obsm['spatial']
X = adata.layers['counts']

if hasattr(X, 'toarray'):
    X_dense = X.toarray()
else:
    X_dense = np.asarray(X)

# PreTSA expects library-size-normalized and log-transformed data
total_counts = X_dense.sum(axis=1, keepdims=True)
total_counts[total_counts == 0] = 1
median_total = np.median(total_counts[total_counts > 0])
X_norm = X_dense / total_counts * median_total
X_log = np.log1p(X_norm)

# PreTSA expects (genes x cells) matrix
expr = X_log.T

coord = {
    'row': coords[:, 0].astype(np.float64),
    'col': coords[:, 1].astype(np.float64),
}

result = spatialTest(expr, coord, knot=0)
# result columns: logpval, pval, fstat
pvals = result[:, 1]
scores = -np.log10(np.clip(pvals, 1e-300, 1))

df = pd.DataFrame({
    'feature_id': list(adata.var_names),
    'pred_spatial_var_score': scores
})

output = ad.AnnData(
    var=df,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': meta['name']
    }
)

print('Write output to file', flush=True)
output.write_h5ad(par['output'])
