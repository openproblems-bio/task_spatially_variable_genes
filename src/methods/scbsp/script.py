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
    'name': 'scbsp'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run scBSP', flush=True)
from scbsp import granp

coords = adata.obsm['spatial']
X = adata.layers['counts']

if hasattr(X, 'toarray'):
    X_dense = X.toarray()
else:
    X_dense = np.asarray(X)

result_df = granp(
    np.ascontiguousarray(coords, dtype=np.float64),
    np.ascontiguousarray(X_dense, dtype=np.float64),
)

# granp returns DataFrame with columns: gene_names, p_values
pvals = result_df['p_values'].values
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
