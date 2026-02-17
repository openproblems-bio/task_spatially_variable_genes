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
    'name': 'hotspot'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run Hotspot', flush=True)
import hotspot

coords = adata.obsm['spatial']
X = adata.layers['counts']

adata_hs = ad.AnnData(X=X.copy())
adata_hs.var_names = list(adata.var_names)
adata_hs.obsm['spatial'] = coords.copy()
adata_hs.obs['total_counts'] = np.array(X.sum(axis=1)).flatten()

hs = hotspot.Hotspot(adata_hs, model='danb', latent_obsm_key='spatial')
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations()

# Score: -log10(FDR) weighted by positive Z-score
logp = -np.log10(np.clip(hs_results['FDR'].values, 1e-300, 1))
z = hs_results['Z'].values
scores_dict = dict(zip(
    hs_results.index.tolist(),
    logp * np.sign(np.clip(z, 0, None))
))

scores = [scores_dict.get(gene, 0.0) for gene in adata.var_names]

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
