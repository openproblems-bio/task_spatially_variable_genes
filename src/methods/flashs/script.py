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
    'name': 'flashs'
}
# VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input_data'])

print('Run FlashS', flush=True)
from flashs import FlashS

coords = adata.obsm['spatial']
X = adata.layers['counts']

# Use raw counts with multi-kernel Cauchy combination.
# adjustment='none' because ranking relies on the composite score.
model = FlashS(adjustment='none', random_state=42)

result = model.fit_test(coords, X, gene_names=list(adata.var_names))

# Composite score: normalized effect size + normalized -log10(p)
es = result.effect_size
logp = -np.log10(np.clip(result.pvalues, 1e-300, 1))

es_norm = (es - es.min()) / (es.max() - es.min() + 1e-10)
logp_norm = (logp - logp.min()) / (logp.max() - logp.min() + 1e-10)
combined = es_norm + logp_norm

# Map back to input gene order (handles filtered genes)
gene_to_score = dict(zip(result.gene_names, combined))
scores = [gene_to_score.get(gene, 0.0) for gene in adata.var_names]

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
