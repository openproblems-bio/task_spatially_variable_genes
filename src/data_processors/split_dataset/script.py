import anndata as ad
import sys 
import openproblems as op

## VIASH START
par = {
    "input": "resources_test/task_spatially_variable_genes/mouse_brain_coronal/simulated_dataset.h5ad",
    "output_dataset": "dataset.h5ad",
    "output_solution": "solution.h5ad",
}
meta = {
    "name": "process_dataset",
    "resources_dir": "src/tasks/task_spatially_variable_genes/process_dataset",
    "config": "target/nextflow/task_spatially_variable_genes/process_dataset/split_dataset/.config.vsh.yaml"
}
## VIASH END

# read viash config
config = op.project.read_viash_config(meta["config"])

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_h5ad_by_format import subset_h5ad_by_format

print(">> Load dataset", flush=True)
adata = ad.read_h5ad(par["input"])

print(">> Create dataset for methods", flush=True)
output_dataset = subset_h5ad_by_format(adata, config, "output_dataset")

print(">> Create solution object for metrics", flush=True)
output_solution = subset_h5ad_by_format(adata, config, "output_solution")

print(">> Write to disk", flush=True)
output_dataset.write_h5ad(par["output_dataset"])
output_solution.write_h5ad(par["output_solution"])
