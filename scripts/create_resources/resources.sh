#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

cat > /tmp/params.yaml << 'HERE'
param_list:
  - id: svg_process_datasets_visium
    input_states: "s3://openproblems-data/resources/datasets/tenx_visium/visium/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200, "coord_type_proc": "grid"}'
  
  - id: svg_process_datasets_zenodo_visium
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial/visium/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200, "coord_type_proc": "grid"}'

  - id: svg_process_datasets_post_xenium
    input_states: "s3://openproblems-data/resources/datasets/tenx_visium/post_xenium/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 100, "coord_type_proc": "grid"}'

  - id: svg_process_datasets_slidetags
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial_slidetags/slidetags/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 50, "coord_type_proc": "grid"}'

  - id: svg_process_datasets_slideseqv2
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial/slideseqv2/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 10, "num_reference_genes": 10, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_dbitseq
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial/dbitseq/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_seqfish
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial/seqfish/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 25, "num_reference_genes": 25, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_starmap
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial/starmap/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 25, "num_reference_genes": 25, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_stereoseq
    input_states: "s3://openproblems-data/resources/datasets/zenodo_spatial/stereoseq/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 50, "coord_type_proc": "generic"}'

rename_keys: 'input:output_dataset'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/task_spatially_variable_genes/datasets
HERE

# cat > /tmp/params.yaml << 'HERE'
# param_list:
#   - id: zenodo_spatial/merfish/human_cortex_1
#     input: "s3://openproblems-data/resources/datasets/zenodo_spatial/merfish/human_cortex_1/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 25
#     num_reference_genes: 25
#     coord_type_proc: generic
  
#   - id: zenodo_spatial/merfish/human_cortex_2
#     input: "s3://openproblems-data/resources/datasets/zenodo_spatial/merfish/human_cortex_2/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 50
#     num_reference_genes: 50
#     coord_type_proc: generic
  
#   - id: zenodo_spatial/merfish/human_cortex_3
#     input: "s3://openproblems-data/resources/datasets/zenodo_spatial/merfish/human_cortex_3/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 50
#     num_reference_genes: 50
#     coord_type_proc: generic
  
#   - id: zenodo_spatial/merfish/human_cortex_4
#     input: "s3://openproblems-data/resources/datasets/zenodo_spatial/merfish/human_cortex_4/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 50
#     num_reference_genes: 50
#     coord_type_proc: generic
  
#   - id: zenodo_spatial/merfish/mouse_cortex
#     input: "s3://openproblems-data/resources/datasets/zenodo_spatial/merfish/mouse_cortex/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 25
#     num_reference_genes: 25
#     coord_type_proc: generic

# output_dataset: "$id/dataset.h5ad"
# output_solution: "$id/solution.h5ad"
# dataset_simulated_normalized: "$id/simulated_dataset.h5ad"
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/task_spatially_variable_genes/datasets
# HERE

tw launch https://github.com/openproblems-bio/task_spatially_variable_genes.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_spatially_variable_genes,process_datasets
