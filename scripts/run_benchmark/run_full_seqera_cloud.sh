#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/task_spatially_variable_genes/results/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
param_list:
  - id: svg_datasets_visium
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/tenx_visium/visium/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'
  
  - id: svg_datasets_zenodo_visium
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/visium/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_post_xenium
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/tenx_visium/post_xenium/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_dbitseq
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/dbitseq/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 4}'

  - id: svg_datasets_merfish
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/merfish/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_seqfish
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/seqfish/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_slidetags
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial_slidetags/slidetags/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_slideseqv2
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/slideseqv2/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_starmap
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/starmap/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_stereoseq
    input_states: "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/stereoseq/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 4}'

rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_spatially_variable_genes.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_template,full