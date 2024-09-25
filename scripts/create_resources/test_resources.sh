#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/task_spatially_variable_genes

mkdir -p $DATASET_DIR

echo "Running process_datasets"
nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -c common/nextflow_helpers/labels_ci.config \
  --id mouse_brain_coronal \
  --input $RAW_DATA/mouse_brain_coronal/dataset.h5ad \
  --output_dataset dataset.h5ad \
  --output_solution solution.h5ad \
  --dataset_simulated_normalized simulated_dataset.h5ad \
  --publish_dir $DATASET_DIR/mouse_brain_coronal \
  --output_state "state.yaml" \
  --gp_k_sim 50 \
  --select_top_variable_genes 50 \
  --num_reference_genes 200 \
  -profile docker 

echo "Running method"
viash run src/methods/moran_i/config.vsh.yaml -- \
  --input_data $DATASET_DIR/mouse_brain_coronal/dataset.h5ad \
  --output $DATASET_DIR/mouse_brain_coronal/output.h5ad

echo "Running metric"
viash run src/metrics/correlation/config.vsh.yaml -- \
    --input_method $DATASET_DIR/mouse_brain_coronal/output.h5ad \
    --input_solution $DATASET_DIR/mouse_brain_coronal/solution.h5ad \
    --output $DATASET_DIR/mouse_brain_coronal/score.h5ad

# only run this if you have access to the openproblems-data bucket
# aws s3 sync --profile op \
#   "$DATASET_DIR" s3://openproblems-data/resources_test/task_spatially_variable_genes \
#   --delete --dryrun