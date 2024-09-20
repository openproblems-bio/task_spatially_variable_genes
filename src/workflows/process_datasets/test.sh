#!/bin/bash

# Run this prior to executing this script:
# viash ns build --setup cb 

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/spatially_variable_genes

if [ ! -d "$DATASET_DIR" ]; then
  mkdir -p "$DATASET_DIR"
fi

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  -c common/nextflow_helpers/labels_ci.config \
  --id mouse_brain_coronal_section1_visium \
  --input $RAW_DATA/mouse_brain_coronal_section1_visium/dataset.h5ad \
  --output_dataset dataset.h5ad \
  --output_solution solution.h5ad \
  --dataset_simulated_normalized simulated_dataset.h5ad \
  --publish_dir $DATASET_DIR/mouse_brain_coronal_section1_visium \
  --output_state "state.yaml" \
  --gp_k_sim 50 \
  --select_top_variable_genes 50 \
  --num_reference_genes 200