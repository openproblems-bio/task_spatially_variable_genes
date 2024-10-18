#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

echo "Running benchmark on test data"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="testrun_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="temp/results/${RUN_ID}"

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels_ci.config \
  --input_states resources_test/task_spatially_variable_genes/mouse_brain_coronal/state.yaml \
  --rename_keys 'input_dataset:output_dataset;input_solution:output_solution' \
  --publish_dir "$publish_dir" \
  --output_state state.yaml \
  --settings '{"coord_type_moran_i":"generic","coord_type_sepal":"grid","max_neighs_sepal":6}'