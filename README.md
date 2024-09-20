# Spatially Variable Genes


<!--
This file is automatically generated from the tasks's api/*.yaml files.
Do not edit this file directly.
-->

Spatially variable genes (SVGs) are genes whose expression levels vary
significantly across different spatial regions within a tissue or across
cells in a spatially structured context.

Repository:
[openproblems-bio/task_spatially_variable_genes](https://github.com/openproblems-bio/task_spatially_variable_genes)

## Description

Recent years have witnessed significant progress in spatially-resolved
transcriptome profiling techniques that simultaneously characterize
cellular gene expression and their physical position, generating spatial
transcriptomic (ST) data. The application of these techniques has
dramatically advanced our understanding of disease and developmental
biology. One common task for all ST profiles, regardless of the employed
protocols, is to identify genes that exhibit spatial patterns. These
genes, defined as spatially variable genes (SVGs), contain additional
information about the spatial structure of the tissues of interest,
compared to highly variable genes (HVGs).

Identification of spatially variable genes is crucial to for studying
spatial domains within tissue microenvironmnets, developmental gradients
and cell signaling pathways. In this task we attempt to evaluate various
methods for detecting SVGs using a number of realistic simulated
datasets with diverse patterns derived from real-world spatial
transcriptomics data using scDesign3. Synthetic data is generated by
mixing a Gaussian Process (GP) model and a non-spatial model (obtained
by shuffling mean parameters of the GP model to remove spatial
correlation between spots) to generate gene expressions with various
spatial variability. For more details, please refer to our
[manuscript](https://www.biorxiv.org/content/10.1101/2023.12.02.569717v1)
and [Github](https://github.com/pinellolab/SVG_Benchmarking).

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Zhijian Li        | author, maintainer |
| Zain M. Patel     | author             |
| Dongyuan Song     | author             |
| Guanao Yan        | author             |
| Jingyi Jessica Li | author             |
| Luca Pinello      | author             |
| Robrecht Cannoodt | contributor        |
| Sai Nirmayi Yasa  | contributor        |

## API

``` mermaid
flowchart LR
  file_common_dataset("Common Dataset")
  file_dataset("Dataset")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  file_output("Output")
  comp_metric[/"Metric"/]
  file_score("Score")
  file_simulated_dataset("Common Dataset")
  file_solution("Solution")
  file_dataset---comp_control_method
  file_dataset---comp_method
  comp_control_method-->file_output
  comp_method-->file_output
  file_output---comp_metric
  comp_metric-->file_score
  file_solution---comp_control_method
  file_solution---comp_metric
```

## File format: Common Dataset

A subset of the common dataset.

Example file: `resources_test/common/mouse_brain_coronal/dataset.h5ad`

Format:

<div class="small">

    AnnData object
     var: 'feature_id', 'feature_name'
     obsm: 'spatial'
     layers: 'counts', 'counts'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `var["feature_id"]` | `string` | (*Optional*) Unique identifier for the feature, usually a ENSEMBL gene id. |
| `var["feature_name"]` | `string` | A human-readable name for the feature, usually a gene symbol. |
| `obsm["spatial"]` | `double` | Spatial coordinates for each spot. |
| `layers["counts"]` | `integer` | Raw counts. |
| `layers["counts"]` | `double` | Normalized expression values. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | (*Optional*) Nicely formatted name. |
| `uns["dataset_url"]` | `string` | (*Optional*) Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | (*Optional*) Short description of the dataset. |
| `uns["dataset_description"]` | `string` | (*Optional*) Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | (*Optional*) The organism of the sample in the dataset. |

</div>

## File format: Dataset

The dataset without spatially variable genes.

Example file:
`resources_test/task_spatially_variable_genes/mouse_brain_coronal/dataset.h5ad`

Format:

<div class="small">

    AnnData object
     var: 'feature_id', 'feature_name'
     obsm: 'spatial'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `var["feature_id"]` | `string` | (*Optional*) Unique identifier for the feature, in this case a ENSEMBL gene id suffixed with alpha value. |
| `var["feature_name"]` | `string` | (*Optional*) A human-readable name for the feature, in this case a gene symbol suffixed with alpha value. |
| `obsm["spatial"]` | `double` | Spatial coordinates for each spot. |
| `layers["counts"]` | `integer` | Raw counts. |
| `layers["normalized"]` | `double` | Normalised expression values. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | (*Optional*) Nicely formatted name. |

</div>

## Component type: Control method

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_data` | `file` | The dataset without spatially variable genes. |
| `--input_solution` | `file` | Anndata with true spatial variability. |
| `--output` | `file` | (*Output*) Anndata with estimate spatial variability. |

</div>

## Component type: Method

A spatially variable gene identification method.

Arguments:

<div class="small">

| Name | Type | Description |
|:---|:---|:---|
| `--input_data` | `file` | The dataset without spatially variable genes. |
| `--output` | `file` | (*Output*) Anndata with estimate spatial variability. |

</div>

## File format: Output

Anndata with estimate spatial variability.

Example file:
`resources_test/task_spatially_variable_genes/mouse_brain_coronal/output.h5ad`

Description:

Anndata with estimated spatial variability score for each gene.

Format:

<div class="small">

    AnnData object
     var: 'feature_id', 'feature_name', 'pred_spatial_var_score'
     uns: 'dataset_id', 'method_id'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `var["feature_id"]` | `string` | Feature ID. |
| `var["feature_name"]` | `string` | (*Optional*) Feature name. |
| `var["pred_spatial_var_score"]` | `double` | Predicted spatial variability score. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["method_id"]` | `string` | A unique identifier for the method. |

</div>

## Component type: Metric

A spatially variable genes identification metric.

Arguments:

<div class="small">

| Name               | Type   | Description                                |
|:-------------------|:-------|:-------------------------------------------|
| `--input_method`   | `file` | Anndata with estimate spatial variability. |
| `--input_solution` | `file` | Anndata with true spatial variability.     |
| `--output`         | `file` | (*Output*) Metric score file.              |

</div>

## File format: Score

Metric score file.

Example file:
`resources_test/task_spatially_variable_genes/mouse_brain_coronal/score.h5ad`

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["method_id"]` | `string` | A unique identifier for the method. |
| `uns["metric_ids"]` | `string` | One or more unique metric identifiers. |
| `uns["metric_values"]` | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>

## File format: Common Dataset

A subset of the common dataset.

Example file:
`resources_test/task_spatially_variable_genes/mouse_brain_coronal/simulated_dataset.h5ad`

Format:

<div class="small">

    AnnData object
     var: 'feature_id', 'feature_name', 'orig_feature_id', 'orig_feature_name', 'true_spatial_var_score'
     obsm: 'spatial'
     layers: 'counts'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `var["feature_id"]` | `string` | (*Optional*) Unique identifier for the feature, in this case a ENSEMBL gene id suffixed with alpha value. |
| `var["feature_name"]` | `string` | A human-readable name for the feature, in this case a gene symbol suffixed with alpha value. |
| `var["orig_feature_id"]` | `string` | (*Optional*) Original unique identifier for the feature, usually a ENSEMBL gene id. |
| `var["orig_feature_name"]` | `string` | Original human-readable name for the feature, usually a gene symbol. |
| `var["true_spatial_var_score"]` | `double` | True spatial variability score. |
| `obsm["spatial"]` | `double` | Spatial coordinates for each spot. |
| `layers["counts"]` | `integer` | Raw counts. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | Nicely formatted name. |
| `uns["dataset_url"]` | `string` | Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | The organism of the sample in the dataset. |

</div>

## File format: Solution

Anndata with true spatial variability.

Example file:
`resources_test/task_spatially_variable_genes/mouse_brain_coronal/solution.h5ad`

Description:

Anndata with true spatial variability score for each gene.

Format:

<div class="small">

    AnnData object
     var: 'feature_id', 'feature_name', 'orig_feature_name', 'true_spatial_var_score'
     obsm: 'spatial'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Data structure:

<div class="small">

| Slot | Type | Description |
|:---|:---|:---|
| `var["feature_id"]` | `string` | (*Optional*) Unique identifier for the feature (e.g., ESEMBL gene id suffixed with alpha value). |
| `var["feature_name"]` | `string` | A human-readable name for the feature, in this case a gene symbol suffixed with alpha value. |
| `var["orig_feature_name"]` | `string` | Original human-readable name for the feature, usually a gene symbol. |
| `var["true_spatial_var_score"]` | `double` | True spatial variability score. |
| `obsm["spatial"]` | `double` | Spatial coordinates for each spot. |
| `uns["dataset_id"]` | `string` | A unique identifier for the dataset. |
| `uns["dataset_name"]` | `string` | Nicely formatted name. |
| `uns["dataset_url"]` | `string` | Link to the original source of the dataset. |
| `uns["dataset_reference"]` | `string` | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]` | `string` | Short description of the dataset. |
| `uns["dataset_description"]` | `string` | Long description of the dataset. |
| `uns["dataset_organism"]` | `string` | The organism of the sample in the dataset. |

</div>

