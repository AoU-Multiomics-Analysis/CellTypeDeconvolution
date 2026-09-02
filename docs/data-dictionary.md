# Workflow data dictionary

All names below use the `cell_type_deconvolution` workflow namespace. All WDL sources use WDL 1.0. The workflow accepts one coordinate-preserving BED or BED.GZ matrix. The matrix contains finite, strictly positive, normalized linear CPM. The workflow applies `log2()` once in each consuming analysis view.

## Inputs

| WDL name | WDL type | Mode | Default | Description |
| --- | --- | --- | --- | --- |
| `expression` | `File` | Both | None | Tab-delimited BED or BED.GZ. Leading columns are `#chr`, `start`, `end`, and `gene_id`, in that order. Remaining columns are samples. Values are finite, strictly positive, normalized linear CPM. |
| `gtf` | `File` | Both | None | GTF or GTF.GZ file with `gene` records and `gene_id` and `gene_name` attributes. No gene-type filter is applied. |
| `lm22` | `File?` | dtangle | None | User-supplied linear LM22 matrix. It has one gene-symbol column and the 22 standard LM22 columns. Values are finite and strictly positive. |
| `precomputed_proportions` | `File?` | precomputed | None | Table with one sample identifier column and the 22 standard LM22 columns. Sample IDs must match `expression`. |
| `covariates` | `File?` | Both | None | Optional mixture-level technical covariates. Sample IDs must match the TCA weights. |
| `docker_image` | `String` | Both | GitHub `latest` | Optional image override. The default is `ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest`. |
| `preemptible_attempts` | `Int` | Both | `2` | Global number of preemptible attempts. Applies to every task. |
| `max_retries` | `Int` | Both | `2` | Global retry limit. Applies to every task. |
| `min_lm22_overlap` | `Float` | dtangle | `0.80` | Minimum fraction of LM22 genes matched to the bulk data. Must be greater than 0 and no greater than 1. |
| `dtangle_marker_fraction` | `Float` | dtangle | `0.10` | Fraction of ranked ratio markers used by dtangle. |
| `dtangle_quantile_normalize` | `Boolean` | dtangle | `false` | Apply joint quantile normalization to dtangle data. |
| `group_mean_threshold` | `Float` | Both | `0.0001` | Retain a group only when its cohort mean is at least this value. |
| `zero_floor` | `Float` | Both | `0.000001` | Replace an exact retained-group zero before row normalization. |
| `tca_max_iters` | `Int` | Both | `10` | Maximum iterations for the cohort-wide TCA fit. |
| `random_seed` | `Int` | Both | `20260901` | Seed for deterministic R operations. |
| `dtangle_cpu` | `Int` | dtangle | `4` | CPU for LM22 validation and dtangle. |
| `dtangle_memory` | `String` | dtangle | `"32 GB"` | Memory for dtangle. |
| `dtangle_disk_gb` | `Int` | dtangle | `100` | Disk for dtangle. |
| `proportions_cpu` | `Int` | Both | `2` | CPU for proportion grouping and filtering. |
| `proportions_memory` | `String` | Both | `"16 GB"` | Memory for proportion processing. |
| `proportions_disk_gb` | `Int` | Both | `50` | Disk for proportion processing. |
| `fit_cpu` | `Int` | Both | `16` | CPU for the full-matrix TCA fit. |
| `fit_memory` | `String` | Both | `"192 GB"` | Memory for the TCA fit. |
| `fit_disk_gb` | `Int` | Both | `750` | Disk for the TCA fit. |
| `export_cpu` | `Int` | Both | `8` | CPU for direct BED export. |
| `export_memory` | `String` | Both | `"128 GB"` | Memory for direct BED export. |
| `export_disk_gb` | `Int` | Both | `500` | Disk for direct BED export. |
| `manifest_cpu` | `Int` | Both | `4` | CPU for manifest creation. |
| `manifest_memory` | `String` | Both | `"32 GB"` | Memory for manifest creation. |
| `manifest_disk_gb` | `Int` | Both | `100` | Disk for manifest creation. |

Set exactly one of `lm22` and `precomputed_proportions`. The workflow records the selected mode before it selects the proportion file. The global retry controls are the only preemptible and retry inputs.

## Outputs

| WDL name | WDL type | Mode | Description |
| --- | --- | --- | --- |
| `proportion_mode_validation_log` | `File` | Both | Log that records the selected proportion mode. |
| `estimated_proportions` | `File?` | dtangle | Raw dtangle proportions. Absent in precomputed mode. |
| `dtangle_markers` | `File?` | dtangle | dtangle marker genes and counts. |
| `dtangle_metadata` | `File?` | dtangle | JSON with dtangle settings and fitted metadata. |
| `dtangle_overlap_report` | `File?` | dtangle | Report of LM22 genes matched to the bulk data. |
| `transformed_lm22` | `File?` | dtangle | LM22 matrix after one `log2()` transform. |
| `dtangle_shared_bulk` | `File?` | dtangle | LM22-shared bulk matrix in LM22 order after one `log2()` transform. |
| `dtangle_log` | `File?` | dtangle | dtangle task log. |
| `proportions_lm22` | `File` | Both | Input proportions in normalized output format. |
| `proportions_combined` | `File` | Both | Ten combined lineage proportions before filtering. |
| `tca_weights` | `File` | Both | Retained, zero-adjusted proportions used by TCA. Each sample row sums to one. |
| `cell_group_filter_report` | `File` | Both | Cohort means and retained status for the ten groups. |
| `proportions_log` | `File` | Both | Proportion-processing task log. |
| `tca_model` | `File` | Both | Serialized cohort-wide TCA model. |
| `tca_model_log` | `File` | Both | TCA model-fit log. |
| `tca_expression` | `File` | Both | Full TCA `log2_cpm` matrix after constant-gene removal. |
| `tca_excluded_genes` | `File` | Both | Genes excluded before TCA fitting because they are constant or invalid. |
| `fit_tca_log` | `File` | Both | TCA-fit task log. |
| `cell_type_beds` | `Array[File]` | Both | One BED.GZ file for each retained group. Files preserve coordinates, modeled-gene order, and sample order. |
| `cell_type_bed_inventory` | `File` | Both | BED output inventory with stable basenames, dimensions, and scale. |
| `reconstruction_by_sample` | `File` | Both | Per-sample TCA reconstruction metrics. |
| `qc_summary` | `File` | Both | Final QC summary with validation, exclusions, reconstruction, and convergence metrics. |
| `qc_plots` | `File` | Both | Reconstruction quality-control PDF. |
| `export_log` | `File` | Both | Direct BED-export task log. |
| `export_detail_log` | `File` | Both | Detailed direct BED-export log. |
| `output_manifest` | `File` | Both | JSON with output checksums and effective parameters. |
| `output_inventory` | `File` | Both | Portable inventory with one stable BED basename per authoritative output. |
| `manifest_log` | `File` | Both | Manifest-creation task log. |
| `effective_parameters_file` | `File` | Both | JSON with effective scientific settings. |

The workflow does not have an expression-preparation task. It does not produce preparation outputs or a pipeline-version output. TCA reads every valid nonconstant gene from the BED and applies `log2()` once. TCA does not use the GTF. dtangle maps gene IDs with the GTF and applies `log2()` once.

## BED output contract

The workflow performs one direct tensor extraction. Each output starts with `#chr`, `start`, `end`, and `gene_id`. The remaining columns are sample IDs in input order. Rows are valid, nonconstant TCA genes in input BED order. Values use the `log2_cpm` model scale. The values are statistical estimates. They are not counts, linear CPM, or absolute expression values.

`cell_type_beds` is the authoritative file array. Public inventories and the manifest use stable BED basenames. The manifest uses localized paths only while it calculates checksums.

The workflow does not produce HDF5 or tensor shards.
