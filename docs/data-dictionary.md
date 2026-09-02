# Workflow data dictionary

All names below use the `cell_type_deconvolution` workflow namespace. The workflow accepts a BED or BED.GZ matrix with positive linear CPM values. It creates `log2(CPM)` without a pseudocount. `Both` means dtangle and precomputed modes.

## Inputs

| WDL name | WDL type | Mode | Default | Scale | Validation | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `expression` | `File` | Both | None | Linear CPM | Tab-delimited BED or BED.GZ. Leading columns are `#chr`, `start`, `end`, and `gene_id`, in that order. Sample values are finite and strictly positive. | Input coordinate-preserving gene-by-sample CPM matrix. |
| `gtf` | `File` | Both | None | N/A | GTF or GTF.GZ `gene` records need `gene_id` and `gene_name`. | Gene annotation. No gene-type filter is applied. |
| `lm22` | `File?` | dtangle | None | Linear expression | One gene-symbol column, 22 standard LM22 columns, unique trimmed genes, and finite positive values. Overlap meets `min_lm22_overlap`. | User-supplied LM22 signature. Omit in precomputed mode. |
| `precomputed_proportions` | `File?` | precomputed | None | Relative proportion | One sample identifier column, 22 standard LM22 columns, and sample IDs that match prepared expression. | Precomputed LM22 proportions. Omit in dtangle mode. |
| `covariates` | `File?` | Both | None | User-defined technical values | Sample IDs match TCA weights. Values are finite and numeric. | Optional mixture-level TCA covariates. |
| `docker_image` | `String` | Both | None | N/A | Production values contain `@sha256:` followed by exactly 64 lowercase hexadecimal characters. The exact `celltype-deconvolution:test` tag is permitted only for local smoke CI. | Image used by every task. Replace the example digest placeholder before a production run. |
| `pipeline_version` | `String` | Both | `"0.1.0"` | N/A | Non-empty version identifier. | Recorded in the manifest and output inventory. |
| `min_lm22_overlap` | `Float` | dtangle | `0.80` | Fraction | Finite value greater than 0 and no greater than 1. | Minimum fraction of LM22 genes matched to bulk data. |
| `dtangle_marker_fraction` | `Float` | dtangle | `0.10` | Fraction | Finite value greater than 0 and no greater than 1. | Fraction of ranked ratio markers used by dtangle. |
| `dtangle_quantile_normalize` | `Boolean` | dtangle | `false` | N/A | Boolean. | Apply joint quantile normalization to dtangle-only data. |
| `group_mean_threshold` | `Float` | Both | `0.0001` | Relative proportion | Finite non-negative value. | Remove a group when its cohort mean is below this value. |
| `zero_floor` | `Float` | Both | `0.000001` | Relative proportion | Finite positive value. | Replace exact retained-group zeros before row normalization. |
| `tca_max_iters` | `Int` | Both | `10` | Iterations | Positive integer. | Maximum iterations for the cohort-wide TCA fit. |
| `random_seed` | `Int` | Both | `20260901` | N/A | Integer. | Seed for deterministic R operations. |
| `prepare_cpu` | `Int` | Both | `4` | CPU | Positive integer runtime value. | CPU for expression preparation. |
| `prepare_memory` | `String` | Both | `"64 GB"` | Memory | Valid Terra memory string. | Memory for expression preparation. |
| `prepare_disk_gb` | `Int` | Both | `400` | GB | Positive integer runtime value. | Disk for expression preparation. |
| `prepare_preemptible_attempts` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value. | Preemptible attempts for expression preparation. |
| `prepare_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value. | Retry limit for expression preparation. |
| `dtangle_cpu` | `Int` | dtangle | `4` | CPU | Positive integer runtime value. | CPU for LM22 validation and dtangle. |
| `dtangle_memory` | `String` | dtangle | `"32 GB"` | Memory | Valid Terra memory string. | Memory for dtangle. |
| `dtangle_disk_gb` | `Int` | dtangle | `100` | GB | Positive integer runtime value. | Disk for dtangle. |
| `dtangle_preemptible_attempts` | `Int` | dtangle | `2` | Attempts | Non-negative integer runtime value. | Preemptible attempts for dtangle. |
| `dtangle_max_retries` | `Int` | dtangle | `2` | Attempts | Non-negative integer runtime value. | Retry limit for dtangle. |
| `proportions_cpu` | `Int` | Both | `2` | CPU | Positive integer runtime value. | CPU for grouping and filtering proportions. |
| `proportions_memory` | `String` | Both | `"16 GB"` | Memory | Valid Terra memory string. | Memory for grouping and filtering proportions. |
| `proportions_disk_gb` | `Int` | Both | `50` | GB | Positive integer runtime value. | Disk for grouping and filtering proportions. |
| `proportions_preemptible_attempts` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value. | Preemptible attempts for proportion processing. |
| `proportions_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value. | Retry limit for proportion processing. |
| `fit_cpu` | `Int` | Both | `16` | CPU | Positive integer runtime value. | CPU for the full-matrix TCA fit. |
| `fit_memory` | `String` | Both | `"192 GB"` | Memory | Valid Terra memory string. | Memory for the full-matrix TCA fit. |
| `fit_disk_gb` | `Int` | Both | `750` | GB | Positive integer runtime value. | Disk for the full-matrix TCA fit. |
| `fit_preemptible_attempts` | `Int` | Both | `0` | Attempts | Non-negative integer runtime value. | Preemptible attempts for TCA fitting. |
| `fit_max_retries` | `Int` | Both | `1` | Attempts | Non-negative integer runtime value. | Retry limit for TCA fitting. |
| `export_cpu` | `Int` | Both | `8` | CPU | Positive integer runtime value. | CPU for direct BED export. |
| `export_memory` | `String` | Both | `"128 GB"` | Memory | Valid Terra memory string. | Memory for direct BED export. |
| `export_disk_gb` | `Int` | Both | `500` | GB | Positive integer runtime value. | Disk for direct BED export. |
| `export_preemptible_attempts` | `Int` | Both | `0` | Attempts | Non-negative integer runtime value. | Preemptible attempts for direct BED export. |
| `export_max_retries` | `Int` | Both | `1` | Attempts | Non-negative integer runtime value. | Retry limit for direct BED export. |
| `manifest_cpu` | `Int` | Both | `4` | CPU | Positive integer runtime value. | CPU for manifest creation. |
| `manifest_memory` | `String` | Both | `"32 GB"` | Memory | Valid Terra memory string. | Memory for manifest creation. |
| `manifest_disk_gb` | `Int` | Both | `100` | GB | Positive integer runtime value. | Disk for manifest creation. |
| `manifest_preemptible_attempts` | `Int` | Both | `1` | Attempts | Non-negative integer runtime value. | Preemptible attempts for manifest creation. |
| `manifest_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value. | Retry limit for manifest creation. |

## Outputs

| WDL name | WDL type | Mode | Default | Scale | Validation | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `prepared_tca_cpm` | `File` | Both | N/A | Linear CPM | Mapped TCA rows are finite and strictly positive. | Prepared CPM matrix for TCA. It keeps mapped `gene_id` rows and order. |
| `prepared_tca_log2_cpm` | `File` | Both | N/A | `log2(CPM)` | Row and sample order match `prepared_tca_cpm`. | Log-scale prepared matrix for TCA. |
| `prepared_dtangle_cpm` | `File` | Both | N/A | Linear CPM | One row per mapped gene symbol after duplicate-symbol summation. | Prepared CPM matrix for dtangle. |
| `prepared_dtangle_log2_cpm` | `File` | Both | N/A | `log2(CPM)` | Row and sample order match `prepared_dtangle_cpm`. | Log-scale prepared matrix for dtangle. |
| `prepared_coordinates` | `File` | Both | N/A | Genomic coordinates | Rows match prepared TCA genes. | Preserved `#chr`, `start`, `end`, and `gene_id` values. |
| `mapping_report` | `File` | Both | N/A | N/A | One record per input expression row. | GTF mapping and duplicate-symbol report. |
| `prepare_excluded_genes` | `File` | Both | N/A | N/A | Lists removed preparation rows and reasons. | Unmapped or invalid rows excluded during preparation. |
| `prepare_log` | `File` | Both | N/A | N/A | Task command log exists. | Expression-preparation log. |
| `proportion_mode_validation_log` | `File` | Both | N/A | N/A | Records exactly one selected mode. | Validation log for the `lm22` versus `precomputed_proportions` selection. |
| `estimated_proportions` | `File?` | dtangle | N/A | Relative proportion | Has sample IDs and 22 LM22 columns. | Raw dtangle estimates. Absent in precomputed mode. |
| `dtangle_markers` | `File?` | dtangle | N/A | N/A | Marker rows are from the dtangle fit. | Marker genes and counts. |
| `dtangle_metadata` | `File?` | dtangle | N/A | N/A | JSON metadata exists. | dtangle settings and fitted gamma information. |
| `dtangle_overlap_report` | `File?` | dtangle | N/A | N/A | Reports matched and missing LM22 genes. | LM22-to-bulk overlap report. |
| `transformed_lm22` | `File?` | dtangle | N/A | `log2(LM22)` | Positive linear LM22 was transformed once. | LM22 matrix used by dtangle. |
| `dtangle_shared_bulk` | `File?` | dtangle | N/A | `log2(CPM)` | Genes are LM22-shared and in LM22 order. | Bulk matrix used only by dtangle. |
| `dtangle_log` | `File?` | dtangle | N/A | N/A | Task command log exists. | dtangle log. |
| `proportions_lm22` | `File` | Both | N/A | Relative proportion | Has the original 22 LM22 columns. | Input proportions in a normalized output format. |
| `proportions_combined` | `File` | Both | N/A | Relative proportion | Has ten combined groups before filtering. | Combined lineage proportions. |
| `tca_weights` | `File` | Both | N/A | Relative proportion | Retained values are positive and each sample row sums to one. | Combined, filtered, zero-adjusted proportions used by TCA. |
| `cell_group_filter_report` | `File` | Both | N/A | Relative proportion | Reports cohort means and retained status for ten groups. | Cohort filter result. |
| `proportions_log` | `File` | Both | N/A | N/A | Task command log exists. | Proportion-processing log. |
| `tca_model` | `File` | Both | N/A | N/A | Serialized cohort-wide model exists. | TCA model fitted with `refit_W = FALSE`. |
| `tca_model_log` | `File` | Both | N/A | N/A | Model log exists. | TCA model-fit log. |
| `tca_expression` | `File` | Both | N/A | `log2(CPM)` | Contains mapped nonconstant genes in stable gene and sample order. | Full TCA matrix after constant-gene removal. |
| `tca_excluded_genes` | `File` | Both | N/A | N/A | Lists constant-gene removal reasons. | Genes excluded before TCA fitting. |
| `fit_tca_log` | `File` | Both | N/A | N/A | Task command log exists. | TCA-fit log. |
| `cell_type_beds` | `Array[File]` | Both | N/A | `log2_cpm` | One unique `.bed.gz` path for each retained group. Each file preserves coordinates, modeled-gene order, and sample order. | Primary cell-group BED.GZ outputs. |
| `cell_type_bed_inventory` | `File` | Both | N/A | `log2_cpm` | One row per retained group with a stable BED basename, dimensions, and scale. | BED output inventory. `cell_type_beds` is authoritative for localized files. |
| `reconstruction_by_sample` | `File` | Both | N/A | `log2_cpm` | Contains per-sample reconstruction metrics. | TCA reconstruction quality control. |
| `qc_summary` | `File` | Both | N/A | N/A | Records reconstruction metrics, LM22 validation when applicable, proportion row-sum errors, normalization adjustment, duplicate and exclusion counts, and TCA convergence. | Final pipeline quality-control summary. TCA convergence is derived from the model log because TCA 1.2.1 has no model convergence field. |
| `qc_plots` | `File` | Both | N/A | N/A | Quality-control PDF exists. | Reconstruction quality-control plots. |
| `export_log` | `File` | Both | N/A | N/A | Task command log exists. | Direct BED-export task log. |
| `export_detail_log` | `File` | Both | N/A | N/A | R export log exists. | Detailed direct BED-export log. |
| `output_manifest` | `File` | Both | N/A | N/A | Machine-readable JSON contains stable BED basenames and SHA-256 checksums. | Output checksums and effective parameters. Localized paths are used only during checksum calculation. |
| `output_inventory` | `File` | Both | N/A | N/A | Contains one stable BED basename for each authoritative `cell_type_beds` file. | Portable final provenance inventory. |
| `manifest_log` | `File` | Both | N/A | N/A | Task command log exists. | Manifest-creation log. |
| `effective_parameters_file` | `File` | Both | N/A | N/A | JSON is written by the top-level workflow. | Effective defaults and user settings. |

## BED output contract

The workflow performs one direct tensor extraction. It writes one BED.GZ file per retained major cell group. Each file starts with `#chr`, `start`, `end`, and `gene_id`. The remaining columns are sample IDs in input order. Rows are mapped, nonconstant TCA genes in input BED order. Values use the `log2_cpm` model scale.

The workflow does not produce HDF5 or tensor shards.
