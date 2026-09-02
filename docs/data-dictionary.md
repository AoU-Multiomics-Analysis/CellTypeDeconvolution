# Workflow data dictionary

All names below use the exact `cell_type_deconvolution` workflow namespace. The workflow accepts gene-by-sample linear CPM. It creates `log2(CPM)` without a pseudocount. "Both" means both proportion modes. "dtangle" means `lm22` is supplied. "Restart" means `precomputed_proportions` is supplied.

## Inputs

| WDL name | Type | Required mode | Default | Scale | Validation | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `expression` | `File` | Both | None | Linear CPM | First column is `gene_id`; values are finite and strictly positive | Gene-by-sample normalized CPM TSV or TSV.GZ. |
| `gtf` | `File` | Both | None | N/A | GTF/GTF.GZ `gene` records need `gene_id` and `gene_name` | Gene annotation for expression mapping. No type filter is applied. |
| `lm22` | `File?` | dtangle | None | Linear expression | Standard 22 LM22 columns; unique trimmed genes; positive finite values; overlap meets `min_lm22_overlap` | User-supplied LM22 signature. Omit in restart mode. |
| `precomputed_proportions` | `File?` | Restart | None | Relative proportion | Sample ID plus all 22 LM22 columns; samples match expression | Precomputed dtangle proportions. Omit in dtangle mode. |
| `covariates` | `File?` | Both | None | User-defined technical values | Sample IDs match TCA weights; finite numeric covariate values | Optional mixture-level TCA technical covariates. |
| `docker_image` | `String` | Both | None | N/A | Valid container image reference accessible to the execution backend | Project image that runs every task. |
| `pipeline_version` | `String` | Both | `"0.1.0"` | N/A | Non-empty version identifier is recommended | Recorded in HDF5 attributes and the manifest. |
| `min_lm22_overlap` | `Float` | dtangle | `0.80` | Fraction | Finite value from 0 to 1 | Minimum fraction of LM22 genes matched to bulk data. |
| `dtangle_marker_fraction` | `Float` | dtangle | `0.10` | Fraction | Finite value greater than 0 and no greater than 1 | Fraction of ranked ratio markers used by dtangle. |
| `dtangle_quantile_normalize` | `Boolean` | dtangle | `false` | N/A | Boolean | Apply joint quantile normalization to dtangle-only inputs. |
| `group_mean_threshold` | `Float` | Both | `0.0001` | Relative proportion | Finite non-negative value | Remove a cell group when its cohort mean is below this value. |
| `zero_floor` | `Float` | Both | `0.000001` | Relative proportion | Finite positive value | Replace exact retained-group zeros before row normalization. |
| `tca_shard_size` | `Int` | Both | `500` | Genes | Positive integer | Number of genes in each tensor-extraction shard. |
| `tca_max_iters` | `Int` | Both | `10` | Iterations | Positive integer | Maximum iterations for the cohort-wide TCA fit. |
| `random_seed` | `Int` | Both | `20260901` | N/A | Integer | Seed for deterministic R operations. |
| `write_tsv` | `Boolean` | Both | `false` | N/A | Boolean | Write optional compressed TSV copies of HDF5 group matrices. |
| `prepare_cpu` | `Int` | Both | `4` | CPU | Positive integer runtime value | CPU for CPM and GTF preparation. |
| `prepare_memory` | `String` | Both | `"64 GB"` | Memory | Valid Terra runtime memory string | Memory for CPM and GTF preparation. |
| `prepare_disk_gb` | `Int` | Both | `400` | GB | Positive integer runtime value | Disk for CPM and GTF preparation. |
| `prepare_preemptible_attempts` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Preemptible attempts for preparation. |
| `prepare_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Retry limit for preparation. |
| `dtangle_cpu` | `Int` | dtangle | `4` | CPU | Positive integer runtime value | CPU for LM22 validation and dtangle. |
| `dtangle_memory` | `String` | dtangle | `"32 GB"` | Memory | Valid Terra runtime memory string | Memory for dtangle. |
| `dtangle_disk_gb` | `Int` | dtangle | `100` | GB | Positive integer runtime value | Disk for dtangle. |
| `dtangle_preemptible_attempts` | `Int` | dtangle | `2` | Attempts | Non-negative integer runtime value | Preemptible attempts for dtangle. |
| `dtangle_max_retries` | `Int` | dtangle | `2` | Attempts | Non-negative integer runtime value | Retry limit for dtangle. |
| `proportions_cpu` | `Int` | Both | `2` | CPU | Positive integer runtime value | CPU for grouping and filtering proportions. |
| `proportions_memory` | `String` | Both | `"16 GB"` | Memory | Valid Terra runtime memory string | Memory for grouping and filtering proportions. |
| `proportions_disk_gb` | `Int` | Both | `50` | GB | Positive integer runtime value | Disk for grouping and filtering proportions. |
| `proportions_preemptible_attempts` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Preemptible attempts for proportions. |
| `proportions_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Retry limit for proportions. |
| `fit_cpu` | `Int` | Both | `16` | CPU | Positive integer runtime value | CPU for the cohort-wide TCA fit. |
| `fit_memory` | `String` | Both | `"192 GB"` | Memory | Valid Terra runtime memory string | Memory for the cohort-wide TCA fit. |
| `fit_disk_gb` | `Int` | Both | `750` | GB | Positive integer runtime value | Disk for the cohort-wide TCA fit. |
| `fit_preemptible_attempts` | `Int` | Both | `0` | Attempts | Non-negative integer runtime value | Preemptible attempts for TCA fitting. |
| `fit_max_retries` | `Int` | Both | `1` | Attempts | Non-negative integer runtime value | Retry limit for TCA fitting. |
| `extract_cpu` | `Int` | Both | `8` | CPU | Positive integer runtime value | CPU for each tensor shard. |
| `extract_memory` | `String` | Both | `"64 GB"` | Memory | Valid Terra runtime memory string | Memory for each tensor shard. |
| `extract_disk_gb` | `Int` | Both | `200` | GB | Positive integer runtime value | Disk for each tensor shard. |
| `extract_preemptible_attempts` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Preemptible attempts for each tensor shard. |
| `extract_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Retry limit for each tensor shard. |
| `assemble_cpu` | `Int` | Both | `8` | CPU | Positive integer runtime value | CPU for HDF5 assembly and quality control. |
| `assemble_memory` | `String` | Both | `"128 GB"` | Memory | Valid Terra runtime memory string | Memory for HDF5 assembly and quality control. |
| `assemble_disk_gb` | `Int` | Both | `500` | GB | Positive integer runtime value | Disk for HDF5 assembly and quality control. |
| `assemble_preemptible_attempts` | `Int` | Both | `0` | Attempts | Non-negative integer runtime value | Preemptible attempts for HDF5 assembly. |
| `assemble_max_retries` | `Int` | Both | `1` | Attempts | Non-negative integer runtime value | Retry limit for HDF5 assembly. |
| `manifest_cpu` | `Int` | Both | `4` | CPU | Positive integer runtime value | CPU for manifest creation. |
| `manifest_memory` | `String` | Both | `"32 GB"` | Memory | Valid Terra runtime memory string | Memory for manifest creation. |
| `manifest_disk_gb` | `Int` | Both | `100` | GB | Positive integer runtime value | Disk for manifest creation. |
| `manifest_preemptible_attempts` | `Int` | Both | `1` | Attempts | Non-negative integer runtime value | Preemptible attempts for manifest creation. |
| `manifest_max_retries` | `Int` | Both | `2` | Attempts | Non-negative integer runtime value | Retry limit for manifest creation. |

## Outputs

| WDL name | Type | Mode | Default | Scale | Validation | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `prepared_cpm` | `File` | Both | N/A | Linear CPM | Finite and strictly positive after GTF mapping and duplicate collapse | Prepared mapped expression matrix. |
| `prepared_log2_cpm` | `File` | Both | N/A | `log2(CPM)` | Exact row and sample order from preparation | Log-scale expression matrix for dtangle and TCA. |
| `mapping_report` | `File` | Both | N/A | N/A | One record per input expression row | GTF mapping and duplicate-resolution report. |
| `prepare_excluded_genes` | `File` | Both | N/A | N/A | Contains removed preparation genes and reasons | Expression rows excluded before preparation output. |
| `prepare_log` | `File` | Both | N/A | N/A | Task command log exists | Preparation stage log. |
| `estimated_proportions` | `File?` | dtangle | N/A | Relative proportion | Has sample IDs and 22 LM22 columns | Raw dtangle output. It is absent in restart mode. |
| `dtangle_markers` | `File?` | dtangle | N/A | N/A | Marker rows are from the dtangle fit | Marker genes and counts. |
| `dtangle_metadata` | `File?` | dtangle | N/A | N/A | JSON metadata exists | dtangle settings and fitted gamma information. |
| `dtangle_overlap_report` | `File?` | dtangle | N/A | N/A | Reports matched and missing LM22 genes | LM22 bulk-overlap report. |
| `transformed_lm22` | `File?` | dtangle | N/A | `log2(LM22)` | Positive linear LM22 was transformed once | LM22 matrix used by dtangle. |
| `dtangle_shared_bulk` | `File?` | dtangle | N/A | `log2(CPM)` | Genes are LM22-shared and ordered as LM22 | Bulk matrix limited to LM22-overlap genes. It is only a dtangle input artifact. |
| `dtangle_log` | `File?` | dtangle | N/A | N/A | Task command log exists | dtangle stage log. |
| `proportions_lm22` | `File` | Both | N/A | Relative proportion | Has the original 22 LM22 columns | Raw input proportions in a normalized output format. |
| `proportions_combined` | `File` | Both | N/A | Relative proportion | Has the ten combined groups before filtering | Combined lineage proportions. |
| `tca_weights` | `File` | Both | N/A | Relative proportion | Retained values are positive and each sample row sums to one | Combined, filtered, zero-adjusted proportions used by TCA. |
| `cell_group_filter_report` | `File` | Both | N/A | Relative proportion | Reports cohort means and retained status | Cohort filter results for all ten groups. |
| `proportions_log` | `File` | Both | N/A | N/A | Task command log exists | Proportion-processing stage log. |
| `tca_model` | `File` | Both | N/A | N/A | Serialized cohort-wide model exists | TCA model fitted with `refit_W = FALSE`. |
| `tca_model_log` | `File` | Both | N/A | N/A | TCA model log exists | Model-fit log from the TCA script. |
| `tca_expression` | `File` | Both | N/A | `log2(CPM)` | Mapped non-constant genes only; gene and sample order is stable | Full mapped matrix after constant-gene removal. It is not LM22-limited. |
| `tca_excluded_genes` | `File` | Both | N/A | N/A | Contains constant-gene removal reasons | Genes excluded before TCA fitting. |
| `gene_shard_manifest` | `File` | Both | N/A | N/A | Deterministic stable shard order | Gene-to-shard mapping. |
| `gene_shards` | `Array[File]` | Both | N/A | N/A | Each listed shard exists | Gene ID lists for tensor extraction. |
| `fit_tca_log` | `File` | Both | N/A | N/A | Task command log exists | TCA-fit stage log. |
| `tensor_shards` | `Array[File]` | Both | N/A | `log2_cpm` | Each shard HDF5 contains the assigned genes | Intermediate TCA tensor HDF5 files. |
| `tensor_shard_logs` | `Array[File]` | Both | N/A | N/A | One task log per tensor shard | Tensor-extraction logs. |
| `group_hdf5` | `Array[File]` | Both | N/A | `log2_cpm` | One HDF5 file per retained group; `expression` is genes by samples | Primary cell-group-specific inferred expression matrices. |
| `group_tsv` | `Array[File]` | Both | N/A | `log2_cpm` | Present only when `write_tsv = true` | Optional compressed TSV copies of group matrices. |
| `reconstruction_by_sample` | `File` | Both | N/A | `log2_cpm` | Contains per-sample reconstruction metrics | TCA reconstruction quality control. |
| `assembly_output_inventory` | `File` | Both | N/A | N/A | Lists outputs before manifest localization | Assembly output inventory. |
| `qc_summary` | `File` | Both | N/A | N/A | Quality-control summary exists | Assembly quality-control summary. |
| `qc_plots` | `File` | Both | N/A | N/A | Quality-control PDF exists | Reconstruction quality-control plots. |
| `assembly_log` | `File` | Both | N/A | N/A | Task command log exists | HDF5-assembly stage log. |
| `output_inventory` | `File` | Both | N/A | N/A | Localized paths resolve during manifest creation | Final provenance inventory. |
| `output_manifest` | `File` | Both | N/A | N/A | Machine-readable JSON exists | Manifest with output checksums and effective parameters. |
| `effective_parameters_file` | `File` | Both | N/A | N/A | JSON is written by the top-level workflow | Effective defaults and user settings. |
| `manifest_log` | `File` | Both | N/A | N/A | Task command log exists | Manifest-creation stage log. |

## HDF5 contract

Every file in `group_hdf5` stores `expression` with genes in rows and samples in columns. It stores `gene_id` and `sample_id` in the same order. Attributes record the cell-group name, `log2_cpm`, the pipeline version, and the TCA version. `tensor_shards` use the same log-scale expression convention, but they are intermediate shard files.
