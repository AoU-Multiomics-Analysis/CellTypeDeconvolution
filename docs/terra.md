# Run the workflow in Terra

Import `workflows/cell_type_deconvolution.wdl` into a Terra workspace. All WDL sources use WDL 1.0. Use [bed.inputs.json](../examples/bed.inputs.json) for dtangle mode. Use [precomputed-proportions.inputs.json](../examples/precomputed-proportions.inputs.json) for precomputed mode. Replace each `gs://YOUR_BUCKET/...` path with a path that your billing project can read.

## Required data

`expression` is a tab-delimited BED or BED.GZ matrix. Its leading columns are `#chr`, `start`, `end`, and `gene_id`, in that order. Remaining columns are samples. The BED is already normalized linear CPM. Values must be finite and strictly positive. The workflow stops on zero or negative values. It does not calculate CPM, renormalize CPM, or add a pseudocount.

`gtf` is a GTF or GTF.GZ file. It needs `gene` records with `gene_id` and `gene_name`. No gene-type filter is applied. The workflow trims surrounding white space before identifier matching. Matching is otherwise exact.

The workflow does not have an expression-preparation task. dtangle maps BED gene IDs with the GTF and applies `log2()` once. Its view uses gene symbols and sums linear CPM for duplicated symbols. TCA reads every valid nonconstant gene from the BED and applies `log2()` once. TCA does not use the GTF. Its view keeps each valid input `gene_id`, with its coordinates and row order. Rows that share a gene symbol stay separate.

## Select a proportion mode

Use dtangle mode when you set `lm22` and omit `precomputed_proportions`. LM22 is user-supplied. It must be a linear LM22 matrix with one gene-symbol column and the 22 standard cell-type columns. Values must be finite and strictly positive. The workflow applies `log2()` once.

Use precomputed mode when you set `precomputed_proportions` and omit `lm22`. The input needs one sample identifier column and all 22 standard LM22 proportion columns. Sample IDs must match the expression matrix.

Provide exactly one of `lm22` and `precomputed_proportions`. A validation task records the selected mode. The workflow stops before proportion-file selection if you provide both files or neither file.

Both modes require `expression` and `gtf`. `docker_image` is optional. The GitHub `latest` image is the default. Set `docker_image` to an immutable digest when your deployment requires one.

## Analysis settings

The workflow creates `log2(CPM)` without a pseudocount. In dtangle mode, it keeps LM22 gene order. dtangle uses RNA-seq mixtures, `marker_method = "ratio"`, a marker fraction of `0.10`, and a minimum LM22 overlap of `0.80`. Joint quantile normalization is off by default. If enabled, it changes only the joined dtangle reference and mixture data. It does not change the TCA matrix.

The workflow combines the 22 LM22 types into ten groups:

| Group | LM22 members |
| --- | --- |
| B cells | B cells naive; B cells memory; Plasma cells |
| CD4 T cells | T cells CD4 naive; T cells CD4 memory resting; T cells CD4 memory activated; T cells follicular helper; T cells regulatory (Tregs) |
| CD8 T cells | T cells CD8 |
| Gamma-delta T cells | T cells gamma delta |
| NK cells | NK cells resting; NK cells activated |
| Monocyte/myeloid | Monocytes; Macrophages M0; Macrophages M1; Macrophages M2 |
| Neutrophils | Neutrophils |
| Eosinophils | Eosinophils |
| Dendritic cells | Dendritic cells resting; Dendritic cells activated |
| Mast cells | Mast cells resting; Mast cells activated |

T follicular helper and regulatory T cells are in CD4 T cells. Gamma-delta T cells stay separate. A group is retained only when its cohort mean is at least `0.0001` (0.01%). No group is forced to remain. Retained exact zeros use `zero_floor`, then each sample row is normalized to one.

TCA fits one full-matrix model with `refit_W = FALSE`. It uses every valid, nonconstant gene from the direct BED in the complete `log2(CPM)` TCA matrix. It does not use only LM22-shared genes. Optional covariates are mixture-level technical covariates.

## Resources and outputs

Default resources support an initial whole-blood run of about 9,000 samples. TCA fitting uses 16 CPU, 192 GB memory, and 750 GB disk by default. Direct BED export uses 8 CPU, 128 GB memory, and 500 GB disk by default. `preemptible_attempts` and `max_retries` are global controls. They apply to every task. Change runtime inputs only after you inspect cohort size and storage needs.

The workflow does one direct tensor extraction. It writes one BED.GZ file per retained major cell group. Each file uses the `log2_cpm` model scale. It preserves `#chr`, `start`, `end`, `gene_id`, modeled-gene order, and sample-column order. It excludes genes that are constant before fitting. The values are statistical TCA estimates. They are not counts, linear CPM, or absolute expression values.

The workflow writes proportions, filter reports, the fitted TCA model, a BED inventory, reconstruction quality control, logs, an output inventory, and a machine-readable manifest. `cell_type_beds` is the authoritative file array. Public inventories and the manifest use stable BED basenames, so task-localized paths do not leak into metadata. The manifest task keeps localized paths only while it calculates checksums.

The final QC summary records LM22 value validation in dtangle mode. It also records maximum proportion row-sum errors, normalization adjustment, adjusted-zero counts, duplicate counts, constant-gene exclusions, reconstruction metrics, and TCA convergence. TCA 1.2.1 does not put a convergence field in its model object. The workflow derives the convergence status and iteration count from the TCA model log.

Migration note: the workflow does not produce HDF5 or tensor shards. Preparation outputs and a pipeline-version output are no longer available.

## Validation and support

GitHub Actions builds the pinned image. It runs R lint, R tests, WDL validation, command-logging checks, and dtangle and precomputed smoke workflows. The smoke workflows check the exact ordered 22 LM22 columns, BED coordinates, and sample order. A local Docker build is not required for a smoke test.

LM22 is not distributed with this repository. Review and comply with the LM22 source license and terms before use.

## Known limitations

- LM22 estimates relative immune proportions. It does not estimate absolute blood-cell counts.
- LM22 was derived from microarray data. RNA-seq and reference platform differences can affect estimates.
- LM22 does not model erythrocyte or platelet expression.
- Retained groups are cohort-specific because the mean filter uses the cohort.
- Very small proportions can be unstable, even after the zero floor is applied.
- TCA outputs are statistical estimates. They are not directly measured cell-type-specific expression values.
