# Cell type deconvolution for human whole blood

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/docker-image.yml)
[![Pipeline CI](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/pipeline-ci.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/pipeline-ci.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/update-readme-badges.yml)
<!-- workflow-badges:end -->

This Terra workflow estimates relative immune-cell proportions with dtangle. It then uses Tensor Composition Analysis (TCA) to estimate expression for each retained major cell group. It supports human whole-blood RNA-seq cohorts. All WDL sources use WDL 1.0. Default resources support about 9,000 samples.

## Input data

Provide `expression` as a tab-delimited BED or BED.GZ matrix. The leading columns must be `#chr`, `start`, `end`, and `gene_id`, in that order. The remaining columns are samples. The BED is already normalized linear CPM. Each value must be finite and strictly positive. The workflow does not calculate CPM, renormalize CPM, or add a pseudocount.

Provide `gtf` as a GTF or GTF.GZ file. It must contain `gene` records with `gene_id` and `gene_name`. The workflow does not apply a gene-type filter.

The workflow trims surrounding white space before it matches identifiers. Matching is otherwise exact. The TCA view keeps one mapped row for each input `gene_id`, with its BED coordinates and order. It does not collapse duplicate gene symbols. The dtangle view uses gene symbols and sums linear CPM for duplicate symbols. This keeps the TCA coordinate identity and gives dtangle one value per symbol.

## Proportion modes

Use dtangle mode when you provide a user-supplied linear LM22 matrix. It must contain one gene-symbol column and the 22 standard LM22 columns. Values must be finite and strictly positive. dtangle maps BED gene IDs with the GTF and applies `log2()` once. dtangle uses RNA-seq mixtures, `marker_method = "ratio"`, and a marker fraction of `0.10`.

Use precomputed mode when you provide `precomputed_proportions` and omit `lm22`. The table must have one sample identifier column and all 22 standard LM22 columns. Sample IDs must match the expression matrix.

Provide exactly one of `lm22` and `precomputed_proportions`. A validation task records the selected mode before the workflow selects the proportion file. The workflow stops if you provide both files or neither file.

Both modes combine the 22 LM22 types into ten groups: B cells, CD4 T cells, CD8 T cells, Gamma-delta T cells, NK cells, Monocyte/myeloid, Neutrophils, Eosinophils, Dendritic cells, and Mast cells. T follicular helper and regulatory T cells are in CD4 T cells. Gamma-delta T cells stay separate. A group is retained only when its cohort mean is at least `0.0001`. Retained exact zeros use `zero_floor` before each sample row is normalized to one.

## TCA and BED outputs

The workflow does not have an expression-preparation task. TCA reads every valid nonconstant gene from the BED and applies `log2()` once. It fits one full-matrix model across the cohort. It does not limit the model to LM22-shared genes. It uses processed weights with `refit_W = FALSE`.

The workflow performs one direct tensor extraction. It writes one BED.GZ file per retained major cell group. Each file uses the `log2_cpm` model scale. These values are inferred log-expression estimates. They are not counts, linear CPM, or absolute cell-type expression values.

Each BED.GZ output preserves `#chr`, `start`, `end`, `gene_id`, modeled-gene order, and sample-column order. It excludes unmapped genes and genes that are constant before fitting. `cell_type_beds` is the authoritative array of files. `cell_type_bed_inventory` and `output_inventory` use stable BED basenames. The output manifest also uses these basenames. The manifest task uses localized paths only for checksum calculation.

The GitHub `latest` image is the default. `docker_image` is optional. `preemptible_attempts` and `max_retries` are global controls that apply to every task. Task resource inputs are optional. See [the data dictionary](docs/data-dictionary.md) for every input and output.

## Run and validation

Import `workflows/cell_type_deconvolution.wdl` into Terra. Start with [bed.inputs.json](examples/bed.inputs.json) for dtangle mode or [precomputed-proportions.inputs.json](examples/precomputed-proportions.inputs.json) for precomputed mode. Replace each example cloud path with a readable path. The GitHub `latest` image is the default. You can set `docker_image` to an immutable digest when your deployment requires one. Set only the global `preemptible_attempts` and `max_retries` retry controls.

The workflow writes a manifest, an output inventory, quality-control files, and task logs. The final QC summary includes LM22 validation when dtangle mode is active, proportion row-sum errors, normalization adjustment, duplicate counts, constant-gene exclusions, reconstruction metrics, and TCA convergence. TCA 1.2.1 does not return a convergence field in the model. The workflow derives the convergence status and iteration count from the TCA model log. GitHub Actions builds the pinned image and runs R tests, R lint, WDL checks, logging checks, and dtangle and precomputed smoke workflows. The smoke workflows check the exact ordered 22 LM22 columns and the BED coordinate and sample-order contract. A local Docker build is not required for a smoke test.

Migration note: the workflow does not produce HDF5 or tensor shards. Preparation outputs and a pipeline-version output are no longer available.

## LM22 license

LM22 is user-supplied. Review and comply with its source license and terms before use. This repository does not redistribute LM22.

## Known limitations

- LM22 estimates relative immune proportions. It does not estimate absolute blood-cell counts.
- LM22 was derived from microarray data. Platform differences can affect RNA-seq estimates.
- LM22 does not model erythrocyte or platelet expression.
- Retained groups are cohort-specific because the filter uses the cohort mean.
- Very small proportions can be unstable, even after the zero floor is applied.
- TCA outputs are statistical estimates. They are not directly measured cell-type-specific expression values.
