# Cell type deconvolution for human whole blood

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/docker-image.yml)
[![Pipeline CI](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/pipeline-ci.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/pipeline-ci.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/CellTypeDeconvolution/actions/workflows/update-readme-badges.yml)
<!-- workflow-badges:end -->

This Terra workflow estimates relative immune-cell proportions with dtangle. It then uses Tensor Composition Analysis (TCA) to estimate expression for each retained major cell group. It supports human whole-blood RNA-seq cohorts. Default resources support about 9,000 samples.

## Input data

Provide `expression` as a tab-delimited BED or BED.GZ matrix. The leading columns must be `#chr`, `start`, `end`, and `gene_id`, in that order. The remaining columns are samples. Each value must be finite, strictly positive, linear CPM. The workflow does not calculate CPM, renormalize CPM, or add a pseudocount.

Provide `gtf` as a GTF or GTF.GZ file. It must contain `gene` records with `gene_id` and `gene_name`. The workflow does not apply a gene-type filter.

The workflow trims surrounding white space before it matches identifiers. Matching is otherwise exact. The TCA view keeps one mapped row for each input `gene_id`, with its BED coordinates and order. It does not collapse duplicate gene symbols. The dtangle view uses gene symbols and sums linear CPM for duplicate symbols. This keeps the TCA coordinate identity and gives dtangle one value per symbol.

## Proportion modes

Use dtangle mode when you provide a user-supplied linear LM22 matrix. It must contain one gene-symbol column and the 22 standard LM22 columns. Values must be finite and strictly positive. The workflow applies `log2()` once. dtangle uses RNA-seq mixtures, `marker_method = "ratio"`, and a marker fraction of `0.10`.

Use precomputed mode when you provide `precomputed_proportions` and omit `lm22`. The table must have one sample identifier column and all 22 standard LM22 columns. Sample IDs must match the prepared expression matrix.

Both modes combine the 22 LM22 types into ten groups: B cells, CD4 T cells, CD8 T cells, Gamma-delta T cells, NK cells, Monocyte/myeloid, Neutrophils, Eosinophils, Dendritic cells, and Mast cells. T follicular helper and regulatory T cells are in CD4 T cells. Gamma-delta T cells stay separate. A group is retained only when its cohort mean is at least `0.0001`. Retained exact zeros use `zero_floor` before each sample row is normalized to one.

## TCA and BED outputs

TCA fits one full-matrix model across the cohort. It uses mapped, nonconstant genes from the TCA `log2(CPM)` view. It does not limit the model to LM22-shared genes. It uses processed weights with `refit_W = FALSE`.

The workflow performs one direct tensor extraction. It writes one BED.GZ file per retained major cell group. Each file uses the `log2_cpm` model scale. These values are inferred log-expression estimates. They are not counts, linear CPM, or absolute cell-type expression values.

Each BED.GZ output preserves `#chr`, `start`, `end`, `gene_id`, modeled-gene order, and sample-column order. It excludes unmapped genes and genes that are constant before fitting. `cell_type_bed_inventory` maps each retained group to its BED.GZ path, dimensions, and scale.

The export task defaults to 8 CPU, 128 GB memory, 500 GB disk, 0 preemptible attempts, and 1 retry. These defaults support an initial cohort of about 9,000 samples. See [the data dictionary](docs/data-dictionary.md) for every input and output.

## Run and validation

Import `workflows/cell_type_deconvolution.wdl` into Terra. Start with [bed.inputs.json](examples/bed.inputs.json) for dtangle mode or [precomputed-proportions.inputs.json](examples/precomputed-proportions.inputs.json) for precomputed mode. Replace each example cloud path with a readable path. Use a published image tag for production.

The workflow writes a manifest, an output inventory, quality-control files, and task logs. GitHub Actions builds the pinned image and runs R tests, R lint, WDL checks, logging checks, and dtangle and precomputed smoke workflows. The smoke workflows check the BED coordinate and sample-order contract. A local Docker build is not required for a smoke test.

Migration note: the workflow does not produce HDF5 or tensor shards.

## LM22 license

LM22 is user-supplied. Review and comply with its source license and terms before use. This repository does not redistribute LM22.

## Known limitations

- LM22 estimates relative immune proportions. It does not estimate absolute blood-cell counts.
- LM22 was derived from microarray data. Platform differences can affect RNA-seq estimates.
- LM22 does not model erythrocyte or platelet expression.
- Retained groups are cohort-specific because the filter uses the cohort mean.
- Very small proportions can be unstable, even after the zero floor is applied.
- TCA outputs are statistical estimates. They are not directly measured cell-type-specific expression values.
