# Direct CPM WDL 1.0 Design

## Goal

Simplify the whole-blood cell-type deconvolution workflow. The workflow will consume one coordinate-preserving CPM BED directly. It will not run a separate expression-preparation task. All WDL files will use WDL 1.0.

## Input contract

The workflow requires an expression BED or BED.GZ file with these columns:

1. `#chr`
2. `start`
3. `end`
4. `gene_id`
5. One or more sample columns

Expression values must be finite, strictly positive, linear CPM values. Gene IDs and sample IDs must be nonempty and unique. The workflow preserves input gene order, sample order, coordinates, and gene IDs.

The GTF remains a required input. In dtangle mode, it maps expression `gene_id` values to LM22 gene symbols. No gene-type filter applies. TCA does not require a successful GTF mapping and uses every valid, nonconstant expression gene.

The proportion-mode contract does not change. The user supplies exactly one of `lm22` or `precomputed_proportions`.

## Workflow structure

The workflow removes `PrepareExpression` and the import of `tasks/expression.wdl`.

### Proportion validation

`ValidateProportionMode` continues to reject runs that supply both proportion inputs or neither input.

### dtangle

`RunDtangle` receives the original CPM BED and GTF. Its R command performs these operations:

1. Validate the BED schema and positive linear CPM values.
2. Read gene mappings from the GTF.
3. Map Ensembl gene IDs to gene symbols.
4. Sum linear CPM rows that map to the same gene symbol.
5. Apply `log2()` once, without a pseudocount.
6. Intersect the bulk matrix with LM22 and run dtangle.

The task retains the current marker, overlap, transformed-reference, shared-bulk, metadata, and proportion outputs. Its overlap output supplies the relevant gene-mapping audit information.

### TCA fitting

`FitTca` receives the original CPM BED. Its R command removes the four BED metadata columns, validates the matrix, applies `log2()` once without a pseudocount, removes constant genes using the existing policy, aligns samples with the processed proportions, and fits one full TCA model.

TCA uses gene IDs as row identifiers. It does not restrict the model to genes that map through the GTF or overlap LM22.

### BED export

`ExportTcaBeds` receives the original CPM BED instead of a separate coordinates file. It reads the four metadata columns, aligns them with the genes retained by TCA, and writes one BED.GZ file for each retained major cell group.

The output BED files remain in log2-CPM model space. They preserve input coordinates, gene IDs, modeled-gene order, and sample order.

### Manifest

`BuildManifest` no longer receives preparation mapping or exclusion files. It continues to include dtangle metadata when dtangle mode runs, processed proportions, TCA information, reconstruction metrics, checksums, parameters, and the container image.

The workflow removes the `pipeline_version` input and the corresponding manifest field. WDL source files declare `version 1.0`; this is not a runtime input.

## Runtime configuration

Every workflow and task WDL file declares `version 1.0`.

The workflow defines this default image:

`ghcr.io/aou-multiomics-analysis/celltypedeconvolution:latest`

The workflow exposes one `preemptible_attempts` input with a default of `2` and one `max_retries` input with a default of `2`. It passes both values to every task. Task-specific preemptible and retry inputs are removed from the public workflow interface.

Task-specific CPU, memory, and disk inputs remain because the TCA tasks need more resources than the validation and proportion tasks.

## Public outputs

The workflow removes these preparation outputs:

- Prepared TCA CPM
- Prepared TCA log2 CPM
- Prepared dtangle CPM
- Prepared dtangle log2 CPM
- Prepared coordinates
- Preparation mapping report
- Preparation exclusion report
- Preparation log

The workflow retains proportion outputs, TCA model outputs, cell-type BED files, quality-control files, inventories, logs, and the output manifest.

## Validation and errors

The workflow stops with a direct error when:

- BED metadata columns are missing or out of order.
- A gene ID or sample ID is missing or duplicated.
- An expression value is missing, nonnumeric, nonfinite, zero, or negative.
- The GTF cannot map enough genes to LM22 in dtangle mode.
- Expression and proportion sample IDs differ.
- A required cell-group or TCA invariant fails.

Each WDL command continues to write start, dimension, output, completion, and failure messages.

## Tests

Tests will cover these requirements:

- Every active WDL file declares `version 1.0`.
- The workflow has no `PrepareExpression` call or preparation outputs.
- The Docker image, preemptible attempts, and retry limit have one workflow-level default each.
- dtangle maps gene IDs through the GTF and transforms CPM once.
- TCA consumes the direct BED and uses all valid nonconstant genes.
- BED export reads coordinates from the direct BED and preserves ordering.
- Both dtangle and precomputed-proportion smoke workflows run through GitHub Actions.
- The public inventory comparison checks schema and values without comparing parser metadata.

The standard smoke test will use GitHub Actions. A local Docker build is not required.

## Documentation migration

The README, Terra guide, data dictionary, example inputs, smoke fixtures, and active contract tests will describe the new direct-input workflow. Historical design and plan documents remain unchanged as records of earlier decisions.
