source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(testthat::test_path("../..", "R", "expression.R"), local = .GlobalEnv)

testthat::test_that("counts convert to TPM by gene length", {
  counts <- matrix(
    c(100, 100, 200, 100),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )

  observed <- counts_to_tpm(counts, c(g1 = 1000, g2 = 2000))

  testthat::expect_equal(colSums(observed), c(s1 = 1e6, s2 = 1e6))
  testthat::expect_equal(
    observed[, "s1"],
    c(g1 = 2 / 3 * 1e6, g2 = 1 / 3 * 1e6)
  )
})

testthat::test_that("count conversion rejects negative values", {
  counts <- matrix(
    c(-1, 2),
    nrow = 2,
    dimnames = list(c("g1", "g2"), "s1")
  )

  testthat::expect_error(
    counts_to_tpm(counts, c(g1 = 1000, g2 = 1000)),
    "nonnegative"
  )
})

testthat::test_that("duplicate symbols are summed and renormalized", {
  values <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("ENSG1", "ENSG2"), c("s1", "s2"))
  )
  annotation <- tibble::tibble(
    gene_id = c("ENSG1", "ENSG2"),
    gene_symbol = c("A", "A"),
    gene_length_bp = c(1000, 1000)
  )

  observed <- collapse_to_symbols(values, annotation)

  testthat::expect_identical(rownames(observed), "A")
  testthat::expect_equal(as.numeric(observed), c(3, 7))
})

testthat::test_that("log expression is log2 TPM plus one", {
  x <- matrix(
    c(0, 1e6),
    nrow = 2,
    dimnames = list(c("A", "B"), "s1")
  )

  result <- prepare_expression(x, "tpm")

  testthat::expect_equal(result$log_expression, log2(result$tpm + 1))
})

testthat::test_that("expression preparation rejects negative values", {
  x <- matrix(
    c(-1, 1),
    nrow = 2,
    dimnames = list(c("A", "B"), "s1")
  )

  testthat::expect_error(prepare_expression(x, "tpm"), "nonnegative")
})

testthat::test_that("count input requires a gene length for every gene", {
  counts <- matrix(
    c(1, 2),
    nrow = 2,
    dimnames = list(c("g1", "g2"), "s1")
  )
  annotation <- tibble::tibble(
    gene_id = "g1",
    gene_symbol = "G1",
    gene_length_bp = 1000
  )

  testthat::expect_error(
    prepare_expression(counts, "counts", annotation),
    "gene_length_bp.*g2"
  )
})

testthat::test_that("count input rejects nonpositive gene lengths", {
  counts <- matrix(1, nrow = 1, dimnames = list("g1", "s1"))
  annotation <- tibble::tibble(
    gene_id = "g1",
    gene_symbol = "G1",
    gene_length_bp = 0
  )

  testthat::expect_error(
    prepare_expression(counts, "counts", annotation),
    "positive.*gene_length_bp"
  )
})

testthat::test_that("expression preparation rejects an invalid expression type", {
  x <- matrix(1, nrow = 1, dimnames = list("A", "s1"))

  testthat::expect_error(prepare_expression(x, "fpkm"), "counts.*tpm")
})

testthat::test_that("expression preparation rejects duplicated sample identifiers", {
  x <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list("A", c("s1", "s1"))
  )

  testthat::expect_error(prepare_expression(x, "tpm"), "sample.*unique")
})

testthat::test_that("whitespace-only symbols are recorded as unmapped", {
  x <- matrix(
    c(1, 2, 2, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2"),
    gene_symbol = c("G1", " ")
  )

  result <- prepare_expression(x, "tpm", annotation)

  testthat::expect_identical(rownames(result$tpm), "G1")
  testthat::expect_identical(
    result$mapping_report$mapping_action,
    c("mapped", "unmapped")
  )
  testthat::expect_true(any(result$excluded_genes$reason == "unmapped_gene_symbol"))
})

testthat::test_that("annotation-missing TPM genes are recorded as unmapped", {
  x <- matrix(
    c(1, 2, 2, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2"))
  )
  annotation <- tibble::tibble(gene_id = "g1", gene_symbol = "G1")

  result <- prepare_expression(x, "tpm", annotation)

  testthat::expect_identical(rownames(result$tpm), "G1")
  testthat::expect_identical(
    result$mapping_report$mapping_action,
    c("mapped", "unmapped")
  )
  testthat::expect_true(any(
    result$excluded_genes$gene_id == "g2" &
      result$excluded_genes$reason == "unmapped_gene_symbol"
  ))
})

testthat::test_that("the expression CLI writes its four declared outputs", {
  input <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()
  dir.create(output_dir)
  x <- matrix(
    c(0, 1e6),
    nrow = 2,
    dimnames = list(c("A", "B"), "s1")
  )
  write_numeric_matrix(x, input, "gene_id")
  original_working_directory <- setwd(pipeline_root)
  on.exit(setwd(original_working_directory), add = TRUE)

  status <- system2(
    "Rscript",
    c(
      "scripts/prepare_expression.R", "--expression", input,
      "--expression-type", "tpm", "--output-dir", output_dir
    )
  )

  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, c(
    "prepared_tpm.tsv.gz", "prepared_log2_tpm_plus_1.tsv.gz",
    "gene_mapping_report.tsv", "excluded_genes.tsv"
  )))))
})
