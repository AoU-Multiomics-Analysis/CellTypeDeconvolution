source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(testthat::test_path("../..", "R", "expression.R"), local = .GlobalEnv)

testthat::test_that("GTF parsing retains all gene types and ignores non-gene records", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(
    "1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\"; gene_type \"protein_coding\";",
    "1\tsrc\tgene\t20\t30\t.\t+\t.\tgene_id \"g2\"; gene_name \"B\"; gene_type \"lncRNA\";",
    "1\tsrc\ttranscript\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\";"
  ), gtf)

  observed <- read_gtf_gene_annotation(gtf, chunk_size = 1L)

  testthat::expect_identical(observed$gene_id, c("g1", "g2"))
  testthat::expect_identical(observed$gene_name, c("A", "B"))
  testthat::expect_identical(observed$gene_type, c("protein_coding", "lncRNA"))
})

testthat::test_that("duplicate gene names are summed without CPM renormalization", {
  cpm <- matrix(c(2, 3, 5, 7), nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1", "s2")))
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2"), gene_name = c("A", "A"), gene_type = c("x", "y")
  )

  result <- prepare_expression(cpm, annotation)

  testthat::expect_identical(rownames(result$cpm), "A")
  testthat::expect_equal(unname(result$cpm), matrix(c(5, 12), nrow = 1))
  testthat::expect_equal(result$log_expression, log2(result$cpm))
  testthat::expect_identical(
    result$mapping_report$mapping_action,
    c("duplicate_gene_name_collapsed", "duplicate_gene_name_collapsed")
  )
})

testthat::test_that("prepared CPM must be strictly positive", {
  cpm <- matrix(c(1, 0), nrow = 1,
    dimnames = list("g1", c("s1", "s2")))
  annotation <- tibble::tibble(
    gene_id = "g1", gene_name = "A", gene_type = NA_character_
  )

  testthat::expect_error(prepare_expression(cpm, annotation), "strictly positive")
})

testthat::test_that("compressed GTF and trimmed IDs are supported", {
  gtf <- tempfile(fileext = ".gtf.gz")
  connection <- gzfile(gtf, "wt")
  writeLines(
    "1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\"; gene_type \"lncRNA\";",
    connection
  )
  close(connection)
  cpm <- matrix(c(2, 4), nrow = 1,
    dimnames = list(" g1 ", c("s1", "s2")))

  result <- prepare_expression(cpm, read_gtf_gene_annotation(gtf))

  testthat::expect_identical(rownames(result$cpm), "A")
  testthat::expect_equal(as.numeric(result$cpm), c(2, 4))
})

testthat::test_that("missing gene names and GTF IDs are reported", {
  cpm <- matrix(1:6, nrow = 3,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2"), gene_name = c("A", NA_character_),
    gene_type = c("protein_coding", "lncRNA")
  )

  result <- prepare_expression(cpm, annotation)

  testthat::expect_identical(
    result$mapping_report$mapping_action,
    c("mapped", "missing_gene_name", "missing_gtf_gene_id")
  )
  testthat::expect_true(any(result$excluded_genes$gene_id == "g2"))
  testthat::expect_true(any(result$excluded_genes$gene_id == "g3"))
})

testthat::test_that("duplicate GTF IDs stop validation", {
  annotation <- tibble::tibble(
    gene_id = c("g1", "g1"), gene_name = c("A", "A"),
    gene_type = c("protein_coding", "protein_coding")
  )

  testthat::expect_error(validate_gtf_gene_annotation(annotation), "gene_id.*unique")
})

testthat::test_that("malformed GTF records stop parsing", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines("1\tsrc\tgene\t1\t10", gtf)

  testthat::expect_error(read_gtf_gene_annotation(gtf), "nine.*fields")
})

testthat::test_that("the CPM CLI writes its four declared outputs", {
  input <- tempfile(fileext = ".tsv")
  gtf <- tempfile(fileext = ".gtf")
  output_dir <- tempfile()
  dir.create(output_dir)
  write_numeric_matrix(
    matrix(c(2, 4), nrow = 1, dimnames = list("g1", c("s1", "s2"))),
    input,
    "gene_id"
  )
  writeLines(
    "1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\"; gene_name \"A\";",
    gtf
  )
  original_working_directory <- setwd(pipeline_root)
  on.exit(setwd(original_working_directory), add = TRUE)

  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("scripts/prepare_expression.R", "--expression", shQuote(input),
      "--gtf", shQuote(gtf), "--output-dir", shQuote(output_dir))
  )
  expected <- c(
    "prepared_cpm.tsv.gz", "prepared_log2_cpm.tsv.gz",
    "gene_mapping_report.tsv", "excluded_genes.tsv"
  )

  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, expected))))
})
