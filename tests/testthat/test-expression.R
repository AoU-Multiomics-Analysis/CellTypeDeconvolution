source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(testthat::test_path("../..", "R", "expression.R"), local = .GlobalEnv)
source(testthat::test_path("../..", "R", "expression_bed.R"), local = .GlobalEnv)

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
    c(
      "duplicate_gene_name_aggregated_for_dtangle",
      "duplicate_gene_name_aggregated_for_dtangle"
    )
  )
})

testthat::test_that("duplicate aggregation preserves first gene-symbol order and exact values", {
  cpm <- matrix(
    c(
      1, 10,
      2, 20,
      4, 40,
      8, 80
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3", "g4"), c("S1", "S2"))
  )
  annotation <- tibble::tibble(
    gene_id = c("g4", "g2", "g1", "g3"),
    gene_name = c("C", "A", "B", "B"),
    gene_type = "protein_coding"
  )

  observed <- collapse_cpm_to_gene_names(cpm, annotation)

  testthat::expect_identical(rownames(observed), c("B", "A", "C"))
  testthat::expect_equal(
    unname(observed),
    matrix(c(5, 50, 2, 20, 8, 80), nrow = 3L, byrow = TRUE)
  )
})

testthat::test_that("direct CPM views separate TCA genes from dtangle symbols", {
  expression <- list(
    coordinates = tibble::tibble(
      `#chr` = c("chr1", "chr1", "chr2"),
      start = c(10L, 20L, 30L),
      end = c(11L, 21L, 31L),
      gene_id = c("ENSG1", "ENSG2", "ENSG3")
    ),
    cpm = matrix(
      c(4, 8, 16, 32, 64, 128),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("ENSG1", "ENSG2", "ENSG3"), c("S1", "S2"))
    )
  )
  annotation <- tibble::tibble(
    gene_id = c("ENSG1", "ENSG2"),
    gene_name = c("MARKER", "MARKER"),
    gene_type = c("protein_coding", "lncRNA")
  )

  tca_expression <- make_tca_expression(expression)
  dtangle_expression <- make_dtangle_expression(expression, annotation)

  testthat::expect_identical(rownames(tca_expression), c("ENSG1", "ENSG2", "ENSG3"))
  testthat::expect_equal(unname(tca_expression[, "S1"]), log2(c(4, 16, 64)))
  testthat::expect_identical(rownames(dtangle_expression$log_expression), "MARKER")
  testthat::expect_equal(unname(dtangle_expression$log_expression[1, ]), log2(c(20, 40)))
})

testthat::test_that("duplicate aggregation does not pivot the full expression matrix", {
  cpm <- matrix(
    c(1, 2, 3, 4),
    nrow = 2L,
    dimnames = list(c("g1", "g2"), c("S1", "S2"))
  )
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2"),
    gene_name = c("A", "A"),
    gene_type = "protein_coding"
  )
  trace(
    "pivot_longer",
    where = asNamespace("tidyr"),
    tracer = quote(stop("full-matrix pivot detected", call. = FALSE)),
    print = FALSE
  )
  on.exit(untrace("pivot_longer", where = asNamespace("tidyr")), add = TRUE)

  testthat::expect_no_error(collapse_cpm_to_gene_names(cpm, annotation))
})

testthat::test_that("BED preparation preserves gene IDs and collapses only dtangle symbols", {
  bed_path <- tempfile(fileext = ".bed.gz")
  input <- tibble::tibble(
    `#chr` = c("chr1", "chr2", "chr3"),
    start = c(0L, 100L, 200L),
    end = c(50L, 150L, 250L),
    gene_id = c("g1", "g2", "g3"),
    S1 = c(2, 3, 5),
    S2 = c(7, 11, 13)
  )
  readr::write_tsv(input, bed_path)
  annotation <- tibble::tibble(
    gene_id = c("g1", "g2", "g3"),
    gene_name = c("A", "A", "B"),
    gene_type = c("protein_coding", "lncRNA", NA_character_)
  )

  expression <- read_expression_bed(bed_path)
  result <- prepare_expression_bed(expression, annotation)

  testthat::expect_identical(rownames(result$tca_cpm), c("g1", "g2", "g3"))
  testthat::expect_identical(result$coordinates, input[1:4])
  testthat::expect_identical(rownames(result$dtangle_cpm), c("A", "B"))
  testthat::expect_equal(unname(result$dtangle_cpm["A", ]), c(5, 18))
  testthat::expect_equal(result$tca_log_expression, log2(result$tca_cpm))
  testthat::expect_equal(result$dtangle_log_expression, log2(result$dtangle_cpm))
})

testthat::test_that("BED validation rejects invalid coordinates and duplicate IDs", {
  invalid <- tibble::tibble(
    `#chr` = c("chr1", "chr1"), start = c(10L, 20L),
    end = c(10L, 30L), gene_id = c("g1", "g1"), S1 = c(1, 2)
  )
  path <- tempfile(fileext = ".bed")
  readr::write_tsv(invalid, path)
  testthat::expect_error(read_expression_bed(path), "start.*less than.*end")
})

testthat::test_that("BED validation rejects non-positive, missing, and nonfinite CPM values", {
  invalid_values <- c("0", "-1", "", "Inf")

  for (invalid_value in invalid_values) {
    path <- tempfile(fileext = ".bed")
    readr::write_tsv(
      tibble::tibble(
        `#chr` = "chr1", start = 0L, end = 10L, gene_id = "g1",
        S1 = invalid_value
      ),
      path,
      na = ""
    )

    testthat::expect_error(read_expression_bed(path), "CPM values")
  }
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

testthat::test_that("the BED CLI writes its seven declared outputs", {
  input <- tempfile(fileext = ".bed.gz")
  gtf <- tempfile(fileext = ".gtf")
  output_dir <- tempfile()
  dir.create(output_dir)
  readr::write_tsv(
    tibble::tibble(
      `#chr` = "chr1", start = 0L, end = 10L, gene_id = "g1", s1 = 2, s2 = 4
    ),
    input
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
    "prepared_tca_cpm.tsv.gz",
    "prepared_tca_log2_cpm.tsv.gz",
    "prepared_dtangle_cpm.tsv.gz",
    "prepared_dtangle_log2_cpm.tsv.gz",
    "prepared_coordinates.tsv",
    "gene_mapping_report.tsv",
    "excluded_genes.tsv"
  )

  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, expected))))
})
