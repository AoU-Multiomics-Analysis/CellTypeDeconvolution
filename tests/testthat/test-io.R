source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

testthat::test_that("the group map partitions all LM22 types once", {
  members <- unlist(lm22_group_map(), use.names = FALSE)
  testthat::expect_setequal(members, lm22_cell_types())
  testthat::expect_length(members, 22L)
  testthat::expect_false(anyDuplicated(members) > 0L)
})

testthat::test_that("numeric matrix round-trip preserves IDs and values", {
  x <- matrix(c(1, 2, 3, 4), nrow = 2,
              dimnames = list(c("G1", "G2"), c("S1", "S2")))
  path <- tempfile(fileext = ".tsv.gz")
  write_numeric_matrix(x, path, "gene_symbol")
  observed <- read_numeric_matrix(path, "gene_symbol")
  testthat::expect_identical(dimnames(observed), dimnames(x))
  testthat::expect_equal(observed, x)
})

testthat::test_that("sample mismatch reports the first difference", {
  testthat::expect_error(
    assert_identical_ids(c("S1", "S2"), c("S2", "S1"), "proportions"),
    "proportions.*S1.*S2"
  )
})
