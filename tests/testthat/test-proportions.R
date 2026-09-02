source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

proportions_path <- testthat::test_path("../..", "R", "proportions.R")
if (file.exists(proportions_path)) {
  source(proportions_path, local = .GlobalEnv)
}

make_lm22_proportions <- function(sample_ids = c("S1", "S2")) {
  proportions <- matrix(
    rep(seq_len(22L), times = length(sample_ids)),
    nrow = length(sample_ids),
    byrow = TRUE,
    dimnames = list(sample_ids, lm22_cell_types())
  )
  proportions / rowSums(proportions)
}

testthat::test_that("the ten groups partition LM22 exactly once", {
  members <- unlist(lm22_group_map(), use.names = FALSE)

  testthat::expect_length(lm22_group_map(), 10L)
  testthat::expect_setequal(members, lm22_cell_types())
  testthat::expect_false(anyDuplicated(members) > 0L)
  testthat::expect_identical(lm22_group_map()[["Gamma-delta T cells"]], "T cells gamma delta")
  testthat::expect_identical(
    lm22_group_map()[["Mast cells"]],
    c("Mast cells resting", "Mast cells activated")
  )
})

testthat::test_that("proportions combine LM22 types in group order", {
  proportions <- make_lm22_proportions()

  combined <- combine_lm22_proportions(proportions)

  testthat::expect_identical(rownames(combined), c("S1", "S2"))
  testthat::expect_identical(colnames(combined), c(
    "B cells", "CD4 T cells", "CD8 T cells", "Gamma-delta T cells",
    "NK cells", "Monocyte/myeloid", "Neutrophils", "Eosinophils",
    "Dendritic cells", "Mast cells"
  ))
  testthat::expect_equal(
    unname(combined[1L, ]),
    c(6, 35, 4, 10, 23, 58, 22, 21, 35, 39) / 253
  )
})

testthat::test_that("filter and zero adjustment follow the cohort rule", {
  combined <- matrix(
    c(0.7, 0.29995, 0.00005, 0.6, 0.4, 0), nrow = 2, byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("B cells", "CD4 T cells", "Eosinophils"))
  )

  result <- filter_and_adjust_groups(combined, 0.0001, 1e-6)

  testthat::expect_false("Eosinophils" %in% colnames(result$weights))
  testthat::expect_true(all(result$weights > 0))
  testthat::expect_equal(rowSums(result$weights), c(S1 = 1, S2 = 1))
  testthat::expect_identical(names(result$report), c(
    "cell_group", "cohort_mean", "threshold", "retained", "filter_reason",
    "zero_count_before", "zero_floor"
  ))
  testthat::expect_identical(
    unname(result$report$filter_reason),
    c("retained", "retained", "below_threshold")
  )
})

testthat::test_that("22-type inputs require exact columns and valid proportions", {
  valid <- matrix(
    1 / 22,
    nrow = 2,
    ncol = 22,
    dimnames = list(c("S1", "S2"), lm22_cell_types())
  )
  testthat::expect_error(
    process_proportions(valid[, -1, drop = FALSE], 0.0001, 1e-6),
    "LM22 columns"
  )
  extra <- cbind(valid, unknown = 0)
  testthat::expect_error(process_proportions(extra, 0.0001, 1e-6), "LM22 columns")
  negative <- valid
  negative[1L, 1L] <- -1
  testthat::expect_error(process_proportions(negative, 0.0001, 1e-6), "nonnegative")
  nonfinite <- valid
  nonfinite[1L, 1L] <- Inf
  testthat::expect_error(process_proportions(nonfinite, 0.0001, 1e-6), "finite")
  bad_sum <- valid
  bad_sum[1L, 1L] <- bad_sum[1L, 1L] + 1e-4
  testthat::expect_error(process_proportions(bad_sum, 0.0001, 1e-6), "sum to one")
})

testthat::test_that("restart row sums are accepted within one millionth", {
  valid <- make_lm22_proportions("S1")
  valid[1L, 1L] <- valid[1L, 1L] + 1e-6

  result <- process_proportions(valid, 0.0001, 1e-6)

  testthat::expect_equal(rowSums(result$tca_weights), c(S1 = 1))
})

testthat::test_that("no lineage is mandatory and TCA still needs two", {
  combined <- matrix(
    c(0.99995, 0.00005),
    nrow = 1,
    dimnames = list("S1", c("CD4 T cells", "Mast cells"))
  )

  testthat::expect_error(filter_and_adjust_groups(combined, 0.01, 1e-6), "two")
})

testthat::test_that("the proportion CLI writes its four declared outputs", {
  proportions <- make_lm22_proportions()
  proportions_path <- tempfile(fileext = ".tsv")
  output_dir <- tempfile()
  dir.create(output_dir)
  write_numeric_matrix(proportions, proportions_path, "sample_id")
  original_working_directory <- setwd(pipeline_root)
  on.exit(setwd(original_working_directory), add = TRUE)

  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "scripts/process_proportions.R", "--proportions", proportions_path,
      "--output-dir", output_dir
    )
  )

  testthat::expect_equal(status, 0L)
  testthat::expect_true(all(file.exists(file.path(output_dir, c(
    "proportions_lm22.tsv", "proportions_combined.tsv",
    "proportions_tca_weights.tsv", "cell_group_filter_report.tsv"
  )))))
})
