source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

tca_stage_path <- testthat::test_path("../..", "R", "tca_stage.R")
if (file.exists(tca_stage_path)) {
  source(tca_stage_path, local = .GlobalEnv)
}

make_tca_inputs <- function() {
  X <- matrix(
    seq_len(24L),
    nrow = 4L,
    dimnames = list(paste0("G", 1:4), paste0("S", 1:6))
  )
  W <- matrix(
    rep(c(0.25, 0.75), times = 6L),
    nrow = 6L,
    byrow = TRUE,
    dimnames = list(colnames(X), c("A", "B"))
  )
  list(X = X, W = W)
}

testthat::test_that("TCA input keeps non-LM22 genes", {
  X <- matrix(
    seq_len(50L),
    nrow = 5L,
    dimnames = list(
      c("LM1", "LM2", "EXTRA1", "EXTRA2", "EXTRA3"),
      paste0("S", 1:10)
    )
  )

  result <- remove_constant_features(X)

  testthat::expect_identical(rownames(result$matrix), rownames(X))
})

testthat::test_that("TCA rejects reordered samples", {
  inputs <- make_tca_inputs()
  rownames(inputs$W) <- rev(rownames(inputs$W))

  testthat::expect_error(
    validate_tca_inputs(inputs$X, inputs$W),
    "sample order"
  )
})

testthat::test_that("TCA validates finite positive normalized inputs", {
  inputs <- make_tca_inputs()
  testthat::expect_true(validate_tca_inputs(inputs$X, inputs$W))

  nonfinite_X <- inputs$X
  nonfinite_X[[1L]] <- Inf
  testthat::expect_error(
    validate_tca_inputs(nonfinite_X, inputs$W),
    "expression.*finite"
  )

  zero_W <- inputs$W
  zero_W[1L, ] <- c(0, 1)
  testthat::expect_error(
    validate_tca_inputs(inputs$X, zero_W),
    "strictly positive"
  )

  bad_sum_W <- inputs$W
  bad_sum_W[1L, 1L] <- bad_sum_W[1L, 1L] + 2e-8
  testthat::expect_error(
    validate_tca_inputs(inputs$X, bad_sum_W),
    "sum to one within 1e-8"
  )

  one_group_W <- matrix(
    1,
    nrow = ncol(inputs$X),
    dimnames = list(colnames(inputs$X), "A")
  )
  testthat::expect_error(
    validate_tca_inputs(inputs$X, one_group_W),
    "at least two groups"
  )
})

testthat::test_that("TCA validates covariate rows, values, and intercepts", {
  inputs <- make_tca_inputs()
  C2 <- matrix(
    seq_len(6L),
    ncol = 1L,
    dimnames = list(colnames(inputs$X), "batch")
  )
  testthat::expect_true(validate_tca_inputs(inputs$X, inputs$W, C2))

  reordered_C2 <- C2[rev(rownames(C2)), , drop = FALSE]
  testthat::expect_error(
    validate_tca_inputs(inputs$X, inputs$W, reordered_C2),
    "covariate sample order"
  )

  missing_C2 <- C2
  missing_C2[[1L]] <- NA_real_
  testthat::expect_error(
    validate_tca_inputs(inputs$X, inputs$W, missing_C2),
    "missing covariates"
  )

  intercept_C2 <- cbind(C2, "(Intercept)" = 1)
  testthat::expect_error(
    validate_tca_inputs(inputs$X, inputs$W, intercept_C2),
    "intercept"
  )
})

testthat::test_that("constant-gene removal reports gene_id and preserves order", {
  X <- matrix(
    c(1, 1, 2, 3, 4, 6),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("S1", "S2"))
  )
  result <- remove_constant_features(X)
  testthat::expect_identical(rownames(result$matrix), c("g2", "g3"))
  testthat::expect_identical(result$report$gene_id, "g1")
  testthat::expect_identical(result$report$reason, "constant_expression")
})

testthat::test_that("the TCA CLI does not create shard artifacts", {
  text <- paste(readLines(
    testthat::test_path("../..", "scripts", "fit_tca.R"),
    warn = FALSE
  ), collapse = "\n")
  testthat::expect_false(grepl("--shard-size", text, fixed = TRUE))
})

testthat::test_that("one model fits all genes without refitting weights", {
  testthat::skip_if_not_installed("TCA")
  testthat::expect_identical(as.character(utils::packageVersion("TCA")), "1.2.1")
  set.seed(20260901)
  data <- TCA::test_data(24, 20, 3, 0, 0, 0.01)

  result <- fit_tca_stage(
    X = data$X,
    W = data$W,
    C2 = NULL,
    num_cores = 1L,
    max_iters = 2L,
    random_seed = 20260901L,
    log_file = tempfile()
  )

  testthat::expect_identical(result$model$W, data$W)
  testthat::expect_equal(dim(result$model$mus_hat), c(20L, 3L))
  testthat::expect_true(is.finite(result$model$tau_hat))
  testthat::expect_identical(rownames(result$X), rownames(data$X))
})
