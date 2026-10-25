test_that("Pearson distances are invariant to extreme finite column scaling", {
  x <- cbind(a = c(-1, -0.5, 0.5, 1), b = c(-1, 0.5, -0.5, 1))
  expected <- -log(abs(cor(x)))
  diag(expected) <- 0
  for (type in c("gaussian", "mixed")) {
    for (scales in list(c(1e200, 1e200), c(1e-200, 1e-200),
                        c(1e200, 1e-200))) {
      expect_equal(information_distance(sweep(x, 2L, scales, "*"), type),
                   expected)
    }
    duplicate <- x[, c(1, 1)] * 1e200
    expect_equal(unname(information_distance(duplicate, type)), matrix(0, 2, 2))
    single <- information_distance(x[, 1, drop = FALSE] * 1e200, type)
    expect_equal(single, matrix(0, 1, 1, dimnames = list("a", "a")))
  }
})

test_that("invalid latent correlations are rejected instead of capped", {
  testthat::local_mocked_bindings(
    latent_correlation = function(...) matrix(c(1, NA_real_, NA_real_, 1), 2),
    .package = "sltree"
  )
  x <- cbind(a = c(0, 1, 0, 1), b = c(0, 0, 1, 1))
  expect_error(information_distance(x), "finite correlations")
})

test_that("data frames with matrix columns are rejected consistently", {
  x <- data.frame(a = I(matrix(1:12, 4, 3)), b = 2:5, c = 3:6)
  for (type in c("gaussian", "mixed", "discrete", "semiparametric", "stm")) {
    expect_error(information_distance(x, type), "one numeric vector per column")
  }
  expect_error(learn_sltree(x), "one numeric vector per column")
  expect_error(infer_variable_types(x), "one numeric vector per column")
})
