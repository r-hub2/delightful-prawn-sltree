test_that("infer_variable_types correctly identifies variable types", {
  set.seed(1)
  x <- data.frame(
    cont = rnorm(50),
    bin  = rep(c(0, 1), 25),
    ord  = rep(0:2, length.out = 50)
  )
  expect_equal(infer_variable_types(x), c("continuous", "binary", "ordinal"))
})

test_that("infer_variable_types validates invalid input", {
  expect_error(infer_variable_types(data.frame(x = c(1, NA, 2))), "NA")
  expect_error(infer_variable_types(data.frame(x = rep(1, 5))), "constant")
})

test_that("ordinal levels must start at 0", {
  # One-based levels fall back to continuous.
  one_based <- rep(1:3, length.out = 50)
  expect_warning(lab <- infer_variable_types(data.frame(x = one_based)),
                 "consecutive integers starting at 0")
  expect_equal(lab, "continuous")

  expect_error(sltree:::prepare_copula_data(data.frame(x = one_based), "ordinal"),
               "starting at 0")

  zero_based <- rep(0:2, length.out = 50)
  expect_silent(out <- sltree:::prepare_copula_data(data.frame(x = zero_based), "ordinal"))
  expect_equal(ncol(out$data), 2L)  # cuts >=1, >=2
})

test_that("ordinal levels must not contain internal gaps", {
  gapped <- rep(c(0, 2, 3), length.out = 50)

  # A missing level would give duplicate binary cuts.
  expect_warning(lab <- infer_variable_types(data.frame(x = gapped)),
                 "consecutive integers starting at 0")
  expect_equal(lab, "continuous")

  expect_error(sltree:::prepare_copula_data(data.frame(x = gapped), "ordinal"),
               "no missing levels")
})
