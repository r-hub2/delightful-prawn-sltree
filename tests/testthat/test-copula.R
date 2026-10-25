test_that("continuous correlations use the compiled Kendall bridge", {
  set.seed(42)
  x <- cbind(a = rnorm(200), b = rnorm(200))
  tau <- pcaPP::cor.fk(x)[1, 2]
  correlation <- latent_correlation(x, c("continuous", "continuous"))

  expect_equal(correlation[1, 2], sin(pi / 2 * tau))
  expect_identical(names(formals(latent_correlation)), c("x", "variable_types"))
})

test_that("an exactly independent pair yields a finite (capped) distance", {
  # tau(a, b) is exactly 0, so -log|rho| is Inf before capping.
  x <- data.frame(
    a = c(0.1, 0.2, 0.3, 0.4),
    b = c(0.2, 0.4, 0.1, 0.3),
    d = c(0.5, 0.6, 0.7, 0.8)
  )
  expect_equal(unname(latent_correlation(x, rep("continuous", 3))[1, 2]), 0)

  dist <- information_distance(x, type = "semiparametric")
  expect_true(all(is.finite(dist)))
  expect_gt(dist["a", "b"], 5)          # independent pair -> large but finite
  expect_equal(dist["a", "b"], -log(.Machine$double.xmin))
})
