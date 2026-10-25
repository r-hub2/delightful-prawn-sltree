# Mixed data from one latent factor.
make_mixed_x <- function(n = 1500, seed = 7) {
  set.seed(seed)
  z <- rnorm(n)
  data.frame(
    c1 = z + rnorm(n, sd = 0.6),
    c2 = z + rnorm(n, sd = 0.6),
    c3 = 0.8 * z + rnorm(n, sd = 0.7),
    b1 = as.integer(z + rnorm(n, sd = 0.6) > 0),
    b2 = as.integer(0.8 * z + rnorm(n, sd = 0.7) > 0.2),
    o1 = as.integer(cut(z + rnorm(n, sd = 0.6),
                        breaks = c(-Inf, -0.7, 0.7, Inf))) - 1L  # 0-based ordinal
  )
}

test_that("learn_sltree default returns a bare ltree", {
  skip_if_not_installed("Matrix")
  x <- make_mixed_x()
  tree <- suppressWarnings(learn_sltree(x, tau1 = 3, tau2 = 5))
  tree_explicit <- suppressWarnings(
    learn_sltree(x, tau1 = 3, tau2 = 5, model = "stm"))
  expect_true(is_ltree(tree))
  expect_identical(tree, tree_explicit)
})

test_that("deprecated distance_type route has been removed", {
  x <- make_mixed_x(n = 100)
  expect_error(learn_sltree(x, distance_type = "mixed"), "unused argument")
})

test_that("return_intermediate hands back artifacts and stage timings", {
  skip_if_not_installed("Matrix")
  x <- make_mixed_x()
  fit <- suppressWarnings(
    learn_sltree(x, tau1 = 3, tau2 = 5, return_intermediate = TRUE))

  expect_type(fit, "list")
  expect_setequal(names(fit), c("tree", "correlation", "distance",
                                "variable_types", "model",
                                "allow_observed_internal", "use_near_pd",
                                "run_time"))
  expect_true(is_ltree(fit$tree))
  expect_identical(fit$model, "stm")
  expect_true(fit$allow_observed_internal)
  expect_true(fit$use_near_pd)

  expect_named(fit$run_time, c("distance", "SLS", "all"))
  expect_type(fit$run_time, "double")
  expect_equal(unname(fit$run_time["all"]),
               unname(fit$run_time["distance"] + fit$run_time["SLS"]))

  p <- ncol(x)
  expect_equal(dim(fit$correlation), c(p, p))
  expect_equal(dim(fit$distance), c(p, p))
  expect_true(all(diag(fit$distance) == 0))
  expect_length(fit$variable_types, p)
})

test_that("return_intermediate$tree matches the default-mode tree", {
  skip_if_not_installed("Matrix")
  x <- make_mixed_x()
  tree <- suppressWarnings(learn_sltree(x, tau1 = 3, tau2 = 5))
  fit  <- suppressWarnings(
    learn_sltree(x, tau1 = 3, tau2 = 5, return_intermediate = TRUE))
  expect_equal(robinson_foulds_distance(tree, fit$tree), 0)
  expect_identical(tree, fit$tree)
})

test_that("gaussian distance uses Pearson correlation", {
  x <- make_mixed_x(n = 600)[, c("c1", "c2", "c3")]
  fit <- learn_sltree(x, model = "gaussian", return_intermediate = TRUE)

  expected_correlation <- stats::cor(x)
  expected_distance <- -log(abs(expected_correlation))
  diag(expected_distance) <- 0

  expect_true(is_ltree(fit$tree))
  expect_equal(fit$correlation, expected_correlation)
  expect_equal(fit$distance, expected_distance)
  expect_equal(fit$distance,
               information_distance(x, type = "gaussian"))
  expect_true(fit$allow_observed_internal)
  expect_false(fit$use_near_pd)
})

test_that("discrete distance uses absolute-determinant distance and has no correlation", {
  x <- make_mixed_x(n = 800)
  x <- data.frame(b1 = x$b1, b2 = x$b2, b3 = as.integer(x$c3 > 0))
  fit <- learn_sltree(x, model = "discrete",
                      return_intermediate = TRUE)

  expect_true(is_ltree(fit$tree))
  expect_null(fit$correlation)
  expect_equal(fit$distance, information_distance(x, type = "discrete"))
  expect_equal(fit$variable_types, c("binary", "binary", "binary"))
  expect_true(fit$allow_observed_internal)
})

test_that("discrete distance requires equal observed state counts", {
  x <- make_mixed_x(n = 200)[, c("b1", "b2", "o1")]
  expect_error(learn_sltree(x, model = "discrete"),
               "same number of observed states")
})

test_that("discrete distance rejects continuous columns", {
  x <- make_mixed_x(n = 200)[, c("c1", "b1", "o1")]
  expect_error(learn_sltree(x, model = "discrete"),
               "requires every column to be binary or ordinal")
})

test_that("parametric mixed model uses Pearson and requires observed leaves", {
  x <- make_mixed_x(n = 600)[, c("c1", "c2", "b1", "b2")]
  fit <- learn_sltree(x, model = "mixed", return_intermediate = TRUE)

  expect_equal(fit$correlation, stats::cor(x))
  expect_equal(fit$distance,
               information_distance(x, type = "mixed"))
  expect_identical(fit$model, "mixed")
  expect_false(fit$allow_observed_internal)
  expect_false(fit$use_near_pd)
})

test_that("nearPD defaults are selected by statistical model", {
  expect_null(formals(learn_sltree)$use_near_pd)
  expect_true(sltree:::.resolve_near_pd_setting("stm"))
  expect_false(sltree:::.resolve_near_pd_setting("gaussian"))
  expect_false(sltree:::.resolve_near_pd_setting("mixed"))
  expect_false(sltree:::.resolve_near_pd_setting("discrete"))
  expect_true(sltree:::.resolve_near_pd_setting("mixed", TRUE))
  expect_false(sltree:::.resolve_near_pd_setting("gaussian", FALSE))
  expect_error(sltree:::.resolve_near_pd_setting("mixed", NA),
               "NULL or one non-missing logical")
})

test_that("parametric mixed model rejects ordinal columns", {
  x <- make_mixed_x(n = 300)[, c("c1", "b1", "o1")]
  expect_error(learn_sltree(x, model = "mixed"), "does not support ordinal")
})

test_that("model-specific tree restriction can be overridden explicitly", {
  x <- make_mixed_x(n = 300)[, c("c1", "c2", "c3")]
  expect_warning(
    fit <- learn_sltree(x, model = "gaussian", allow_observed_internal = FALSE,
                        return_intermediate = TRUE),
    "overrides the structural restriction"
  )
  expect_false(fit$allow_observed_internal)
})

test_that("semiparametric aliases stm", {
  skip_if_not_installed("Matrix")
  x <- make_mixed_x(n = 300)
  fit <- suppressWarnings(
    learn_sltree(x, model = "semiparametric", return_intermediate = TRUE))
  expect_identical(fit$model, "stm")
})

test_that("public learning APIs expose one correlation implementation", {
  expect_false("engine" %in% names(formals(latent_correlation)))
  expect_false("engine" %in% names(formals(information_distance)))
  expect_false("engine" %in% names(formals(learn_sltree)))
  expect_false("exact" %in% names(formals(latent_correlation)))
  expect_false("exact" %in% names(formals(information_distance)))
  expect_false("exact" %in% names(formals(learn_sltree)))
})

test_that("stm and information_distance share the semiparametric-distance pipeline", {
  x <- make_mixed_x(n = 500)
  fit <- learn_sltree(x, model = "stm", use_near_pd = FALSE,
                      return_intermediate = TRUE)

  expected <- information_distance(x, type = "semiparametric", use_near_pd = FALSE)
  expect_equal(fit$distance, expected)
})

test_that("nearPD is shared by learn_sltree and information_distance", {
  skip_if_not_installed("Matrix")
  x <- make_mixed_x(n = 400)
  fit <- learn_sltree(x, model = "stm", use_near_pd = TRUE,
                      return_intermediate = TRUE)

  expected <- information_distance(x, type = "semiparametric", use_near_pd = TRUE)
  expect_equal(fit$distance, expected)
})
