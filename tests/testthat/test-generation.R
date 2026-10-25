test_that("tree generation requires a subcritical final branching rate", {
  expect_error(generate_latent_tree(branching_rates = 1),
               "last value.*less than 1")
  expect_error(generate_latent_tree(branching_rates = c(2, 1)),
               "last value.*less than 1")
  expect_error(generate_latent_tree(branching_rates = c(0, 2)),
               "last value.*less than 1")
})

test_that("tree generation rejects excessive growth before allocating nodes", {
  set.seed(1)
  expect_error(generate_latent_tree(branching_rates = c(1e6, 0)),
               "100,000 nodes")
})

test_that("data simulation rejects a sample size of one at the public boundary", {
  model <- data.frame(nodes = "V1", type = "observed_Gaussian", parent = "")
  expect_error(simulate_ltree_data(1, model), "positive integer of at least 2")
})

test_that("data simulation supports a single observed node", {
  skip_if_not_installed("pcalg")
  skip_if_not_installed("graph")
  model <- data.frame(nodes = "V1", type = "observed_Gaussian", parent = "")
  set.seed(1)
  result <- simulate_ltree_data(10, model)
  expect_s3_class(result$data, "data.frame")
  expect_identical(dim(result$data), c(10L, 1L))
  expect_identical(names(result$data), "V1")
  expect_length(result$edge_correlations, 0L)
  expect_identical(result$gaussian_correlation,
                   matrix(1, 1, 1, dimnames = list("V1", "V1")))
  expect_true(all(is.finite(result$data$V1)))
  expect_gt(stats::var(result$data$V1), 0)
})
