test_that("robinson_foulds_distance of a tree with itself is zero", {
  model <- latent_tree_models(1)
  node_types <- stats::setNames(model[, 2], model[, 1])
  true_tree  <- ltree_from_edge(model[, -2], node_types)
  expect_equal(robinson_foulds_distance(true_tree, true_tree), 0)
})

test_that("latent_tree_models returns expected structure", {
  expect_length(latent_tree_models(), 6L)
  m <- latent_tree_models(4)
  expect_identical(names(m), c("nodes", "type", "parent"))
  expect_error(latent_tree_models(7), "1 to 6")
  expect_error(latent_tree_models(99), "1 to 6")
})

test_that("learn_sltree returns a valid ltree on large-sample mixed data", {
  skip_on_cran()
  skip_if_not_installed("pcalg")
  skip_if_not_installed("graph")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("Matrix")
  set.seed(1)
  model <- latent_tree_models(1)
  dat   <- simulate_ltree_data(2e4, model, rho_lower = 0.3, rho_upper = 0.7)
  tree  <- learn_sltree(dat$data, tau1 = 3, tau2 = 5)
  expect_true(is_ltree(tree))
  expect_setequal(tree$observed_nodes, paste0("V", 1:16))
})
