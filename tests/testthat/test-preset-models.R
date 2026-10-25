test_that("every preset model is a well-formed rooted tree", {
  for (k in seq_along(latent_tree_models())) {
    model <- latent_tree_models(k)

    expect_identical(names(model), c("nodes", "type", "parent"))
    expect_false(anyDuplicated(model$nodes) > 0)
    expect_equal(sum(model$parent == ""), 1L)
    expect_true(all(model$parent[model$parent != ""] %in% model$nodes))

    tree <- ltree_from_edge(model[, -2], stats::setNames(model[, 2], model[, 1]))
    expect_true(is_ltree(tree))
    expect_equal(robinson_foulds_distance(tree, tree), 0)
  }
})


test_that("latent node identifiers are the H-prefixed ones", {
  model <- latent_tree_models(1)
  tree <- ltree_from_edge(model[, -2], stats::setNames(model[, 2], model[, 1]))
  expect_setequal(tree$latent_nodes, grep("^H", model$nodes, value = TRUE))
  expect_setequal(tree$observed_nodes, grep("^H", model$nodes, value = TRUE,
                                            invert = TRUE))
})
