test_that("discrete distance handles multi-level duplicates correctly", {
  x <- data.frame(
    a = rep(0:2, each = 20),
    b = rep(0:2, each = 20),
    z = rep(0:2, length.out = 60)
  )
  d <- information_distance(x, type = "discrete")
  expect_equal(d[1, 2], 0, tolerance = 1e-12)
  expect_true(all(d >= 0))
  expect_error(
    information_distance(data.frame(a = rep(0:2, 20), b = rep(0:1, 30)),
                         type = "discrete"),
    "same number of observed states"
  )
  expect_error(
    information_distance(data.frame(a = rep(1, 20), b = rep(0:1, 10)),
                         type = "discrete"),
    "constant"
  )
})

test_that("sls validates distances and cannot hang on a zero step", {
  d <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
  expect_error(sls(d, step_size = 0), "step_size")
  expect_error(sls(d, step_size = -0.1), "step_size")
  expect_error(sls(d, step_size = 1e-300, epsilon_max = 1),
               "1,000,000 steps")

  asymmetric <- d
  asymmetric[1, 2] <- 2
  expect_error(sls(asymmetric), "symmetric")
  negative <- d
  negative[1, 2] <- negative[2, 1] <- -1
  expect_error(sls(negative), "non-negative")
  nonfinite <- d
  nonfinite[1, 2] <- nonfinite[2, 1] <- Inf
  expect_error(sls(nonfinite), "finite")
})

test_that("learning and direct mixed helpers reject malformed inputs", {
  x2 <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_error(learn_sltree(x2), "at least 3 columns")
  expect_error(learn_sltree(x2[, 1, drop = FALSE]), "at least 3 columns")

  bad_binary <- data.frame(a = rep(c(0, 2), 10), b = rnorm(20))
  expect_error(latent_correlation(bad_binary, c("binary", "continuous")),
               "exactly 0 and 1")
  expect_error(sltree:::prepare_copula_data(data.frame(a = rnorm(20)), character()),
               "one value per column")
  expect_error(sltree:::kendall_tau_a(cbind(1:10, rep(1, 10))), "constant")
})

test_that("is_ltree returns a scalar logical for arbitrary objects", {
  expect_false(is_ltree(NULL))
  expect_false(is_ltree(list()))
  expect_false(is_ltree(structure(list(), class = "not_ltree")))
})

test_that("random tree generation is robust on boundary probabilities", {
  set.seed(1)
  no_latent <- generate_latent_tree(latent_probability = 0)
  expect_identical(names(no_latent), c("nodes", "type", "parent"))
  expect_false(any(startsWith(no_latent$nodes, "H")))
  expect_true(all(startsWith(no_latent$type, "observed_")))
  expect_equal(sum(no_latent$parent == ""), 1L)

  for (seed in seq_len(50)) {
    set.seed(seed)
    model <- generate_latent_tree()
    expect_true(is.data.frame(model))
    expect_false(anyNA(model))
    expect_equal(anyDuplicated(model$nodes), 0L)
    expect_equal(sum(model$parent == ""), 1L)
    expect_true(all(model$parent[model$parent != ""] %in% model$nodes))
  }
})

test_that("simulate_ltree_data rejects invalid models before optional dependencies", {
  bad <- data.frame(
    node = c("V1", "V2"),
    type = c("observed_Gaussian", "observed_Gaussian"),
    parent = c("V2", "V1"),
    stringsAsFactors = FALSE
  )
  expect_error(simulate_ltree_data(10, bad), "exactly one root")
  good <- data.frame(
    node = c("V1", "V2"),
    type = c("observed_Gaussian", "observed_Gaussian"),
    parent = c("", "V1"),
    stringsAsFactors = FALSE
  )
  expect_error(simulate_ltree_data(0, good), "positive integer")
  expect_error(simulate_ltree_data(10, good, rho_lower = 0.9, rho_upper = 0.2),
               "rho_lower")
})

test_that("ltree_from_edge validates and aligns tree metadata", {
  edge <- data.frame(
    node = factor(c("V1", "V2", "H2", "H1")),
    parent = factor(c("H1", "H1", "H1", "")),
    stringsAsFactors = TRUE
  )
  node_type <- c(H2 = "latent_Gaussian", V2 = "observed_binary",
                 H1 = "latent_Gaussian", V1 = "observed_ordinal")
  tree <- ltree_from_edge(edge, node_type)
  expect_identical(unname(tree$node_types[c("V1", "V2", "H2", "H1")]),
                   c("ordinal", "binary", "latent", "latent"))
  expect_equal(tree$node_metadata["H1", "observed_child_count"], 2)

  reordered <- ltree_from_edge(edge[c(4, 2, 1, 3), ], node_type)
  expect_equal(robinson_foulds_distance(tree, reordered), 0)

  cyclic <- data.frame(node = c("V1", "H1"), parent = c("H1", "V1"))
  expect_error(ltree_from_edge(cyclic), "root")
  expect_error(ltree_from_edge(edge, c("continuous", "continuous")),
               "one value per node")
})

test_that("internal structure-learning helpers are not exported", {
  exported <- getNamespaceExports("sltree")
  expect_false(any(c(".initialize_ltree_state", ".recover_tree_at_epsilon", ".judge_observed_pairs", ".find_parent_child_relations",
                     ".find_sibling_pairs", ".add_sibling_groups", ".attach_remaining_children") %in%
                   exported))
})

test_that("preset model cache is stable and protected by copy-on-modify", {
  models1 <- latent_tree_models()
  models2 <- latent_tree_models()
  expect_identical(models1, models2)

  original <- latent_tree_models(1)[1, 1]
  models1[[1]][1, 1] <- "changed-locally"
  expect_identical(latent_tree_models(1)[1, 1], original)
  expect_false(".sltree_model_cache" %in% getNamespaceExports("sltree"))
})
