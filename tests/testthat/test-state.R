test_that("tree state transitions keep all cached topology views aligned", {
  distance <- matrix(
    c(0, 1, 1,
      1, 0, 1,
      1, 1, 0),
    nrow = 3L,
    byrow = TRUE
  )

  tree <- sltree:::.initialize_ltree_state(distance)
  expect_silent(assert_ltree_state(tree))

  siblings <- data.frame(node1 = "V1", node2 = "V2")
  tree <- sltree:::.add_sibling_groups(siblings, tree)
  expect_silent(assert_ltree_state(tree))

  remaining <- data.frame(child = "V3", parent = "H1")
  tree <- sltree:::.attach_remaining_children(remaining, tree)
  expect_silent(assert_ltree_state(tree, allow_forest = FALSE))
  expect_identical(unname(tree$parent_by_node[c("V1", "V2", "V3", "H1")]),
                   c("H1", "H1", "H1", ""))
  expect_setequal(tree$children_by_node$H1, c("V1", "V2", "V3"))
  expect_setequal(tree$observed_descendants$H1, c("V1", "V2", "V3"))
})


test_that("state validation detects stale cached topology", {
  distance <- matrix(c(0, 1, 1, 0), nrow = 2L)
  tree <- sltree:::.initialize_ltree_state(distance)
  tree$children_by_node$V1 <- "V2"

  expect_error(
    assert_ltree_state(tree),
    "children of `V1` disagree"
  )
})


test_that("sls retains the exact three-node latent-class result", {
  distance <- matrix(
    c(0, 1, 1,
      1, 0, 1,
      1, 1, 0),
    nrow = 3L,
    byrow = TRUE
  )

  tree <- sls(distance, epsilon_init = 0.1, epsilon_max = 0.1)

  expect_silent(assert_ltree_state(tree, allow_forest = FALSE))
  expect_identical(tree$latent_nodes, "H1")
  expect_identical(unname(tree$parent_by_node[c("V1", "V2", "V3", "H1")]),
                   c("H1", "H1", "H1", ""))
})
