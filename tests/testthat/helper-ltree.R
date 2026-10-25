# Check that all cached views of an ltree describe the same tree.
# Intermediate SLS states are forests, so a single root is optional.
assert_ltree_state <- function(tree, allow_forest = TRUE) {
  if (!inherits(tree, "ltree") || !identical(names(tree), sltree:::.LTREE_FIELDS)) {
    stop("Malformed 'ltree': unexpected class or components.", call. = FALSE)
  }
  if (!is.logical(allow_forest) || length(allow_forest) != 1L || is.na(allow_forest)) {
    stop("`allow_forest` must be one non-missing logical value.", call. = FALSE)
  }

  ids <- names(tree$parent_by_node)
  if (is.null(ids) || any(ids == "") || anyDuplicated(ids)) {
    stop("Malformed 'ltree': parent vector must have unique node names.", call. = FALSE)
  }
  obs <- tree$observed_nodes
  lat <- tree$latent_nodes
  if (anyDuplicated(c(obs, lat)) || !setequal(c(obs, lat), ids)) {
    stop("Malformed 'ltree': observed and latent nodes do not partition the nodes.",
         call. = FALSE)
  }
  if (length(tree$parent_by_node) != length(ids) || anyNA(tree$parent_by_node) ||
      any(tree$parent_by_node != "" & !tree$parent_by_node %in% ids) || any(tree$parent_by_node == ids)) {
    stop("Malformed 'ltree': invalid parent vector.", call. = FALSE)
  }
  roots <- ids[tree$parent_by_node == ""]
  if (!allow_forest && length(roots) != 1L) {
    stop("Malformed 'ltree': a completed tree must have exactly one root.",
         call. = FALSE)
  }

  named_components <- c(
    "children_by_node", "observed_descendants", "node_types", "node_names"
  )
  for (component in named_components) {
    value <- tree[[component]]
    if (is.null(names(value)) || anyDuplicated(names(value)) ||
        !setequal(names(value), ids)) {
      stop(sprintf("Malformed 'ltree': `%s` is not aligned to the nodes.", component),
           call. = FALSE)
    }
  }
  if (is.null(rownames(tree$node_metadata)) || !setequal(rownames(tree$node_metadata), ids) ||
      !all(c("node_id", "node_kind", "degree", "observed_child_count") %in%
           names(tree$node_metadata))) {
    stop("Malformed 'ltree': `node_metadata` is not aligned to the nodes.", call. = FALSE)
  }
  if (!is.matrix(tree$distance) || nrow(tree$distance) != length(obs) ||
      ncol(tree$distance) != length(obs) ||
      !identical(rownames(tree$distance), obs) ||
      !identical(colnames(tree$distance), obs)) {
    stop("Malformed 'ltree': distance matrix is not aligned to observed nodes.",
         call. = FALSE)
  }

  expected_active <- ids[tree$parent_by_node == ""]
  expected_discovered <- ids[tree$parent_by_node != ""]
  if (!setequal(tree$active_nodes, expected_active) || !setequal(tree$discovered_nodes, expected_discovered)) {
    stop("Malformed 'ltree': active/discovered sets disagree with parents.",
         call. = FALSE)
  }

  expected_internal <- unique(tree$parent_by_node[tree$parent_by_node != ""])
  missing_internal <- setdiff(expected_internal, tree$internal_nodes)
  omitted_unary_roots <- missing_internal %in% roots &
    vapply(missing_internal, function(id) sum(tree$parent_by_node == id) == 1L, logical(1))
  if (length(setdiff(tree$internal_nodes, expected_internal)) ||
      any(!omitted_unary_roots)) {
    stop("Malformed 'ltree': internal nodes disagree with parents.", call. = FALSE)
  }

  for (id in ids) {
    expected_children <- ids[tree$parent_by_node == id]
    if (!setequal(tree$children_by_node[[id]], expected_children)) {
      stop(sprintf("Malformed 'ltree': children of `%s` disagree with parents.", id),
           call. = FALSE)
    }
  }

  expected_degree <- as.integer(tree$parent_by_node != "") +
    vapply(ids, function(id) sum(tree$parent_by_node == id), integer(1))
  expected_observed_children <- vapply(
    ids, function(id) sum(tree$parent_by_node[obs] == id), integer(1)
  )
  if (!identical(as.numeric(tree$node_metadata[ids, "degree"]),
                 as.numeric(expected_degree)) ||
      !identical(as.numeric(tree$node_metadata[ids, "observed_child_count"]),
                 as.numeric(expected_observed_children))) {
    stop("Malformed 'ltree': cached node degrees disagree with parents.",
         call. = FALSE)
  }

  for (id in ids) {
    seen <- character(0)
    frontier <- tree$children_by_node[[id]]
    while (length(frontier)) {
      if (any(frontier %in% seen)) {
        stop("Malformed 'ltree': parent/child relations contain a cycle.",
             call. = FALSE)
      }
      seen <- union(seen, frontier)
      frontier <- unique(unlist(tree$children_by_node[frontier], use.names = FALSE))
    }
    expected_descendants <- intersect(obs, seen)
    if (!setequal(tree$observed_descendants[[id]], expected_descendants)) {
      stop(sprintf("Malformed 'ltree': descendants of `%s` are stale.", id),
           call. = FALSE)
    }
  }

  invisible(tree)
}
