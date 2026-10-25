# Standardize node types; match latent first so "latent_Gaussian" is not Gaussian.
.normalize_node_types <- function(node_types) {
  index <- grep("latent", node_types)
  if (length(index) != 0) node_types[index] <- "latent"
  index <- grep("hidden", node_types)
  if (length(index) != 0) node_types[index] <- "latent"
  index <- grep("continuous", node_types)
  if (length(index) != 0) node_types[index] <- "continuous"
  index <- grep("Gauss", node_types)
  if (length(index) != 0) node_types[index] <- "Gaussian"
  index <- grep("Normal", node_types)
  if (length(index) != 0) node_types[index] <- "Gaussian"
  index <- grep("ordinal", node_types)
  if (length(index) != 0) node_types[index] <- "ordinal"
  index <- grep("binary", node_types)
  if (length(index) != 0) node_types[index] <- "binary"
  index <- grep("count", node_types)
  if (length(index) != 0) node_types[index] <- "count"
  index <- grep("categorical", node_types)
  if (length(index) != 0) node_types[index] <- "categorical"
  index <- grep("nominal", node_types)
  if (length(index) != 0) node_types[index] <- "nominal"
  index <- grep("discrete", node_types)
  if (length(index) != 0) node_types[index] <- "discrete"
  valid <- c("ordinal", "binary", "count", "categorical", "nominal", "discrete",
             "continuous", "Gaussian", "latent")
  node_types[!(node_types %in% valid)] <- "continuous"
  node_types
}

# Initialize the internal tree state.
.initialize_ltree_state <- function(distance, node_types = rep("continuous", nrow(distance)),
                 node_names = NULL) {
  node_count <- nrow(distance)
  observed_nodes <- paste0("V", seq_len(node_count))
  if (is.null(node_types)) {
    node_types <- rep("continuous", node_count)
  } else if (length(node_types) < node_count) {
    warning("The length of node_types is not equal to the number of nodes!")
    node_types <- c(node_types, rep("continuous", node_count - length(node_types)))
  } else if (length(node_types) > node_count) {
    warning("The length of node_types is not equal to the number of nodes!")
    node_types <- node_types[seq_len(node_count)]
  }
  node_types <- .normalize_node_types(node_types)

  if (is.null(node_names)) {
    distance_names <- colnames(distance)
    node_names <- if (!is.null(distance_names)) distance_names else observed_nodes
  } else if (length(node_names) < node_count) {
    warning("The length of node_names is not equal to the number of nodes!")
    node_names <- c(node_names,
                    observed_nodes[(length(node_names) + 1):node_count])
  } else if (length(node_names) > node_count) {
    warning("The length of node_names is not equal to the number of nodes!")
    node_names <- node_names[seq_len(node_count)]
  }

  colnames(distance) <- observed_nodes; rownames(distance) <- observed_nodes
  names(node_types) <- observed_nodes

  tree <- .new_ltree()
  tree$observed_nodes <- observed_nodes
  node_id <- observed_nodes
  node_kind <- paste('observed', node_types, sep = "_")
  degree <- rep(0, node_count)
  observed_child_count <- degree
  node_metadata <- data.frame(
    node_id, node_kind, degree, observed_child_count,
    stringsAsFactors = FALSE
  )
  row.names(node_metadata) <- observed_nodes
  tree$node_metadata <- node_metadata
  parent_by_node <- rep("", node_count); names(parent_by_node) <- observed_nodes
  tree$parent_by_node <- parent_by_node
  children_by_node <- vector(mode = "list", length = node_count)
  names(children_by_node) <- observed_nodes
  tree$children_by_node <- children_by_node
  tree$observed_descendants <- children_by_node
  tree$active_nodes <- observed_nodes
  tree$distance <- distance
  tree$node_types <- node_types
  names(node_names) <- observed_nodes
  tree$node_names <- node_names

  return(tree)
}

# Recompute the active and discovered node sets from the parent vector.
.refresh_active_sets <- function(tree) {
  nodes <- names(tree$parent_by_node)
  has_parent <- tree$parent_by_node != ""
  tree$active_nodes <- nodes[!has_parent]
  tree$discovered_nodes <- nodes[has_parent]
  tree
}

# Add maximal sibling groups to the internal tree state.
#' @importFrom igraph graph_from_data_frame components groups
#' @noRd
.add_sibling_groups <- function(sibling_pairs, tree) {
  if (is.null(sibling_pairs)) return(tree)
  if (!inherits(tree, "ltree")) {
    stop("`tree` must be an 'ltree' object.")
  }

  sibling_graph <- graph_from_data_frame(sibling_pairs)
  graph_components <- components(sibling_graph)
  sibling_groups <- groups(graph_components)
  sibling_nodes <- unlist(sibling_groups)
  if (is.null(sibling_nodes)) return(tree)

  latent_nodes <- tree$latent_nodes
  internal_nodes <- tree$internal_nodes
  existing_latent_count <- length(latent_nodes)
  new_latent_count <- length(sibling_groups)
  new_latent_nodes <- paste0(
    "H", existing_latent_count + seq_len(new_latent_count)
  )

  tree$latent_nodes <- c(latent_nodes, new_latent_nodes)
  tree$internal_nodes <- c(internal_nodes, new_latent_nodes)
  node_id <- new_latent_nodes
  node_kind <- rep("latent", new_latent_count)
  degree <- rep(0, new_latent_count)
  observed_child_count <- degree
  for (i in seq_len(new_latent_count)) {
    degree[i] <- length(sibling_groups[[i]])
    is_observed_child <- sibling_groups[[i]] %in% tree$observed_nodes
    observed_child_count[i] <- sum(is_observed_child)
  }
  new_metadata <- data.frame(
    node_id, node_kind, degree, observed_child_count,
    stringsAsFactors = FALSE
  )
  row.names(new_metadata) <- new_latent_nodes
  existing_metadata <- tree$node_metadata
  existing_sibling_nodes <- intersect(sibling_nodes, existing_metadata[, 1])
  existing_metadata[existing_sibling_nodes, 3] <-
    existing_metadata[existing_sibling_nodes, 3] + 1
  tree$node_metadata <- rbind(existing_metadata, new_metadata)
  existing_parents <- tree$parent_by_node
  new_parents <- rep("", new_latent_count)
  names(new_parents) <- new_latent_nodes
  for (i in seq_len(new_latent_count)) {
    assignable_children <- intersect(names(existing_parents), sibling_groups[[i]])
    if (sum(existing_parents[assignable_children] != "") != 0)
      warning("The nodes in the sibling group already have parent nodes!")
    existing_parents[assignable_children] <- new_latent_nodes[i]
  }
  tree$parent_by_node <- c(existing_parents, new_parents)
  names(sibling_groups) <- new_latent_nodes
  children_by_node <- c(tree$children_by_node, sibling_groups)
  tree$children_by_node <- children_by_node
  existing_descendants <- tree$observed_descendants
  new_descendants <- vector(mode = "list", length = new_latent_count)
  names(new_descendants) <- new_latent_nodes
  for (i in seq_len(new_latent_count)) {
    internal_children <- intersect(sibling_groups[[i]], tree$internal_nodes)
    observed_children <- intersect(sibling_groups[[i]], tree$observed_nodes)
    inherited_descendants <- NULL
    if (length(internal_children) != 0L) {
      inherited_descendants <- do.call("c", existing_descendants[internal_children])
    }
    new_descendants[[i]] <- c(inherited_descendants, observed_children)
  }
  observed_descendants <- c(existing_descendants, new_descendants)
  tree$observed_descendants <- observed_descendants
  tree <- .refresh_active_sets(tree)
  node_types <- tree$node_types
  new_node_types <- rep("latent", new_latent_count)
  names(new_node_types) <- new_latent_nodes
  node_types <- c(node_types, new_node_types)
  tree$node_types <- node_types
  new_node_names <- new_latent_nodes
  names(new_node_names) <- new_latent_nodes
  node_names <- c(tree$node_names, new_node_names)
  tree$node_names <- node_names

  return(tree)
}

# Add remaining-child relationships to the internal tree state.
.attach_remaining_children <- function(parent_child_relations, tree) {
  if (is.null(parent_child_relations)) return(tree)
  if (!inherits(tree, "ltree")) {
    stop("`tree` must be an 'ltree' object.")
  }

  # Character matrix, so factor columns do not become integer codes.
  parent_child_relations <- as.matrix(parent_child_relations)
  storage.mode(parent_child_relations) <- "character"

  relation_count <- nrow(parent_child_relations)
  internal_nodes <- tree$internal_nodes
  internal_nodes <- union(internal_nodes, parent_child_relations[, 2])
  tree$internal_nodes <- internal_nodes
  # Update child and parent degrees.
  node_metadata <- tree$node_metadata
  child_nodes <- parent_child_relations[, 1]
  parent_nodes <- parent_child_relations[, 2]
  node_metadata[child_nodes, 3] <- node_metadata[child_nodes, 3] + 1
  children_per_parent <- table(parent_nodes)
  node_metadata[names(children_per_parent), 3] <-
    node_metadata[names(children_per_parent), 3] + as.numeric(children_per_parent)
  observed_relations <- child_nodes %in% tree$observed_nodes
  if (any(observed_relations)) {
    observed_per_parent <- table(parent_nodes[observed_relations])
    node_metadata[names(observed_per_parent), 4] <-
      node_metadata[names(observed_per_parent), 4] + as.numeric(observed_per_parent)
  }
  tree$node_metadata <- node_metadata
  # Assign parents to child nodes.
  parent_by_node <- tree$parent_by_node
  already_assigned <- parent_by_node[child_nodes] != ""
  if (any(already_assigned)) {
    warning(paste("The child nodes",
                  paste(child_nodes[already_assigned], sep = "", collapse = ","),
                  "already have parent!", sep = ""))
  }
  parent_by_node[child_nodes] <- parent_nodes
  tree$parent_by_node <- parent_by_node
  # Add children to parent nodes.
  children_by_node <- tree$children_by_node
  for (i in seq_len(relation_count)) {
    children_by_node[[parent_nodes[i]]] <- c(
      children_by_node[[parent_nodes[i]]], child_nodes[i]
    )
  }
  tree$children_by_node <- children_by_node
  # Add observed descendants to parent nodes.
  observed_descendants <- tree$observed_descendants
  for (i in seq_len(relation_count)) {
    inherited_descendants <- unlist(observed_descendants[child_nodes[i]])
    if (child_nodes[i] %in% tree$observed_nodes) {# if child is observed
      observed_descendants[[parent_nodes[i]]] <- c(
        inherited_descendants, child_nodes[i],
        observed_descendants[[parent_nodes[i]]]
      )
    } else {# if child is latent
      observed_descendants[[parent_nodes[i]]] <- c(
        inherited_descendants, observed_descendants[[parent_nodes[i]]]
      )
    }
  }
  tree$observed_descendants <- observed_descendants
  tree <- .refresh_active_sets(tree)

  return(tree)
}

.is_phi_constant <- function(phi_vec, epsilon) {
  if (length(phi_vec) <= 1) return(TRUE)
  return(diff(range(phi_vec)) <= epsilon)
}

.is_phi_near <- function(phi_vec, phi, epsilon) {
  if (length(phi_vec) == 0) return(TRUE)
  return(all(abs(phi_vec-phi) <= epsilon))
}

.is_phi_above <- function(phi_vec, phi, relation_tolerance = 0.1) {
  if (length(phi_vec) == 0) return(TRUE)
  return(all(phi_vec > phi + relation_tolerance))
}

.is_phi_below <- function(phi_vec, phi, relation_tolerance = 0.1) {
  if (length(phi_vec) == 0) return(TRUE)
  return(all(phi_vec < phi - relation_tolerance))
}

# Select two observed representatives below a latent node.
.select_bifurcation_pair <- function(u, distance, tree) {
  children_by_node <- tree$children_by_node
  ch_u <- children_by_node[[u]]
  ch_obs <- intersect(ch_u, tree$observed_nodes)
  ch_lat <- setdiff(ch_u, tree$observed_nodes)
  num <- length(ch_obs)
  observed_descendants <- tree$observed_descendants
  diag(distance) <- -1

  if (num >= 2) {
    ij <- order(distance[ch_obs, ch_obs]) [num+1]
    i <- 1 + (ij - 1) %% num; j <- 1 + (ij - 1) %/% num
    a <- ch_obs[c(i, j)]
    a <- intersect(tree$observed_nodes, a)# sort a
  } else if (num == 1) {
    a <- setdiff(tree$observed_nodes, c(ch_obs, observed_descendants[[ch_obs]]))
    b <- intersect(observed_descendants[[u]], a)
    d_ij <- distance[ch_obs, b]
    i <- (order(d_ij))[1]
    a <- c(ch_obs, b[i])
  } else {
    if (length(ch_lat) < 2)
      stop("A latent node should have two child-nodes at least!")
    des_lat1 <- intersect(observed_descendants[[ch_lat[1]]], tree$observed_nodes)
    des_lat2 <- intersect(observed_descendants[[ch_lat[2]]], tree$observed_nodes)
    temp <- length(intersect(des_lat1, des_lat2))
    ij <- order(distance[des_lat1, des_lat2]) [temp+1]
    i <- 1 + (ij - 1) %% length(des_lat1)
    j <- 1 + (ij - 1) %/% length(des_lat1)
    a <- c(des_lat1[i], des_lat2[j])
    a <- intersect(tree$observed_nodes, a)# sort a
  }

  return(a)
}

# Judge an observed-observed active-node pair.
.judge_observed_pair <- function(nodes, allow_observed_internal, is_sibling, tree, epsilon, tau1, tau2, relation_tolerance) {
  u <- nodes[1]
  v <- nodes[2]

  distance <- tree$distance
  observed_nodes <- tree$observed_nodes
  de_u <- tree$observed_descendants[[u]]
  de_v <- tree$observed_descendants[[v]]
  tree <- NULL

  temp <- union(de_u, de_v)
  temp <- union(c(u, v), temp)
  z <- setdiff(observed_nodes, temp)

  if (length(z) == 0) {
    if (!is_sibling) {
      return(ifelse(distance[u, v] < tau1, -1L, 0L))
    } else {
      if (distance[u, v] < tau1)
        warning(paste(u, ",", v, "should not be judged as sibling pair!", sep = ""))
      return(FALSE)
    }
  }

  if (distance[u, v] > tau1) {
    return(if (!is_sibling) 0L else FALSE)
  }

  logi1 <- observed_nodes %in% z
  logi2 <- distance[u, ] <= tau2 & distance[v, ] <= tau2
  logi <- logi1 & logi2

  d_uv <- distance[u, v]
  relation_tolerance <- relation_tolerance * (d_uv + (1 - d_uv) * (d_uv >= 1))

  flag <- flag11 <- flag12 <- flag21 <- flag22 <- FALSE
  if (sum(logi) != 0) {
    phi_uv <- distance[u, logi] - distance[v, logi]
    flag <- .is_phi_constant(phi_vec = phi_uv, epsilon = epsilon)
    if (allow_observed_internal) {
      flag11 <- .is_phi_near(phi_vec = phi_uv, phi = d_uv, epsilon = epsilon)
      flag12 <- .is_phi_above(phi_vec = phi_uv, phi = -d_uv, relation_tolerance = relation_tolerance)
      flag21 <- .is_phi_below(phi_vec = phi_uv, phi = d_uv, relation_tolerance = relation_tolerance)
      flag22 <- .is_phi_near(phi_vec = phi_uv, phi = -d_uv, epsilon = epsilon)
    } else {
      flag11 <- flag22 <- FALSE
      flag12 <- flag21 <- TRUE
    }

  }

  return(.classify_pair_relation(TRUE, is_sibling, logi, flag, flag11, flag12, flag21, flag22))
}

# Judge an observed-latent active-node pair.
.judge_observed_latent_pair <- function(nodes, allow_observed_internal, is_sibling, tree, epsilon, tau1, tau2, relation_tolerance) {
  u <- nodes[1]
  h <- nodes[2]

  distance <- tree$distance
  observed_nodes <- tree$observed_nodes
  ch_h <- tree$children_by_node[[h]]
  de_u <- tree$observed_descendants[[u]]
  de_h <- tree$observed_descendants[[h]]

  ij <- .select_bifurcation_pair(h, distance, tree)
  i <- ij[1]
  j <- ij[2]
  tree <- NULL

  if (!is_sibling) {
    ch_h_obs <- intersect(observed_nodes, ch_h)
    if (length(ch_h_obs) > 0) {
      if (min(distance[ch_h_obs, u]) > tau1) {
        return(0L)
      }
    }
  }

  temp <- union(de_u, de_h)
  temp <- union(u, temp)
  z <- setdiff(observed_nodes, temp)

  if (length(z) == 0) {
    if (!is_sibling) {
      ch_h_obs <- intersect(observed_nodes, ch_h)
      if (length(ch_h_obs) > 0) {
        if (min(distance[ch_h_obs, u]) < tau1) {
          return(1L)
        } else {
          return(0L)
        }
      } else {
        return(1L)
      }
    } else {
      return(FALSE)
    }
  }

  logi1 <- observed_nodes %in% z
  logi2 <- distance[u, ] <= tau2 & distance[i, ] <= tau2
  logi <- logi1 & logi2

  d_ui <- distance[u, i]
  phi_uij <- distance[u, j] - distance[i, j]
  temp <- d_ui + (1 - d_ui) * (d_ui >= 1)

  flag <- flag11 <- flag12 <- flag21 <- flag22 <- FALSE
  if (sum(logi) != 0) {
    phi_uv <- distance[u, logi] - distance[i, logi]
    flag <- .is_phi_constant(phi_vec = phi_uv, epsilon = epsilon)
    flag21 <- .is_phi_below(phi_vec = phi_uv, phi = phi_uij, relation_tolerance = relation_tolerance * temp)
    if (allow_observed_internal) {
      flag11 <- .is_phi_near(phi_vec = phi_uv, phi = phi_uij, epsilon = epsilon)
      flag12 <- .is_phi_above(phi_vec = phi_uv, phi = -d_ui, relation_tolerance = relation_tolerance * temp)
      flag22 <- .is_phi_near(phi_vec = phi_uv, phi = -d_ui, epsilon = epsilon)
    } else {
      flag11 <- flag22 <- FALSE
      flag12 <- TRUE
    }

  }
  distance <- NULL

  return(.classify_pair_relation(FALSE, is_sibling, logi, flag, flag11, flag12, flag21, flag22))
}

# Judge a latent-latent active-node pair.
.judge_latent_pair <- function(nodes, is_sibling, tree, epsilon, tau2, relation_tolerance) {
  h1 <- nodes[1]
  h2 <- nodes[2]

  distance <- tree$distance
  observed_nodes <- tree$observed_nodes
  de_h1 <- tree$observed_descendants[[h1]]
  de_h2 <- tree$observed_descendants[[h2]]

  ij <- .select_bifurcation_pair(h1, distance, tree)
  i <- ij[1]
  j <- ij[2]
  kl <- .select_bifurcation_pair(h2, distance, tree)
  k <- kl[1]
  l <- kl[2]
  tree <- NULL

  temp <- union(de_h1, de_h2)
  z <- setdiff(observed_nodes, temp)

  if (length(z) == 0) {
    if (!is_sibling) {
      return(1L)
    } else {
      warning(paste(h1, ",", h2, "should not be judged as sibling pair!", sep = ""))
      return(FALSE)
    }
  }

  logi1 <- observed_nodes %in% z
  logi2 <- distance[i, ] <= tau2 & distance[k, ] <= tau2
  logi <- logi1 & logi2

  phi_ikj <- distance[i, j] - distance[k, j]
  phi_ikl <- distance[i, l] - distance[k, l]

  flag <- flag11 <- flag12 <- flag21 <- flag22 <- FALSE
  if (sum(logi) != 0) {
    phi_ik <- distance[i, logi] - distance[k, logi]
    flag <- .is_phi_constant(phi_vec = phi_ik, epsilon = epsilon)
    flag11 <- .is_phi_near(phi_vec = phi_ik, phi = phi_ikl, epsilon = epsilon)
    flag12 <- .is_phi_above(phi_vec = phi_ik, phi = phi_ikj, relation_tolerance = relation_tolerance)
    flag21 <- .is_phi_below(phi_vec = phi_ik, phi = phi_ikl, relation_tolerance = relation_tolerance)
    flag22 <- .is_phi_near(phi_vec = phi_ik, phi = phi_ikj, epsilon = epsilon)
  }

  return(.classify_pair_relation(FALSE, is_sibling, logi, flag, flag11, flag12, flag21, flag22))
}


# Judge observed-observed pairs in the active set.
.judge_observed_pairs <- function(is_sibling, tree, allow_observed_internal, epsilon, tau1, tau2, relation_tolerance = 0.1) {

  active_nodes <- tree$active_nodes
  observed_nodes <- tree$observed_nodes
  upper_distances <- tree$distance[upper.tri(tree$distance)]

  if (length(active_nodes) < 2L) return(matrix("", 0, 2))

  candidate_pairs <- t(as.matrix(utils::combn(active_nodes, 2)))

  both_observed <- candidate_pairs[, 1] %in% observed_nodes &
    candidate_pairs[, 2] %in% observed_nodes
  candidate_pairs <- candidate_pairs[both_observed, , drop = FALSE]

  if (!any(both_observed)) return(matrix("", 0, 2))

  distance_index <- .upper_triangle_index(candidate_pairs, observed_nodes)
  within_tau1 <- upper_distances[distance_index] <= tau1

  if (!any(within_tau1)) return(matrix("", 0, 2))

  candidate_pairs <- candidate_pairs[within_tau1, , drop = FALSE]

  relation_results <- apply(
    candidate_pairs, MARGIN = 1, FUN = .judge_observed_pair,
    allow_observed_internal = allow_observed_internal,
    is_sibling = is_sibling, tree = tree,
    epsilon = epsilon, tau1 = tau1, tau2 = tau2,
    relation_tolerance = relation_tolerance
  )

  if (!is_sibling) {
    candidate_pairs[relation_results == -1, ] <-
      candidate_pairs[relation_results == -1, 2:1, drop = FALSE]
    parent_child_relations <-
      candidate_pairs[relation_results != 0, , drop = FALSE]

    if (!nrow(parent_child_relations)) return(parent_child_relations)

    .remove_ambiguous_relations(parent_child_relations)
  } else {
    candidate_pairs[relation_results, , drop = FALSE]
  }
}

# Find remaining-child relationships in the active set.
.find_parent_child_relations <- function(tree, allow_observed_internal, epsilon, tau1, tau2, relation_tolerance = 0.1) {
  judgments <- .judge_active_pairs(
    tree = tree, is_sibling = FALSE, allow_observed_internal = allow_observed_internal,
    epsilon = epsilon, tau1 = tau1, tau2 = tau2, relation_tolerance = relation_tolerance
  )
  if (is.null(judgments)) return(matrix("", 0, 2))
  candidate_pairs <- judgments$pairs
  relation_results <- judgments$results

  # Orient all remaining-child pairs consistently.
  candidate_pairs[relation_results == -1, ] <-
    candidate_pairs[relation_results == -1, 2:1, drop = FALSE]

  parent_child_relations <-
    candidate_pairs[relation_results != 0, , drop = FALSE]

  if (!nrow(parent_child_relations)) return(parent_child_relations)

  # Remove ambiguous or cyclic relationships.
  .remove_ambiguous_relations(parent_child_relations)
}

# Find sibling relationships in the active set.
.find_sibling_pairs <- function(tree, allow_observed_internal, epsilon, tau1, tau2, relation_tolerance = 0.1) {
  judgments <- .judge_active_pairs(
    tree = tree, is_sibling = TRUE, allow_observed_internal = allow_observed_internal,
    epsilon = epsilon, tau1 = tau1, tau2 = tau2, relation_tolerance = relation_tolerance
  )
  if (is.null(judgments)) return(matrix("", 0, 2))
  judgments$pairs[judgments$results, , drop = FALSE]
}



# Combine comparison flags into a relationship.
.classify_pair_relation <- function(non_latent, is_sibling, has_comparisons,
                                    flag, flag11, flag12, flag21, flag22) {
  if (!is_sibling) {
    if (!any(has_comparisons)) {
      return(0L)
    } else {
      if (flag && flag11 && flag12) {
        return(1L)
      } else if (flag && flag21 && flag22) {
        return(-1L)
      } else {
        return(0L)
      }
    }
  } else {
    if (!any(has_comparisons)) {
      return(non_latent)
    } else {
      if (flag && flag12 && flag21) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
}

.remove_ambiguous_relations <- function(relations) {
  child_nodes <- relations[, 1]
  if (anyDuplicated(child_nodes)) {
    first_index <- seq_along(child_nodes)
    names(first_index) <- child_nodes
    relations <- relations[first_index[unique(child_nodes)], , drop = FALSE]
  }

  forms_cycle <- relations[, 1] %in% relations[, 2]
  if (all(forms_cycle)) {
    relations <- relations[1, , drop = FALSE]
  } else {
    relations <- relations[!forms_cycle, , drop = FALSE]
  }

  relations
}


# Map observed-node pairs to positions in distance[upper.tri(distance)].
.upper_triangle_index <- function(pairs, observed_nodes) {
  node_index <- seq_along(observed_nodes)
  names(node_index) <- observed_nodes
  i <- node_index[pairs[, 1]]
  j <- node_index[pairs[, 2]]
  (j - 2) * (j - 1) / 2 + i
}

# Enumerate and judge candidate pairs in the active set.
.judge_active_pairs <- function(tree, is_sibling, allow_observed_internal, epsilon, tau1, tau2,
                               relation_tolerance = 0.1) {
  active_nodes <- tree$active_nodes
  observed_nodes <- tree$observed_nodes
  distance <- tree$distance
  latent_nodes <- tree$latent_nodes
  upper_distances <- distance[upper.tri(distance)]
  distance <- NULL

  candidate_pairs <- t(as.matrix(utils::combn(active_nodes, 2)))

  observed_latent_pairs <- candidate_pairs[, 1] %in% observed_nodes &
    candidate_pairs[, 2] %in% latent_nodes
  latent_pairs <- candidate_pairs[, 1] %in% latent_nodes &
    candidate_pairs[, 2] %in% latent_nodes
  if (allow_observed_internal) {
    observed_pairs <- candidate_pairs[, 1] %in% observed_nodes &
      candidate_pairs[, 2] %in% observed_nodes
    supported_pairs <- observed_pairs | observed_latent_pairs | latent_pairs
  } else {
    supported_pairs <- observed_latent_pairs | latent_pairs
  }

  if (!any(supported_pairs)) return(NULL)

  observed_latent_candidates <-
    candidate_pairs[observed_latent_pairs, , drop = FALSE]
  latent_candidates <- candidate_pairs[latent_pairs, , drop = FALSE]

  observed_results <- NULL
  if (allow_observed_internal) {
    observed_candidates <- candidate_pairs[observed_pairs, , drop = FALSE]
    observed_within_tau1 <- observed_pairs
    if (any(observed_pairs)) {
      observed_within_tau1 <- upper_distances[
        .upper_triangle_index(observed_candidates, observed_nodes)
      ] <= tau1
      observed_candidates <-
        observed_candidates[observed_within_tau1, , drop = FALSE]
    }
    if (any(observed_within_tau1)) {
      observed_results <- apply(
        observed_candidates, MARGIN = 1, FUN = .judge_observed_pair,
        allow_observed_internal = allow_observed_internal,
        is_sibling = is_sibling, tree = tree,
        epsilon = epsilon, tau1 = tau1, tau2 = tau2,
        relation_tolerance = relation_tolerance
      )
    }
  }

  observed_latent_results <- NULL
  latent_results <- NULL
  if (any(observed_latent_pairs)) {
    observed_latent_results <- apply(observed_latent_candidates, MARGIN = 1,
                     FUN = .judge_observed_latent_pair,
                     allow_observed_internal = allow_observed_internal,
                     is_sibling = is_sibling, tree = tree,
                     epsilon = epsilon, tau1 = tau1, tau2 = tau2,
                     relation_tolerance = relation_tolerance)
  }
  if (any(latent_pairs)) {
    latent_results <- apply(latent_candidates, MARGIN = 1,
                     FUN = .judge_latent_pair,
                     is_sibling = is_sibling, tree = tree,
                     epsilon = epsilon, tau2 = tau2, relation_tolerance = relation_tolerance)
  }

  if (allow_observed_internal) {
    list(
      pairs = rbind(observed_candidates, observed_latent_candidates,
                    latent_candidates),
      results = c(observed_results, observed_latent_results, latent_results)
    )
  } else {
    list(
      pairs = rbind(observed_latent_candidates, latent_candidates),
      results = c(observed_latent_results, latent_results)
    )
  }
}


