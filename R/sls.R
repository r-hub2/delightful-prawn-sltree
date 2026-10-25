#' @include ltree_state.R sls_internal.R
NULL


#' Structure learning for semi-parametric latent trees
#'
#' Input all information distances between the observed variables,
#'   \code{sls} outputs a latent tree structure, where leaf nodes are observed
#'   variables and latent nodes are latent or observed variables.
#'   The result has one root, leaf degree 1, and latent-node degree at least 3.
#'   The root has an empty parent identifier.
#'
#' The \emph{SLS} algorithm is a bottom-up structure learning algorithm
#'   based on an additive information distance.
#'   The \emph{SLS} algorithm determines the basic sibling pair by using
#'   a prescribed threshold (i.e. argument \code{epsilon}).
#'   For fixed values of the thresholds \code{epsilon} and \code{tau2}, the
#'   algorithm may not produce a complete tree. To obtain a complete tree
#'   structure, \code{sls} starts the value of the threshold
#'   \code{epsilon} from \code{epsilon_init}. Then increase \code{epsilon}
#'   (if \code{epsilon} is less than \code{epsilon_max}) with the step size
#'   \code{step_size} until the \emph{SLS} algorithm obtains a tree.
#'   If no tree is obtained for values up to \code{epsilon_max}, \code{sls}
#'   returns a latent class model.
#'
#' Let \eqn{V} be the set of observed nodes. For \eqn{u,v \in V}, let
#'   \eqn{d_{uv}} be the information distance between \eqn{u} and \eqn{v}.
#'   Let \eqn{\Phi_{uvw} = d_{uw} - d_{vw}} be the difference of information
#'   distance.
#'   Particularly, if the difference \deqn{|\Phi_{uvw} - \Phi_{uvz}| < epsilon}
#'   for \eqn{w,z \in V \ {u,v}}, we judge that \eqn{{u,v}} is a sibling
#'   pair.
#'
#' Since a longer distance estimate is less accurate for a given number of
#'   samples, not all estimated distances are reliable for structure learning.
#'   The algorithm considers possible sibling pairs \eqn{{u,v}} whose estimated
#'   distances \eqn{d_{uv}, d_{uw}, d_{vw}} are controlled by two thresholds
#'   \code{tau1}, \code{tau2}.
#'   In particular, for each pair of nodes \eqn{{u,v}} satisfied \eqn{d_{uv}
#'   < tau1}, the estimated difference \eqn{\Phi_{uvw}} was computed only for
#'   node \deqn{w in K_{uv} = {w \in V \ {u,v} | max(d_{uw}, d_{vw}) < tau2}}.
#'   Furthermore, if the difference \deqn{|\Phi_{uvw} - \Phi_{uvz}| < epsilon}
#'   for any \eqn{w,z \in K_{uv}} and \eqn{d_{uv} < tau1}, we consider
#'   \eqn{{u,v}} to be a sibling pair.
#'
#' @return \code{sls} returns a \code{ltree} object containing
#'   the following components: \code{observed_nodes}, \code{latent_nodes},
#'   \code{internal_nodes}, \code{node_metadata}, \code{parent_by_node},
#'   \code{children_by_node}, \code{observed_descendants}, \code{active_nodes},
#'   \code{discovered_nodes}, \code{distance}, \code{node_types}, and
#'   \code{node_names}.
#'   See \code{\link{ltree_from_edge}} for more details about tree structure.
#'
#'
#' @param distance A finite, symmetric, non-negative numeric matrix with a zero
#'   diagonal and at least three rows/columns.
#' @param allow_observed_internal A logical scalar indicates whether the observed
#'   variables can be internal nodes. The default is \code{TRUE}.
#' @param tau1 A numeric scalar. \code{tau1} is the upper bound on the
#'   information distance between sibling pair. If the information distance
#'   between two nodes is greater than \code{tau1}, they are considered to
#'   be not a sibling pair, and the judgment step is skipped to save time.
#' @param tau2 A numeric scalar. The greater the information distance, the
#'   greater the error. All information distances greater than \code{tau2}
#'   are considered unreliable. But \code{tau2} should not be too small,
#'   otherwise it will result in too little information available.
#' @param relation_tolerance A numeric scalar. The threshold for judging whether the
#'   differences of information distances are different.
#' @param epsilon_init,epsilon_max,step_size Finite numeric scalars. Thresholds
#'   must be non-negative and \code{step_size} must be strictly positive. Search
#'   grids requiring more than one million steps are rejected.
#' @param node_types An optional character vector of node display types.
#'   Recognized continuous, Gaussian, and discrete-type labels are preserved;
#'   unrecognized values are treated as \code{"continuous"}.
#' @param node_names A character vector with one display name per node.
#' @seealso \code{\link{ltree_from_edge}} for more details about tree structure.
#'   \code{\link{information_distance}} for getting information distance.
#'
#' @examples
#' \dontrun{
#' require(igraph)
#'
#' set.seed(1234)
#' # needs pcalg and graph
#' model <- latent_tree_models(1)
#' data  <- simulate_ltree_data(n = 1e5, model = model)
#' x     <- data$data
#'
#' variable_types <- infer_variable_types(x)
#' correlation <- latent_correlation(x, variable_types)
#' distance <- -log(abs(correlation)); diag(distance) <- 0
#'
#' est_tree <- sls(distance, tau1 = 3, tau2 = 5,
#'                  node_types = variable_types, node_names = colnames(x))
#' edges_hat    <- edges_from_ltree(est_tree)
#' vertices_hat <- vertices_from_ltree(est_tree)
#' root <- setdiff(vertices_hat, edges_hat[, 1])
#'
#' g <- graph_from_data_frame(edges_hat, directed = FALSE, vertices = vertices_hat)
#' plot(g, layout = layout_as_tree(g, root = root), vertex.label = vertices_hat)
#'
#' true_tree <- ltree_from_edge(model[, -2], setNames(model[, 2], model[, 1]))
#' robinson_foulds_distance(true_tree, est_tree)
#' }
#' @name sls
NULL

# A node is latent-like if its name starts with "H" or it is an internal node.
.is_latent_like <- function(nodes, internal_nodes) {
  substr(nodes, 1, 1) == "H" | nodes %in% internal_nodes
}

# Run one SLS pass at a fixed epsilon value.
.recover_tree_at_epsilon <- function(distance, allow_observed_internal = TRUE, epsilon = 0.1, tau1 = 3, tau2 = 5,
                    relation_tolerance = 0.1, node_types = rep("continuous", nrow(distance)),
                    node_names = NULL) {
  tree <- .initialize_ltree_state(distance, node_types, node_names)

  # Step 1: observed-node child/sibling pairs
  repeat {
    if (allow_observed_internal) {
      parent_child_relations <- .judge_observed_pairs(
        is_sibling = FALSE, tree = tree, allow_observed_internal = allow_observed_internal,
        epsilon = epsilon, tau1 = tau1, tau2 = tau2, relation_tolerance = relation_tolerance
      )
      if (nrow(parent_child_relations) != 0L) {
        tree <- .attach_remaining_children(parent_child_relations, tree = tree)
        next()
      }
    }
    sibling_pairs <- .judge_observed_pairs(
      is_sibling = TRUE, tree = tree, allow_observed_internal = allow_observed_internal,
      epsilon = epsilon, tau1 = tau1, tau2 = tau2, relation_tolerance = relation_tolerance
    )
    tree <- .add_sibling_groups(sibling_pairs = sibling_pairs, tree = tree)
    break
  }

  # Step 2: latent-node child/sibling pairs
  if (!is.null(tree$latent_nodes)) {
    repeat {
      active_nodes <- tree$active_nodes
      internal_nodes <- tree$internal_nodes
      latent_like <- .is_latent_like(active_nodes, internal_nodes)
      ready_to_finish <- sum(latent_like) <= 2L && !any(!latent_like)
      if (ready_to_finish) break

      parent_child_relations <- .find_parent_child_relations(
        tree = tree, allow_observed_internal = allow_observed_internal, epsilon = epsilon,
        tau1 = tau1, tau2 = tau2, relation_tolerance = relation_tolerance
      )
      if (nrow(parent_child_relations) != 0L) {
        tree <- .attach_remaining_children(parent_child_relations, tree = tree)
        next()
      }

      sibling_pairs <- .find_sibling_pairs(
        tree = tree, allow_observed_internal = allow_observed_internal, epsilon = epsilon,
        tau1 = tau1, tau2 = tau2, relation_tolerance = relation_tolerance
      )
      if (nrow(sibling_pairs) != 0L) {
        tree <- .add_sibling_groups(sibling_pairs = sibling_pairs, tree = tree)
        next()
      }
      break
    }
  }

  # Step 3: two-node case
  active_nodes <- tree$active_nodes
  internal_nodes <- tree$internal_nodes
  if (length(active_nodes) == 2L) {
    latent_like <- .is_latent_like(active_nodes, internal_nodes)
    if (!any(latent_like)) {
      pair <- data.frame(node1 = active_nodes[1], node2 = active_nodes[2],
                         stringsAsFactors = FALSE)
      tree <- .add_sibling_groups(pair, tree)
    } else {
      parent_node <- active_nodes[latent_like][sum(latent_like)]
      child_node <- setdiff(active_nodes, parent_node)
      parent_child_relation <- data.frame(
        child = child_node, parent = parent_node, stringsAsFactors = FALSE
      )
      tree <- .attach_remaining_children(parent_child_relation, tree = tree)
    }
  }

  # Step 4: update internal nodes
  tree$internal_nodes <- .derive_internal_nodes(tree$parent_by_node)

  return(tree)
}

#' @rdname sls
#' @export
sls <- function(distance, allow_observed_internal = TRUE, epsilon_init = 0.1,
                 epsilon_max = max(distance), tau1 = 3, tau2 = 5, step_size = 0.1,
                 relation_tolerance = 0.1, node_types = NULL, node_names = NULL) {

  distance <- .validate_distance_matrix(distance, min_nodes = 3L)
  .validate_sls_scalar(epsilon_init, "epsilon_init", lower = 0)
  .validate_sls_scalar(epsilon_max, "epsilon_max", lower = 0)
  .validate_sls_scalar(step_size, "step_size", lower = 0, strict = TRUE)
  .validate_sls_scalar(relation_tolerance, "relation_tolerance", lower = 0)
  .validate_sls_scalar(tau1, "tau1", lower = 0)
  .validate_sls_scalar(tau2, "tau2", lower = 0)
  if (!is.logical(allow_observed_internal) || length(allow_observed_internal) != 1L ||
      is.na(allow_observed_internal)) {
    stop("`allow_observed_internal` must be one non-missing logical value.")
  }
  node_count <- nrow(distance)
  if (!is.null(node_types) &&
      (!is.character(node_types) || length(node_types) != node_count ||
       anyNA(node_types))) {
    stop("`node_types` must be NULL or a non-missing character vector with one value per node.")
  }
  if (!is.null(node_names) &&
      (!is.character(node_names) || length(node_names) != node_count ||
       anyNA(node_names))) {
    stop("`node_names` must be NULL or a non-missing character vector with one value per node.")
  }
  epsilon_max <- max(epsilon_init, epsilon_max)
  required_steps <- ceiling((epsilon_max - epsilon_init) / step_size)
  if (!is.finite(required_steps) || required_steps > 1e6) {
    stop("The epsilon search would require more than 1,000,000 steps; increase `step_size` or reduce `epsilon_max`.")
  }

  epsilon <- epsilon_init
  repeat {
    tree <- .recover_tree_at_epsilon(distance, allow_observed_internal, epsilon, tau1, tau2, relation_tolerance,
                    node_types, node_names)

    parent_by_node <- tree$parent_by_node
    node_metadata <- tree$node_metadata
    latent_nodes <- tree$latent_nodes
    observed_nodes <- tree$observed_nodes
    internal_nodes <- tree$internal_nodes
    leaf_nodes <- setdiff(observed_nodes, internal_nodes)

    is_valid_tree <- .is_valid_recovered_tree(
      parent_by_node, node_metadata, latent_nodes, leaf_nodes, internal_nodes
    )

    if (is_valid_tree) {
      break
    }

    if (epsilon < epsilon_max - 1e-7) {
      epsilon_next <- min(epsilon + step_size, epsilon_max)
      if (!is.finite(epsilon_next) || epsilon_next <= epsilon) {
        stop("`step_size` is too small to advance the epsilon search numerically.")
      }
      epsilon <- epsilon_next
    } else {
      # No tree found: fall back to a latent class model.
      all_pairs <- .all_observed_pairs(tree$observed_nodes)
      tree <- .initialize_ltree_state(distance, node_types, node_names)
      tree <- .add_sibling_groups(sibling_pairs = all_pairs, tree = tree)
      break
    }

    tree <- NULL
  }

  return(tree)
}


.validate_sls_scalar <- function(x, name, lower = -Inf, strict = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  valid <- valid && if (strict) x > lower else x >= lower
  if (!valid) {
    relation <- if (strict) "greater than" else "at least"
    stop(sprintf("`%s` must be one finite numeric value %s %s.",
                 name, relation, format(lower)))
  }
  invisible(x)
}


.validate_distance_matrix <- function(distance, min_nodes = 3L) {
  if (!is.matrix(distance) || !is.numeric(distance) ||
      nrow(distance) != ncol(distance)) {
    stop("`distance` must be a square numeric matrix.")
  }
  if (nrow(distance) < min_nodes) {
    stop(sprintf("`distance` must contain at least %d nodes.", min_nodes))
  }
  if (anyNA(distance) || !all(is.finite(distance))) {
    stop("`distance` must contain only finite values.")
  }

  tolerance <- sqrt(.Machine$double.eps)
  if (max(abs(distance - t(distance))) > tolerance) {
    stop("`distance` must be symmetric.")
  }
  if (max(abs(diag(distance))) > tolerance) {
    stop("The diagonal of `distance` must be zero.")
  }
  if (any(distance < -tolerance)) {
    stop("`distance` must contain only non-negative values.")
  }

  distance <- (distance + t(distance)) / 2
  distance[distance < 0] <- 0
  diag(distance) <- 0
  if (!is.null(rownames(distance)) && !is.null(colnames(distance)) &&
      !identical(rownames(distance), colnames(distance))) {
    stop("The row and column names of `distance` must match in the same order.")
  }
  distance
}

.is_valid_recovered_tree <- function(parent_by_node, node_metadata, latent_nodes,
                                     leaf_nodes, internal_nodes) {
  has_one_root <- sum(parent_by_node == "") == 1L
  latent_degrees_valid <- sum(node_metadata[latent_nodes, 3] < 3, na.rm = TRUE) == 0L
  leaf_degrees_valid <- sum(node_metadata[leaf_nodes, 3] != 1, na.rm = TRUE) == 0L
  internal_degrees_valid <-
    sum(node_metadata[internal_nodes, 3] < 2, na.rm = TRUE) == 0L

  has_one_root && latent_degrees_valid && leaf_degrees_valid &&
    internal_degrees_valid
}

.all_observed_pairs <- function(observed_nodes) {
  return(as.data.frame(t(as.matrix(utils::combn(observed_nodes, 2))),
                       stringsAsFactors = FALSE))
}

.derive_internal_nodes <- function(parent_by_node) {
  edge1 <- names(parent_by_node)
  edge2 <- parent_by_node
  names(edge2) <- NULL
  has_parent <- edge2 != ""
  nodes <- c(edge1[has_parent], edge2[has_parent])
  nodes_tab <- table(nodes)
  internal_nodes <- names(nodes_tab)[nodes_tab >= 2]
  return(internal_nodes)
}

