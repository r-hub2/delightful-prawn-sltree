#' Randomly generate a mixed latent tree structure
#'
#' Generate a tree by Poisson branching and optionally mark nodes of degree at
#' least 3 as latent.
#'
#' @param variable_type Observed variable type: "Gaussian", "continuous", "binary",
#'   "ordinal", or "discrete"; use "mixed" to sample discrete, ordinal, or
#'   binary type independently for each observed node.
#' @param latent_probability Probability that a node with degree >= 3 is set as a latent node.
#' @param branching_rates Non-negative, finite Poisson branching means for
#'   successive levels. The last value is reused for deeper levels and must be
#'   less than 1. Generation stops with an error before exceeding 100,000 nodes.
#' @return A three-column data frame \code{nodes}/\code{type}/\code{parent}
#'   (\code{""} for the root), as used by \code{\link{simulate_ltree_data}}.
#' @seealso \code{\link{simulate_ltree_data}}, \code{\link{latent_tree_models}}
#' @importFrom stats rpois runif
#' @export
generate_latent_tree <- function(variable_type = "Gaussian", latent_probability = 0.5,
                              branching_rates = c(2, 2, 1, 0.5, 0.05)) {
  valid_types <- c("Gaussian", "continuous", "binary", "ordinal",
                   "discrete", "mixed")
  if (!is.character(variable_type) || length(variable_type) != 1L ||
      is.na(variable_type) || !variable_type %in% valid_types) {
    stop("`variable_type` must be one of: ", paste(valid_types, collapse = ", "), ".")
  }
  if (!is.numeric(latent_probability) || length(latent_probability) != 1L ||
      !is.finite(latent_probability) || latent_probability < 0 || latent_probability > 1) {
    stop("`latent_probability` must be one finite number between 0 and 1.")
  }
  if (!is.numeric(branching_rates) || length(branching_rates) < 1L ||
      any(!is.finite(branching_rates)) || any(branching_rates < 0)) {
    stop("`branching_rates` must be a non-empty vector of finite, non-negative numbers.")
  }
  if (branching_rates[length(branching_rates)] >= 1) {
    stop("The last value of `branching_rates` must be less than 1.")
  }

  max_nodes <- 100000L
  node <- "n1"
  parent <- ""
  frontier <- node
  next_id <- 0L
  level <- 0L
  while (length(frontier)) {
    if (level == 0L) {
      child_counts <- stats::rpois(1L, 1) + 3L
    } else {
      lambda <- branching_rates[min(level, length(branching_rates))]
      child_counts <- stats::rpois(length(frontier), lambda)
    }
    total <- sum(child_counts)
    if (!is.finite(total) || total > max_nodes - length(node)) {
      stop("Tree generation would exceed 100,000 nodes; reduce `branching_rates`.")
    }
    if (total == 0L) break

    new_nodes <- paste0("N", next_id + seq_len(total))
    next_id <- next_id + total
    node <- c(node, new_nodes)
    parent <- c(parent, rep(frontier, child_counts))
    frontier <- new_nodes
    level <- level + 1L
  }

  edges <- parent != ""
  degree <- stats::setNames(integer(length(node)), node)
  degree[] <- as.integer(table(factor(c(node[edges], parent[edges]), levels = node)))
  candidates <- node[degree >= 3L]
  nodes_to_lat <- candidates[stats::runif(length(candidates)) < latent_probability]
  is_latent <- node %in% nodes_to_lat

  type <- character(length(node))
  type[is_latent] <- "latent_Gaussian"
  if (variable_type == "mixed") {
    observed_types <- sample(c("discrete", "ordinal", "binary"),
                             sum(!is_latent), replace = TRUE)
  } else {
    observed_types <- rep(variable_type, sum(!is_latent))
  }
  type[!is_latent] <- paste0("observed_", observed_types)

  translated <- character(length(node))
  translated[is_latent] <- paste0("H", seq_len(sum(is_latent)))
  translated[!is_latent] <- paste0("V", seq_len(sum(!is_latent)))
  names(translated) <- node
  model <- data.frame(
    nodes = unname(translated[node]),
    type = type,
    parent = ifelse(parent == "", "", unname(translated[parent])),
    stringsAsFactors = FALSE
  )
  model <- rbind(model[!is_latent, , drop = FALSE],
                 model[is_latent, , drop = FALSE])
  rownames(model) <- NULL
  model
}

#' Generate mixed observed data from a given latent tree
#'
#' Assign edge correlations, generate latent Gaussian data along the tree, and
#' transform observed variables according to their declared types.
#'
#' @param n Observed sample size (number of rows), an integer of at least 2.
#' @param model Valid rooted latent tree model (three-column data frame, see
#'   \code{\link{generate_latent_tree}} or \code{\link{latent_tree_models}}). Node
#'   identifiers must be unique; latent nodes use type
#'   \code{"latent_Gaussian"}.
#' @param rho_lower,rho_upper Bounds of the uniform distribution for edge correlations.
#' @return A list with \code{data}, \code{edge_correlations}, and
#'   \code{gaussian_correlation}.
#' @seealso \code{\link{generate_latent_tree}}, \code{\link{learn_sltree}}
#' @importFrom stats cor quantile runif
#' @importFrom mvtnorm rmvnorm
#' @importFrom methods as
#' @export
simulate_ltree_data <- function(n, model, rho_lower = 0.2, rho_upper = 0.8) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 2 || n > .Machine$integer.max || n != floor(n)) {
    stop("`n` must be one positive integer of at least 2.")
  }
  n <- as.integer(n)
  if (!(is.matrix(model) || is.data.frame(model)) || ncol(model) != 3L ||
      nrow(model) < 1L) {
    stop("`model` must be a three-column matrix or data frame with at least one row.")
  }
  model <- as.data.frame(lapply(as.data.frame(model), as.character),
                         stringsAsFactors = FALSE)
  colnames(model) <- c("nodes", "type", "parent")
  if (anyNA(model) || any(model$nodes == "") || anyDuplicated(model$nodes)) {
    stop("`model` must contain unique, non-empty node identifiers and no missing values.")
  }
  if (sum(model$parent == "") != 1L ||
      any(!model$parent[model$parent != ""] %in% model$nodes) ||
      any(model$nodes == model$parent)) {
    stop("`model` must describe a rooted tree with exactly one root and valid parents.")
  }
  valid_node_type <- grepl(
    "^(latent_Gaussian|observed_(Gaussian|continuous|binary|ordinal|discrete))$",
    model$type
  )
  if (any(!valid_node_type)) {
    stop("`model` contains an unsupported node type.")
  }
  if (!is.numeric(rho_lower) || length(rho_lower) != 1L || !is.finite(rho_lower) ||
      !is.numeric(rho_upper) || length(rho_upper) != 1L || !is.finite(rho_upper) ||
      rho_lower < 0 || rho_lower > rho_upper || rho_upper >= 1) {
    stop("`rho_lower` and `rho_upper` must satisfy 0 <= rho_lower <= rho_upper < 1.")
  }

  node_count <- nrow(model)
  root <- model$nodes[model$parent == ""]
  topological_order <- root
  frontier <- root
  while (length(frontier)) {
    frontier <- model$nodes[model$parent %in% frontier]
    topological_order <- c(topological_order, frontier)
  }
  if (length(topological_order) != node_count || anyDuplicated(topological_order)) {
    stop("`model` must be connected and acyclic.")
  }
  if (!requireNamespace("pcalg", quietly = TRUE) ||
      !requireNamespace("graph", quietly = TRUE)) {
    stop("simulate_ltree_data() requires the 'pcalg' and 'graph' packages (see Suggests in DESCRIPTION).")
  }
  has_parent <- model$parent != ""
  edge_count <- sum(has_parent)
  edge_correlations <- runif(edge_count, min = rho_lower, max = rho_upper)
  names(edge_correlations) <- model$nodes[has_parent]

  weighted_adjacency <- matrix(0, node_count, node_count)
  rownames(weighted_adjacency) <- colnames(weighted_adjacency) <- topological_order
  for (i in seq_len(node_count)) {
    if (model$parent[i] == "") next
    weighted_adjacency[model$nodes[i], model$parent[i]] <-
      edge_correlations[model$nodes[i]]
  }
  graph <- methods::as(t(weighted_adjacency), "graphNEL")
  residual_variances <- rep(1, node_count)
  names(residual_variances) <- model$nodes
  residual_variances[names(edge_correlations)] <- 1 - edge_correlations^2
  residual_mean <- rep(0, node_count)
  residual_covariance <- diag(residual_variances)
  residuals <- mvtnorm::rmvnorm(
    n, mean = residual_mean, sigma = residual_covariance
  )
  colnames(residuals) <- model$nodes
  residuals <- residuals[, topological_order, drop = FALSE]
  gaussian_data <- if (node_count == 1L) {
    residuals
  } else {
    pcalg::rmvDAG(n, graph, errMat = residuals)
  }
  colnames(gaussian_data) <- topological_order

  node_type_parts <- do.call("rbind", strsplit(model$type, split = "_"))
  rownames(node_type_parts) <- model$nodes

  gaussian_correlation <- cor(gaussian_data)
  gaussian_correlation <- gaussian_correlation[model$nodes, model$nodes,
                                               drop = FALSE]
  for (j in seq_len(node_count)) {
    current_type <- node_type_parts[topological_order[j], 2]
    if (current_type == "continuous") {
      gaussian_data[, j] <- gaussian_data[, j]^3
    } else if (current_type == "binary") {
      quantile_bounds <- quantile(gaussian_data[, j], c(0.25, 0.75))
      threshold <- runif(1, min = quantile_bounds[1], max = quantile_bounds[2])
      gaussian_data[, j] <- as.numeric(gaussian_data[, j] > threshold)
    } else if (current_type == "ordinal") {
      quantile_bounds <- quantile(gaussian_data[, j], c(0.2, 0.4, 0.6, 0.8))
      lower_threshold <- runif(1, min = quantile_bounds[1], max = quantile_bounds[2])
      upper_threshold <- runif(1, min = quantile_bounds[3], max = quantile_bounds[4])
      gaussian_data[, j] <- (gaussian_data[, j] > lower_threshold) +
        (gaussian_data[, j] > upper_threshold)
    } else if (current_type == "discrete") {
      quantile_bounds <- quantile(
        gaussian_data[, j], c(0.1, 0.3, 0.4, 0.6, 0.7, 0.9)
      )
      lower_threshold <- runif(1, min = quantile_bounds[1], max = quantile_bounds[2])
      middle_threshold <- runif(1, min = quantile_bounds[3], max = quantile_bounds[4])
      upper_threshold <- runif(1, min = quantile_bounds[5], max = quantile_bounds[6])
      gaussian_data[, j] <- (gaussian_data[, j] > lower_threshold) +
        (gaussian_data[, j] > middle_threshold) +
        (gaussian_data[, j] > upper_threshold)
    }
  }
  is_observed <- node_type_parts[, 1] == "observed"
  data <- gaussian_data[, model$nodes[is_observed], drop = FALSE]
  data <- as.data.frame(data)

  list(
    data = data,
    edge_correlations = edge_correlations,
    gaussian_correlation = gaussian_correlation
  )
}
