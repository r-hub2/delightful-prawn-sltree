#' @include sls.R
NULL

#' Create a Latent Tree Object from Edge Set
#'
#' Creates an \code{ltree} object from the edge set of a directed latent tree.
#'
#' @param edge A two-column matrix or data frame, representing the edge set of a directed
#'   latent tree structure. The first column contains the node labels (with observed nodes
#'   labeled as "V1", "V2", ..., and latent nodes as "H1", "H2", ...). The second column
#'   contains the parent nodes, with the root node having an empty string as its parent.
#' @param node_types A character vector of length equal to the number of nodes,
#'   indicating node types. If named, names must be exactly the node identifiers
#'   and values are reordered to match \code{edge}; otherwise values follow row
#'   order. Identifiers beginning with \code{"H"} must have latent type.
#'
#' @return A list of class \code{ltree} with components
#' \tabular{ll}{
#'   \code{observed_nodes} \tab A character vector containing all observed nodes
#'                         of the latent tree structure. \cr
#'   \code{latent_nodes} \tab A character vector containing all latent nodes of
#'                         the latent tree structure. \cr
#'   \code{internal_nodes} \tab A character vector containing all internal nodes
#'                           of the latent tree structure. \cr
#'   \code{node_metadata} \tab A data frame with 4 columns, where:
#'                   The first column is a character vector containing all
#'                   nodes of the latent tree model;
#'                   The second column records \code{"latent"} for latent
#'                   nodes and \code{"observed_<type>"} for observed nodes,
#'                   such as \code{"observed_continuous"},
#'                   \code{"observed_binary"}, or \code{"observed_ordinal"};
#'                   The third column is an integer vector indicating the
#'                   degree of each node;
#'                   The fourth column is an integer vector indicating the
#'                   number of observed nodes among the child nodes of each
#'                   node. \cr
#'   \code{parent_by_node} \tab A character vector indicating the parent node of each
#'                  node. This is a named vector, where the \code{names}
#'                  attribute is all nodes.
#'                  \code{""} means that the node has no parent node. \cr
#'   \code{children_by_node} \tab A named list in which each element is a set of child
#'                     nodes of this node. \cr
#'   \code{observed_descendants} \tab A named list in which each element is a set of
#'                         observed descendant nodes of this node. \cr
#'   \code{active_nodes} \tab A character vector indicating the active set in the
#'                 current step of structure learning algorithm, that is, the
#'                 nodes whose structure need to be explored in the current
#'                 step. \cr
#'   \code{discovered_nodes} \tab A character vector containing the nodes whose
#'                 structure has been found. \cr
#'   \code{distance} \tab A placeholder matrix for the observed nodes, with
#'                        zero diagonal and all off-diagonal entries equal to
#'                        one. These values are not information distances
#'                        estimated from the edge list. \cr
#'   \code{node_types} \tab A named character vector with node identifiers
#'                          as names. Types are \code{"continuous"},
#'                          \code{"Gaussian"}, \code{"binary"},
#'                          \code{"ordinal"}, \code{"count"},
#'                          \code{"categorical"}, \code{"nominal"},
#'                          \code{"discrete"}, or \code{"latent"}. \cr
#'   \code{node_names} \tab A named character vector indicate nodes name,
#'                          where its \code{names} attribute vector is nodes
#'                          labels. \cr
#' }
#'
#' The components \code{active_nodes}, \code{discovered_nodes}, and
#' \code{distance} are
#'   used in structure learning algorithm.
#'
#' @examples
#' edges <- data.frame(
#'   node = c("V1", "V2", "V3", "V4", "V5", "H1", "H2"),
#'   parent = c("H1", "V1", "H1", "H2", "H2", "", "H1"),
#'   stringsAsFactors = FALSE
#' )
#' node_types <- rep(c("ordinal", "continuous", "latent"), c(2, 3, 2))
#' tree <- ltree_from_edge(edges, node_types)
#' @export
ltree_from_edge <- function(edge, node_types = NULL) {
  if (!(is.matrix(edge) || is.data.frame(edge)) || ncol(edge) != 2L ||
      nrow(edge) < 1L) {
    stop("`edge` must be a two-column matrix or data frame with at least one row.")
  }
  edge <- as.data.frame(lapply(as.data.frame(edge), as.character),
                        stringsAsFactors = FALSE)
  colnames(edge) <- c("node", "parent")
  if (anyNA(edge) || any(edge$node == "") || anyDuplicated(edge$node)) {
    stop("`edge` must contain unique, non-empty node identifiers and no missing values.")
  }

  # Add the root row if the edge list omits it.
  if (!any(edge$parent == "")) {
    root <- setdiff(edge$parent, edge$node)
    root <- root[root != ""]
    if (length(root) != 1L) {
      stop("`edge` must identify exactly one root.")
    }
    edge <- rbind(edge, data.frame(node = root, parent = "",
                                   stringsAsFactors = FALSE))
  }
  root <- edge$node[edge$parent == ""]
  if (length(root) != 1L ||
      any(!edge$parent[edge$parent != ""] %in% edge$node) ||
      any(edge$node == edge$parent)) {
    stop("`edge` must describe a rooted tree with exactly one root and valid parents.")
  }

  reached <- root
  frontier <- root
  while (length(frontier)) {
    frontier <- edge$node[edge$parent %in% frontier]
    reached <- c(reached, frontier)
  }
  if (length(reached) != nrow(edge) || anyDuplicated(reached)) {
    stop("`edge` must describe a connected, acyclic tree.")
  }

  tree <- .new_ltree()

  # Observed nodes
  is_latent <- startsWith(edge$node, "H")
  observed_nodes <- edge$node[!is_latent]
  tree$observed_nodes <- observed_nodes

  # Latent nodes
  latent_nodes <- edge$node[is_latent]
  tree$latent_nodes <- latent_nodes

  # Node type, degree, and observed child count
  node_ids <- edge$node
  node_count <- nrow(edge)
  observed_count <- sum(!is_latent)
  if (observed_count < 1L) {
    stop("`edge` must contain at least one observed node (an identifier not beginning with 'H').")
  }

  if (is.null(node_types)) {
    node_types <- c(rep("continuous", observed_count),
                    rep("latent", node_count - observed_count))
    names(node_types) <- c(observed_nodes, latent_nodes)
  } else {
    if (!is.character(node_types) || length(node_types) != node_count || anyNA(node_types)) {
      stop("`node_types` must be NULL or a non-missing character vector with one value per node.")
    }
    if (!is.null(names(node_types))) {
      if (any(names(node_types) == "") || anyDuplicated(names(node_types)) ||
          !setequal(names(node_types), node_ids)) {
        stop("Named `node_types` must have exactly the node identifiers in `edge`.")
      }
      node_types <- node_types[node_ids]
    } else {
      names(node_types) <- node_ids
    }
    node_types <- .normalize_node_types(node_types)
    names(node_types) <- node_ids
    if (any(node_types[latent_nodes] != "latent") ||
        any(node_types[observed_nodes] == "latent")) {
      stop("Node identifiers beginning with 'H' must be latent, and all other identifiers observed.")
    }
  }
  node_types <- node_types[node_ids]

  # Degree of each node
  degree <- integer(node_count)
  names(degree) <- node_ids
  has_parent <- edge$parent != ""
  degree_counts <- table(c(edge$node[has_parent], edge$parent[has_parent]))
  degree[names(degree_counts)] <- as.numeric(degree_counts)

  # Number of observed children for each node
  observed_child_count <- integer(node_count)
  names(observed_child_count) <- node_ids
  observed_child_counts <- table(
    edge$parent[has_parent & edge$node %in% observed_nodes]
  )
  observed_child_count[names(observed_child_counts)] <-
    as.numeric(observed_child_counts)

  node_kind <- node_types
  node_kind[node_kind != "latent"] <- paste0("observed_", node_kind[node_kind != "latent"])
  names(node_kind) <- names(degree) <- names(observed_child_count) <- NULL
  node_metadata <- data.frame(
    node_id = node_ids,
    node_kind = node_kind,
    degree = degree,
    observed_child_count = observed_child_count,
    stringsAsFactors = FALSE
  )
  row.names(node_metadata) <- node_ids
  tree$node_metadata <- node_metadata

  # Internal nodes
  internal_nodes <- node_ids[node_ids %in% unique(edge$parent[edge$parent != ""])]
  tree$internal_nodes <- internal_nodes

  # Parent nodes of each node
  parent_by_node <- edge$parent
  names(parent_by_node) <- edge$node
  tree$parent_by_node <- parent_by_node

  # Child nodes of each node
  children_by_node <- vector(mode = "list", length = node_count)
  names(children_by_node) <- node_ids
  for (h in internal_nodes) {
    children_by_node[[h]] <- edge$node[edge$parent == h]
  }
  tree$children_by_node <- children_by_node

  # Descendant nodes of each node
  observed_descendants <- vector(mode = "list", length = node_count)
  names(observed_descendants) <- node_ids
  collect_descendants <- function(node, children_by_node,
                                  observed_descendants = NULL) {
    child_nodes <- children_by_node[[node]]
    if (is.null(child_nodes)) return(NULL)

    internal_children <- child_nodes[child_nodes %in% internal_nodes]
    if (!length(internal_children)) return(child_nodes)

    descendants <- child_nodes
    for (internal_child in internal_children) {
      child_descendants <- observed_descendants[[internal_child]]
      if (is.null(child_descendants)) {
        child_descendants <- collect_descendants(
          internal_child,
          children_by_node = children_by_node,
          observed_descendants = observed_descendants
        )
      }
      descendants <- union(descendants, child_descendants)
    }
    descendants
  }

  for (internal_node in internal_nodes[rev(seq_along(internal_nodes))]) {
    descendants <- collect_descendants(
      internal_node,
      children_by_node = children_by_node,
      observed_descendants = observed_descendants
    )
    observed_descendants[[internal_node]] <-
      intersect(observed_nodes, descendants)
  }
  tree$observed_descendants <- observed_descendants

  # Active nodes (initially root node)
  tree$active_nodes <- root

  # Discovered nodes (everything except root)
  tree$discovered_nodes <- setdiff(node_ids, root)

  # Information distance matrix
  distance <- 1 - diag(observed_count)
  colnames(distance) <- rownames(distance) <- observed_nodes
  tree$distance <- distance

  tree$node_types <- node_types

  node_names <- c(observed_nodes, latent_nodes)
  names(node_names) <- node_names
  tree$node_names <- node_names

  return(tree)
}


#' Compute Robinson-Foulds (rf) Distance Between Two Latent Trees
#'
#' \code{robinson_foulds_distance} computes the Robinson-Foulds distance between two
#' undirected latent tree structures. The rf distance measures the
#' topological difference between trees. A distance of 0 indicates
#' identical tree structures (ignoring root position and hidden node names).
#'
#' @param reference_tree An `ltree` object representing the reference tree.
#' Must include `observed_nodes`, `internal_nodes`, and `observed_descendants`.
#' @param estimated_tree An `ltree` object representing the estimated tree.
#' Must include `observed_nodes`, `internal_nodes`, and `observed_descendants`.
#'
#' @return A non-negative numeric value giving the RF distance.
#'
#' @details This function compares two trees based on their topologies.
#' The names of latent nodes and the root node do not affect the rf distance.
#'
#' @seealso \code{\link{ltree_from_edge}} for constructing latent tree structures.
#' @export
robinson_foulds_distance <- function(reference_tree, estimated_tree) {
  if (!is_ltree(reference_tree) || !is_ltree(estimated_tree)) {
    stop("`reference_tree` and `estimated_tree` must both be 'ltree' objects.")
  }
  if (is.null(reference_tree$observed_nodes) || is.null(estimated_tree$observed_nodes) ||
      is.null(reference_tree$internal_nodes) || is.null(estimated_tree$internal_nodes) ||
      is.null(reference_tree$observed_descendants) ||
      is.null(estimated_tree$observed_descendants)) {
    stop("The two trees must contain `observed_nodes`, `internal_nodes` and `observed_descendants`.")
  }

  reference_observed <- reference_tree$observed_nodes
  estimated_observed <- estimated_tree$observed_nodes
  if (!setequal(reference_observed, estimated_observed)) {
    stop("The two trees must have the same observed nodes.")
  }
  observed_nodes <- sort(reference_observed)

  get_signatures <- function(internal_nodes, observed_descendants, observed_nodes) {
    leaf_nodes <- setdiff(observed_nodes, internal_nodes)
    signatures <- character(length = length(leaf_nodes) + length(internal_nodes))
    names(signatures) <- c(leaf_nodes, internal_nodes)
    complement_signatures <- signatures

    for (node in leaf_nodes) {
      signatures[node] <- node
      complement_signatures[node] <-
        paste(setdiff(observed_nodes, node), sep = "", collapse = ",")
    }
    for (node in internal_nodes) {
      descendants <- sort(intersect(
        observed_nodes, c(observed_descendants[[node]], node)
      ))
      signatures[node] <- paste(descendants, collapse = ",")
      complement <- setdiff(observed_nodes, descendants)
      complement_signatures[node] <- paste(complement, collapse = ",")
    }
    list(signatures = signatures, complements = complement_signatures)
  }

  reference_splits <- get_signatures(
    reference_tree$internal_nodes,
    reference_tree$observed_descendants,
    observed_nodes
  )
  estimated_splits <- get_signatures(
    estimated_tree$internal_nodes,
    estimated_tree$observed_descendants,
    observed_nodes
  )

  reference_matched <-
    reference_splits$signatures %in% estimated_splits$signatures |
    reference_splits$complements %in% estimated_splits$signatures
  estimated_matched <-
    estimated_splits$signatures %in% reference_splits$signatures |
    estimated_splits$complements %in% reference_splits$signatures

  (sum(!reference_matched) + sum(!estimated_matched)) / 2
}

