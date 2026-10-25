#' Is this object an ltree latent tree?
#'
#' Checks whether an object has class \code{ltree}.
#'
#' @param tree An R object.
#' @return \code{TRUE} or \code{FALSE}.
#' @export
is_ltree <- function(tree) {
  inherits(tree, "ltree")
}


# Resolve an ltree object or preset-model number.
.resolve_ltree <- function(value, arg = "tree") {
  if (inherits(value, "ltree")) return(value)

  is_preset <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value == as.integer(value) && value >= 1L && value <= 6L
  if (is_preset) {
    model <- latent_tree_models(as.integer(value))
    node_types <- stats::setNames(model$type, model$nodes)
    return(ltree_from_edge(
      model[, c("nodes", "parent"), drop = FALSE],
      node_types
    ))
  }

  stop(sprintf("`%s` must be an 'ltree' object or a preset model number from 1 to 6.",
               arg))
}


#' View a latent tree
#'
#' Shows an \code{ltree} object or one of the six preset trees as an
#' interactive HTML widget, a static plot, or a printed summary.
#'
#' @param tree An \code{ltree} object, or a single integer from 1 to 6 selecting a
#'   tree from \code{\link{latent_tree_models}}.
#' @param type Display type: \code{"interactive"} (default) returns the
#'   interactive \pkg{visNetwork} widget, \code{"plot"} draws a static image,
#'   and \code{"summary"} prints the tree structure. When \code{type} is not
#'   given and \pkg{visNetwork} is not installed, a static plot is drawn
#'   instead, with a message.
#' @param ... Display options: see the section \emph{Interactive options} for
#'   \code{type = "interactive"}, and \code{\link{plot.ltree}} for
#'   \code{type = "plot"}. When the default falls back to a static plot, only
#'   the options shared by both views (\code{direction}, \code{legend},
#'   \code{vertex.color}, \code{vertex.frame.color}, \code{main}, and
#'   \code{sub}) are passed on.
#'
#' @section Interactive options:
#' \describe{
#'   \item{\code{direction}}{\code{"UD"} (top down, default), \code{"DU"},
#'     \code{"LR"}, or \code{"RL"}.}
#'   \item{\code{legend}, \code{legend.args}}{Whether to show the node-type
#'     legend (default \code{TRUE}), and a named list of extra arguments for
#'     \code{visNetwork::visLegend} (e.g. \code{list(position = "left")}).}
#'   \item{\code{vertex.color}, \code{vertex.frame.color}}{Fill and border
#'     colours, in the same forms as for \code{\link{plot.ltree}}.}
#'   \item{\code{vertex.args}, \code{edge.args}}{Named lists of extra
#'     arguments for \code{visNetwork::visNodes} and
#'     \code{visNetwork::visEdges} (e.g. \code{list(font = list(size = 20))}
#'     or \code{list(width = 3, color = list(color = "black"))}). Nested lists
#'     are merged with the defaults.}
#'   \item{\code{layout}}{\code{"tree"} (default; fixed Reingold-Tilford
#'     coordinates) or \code{"hierarchical"} (visNetwork's own layout).}
#'   \item{\code{edge_style}}{\code{"curved"} (default) or \code{"straight"}.}
#'   \item{\code{level_sep}, \code{node_spacing}, \code{tree_spacing}}{Spacing
#'     in pixels between depths, between siblings, and between subtrees
#'     (\code{"hierarchical"} only). Defaults 150, 130, and 60.}
#'   \item{\code{height}, \code{width}}{Widget size; defaults \code{"97vh"}
#'     and \code{"100\%"}.}
#'   \item{other arguments}{Passed to \code{visNetwork::visNetwork}, e.g.
#'     \code{main}, \code{submain}, \code{footer}, \code{background}.}
#' }
#'
#' @return For the interactive view, a \code{visNetwork} htmlwidget, which
#'   is displayed when printed. Otherwise, the displayed \code{ltree} object
#'   invisibly. A preset number is converted to an \code{ltree} before display.
#' @seealso \code{\link{save_ltree}}, \code{\link{learn_sltree}},
#'   \code{\link{latent_tree_models}}
#' @examples
#' \dontrun{
#' view_ltree(1)
#' view_ltree(1, "plot")
#' view_ltree(1, direction = "LR", vertex.args = list(font = list(size = 20)))
#' view_ltree(1, "plot", direction = "UD", main = "Preset tree 1")
#' fit <- learn_sltree(x)
#' view_ltree(fit, "summary")
#' view_ltree(fit)
#' }
#' @export
view_ltree <- function(tree, type = c("interactive", "plot", "summary"), ...) {
  use_default <- missing(type)
  type <- match.arg(type)
  tree <- .resolve_ltree(tree)

  if (type == "interactive" && use_default && !.has_visnetwork()) {
    message("The interactive view needs the 'visNetwork' package; ",
            "drawing a static plot instead. ",
            "Install it with install.packages('visNetwork').")
    plot_args <- list(...)
    shared <- c("direction", "legend", "vertex.color", "vertex.frame.color",
                "main", "sub")
    do.call(graphics::plot,
            c(list(tree), plot_args[intersect(names(plot_args), shared)]))
    return(invisible(tree))
  }

  if (type == "summary") {
    print(tree)
  } else if (type == "interactive") {
    return(plot_ltree_interactive(tree, ...))
  } else {
    plot(tree, ...)
  }

  invisible(tree)
}

.has_visnetwork <- function() {
  requireNamespace("visNetwork", quietly = TRUE)
}

#' Extract Edges and Vertices from a Latent Tree
#'
#' \code{edges_from_ltree} extracts the edges from a latent tree structure
#'   ('ltree' object). \code{vertices_from_ltree} extracts the vertices from a
#'   latent tree. The vertices are either node names or node IDs, depending on
#'   the `use_node_names` parameter.
#' @return \code{edges_from_ltree} returns a data frame with two columns:
#'   "nodes" and "parent", representing the edges. \code{vertices_from_ltree}
#'   returns a character vector containing the vertices (nodes) of the latent tree.
#'
#' @param tree A 'ltree' object representing a latent tree structure.
#' @param use_node_names Logical; return node names instead of node IDs.
#' @export
edges_from_ltree <- function(tree, use_node_names = FALSE) {
  if (!is_ltree(tree)) stop("`tree` must be an 'ltree' object.")
  node_ids <- names(tree$parent_by_node)
  parent_by_node <- unname(tree$parent_by_node)
  if (use_node_names) {
    node_names <- tree$node_names
    if (is.null(node_names)) {
      warning("The nodes name of `tree' is NULL!")
    } else {
      node_ids <- unname(node_names[node_ids])
      has_parent <- parent_by_node != ""
      parent_by_node[has_parent] <- unname(node_names[parent_by_node[has_parent]])
    }
  }
  edges <- data.frame(nodes = node_ids, parent = parent_by_node,
                      stringsAsFactors = FALSE)
  edges <- edges[edges$parent != "", , drop = FALSE]

  return(edges)
}

#' @rdname edges_from_ltree
#' @export
vertices_from_ltree <- function(tree, use_node_names = FALSE) {
  if (!is_ltree(tree)) {
    stop("`tree` must be an 'ltree' object.")
  }

  if (use_node_names) {
    node_names <- tree$node_names
    if (is.null(node_names)) {
      warning("The nodes name of `tree` is NULL!")
      return(names(tree$parent_by_node))
    } else {
      return(node_names[names(tree$parent_by_node)])
    }
  }

  return(names(tree$parent_by_node))
}


#' Print a latent tree
#'
#' Prints the nodes and edges of an \code{ltree} object.
#'
#' @param x An \code{ltree} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.ltree <- function(x, ...) {
  if (!is_ltree(x)) {
    stop("`x` must be an 'ltree' object.")
  }

  if (is.null(x$observed_nodes)) {
    stop("`x` has no observed nodes.")
  }

  observed_nodes <- x$observed_nodes
  has_latent_nodes <- length(x$latent_nodes) > 0L
  latent_nodes <- x$latent_nodes

  if (is.null(x$node_metadata)) {
    stop("`x` has no 'node_metadata' component.")
  }
  node_metadata <- x$node_metadata

  is_redundant <- if (has_latent_nodes) {
    any(node_metadata[latent_nodes, "degree"] < 3)
  } else FALSE

  if (is.null(x$parent_by_node)) {
    stop("`x` has no 'parent_by_node' component.")
  }
  parent_by_node <- x$parent_by_node

  has_single_root <- sum(parent_by_node == "") == 1

  if (is.null(x$node_names) || length(x$node_names) != length(latent_nodes) + length(observed_nodes)) {
    node_names <- c(observed_nodes, if (has_latent_nodes) latent_nodes)
    names(node_names) <- node_names
  } else {
    node_names <- x$node_names
  }

  txt <- "'ltree' object: \n"

  if (!has_single_root) {
    txt <- paste0(txt, "The object is not a tree! \n")
  } else if (!has_latent_nodes) {
    txt <- paste0(txt, "A tree without latent nodes! \n")
  } else {
    if (is_redundant) {
      redundant_nodes <- latent_nodes[node_metadata[latent_nodes, "degree"] < 3]
      txt <- paste0(txt, "A latent tree with redundant nodes (",
                    paste0(redundant_nodes, collapse = ","), ")! \n")
    } else {
      txt <- paste0(txt, "A latent tree with ",
                    length(latent_nodes) + length(observed_nodes), " nodes. \n")
    }
  }

  cat(txt)
  cat("\nAvailable components: \n")
  print(names(x))

  is_obs <- grepl("^observed_", node_metadata[, 2])
  is_con <- node_metadata[, 2] %in% c("observed_continuous", "observed_Gaussian")
  is_dis <- is_obs & !is_con
  obs_dis <- sum(is_dis)
  obs_con <- sum(is_con)
  txt <- "\nVertices: "
  txt <- paste0(txt, length(latent_nodes) + length(observed_nodes), " \n")
  txt <- paste0(txt, "    observed nodes: ", length(observed_nodes), " \n")
  txt <- paste0(txt, "        observed discrete nodes: ", obs_dis)

  if (obs_dis > 0 && obs_dis < 20) {
    temp_char <- node_names[is_dis]
    names(temp_char) <- NULL
    temp_char <- paste0(temp_char, collapse = ", ")
    txt <- paste0(txt, "(", temp_char, ") \n")
  } else {
    txt <- paste0(txt, " \n")
  }

  txt <- paste0(txt, "        observed continuous nodes: ", obs_con)
  if (obs_con > 0 && obs_con < 20) {
    temp_char <- node_names[is_con]
    names(temp_char) <- NULL
    temp_char <- paste0(temp_char, collapse = ", ")
    txt <- paste0(txt, "(", temp_char, ") \n")
  } else {
    txt <- paste0(txt, " \n")
  }

  if (has_latent_nodes) {
    txt <- paste0(txt, "    latent nodes: ", length(latent_nodes), " \n")
  }
  cat(txt)

  txt <- "\nEdges: \n"
  has_parent <- parent_by_node != ""
  if (!any(has_parent)) {
    txt <- paste0(txt, "empty! ")
  } else {
    child_nodes <- names(parent_by_node[has_parent])
    parent_nodes <- unname(parent_by_node[has_parent])
    edge_text <- paste0("    ", node_names[child_nodes], " -- ",
                        node_names[parent_nodes], collapse = "\n")
    txt <- paste0(txt, edge_text)
  }
  cat(txt)

  invisible(x)
}



