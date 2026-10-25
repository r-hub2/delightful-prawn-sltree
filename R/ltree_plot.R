# Base-graphics shape helpers for plot.ltree().
.rrect_path <- function(cx, cy, hw, hh, r, k = 12) {
  r  <- min(r, hw, hh)
  t  <- seq(0, pi/2, length.out = k)
  xr <- cx + hw - r; xl <- cx - hw + r
  yt <- cy + hh - r; yb <- cy - hh + r
  list(x = c(xr + r*cos(t), xl + r*cos(t + pi/2),
             xl + r*cos(t + pi), xr + r*cos(t + 3*pi/2)),
       y = c(yt + r*sin(t), yt + r*sin(t + pi/2),
             yb + r*sin(t + pi), yb + r*sin(t + 3*pi/2)))
}
.circle_path <- function(cx, cy, r, n = 64) {
  t <- seq(0, 2*pi, length.out = n)
  list(x = cx + r*cos(t), y = cy + r*sin(t))
}

.merge_plot_args <- function(defaults, overrides, arg) {
  if (!is.list(overrides) || (length(overrides) && is.null(names(overrides))) ||
      any(!nzchar(names(overrides)))) {
    stop("`", arg, "` must be a named list.", call. = FALSE)
  }
  utils::modifyList(defaults, overrides)
}

.tidy_tree_layout <- function(parent_by_node, heights, gap) {
  node_ids <- names(parent_by_node)
  children <- split(node_ids, factor(unname(parent_by_node), levels = node_ids))
  roots <- node_ids[parent_by_node == ""]
  place_subtree <- function(node) {
    is_virtual <- node == ""
    height <- if (is_virtual) 0 else heights[[node]]
    y <- stats::setNames(0, node)
    depth <- stats::setNames(0L, node)
    lower <- numeric(0)
    upper <- numeric(0)
    offsets <- numeric(0)
    for (child in if (is_virtual) roots else children[[node]]) {
      subtree <- place_subtree(child)
      offset <- 0
      shared <- seq_len(min(length(upper), length(subtree$lower)))
      if (length(shared)) {
        offset <- max(upper[shared] - subtree$lower[shared]) + gap
      }
      offsets <- c(offsets, offset)
      y <- c(y, subtree$y + offset)
      depth <- c(depth, subtree$depth + 1L)
      levels_needed <- length(subtree$upper)
      if (levels_needed > length(upper)) {
        added <- levels_needed - length(upper)
        upper <- c(upper, rep(-Inf, added))
        lower <- c(lower, rep(Inf, added))
      }
      level <- seq_len(levels_needed)
      upper[level] <- pmax(upper[level], subtree$upper + offset)
      lower[level] <- pmin(lower[level], subtree$lower + offset)
    }
    centre <- if (length(offsets)) mean(offsets[c(1L, length(offsets))]) else 0
    y[-1] <- y[-1] - centre
    list(y = y, depth = depth,
         lower = c(-height / 2, lower - centre),
         upper = c(height / 2, upper - centre))
  }
  placed <- place_subtree("")
  list(y = placed$y[node_ids], depth = placed$depth[node_ids] - 1L)
}

.ltree_palette <- function() {
  list(
    fill = c(latent = "#6F7CA8", continuous = "#86BE92", discrete = "#E3B36B"),
    border = c(latent = "#4C577D", continuous = "#558662", discrete = "#B5853A")
  )
}

.ltree_node_groups <- function(node_types) {
  node_groups <- rep("discrete", length(node_types))
  node_groups[node_types %in% c("continuous", "Gaussian")] <- "continuous"
  node_groups[node_types %in% "latent"] <- "latent"
  node_groups
}

.resolve_vertex_colors <- function(value, defaults, node_groups, node_names, arg) {
  group_colors <- defaults
  node_colors <- unname(group_colors[node_groups])
  if (is.null(value)) return(list(node = node_colors, group = group_colors))
  if (!is.character(value) || !length(value) || anyNA(value)) {
    stop("`", arg, "` must be a character vector of colours.", call. = FALSE)
  }
  value_names <- names(value)
  value <- unname(value)
  if (is.null(value_names)) {
    if (length(value) == 1L) {
      group_colors[] <- value
      return(list(node = rep(value, length(node_groups)), group = group_colors))
    }
    if (length(value) == length(node_groups)) {
      return(list(node = value, group = group_colors))
    }
    stop("`", arg, "` must be one colour, one colour per node, or a named ",
         "vector.", call. = FALSE)
  }
  if (all(value_names %in% names(group_colors))) {
    group_colors[value_names] <- value
    return(list(node = unname(group_colors[node_groups]), group = group_colors))
  }
  index <- match(value_names, unname(node_names))
  node_colors[index[!is.na(index)]] <- value[!is.na(index)]
  list(node = node_colors, group = group_colors)
}

.static_tree_coordinates <- function(parent_by_node, is_latent, half_width,
                                     half_height, direction,
                                     available = c(2, 2)) {
  horizontal <- direction %in% c("LR", "RL")
  depth_room <- if (horizontal) available[1] else available[2]
  breadth_size <- if (horizontal) {
    ifelse(is_latent, 2 * half_width, 2 * half_height)
  } else rep(2 * half_width, length(is_latent))
  tidy <- .tidy_tree_layout(
    parent_by_node, stats::setNames(breadth_size, names(parent_by_node)),
    gap = 0.6 * half_height
  )
  breadth <- unname(tidy$y)
  depth <- unname(tidy$depth)
  span <- diff(range(breadth))
  if (span > 0 && span < 2) breadth <- breadth * 2 / span
  breadth <- breadth - mean(range(breadth))
  max_depth <- max(depth)
  step <- if (max_depth > 0) {
    max(min(max(2, span), depth_room) / max_depth,
        2 * half_width + 1.5 * half_height)
  } else 0
  level <- depth * step
  level <- level - mean(range(level))
  switch(direction,
         LR = list(x = level, y = -breadth),
         RL = list(x = -level, y = -breadth),
         UD = list(x = breadth, y = -level),
         DU = list(x = breadth, y = level))
}

.title_margins <- function(title_args) {
  c(if (!is.null(title_args$sub)) 5.1 else if (!is.null(title_args$xlab)) 4.1 else 0.4,
    if (!is.null(title_args$ylab)) 4.1 else 0.4,
    if (!is.null(title_args$main)) 2.6 else 0.4,
    0.4)
}

.ltree_legend_spec <- function(node_groups, fill_palette, border_palette) {
  groups <- intersect(c("latent", "continuous", "discrete"), node_groups)
  labels <- c(latent = "Latent", continuous = "Observed (continuous)",
              discrete = "Observed (discrete)")
  shapes <- c(latent = 21, continuous = 22, discrete = 22)
  list(legend = unname(labels[groups]), pch = unname(shapes[groups]),
       pt.bg = unname(fill_palette[groups]),
       col = unname(border_palette[groups]))
}

#' Plot a latent tree
#'
#' Draws a latent tree on the current graphics device. Latent nodes are
#'   circles and observed nodes are boxes coloured by type (continuous green,
#'   discrete amber). A legend of node types is drawn below the tree.
#'
#' Argument names follow \pkg{igraph} (\code{vertex.*}, \code{edge.*}) and base
#'   graphics (\code{main}, \code{xlab}, \code{col.main}, ...). Other arguments
#'   of the underlying graphics functions go in the \code{*.args} lists, which
#'   replace the defaults; a \code{NULL} entry removes a default.
#'
#' @param x An \code{ltree} object.
#' @param vertex.size,edge.width,vertex.label.cex Vertex size, edge width and
#'   label size. Large trees are scaled down together with their labels, so
#'   set label size with \code{vertex.label.cex} rather than \code{label.args}.
#' @param layout Optional two-column matrix of vertex coordinates. By default
#'   nodes are placed by a tree layout that avoids overlaps.
#' @param direction Direction of the default tree layout: \code{"LR"}
#'   (left to right, default), \code{"RL"}, \code{"UD"} (top down), or
#'   \code{"DU"}. Ignored when \code{layout} is supplied.
#' @param vertex.color,vertex.frame.color Fill and border colours of the
#'   vertices, given as one colour for all vertices, one colour per vertex (in
#'   the order of \code{vertices_from_ltree(x)}), a vector named by node name
#'   (\code{x$node_names}) to recolour selected vertices, or a vector named by
#'   node type (\code{"latent"}, \code{"continuous"}, \code{"discrete"}) to
#'   change the colour of a whole type. The legend follows single-colour and
#'   per-type settings.
#' @param legend Logical; whether to draw the node-type legend. Default
#'   \code{TRUE}.
#' @param edge.args Named list of extra arguments for
#'   \code{\link[graphics]{segments}}, which draws the edges (e.g.
#'   \code{list(lty = 2, col = "black")}).
#' @param vertex.args Named list of extra arguments for
#'   \code{\link[graphics]{polygon}}, which draws every vertex shape (e.g.
#'   \code{list(lwd = 2)}).
#' @param label.args Named list of extra arguments for
#'   \code{\link[graphics]{text}}, which draws the vertex labels (e.g.
#'   \code{list(font = 2, family = "serif")}).
#' @param legend.args Named list of extra arguments for
#'   \code{\link[graphics]{legend}} (e.g. \code{list(x = "topright",
#'   title = "Node type")}). By default the legend is a single row centred
#'   below the tree.
#' @param par.args Named list of graphical parameters set with
#'   \code{\link[graphics]{par}} while drawing and restored afterwards (e.g.
#'   \code{list(bg = "white", mar = c(1, 1, 3, 1))}).
#' @param ... Further arguments passed to \code{\link[graphics]{title}}, such
#'   as \code{main}, \code{sub}, \code{xlab}, \code{ylab}, \code{cex.main}, or
#'   \code{col.main}. Margins are enlarged automatically for \code{main},
#'   \code{sub}, \code{xlab}, and \code{ylab}.
#'
#' @return \code{x}, invisibly.
#' @examples
#' \dontrun{
#' tree <- view_ltree(1, type = "summary")
#' plot(tree)
#' plot(tree, main = "Preset tree 1", direction = "UD", legend = FALSE)
#' plot(tree, vertex.color = c(latent = "grey60"), vertex.frame.color = "black")
#' plot(tree, edge.args = list(lty = 2),
#'      legend.args = list(x = "topright", title = "Node type"))
#' }
#'
#' @export
plot.ltree <- function(x, vertex.size = 24, edge.width = 2, vertex.label.cex = 0.9,
                       layout = NULL, direction = c("LR", "RL", "UD", "DU"),
                       vertex.color = NULL, vertex.frame.color = NULL,
                       legend = TRUE, edge.args = list(), vertex.args = list(),
                       label.args = list(), legend.args = list(),
                       par.args = list(), ...) {
  direction <- match.arg(direction)
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("`legend` must be one non-missing logical value.")
  }
  title_args <- list(...)
  edges <- edges_from_ltree(x, use_node_names = TRUE)
  vertex_names <- vertices_from_ltree(x, use_node_names = TRUE)
  node_ids <- names(x$parent_by_node)
  node_types <- x$node_types
  if (!is.null(names(node_types))) node_types <- node_types[node_ids]
  graph <- igraph::graph_from_data_frame(
    edges, directed = FALSE, vertices = vertex_names
  )

  palette <- .ltree_palette()
  node_groups <- .ltree_node_groups(node_types)
  fill <- .resolve_vertex_colors(vertex.color, palette$fill, node_groups,
                                 vertex_names, "vertex.color")
  border <- .resolve_vertex_colors(vertex.frame.color, palette$border,
                                   node_groups, vertex_names,
                                   "vertex.frame.color")
  if (length(vertex.args$col) == 1L) fill$group[] <- vertex.args$col
  if (length(vertex.args$border) == 1L) border$group[] <- vertex.args$border
  fill_colors <- fill$node
  border_colors <- border$node
  label_colors <- ifelse(node_groups == "latent", "#ffffff", "#1f2a17")
  is_latent <- node_groups == "latent"

  half_width <- vertex.size / 200
  half_height <- vertex.size * 0.6 / 200
  rescale_coordinate <- function(values) {
    value_range <- range(values)
    if (diff(value_range) == 0) return(rep(0, length(values)))
    2 * (values - value_range[1]) / diff(value_range) - 1
  }

  par_args <- .merge_plot_args(list(mar = .title_margins(title_args)),
                               par.args, "par.args")
  previous_parameters <- do.call(graphics::par, par_args)
  on.exit(graphics::par(previous_parameters))

  legend_spec <- .ltree_legend_spec(node_groups, fill$group, border$group)
  show_legend <- legend && length(legend_spec$legend) > 0L
  graphics::plot.new()
  plot_size <- graphics::par("pin")
  legend_band <- if (show_legend) min(0.5, plot_size[2] / 2) else 0
  usable_height <- plot_size[2] - legend_band
  base_scale <- max(2.3 / plot_size[1], 2.2 / usable_height)

  if (is.null(layout)) {
    coordinates <- .static_tree_coordinates(
      x$parent_by_node, is_latent, half_width, half_height, direction,
      available = c(plot_size[1], usable_height) * base_scale - c(0.3, 0.2)
    )
    x_coordinates <- coordinates$x
    y_coordinates <- coordinates$y
  } else {
    x_coordinates <- rescale_coordinate(layout[, 1])
    y_coordinates <- rescale_coordinate(layout[, 2])
  }
  x_half <- max(1.15, max(abs(x_coordinates)) + half_width + 0.03)
  y_half <- max(1.1, max(abs(y_coordinates)) + half_width + 0.03)
  scale <- max(2 * x_half / plot_size[1], 2 * y_half / usable_height)
  zoom <- scale / base_scale
  graphics::plot.window(xlim = c(-x_half, x_half),
                        ylim = c(-y_half - legend_band * scale, y_half),
                        asp = 1)

  edge_indices <- igraph::as_edgelist(graph, names = FALSE)
  do.call(graphics::segments, .merge_plot_args(list(
    x0 = x_coordinates[edge_indices[, 1]],
    y0 = y_coordinates[edge_indices[, 1]],
    x1 = x_coordinates[edge_indices[, 2]],
    y1 = y_coordinates[edge_indices[, 2]],
    col = "#b7bdc6", lwd = edge.width
  ), edge.args, "edge.args"))
  for (i in seq_along(vertex_names)) {
    node_path <- if (is_latent[i]) {
      .circle_path(x_coordinates[i], y_coordinates[i], half_width)
    } else {
      .rrect_path(
        x_coordinates[i], y_coordinates[i], half_width, half_height,
        r = min(half_width, half_height) * 0.5
      )
    }
    do.call(graphics::polygon, .merge_plot_args(list(
      x = node_path$x, y = node_path$y, col = fill_colors[[i]],
      border = border_colors[[i]], lwd = 1.5
    ), vertex.args, "vertex.args"))
  }
  do.call(graphics::text, .merge_plot_args(list(
    x = x_coordinates, y = y_coordinates, labels = vertex_names,
    col = label_colors, cex = vertex.label.cex / zoom, family = "sans"
  ), label.args, "label.args"))

  if (show_legend) {
    legend_args <- .merge_plot_args(c(list(
      x = "bottom", horiz = TRUE, bty = "n", pt.cex = 2, pt.lwd = 1.5,
      cex = vertex.label.cex, text.col = "#333333", x.intersp = 0.8
    ), legend_spec), legend.args, "legend.args")
    if (isTRUE(legend_args$horiz) && is.null(legend_args$text.width) &&
        getRversion() >= "4.2.0") {
      legend_cex <- if (is.null(legend_args$cex)) 1 else legend_args$cex
      legend_args$text.width <- graphics::strwidth(
        paste0(legend_args$legend, "m"), cex = legend_cex
      )
    }
    do.call(graphics::legend, legend_args)
  }
  if (length(title_args)) do.call(graphics::title, title_args)
  invisible(x)
}

#' Interactive latent tree plot
#'
#' Draws an \code{ltree} object as a \pkg{visNetwork} HTML widget.
#'
#' @param tree An \code{ltree} object.
#' @param direction \code{"UD"} (default), \code{"DU"}, \code{"LR"} or
#'   \code{"RL"}.
#' @param layout \code{"tree"} (default) uses igraph's Reingold-Tilford
#'   layout; \code{"hierarchical"} uses visNetwork's hierarchical layout.
#' @param edge_style \code{"curved"} (default) or \code{"straight"}.
#' @param legend Logical; whether to add a node-type legend. Default
#'   \code{TRUE}.
#' @param legend.args Named list of extra arguments for
#'   \code{visNetwork::visLegend} (e.g. \code{list(position = "left")}).
#' @param vertex.color,vertex.frame.color Fill and border colours, in the same
#'   forms as for \code{\link{plot.ltree}}.
#' @param vertex.args Named list of extra arguments for
#'   \code{visNetwork::visNodes} (e.g. \code{list(font = list(size = 20))}).
#' @param edge.args Named list of extra arguments for
#'   \code{visNetwork::visEdges} (e.g. \code{list(width = 3)}).
#' @param level_sep,node_spacing,tree_spacing Gaps between depths, between
#'   siblings, and between subtrees (\code{"hierarchical"} only).
#' @param height,width Widget size passed to \code{visNetwork::visNetwork}.
#' @param ... Passed to \code{visNetwork::visNetwork}.
#'
#' @return A \code{visNetwork} htmlwidget.
#' @noRd
plot_ltree_interactive <- function(tree, direction = c("UD", "DU", "LR", "RL"),
                      layout = c("tree", "hierarchical"),
                      edge_style = c("curved", "straight"), legend = TRUE,
                      legend.args = list(), vertex.color = NULL,
                      vertex.frame.color = NULL, vertex.args = list(),
                      edge.args = list(), level_sep = 150,
                      node_spacing = 130, tree_spacing = 60,
                      height = "97vh", width = "100%", ...) {
  if (!is_ltree(tree)) stop("`tree` must be an 'ltree' object.")
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("plot_ltree_interactive() requires the 'visNetwork' package (in Suggests); ",
         "install it with install.packages('visNetwork').")
  }
  direction  <- match.arg(direction)
  layout     <- match.arg(layout)
  edge_style <- match.arg(edge_style)

  node_ids <- names(tree$parent_by_node)
  node_types <- tree$node_types[node_ids]
  node_labels <- tree$node_names[node_ids]
  if (is.null(node_labels)) node_labels <- node_ids
  node_labels <- unname(node_labels)

  palette <- .ltree_palette()
  node_groups <- .ltree_node_groups(node_types)
  fill <- .resolve_vertex_colors(vertex.color, palette$fill, node_groups,
                                 node_labels, "vertex.color")
  border <- .resolve_vertex_colors(vertex.frame.color, palette$border,
                                   node_groups, node_labels,
                                   "vertex.frame.color")
  fill_colors <- fill$node

  nodes <- data.frame(
    id = node_ids, label = node_labels, group = node_groups,
    shape = ifelse(node_groups == "latent", "circle", "box"),
    color.background = fill_colors,
    color.border = border$node,
    color.highlight.background = fill_colors,
    color.highlight.border = "#2f3340",
    font.color = ifelse(node_groups == "latent", "#ffffff", "#1f2a17"),
    stringsAsFactors = FALSE)
  edge_df <- edges_from_ltree(tree, use_node_names = FALSE)
  edges <- data.frame(from = edge_df$parent, to = edge_df$nodes,
                      stringsAsFactors = FALSE)

  if (layout == "tree") {
    graph <- igraph::graph_from_data_frame(
      edges, directed = TRUE, vertices = data.frame(name = nodes$id)
    )
    root <- setdiff(nodes$id, edges$to)
    root_index <- match(root, nodes$id)
    horizontal_position <- igraph::layout_as_tree(
      graph, root = root_index
    )[, 1]
    depth <- apply(
      igraph::distances(graph, v = root_index, mode = "out"), 2, min
    )
    depth[!is.finite(depth)] <- max(c(depth[is.finite(depth)], 0)) + 1
    if (direction %in% c("UD", "DU")) {
      nodes$x <- horizontal_position * node_spacing
      nodes$y <- depth * level_sep * (if (direction == "UD") 1 else -1)
    } else {
      nodes$x <- depth * level_sep * (if (direction == "LR") 1 else -1)
      nodes$y <- horizontal_position * node_spacing
    }
  }

  force_direction <- if (direction %in% c("UD", "DU")) {
    "vertical"
  } else "horizontal"
  edge_smoothing <- if (edge_style == "curved")
    list(enabled = TRUE, type = "cubicBezier", forceDirection = force_direction,
         roundness = 0.5) else FALSE

  widget <- visNetwork::visNetwork(nodes, edges, height = height, width = width, ...)
  widget <- do.call(visNetwork::visNodes, .merge_plot_args(list(
    graph = widget, borderWidth = 1.5, margin = 10,
    shadow = list(enabled = TRUE, size = 10, x = 2, y = 4,
                  color = "rgba(0,0,0,0.18)"),
    font = list(size = 15, face = "Helvetica")
  ), vertex.args, "vertex.args"))
  if (layout == "hierarchical") {
    widget <- visNetwork::visHierarchicalLayout(
      widget, direction = direction, sortMethod = "directed",
      levelSeparation = level_sep, nodeSpacing = node_spacing,
      treeSpacing = tree_spacing, blockShifting = TRUE,
      edgeMinimization = TRUE, parentCentralization = TRUE)
  }
  widget <- do.call(visNetwork::visEdges, .merge_plot_args(list(
    graph = widget, width = 1.4,
    color = list(color = "#b7bdc6", highlight = "#7c8694"),
    smooth = edge_smoothing
  ), edge.args, "edge.args"))
  widget <- visNetwork::visPhysics(widget, enabled = FALSE)
  widget <- visNetwork::visInteraction(widget, hover = TRUE, tooltipDelay = 120,
                                   dragView = TRUE, zoomView = TRUE,
                                   navigationButtons = TRUE, keyboard = TRUE)
  widget <- visNetwork::visOptions(
    widget, highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE))
  if (legend) {
    legend_groups <- intersect(
      c("latent", "continuous", "discrete"), unique(node_groups)
    )
    legend_labels <- c(
      latent = "latent", continuous = "observed (cont.)",
      discrete = "observed (disc.)"
    )
    legend_shapes <- c(latent = "circle", continuous = "box", discrete = "box")
    legend_nodes <- data.frame(
      id = seq_along(legend_groups),
      label = unname(legend_labels[legend_groups]),
      shape = unname(legend_shapes[legend_groups]),
      color.background = unname(fill$group[legend_groups]),
      color.border = unname(border$group[legend_groups]),
      font.color = "#333333",
      stringsAsFactors = FALSE
    )
    widget <- do.call(visNetwork::visLegend, .merge_plot_args(list(
      graph = widget, addNodes = legend_nodes, useGroups = FALSE,
      position = "right", width = 0.16, zoom = FALSE
    ), legend.args, "legend.args"))
  }
  widget
}

