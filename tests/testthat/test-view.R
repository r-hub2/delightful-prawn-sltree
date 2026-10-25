test_that("view_ltree resolves preset numbers and recovered trees", {
  capture.output(preset <- view_ltree(1, type = "summary"))
  expect_true(is_ltree(preset))

  capture.output(shown <- view_ltree(preset, type = "summary"))
  expect_identical(shown, preset)
})

test_that("view_ltree rejects invalid inputs", {
  expect_error(view_ltree(0, type = "summary"), "1 to 6")
  expect_error(view_ltree(7, type = "summary"), "1 to 6")
  expect_error(view_ltree(data.frame(x = 1), type = "summary"),
               "'ltree' object")
})

test_that("save_ltree defaults to PNG, EPS, and HTML", {
  expect_identical(formals(save_ltree)$html, TRUE)
  expect_identical(formals(save_ltree)$eps, TRUE)
  expect_identical(formals(save_ltree)$png, TRUE)
})

test_that("save_ltree strips a known extension from the output base name", {
  tree <- latent_tree_models(1)
  node_types <- stats::setNames(tree[, 2], tree[, 1])
  tree <- ltree_from_edge(tree[, c(1, 3)], node_types)
  prefix <- tempfile("sltree-save-")

  written <- save_ltree(tree, paste0(prefix, ".png"),
                        html = FALSE, eps = FALSE)
  expect_identical(unname(written["png"]), paste0(prefix, ".png"))
  expect_true(file.exists(paste0(prefix, ".png")))
  unlink(paste0(prefix, ".png"))
})

preset_ltree <- function(k = 1) {
  model <- latent_tree_models(k)
  ltree_from_edge(model[, c(1, 3)], stats::setNames(model[, 2], model[, 1]))
}

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  code
}

test_that("plot legend lists only the node groups present", {
  palette <- c(latent = "#111111", continuous = "#222222", discrete = "#333333")
  spec <- sltree:::.ltree_legend_spec(c("discrete", "latent", "discrete"),
                                      palette, palette)
  expect_identical(spec$legend, c("Latent", "Observed (discrete)"))
  expect_identical(spec$pch, c(21, 22))
  expect_identical(spec$pt.bg, c("#111111", "#333333"))
  expect_identical(formals(plot.ltree)$legend, TRUE)
})

test_that("plot.ltree forwards extra arguments to the drawing functions", {
  tree <- preset_ltree(1)
  with_null_device({
    old_mar <- graphics::par("mar")
    expect_silent(plot(tree))
    expect_silent(plot(tree, legend = FALSE))
    expect_silent(plot(tree, main = "Preset 1", sub = "subtitle",
                       cex.main = 1.2,
                       edge.args = list(lty = 2, col = "black"),
                       vertex.args = list(lwd = 2),
                       label.args = list(font = 2),
                       legend.args = list(x = "topright", title = "Node type"),
                       par.args = list(bg = "white")))
    expect_equal(graphics::par("mar"), old_mar)
    expect_identical(view_ltree(tree, "plot", main = "via view_ltree"), tree)
  })
})

test_that("plot.ltree validates legend and argument lists", {
  tree <- preset_ltree(1)
  with_null_device({
    expect_error(plot(tree, legend = NA), "`legend` must be")
    expect_error(plot(tree, edge.args = c(lty = 2)), "`edge.args` must be a named list")
    expect_error(plot(tree, legend.args = list("topright")),
                 "`legend.args` must be a named list")
  })
})

test_that("save_ltree forwards plot_args to the static images", {
  prefix <- tempfile("sltree-args-")
  expect_message(
    written <- save_ltree(preset_ltree(1), prefix, html = FALSE, eps = FALSE,
                          plot_args = list(main = "Saved", legend = FALSE)),
    "Wrote"
  )
  expect_true(file.exists(written[["png"]]))
  unlink(written)
  expect_error(save_ltree(preset_ltree(1), prefix, html = FALSE, eps = FALSE,
                          plot_args = "main"),
               "`plot_args` must be a named list")
})

test_that("interactive view shows a legend by default and forwards legend.args", {
  skip_if_not_installed("visNetwork")
  widget <- view_ltree(1, "interactive")
  expect_identical(widget$x$legend$position, "right")
  expect_identical(widget$x$legend$nodes$label,
                   c("latent", "observed (cont.)", "observed (disc.)"))

  moved <- view_ltree(1, "interactive", legend.args = list(position = "left"))
  expect_identical(moved$x$legend$position, "left")
  expect_null(view_ltree(1, "interactive", legend = FALSE)$x$legend)
})

test_that("tidy layout keeps every column free of overlapping nodes", {
  gap <- 0.05
  for (k in seq_along(latent_tree_models())) {
    tree <- preset_ltree(k)
    is_latent <- tree$node_types[names(tree$parent_by_node)] == "latent"
    heights <- stats::setNames(ifelse(is_latent, 0.24, 0.144),
                               names(tree$parent_by_node))
    tidy <- sltree:::.tidy_tree_layout(tree$parent_by_node, heights, gap)
    expect_identical(names(tidy$y), names(tree$parent_by_node))
    for (column in split(names(tidy$y), tidy$depth)) {
      if (length(column) < 2L) next
      column <- column[order(tidy$y[column])]
      lower <- column[-length(column)]
      upper <- column[-1L]
      clearance <- tidy$y[upper] - tidy$y[lower] -
        (heights[upper] + heights[lower]) / 2
      expect_true(all(clearance >= gap - 1e-9), info = paste("preset", k))
    }
  }
})

test_that("tidy layout places roots at depth zero and centres parents", {
  parent <- c(R = "", A = "R", B = "R", C = "A")
  tidy <- sltree:::.tidy_tree_layout(parent, c(R = 1, A = 1, B = 1, C = 1), 0)
  expect_identical(unname(tidy$depth), c(0L, 1L, 1L, 2L))
  expect_equal(unname(tidy$y["R"]), mean(tidy$y[c("A", "B")]))
  expect_equal(unname(tidy$y["C"]), unname(tidy$y["A"]))
})

test_that("view_ltree defaults to the interactive view", {
  expect_identical(eval(formals(view_ltree)$type)[1], "interactive")
  skip_if_not_installed("visNetwork")
  expect_s3_class(view_ltree(1), "visNetwork")
})

test_that("default view falls back to a static plot without visNetwork", {
  local_mocked_bindings(.has_visnetwork = function() FALSE)
  tree <- preset_ltree(1)
  with_null_device({
    expect_message(shown <- view_ltree(tree, main = "fallback",
                                       direction = "LR"),
                   "drawing a static plot instead")
    expect_identical(shown, tree)
  })
})

test_that("vertex colours accept igraph forms and per-type names", {
  palette <- sltree:::.ltree_palette()$fill
  groups <- c("latent", "continuous", "discrete", "latent")
  names <- c("H1", "V1", "V2", "H2")
  resolve <- function(value) {
    sltree:::.resolve_vertex_colors(value, palette, groups, names, "vertex.color")
  }
  expect_identical(resolve(NULL)$node, unname(palette[groups]))
  expect_identical(resolve("red")$node, rep("red", 4))
  expect_identical(unname(resolve("red")$group), rep("red", 3))
  expect_identical(resolve(c("a", "b", "c", "d"))$node, c("a", "b", "c", "d"))
  by_type <- resolve(c(latent = "grey"))
  expect_identical(by_type$node, c("grey", palette[["continuous"]],
                                   palette[["discrete"]], "grey"))
  expect_identical(by_type$group[["latent"]], "grey")
  by_node <- resolve(c(V2 = "pink"))
  expect_identical(by_node$node[3], "pink")
  expect_identical(by_node$group, palette)
  expect_error(resolve(c("a", "b")), "one colour per node")
  expect_error(resolve(1), "character vector of colours")
})

test_that("static layout honours every direction", {
  tree <- preset_ltree(1)
  is_latent <- tree$node_types[names(tree$parent_by_node)] == "latent"
  root <- names(tree$parent_by_node)[tree$parent_by_node == ""]
  coordinates <- function(direction) {
    xy <- sltree:::.static_tree_coordinates(tree$parent_by_node, is_latent,
                                            0.12, 0.072, direction)
    lapply(xy, stats::setNames, names(tree$parent_by_node))
  }
  expect_equal(coordinates("LR")$x[[root]], min(coordinates("LR")$x))
  expect_equal(coordinates("RL")$x[[root]], max(coordinates("RL")$x))
  expect_equal(coordinates("UD")$y[[root]], max(coordinates("UD")$y))
  expect_equal(coordinates("DU")$y[[root]], min(coordinates("DU")$y))
  with_null_device({
    for (direction in c("LR", "RL", "UD", "DU")) {
      expect_silent(plot(tree, direction = direction))
    }
    expect_error(plot(tree, direction = "up"), "should be one of")
  })
})

test_that("axis titles get margins and legend follows type colours", {
  expect_identical(sltree:::.title_margins(list()), c(0.4, 0.4, 0.4, 0.4))
  expect_identical(sltree:::.title_margins(list(xlab = "x", ylab = "y")),
                   c(4.1, 4.1, 0.4, 0.4))
  expect_identical(sltree:::.title_margins(list(main = "m", sub = "s")),
                   c(5.1, 0.4, 2.6, 0.4))
  with_null_device(
    expect_silent(plot(preset_ltree(1), xlab = "x", ylab = "y",
                       vertex.color = c(latent = "grey60"),
                       vertex.frame.color = "black"))
  )
})

test_that("save_ltree forwards device arguments", {
  prefix <- tempfile("sltree-device-")
  written <- suppressMessages(save_ltree(
    preset_ltree(1), prefix, html = FALSE, width = 4, height = 3,
    png_args = list(res = 100), eps_args = list(family = "serif")
  ))
  header <- readBin(written[["png"]], "raw", 24L)
  png_width <- sum(as.integer(header[17:20]) * 256^(3:0))
  expect_equal(png_width, 400)
  expect_true(file.exists(written[["eps"]]))
  unlink(written)
  expect_error(save_ltree(preset_ltree(1), prefix, html = FALSE,
                          png_args = list(300)),
               "`png_args` must be a named list")
})

test_that("interactive view forwards vertex, edge, and html arguments", {
  skip_if_not_installed("visNetwork")
  widget <- view_ltree(1, vertex.args = list(font = list(size = 22)),
                       edge.args = list(width = 3, color = list(color = "black")),
                       vertex.color = c(latent = "grey60"))
  nodes <- widget$x$options$nodes
  edges <- widget$x$options$edges
  expect_identical(nodes$font$size, 22)
  expect_identical(nodes$font$face, "Helvetica")
  expect_identical(edges$width, 3)
  expect_identical(edges$color$color, "black")
  expect_identical(edges$color$highlight, "#7c8694")
  expect_identical(widget$x$legend$nodes$color.background[1], "grey60")

  prefix <- tempfile("sltree-html-")
  written <- suppressMessages(save_ltree(
    preset_ltree(1), prefix, eps = FALSE, png = FALSE,
    html_args = list(selfcontained = FALSE)
  ))
  expect_true(file.exists(written[["html"]]))
  unlink(c(written, paste0(prefix, "_files")), recursive = TRUE)
})
