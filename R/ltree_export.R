#' Save a latent tree to HTML, EPS and/or PNG
#'
#' \code{save_ltree} writes an interactive HTML file with \pkg{visNetwork} and
#'   EPS and PNG images with \code{\link{plot.ltree}}.
#'
#' @param tree A \code{ltree} object, or a single integer from 1 to 6 selecting
#'   a tree from \code{\link{latent_tree_models}}.
#' @param file Output path, preferably without an extension. By default
#'   \code{"<file>.png"}, \code{"<file>.eps"}, and \code{"<file>.html"} are
#'   all written. A trailing \code{.png}, \code{.eps}, or \code{.html} is
#'   removed before the three output names are constructed.
#' @param html,eps,png Logical; which formats to write. All three default to
#'   \code{TRUE}. Set a format to \code{FALSE} to suppress it.
#' @param width,height EPS/PNG size in inches (PNG is rendered at 150 dpi
#'   unless \code{png_args} sets \code{res}).
#' @param plot_args Named list of extra arguments forwarded to
#'   \code{\link{plot.ltree}} for the EPS and PNG images (e.g.
#'   \code{list(main = "True tree", direction = "UD")}).
#' @param png_args Named list of extra arguments for
#'   \code{\link[grDevices]{png}} (e.g. \code{list(res = 300, bg = "white")}).
#' @param eps_args Named list of extra arguments for
#'   \code{\link[grDevices]{postscript}} (e.g. \code{list(family = "serif")}).
#' @param html_args Named list of extra arguments for
#'   \code{visNetwork::visSave} (e.g. \code{list(selfcontained = FALSE)}, which
#'   also avoids the need for pandoc).
#' @param ... Extra arguments forwarded to the interactive \pkg{visNetwork}
#'   viewer for the HTML; see the interactive options in
#'   \code{\link{view_ltree}} (e.g. \code{direction}, \code{legend}).
#' @return Invisibly, a named character vector of the files written.
#' @seealso \code{\link{plot.ltree}}, \code{\link{view_ltree}}
#' @examples
#' \dontrun{
#' model <- latent_tree_models(1)
#' true_tree <- ltree_from_edge(model[, -2], setNames(model[, 2], model[, 1]))
#' save_ltree(true_tree, "true_tree")
#' save_ltree(1, "preset1")
#' save_ltree(true_tree, "png_only", html = FALSE, eps = FALSE)
#' save_ltree(true_tree, "titled", plot_args = list(main = "True tree"))
#' save_ltree(true_tree, "print", html = FALSE, png_args = list(res = 300))
#' }
#' @importFrom grDevices setEPS postscript png dev.off
#' @importFrom graphics plot
#' @export
save_ltree <- function(tree, file = "tree", html = TRUE, eps = TRUE, png = TRUE,
                       width = 12, height = 7, plot_args = list(),
                       png_args = list(), eps_args = list(),
                       html_args = list(), ...) {
  tree <- .resolve_ltree(tree, arg = "tree")
  if (!any(c(html, eps, png))) {
    stop("At least one of `html`, `eps`, or `png` must be TRUE.")
  }
  plot_args <- .merge_plot_args(list(x = tree), plot_args, "plot_args")
  if (!is.character(file) || length(file) != 1L || is.na(file) ||
      !nzchar(file)) {
    stop("`file` must be one non-empty character string.")
  }
  file <- sub("\\.(png|eps|html)$", "", file, ignore.case = TRUE)
  png_args <- .merge_plot_args(list(
    filename = paste0(file, ".png"), width = width, height = height,
    units = "in", res = 150
  ), png_args, "png_args")
  eps_args <- .merge_plot_args(list(
    file = paste0(file, ".eps"), width = width, height = height,
    horizontal = FALSE, paper = "special",
    family = "sans", fonts = c("serif", "sans", "mono")
  ), eps_args, "eps_args")
  html_args <- .merge_plot_args(list(file = paste0(file, ".html")),
                                html_args, "html_args")
  written <- character(0)
  if (html) {
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
      warning("HTML output requires the 'visNetwork' package; skipping .html.")
    } else {
      do.call(visNetwork::visSave,
              c(list(graph = plot_ltree_interactive(tree, ...)), html_args))
      written["html"] <- html_args$file
    }
  }
  if (eps) {
    grDevices::setEPS()
    do.call(grDevices::postscript, eps_args)
    tryCatch(do.call(plot, plot_args), finally = grDevices::dev.off())
    written["eps"] <- eps_args$file
  }
  if (png) {
    do.call(grDevices::png, png_args)
    tryCatch(do.call(plot, plot_args), finally = grDevices::dev.off())
    written["png"] <- png_args$filename
  }
  message("Wrote: ", paste(written, collapse = ", "))
  invisible(written)
}

