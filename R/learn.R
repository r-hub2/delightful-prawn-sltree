#' Learn a latent tree under one of four statistical models
#'
#' Estimates the information distance for the chosen model and recovers the
#' tree with \code{\link{sls}}. Semi-parametric, Gaussian and discrete trees
#' may have observed internal nodes; in the parametric mixed model every
#' observed node is a leaf.
#'
#' @param x Numeric data frame or matrix with at least three columns and two
#'   rows; it must contain no NA, infinite, or constant columns.
#' @param tau1,tau2 distance thresholds for the SLS algorithm.
#' @param allow_observed_internal Whether observed variables may be internal
#'   nodes. \code{NULL} (default) uses \code{TRUE} for \code{"stm"},
#'   \code{"gaussian"} and \code{"discrete"}, and \code{FALSE} for
#'   \code{"mixed"}. A value that conflicts with the model gives a warning.
#' @param epsilon_init,step_size,relation_tolerance SLS tolerance-search parameters.
#' @param use_near_pd Whether to project the correlation matrix with
#'   \code{Matrix::nearPD}. \code{NULL} (default) uses \code{TRUE} for
#'   \code{"stm"} and \code{FALSE} for the other models. Requires \pkg{Matrix}.
#' @param maxit maximum number of iterations for \code{Matrix::nearPD}.
#' @param return_intermediate logical; if \code{TRUE}, also return the
#'   estimated correlation and distance matrices, variable types, settings and
#'   timings.
#' @param model One of
#'   \itemize{
#'     \item \code{"stm"} (default): semi-parametric tree model for continuous,
#'       binary and ordinal data, with latent correlations estimated from
#'       Kendall's tau (Zhou et al., 2026).
#'     \item \code{"gaussian"}: Gaussian latent tree with Pearson correlations;
#'       all columns continuous.
#'     \item \code{"discrete"}: discrete latent tree with the determinant-based
#'       distance; all columns coded 0, 1, ... with the same number of states.
#'     \item \code{"mixed"}: parametric mixed latent tree of Zhou et al. (2020)
#'       for binary and conditionally Gaussian variables; ordinal variables are
#'       not supported.
#'     \item \code{"semiparametric"}: same as \code{"stm"}.
#'   }
#' @return If \code{return_intermediate = FALSE}, an \code{ltree}. Otherwise a
#'   list with \code{tree}, \code{correlation}, \code{distance},
#'   \code{variable_types},
#'   \code{model}, \code{allow_observed_internal}, \code{use_near_pd}, and
#'   \code{run_time}.
#' @seealso \code{\link{information_distance}}, \code{\link{sls}},
#'   \code{\link{latent_correlation}}, \code{\link{robinson_foulds_distance}}
#' @examples
#' set.seed(1234)
#' n <- 400
#' h <- rnorm(n)
#' gaussian_x <- as.data.frame(replicate(4, 0.8 * h + 0.6 * rnorm(n)))
#' names(gaussian_x) <- paste0("V", 1:4)
#' stm_x <- data.frame(
#'   V1 = exp(gaussian_x$V1),
#'   V2 = gaussian_x$V2^3,
#'   V3 = as.integer(gaussian_x$V3 > 0),
#'   V4 = findInterval(gaussian_x$V4, c(-0.5, 0.5))
#' )
#' est_stm <- learn_sltree(stm_x, model = "stm", use_near_pd = FALSE)
#'
#' est_gaussian <- learn_sltree(gaussian_x, model = "gaussian")
#'
#' h_binary <- rbinom(n, 1, 0.5)
#' discrete_x <- as.data.frame(replicate(4, rbinom(n, 1, 0.1 + 0.8 * h_binary)))
#' names(discrete_x) <- paste0("V", 1:4)
#' est_discrete <- learn_sltree(discrete_x, model = "discrete")
#'
#' binary_and_gaussian_x <- data.frame(
#'   V1 = rnorm(n, mean = 2 * h_binary - 1, sd = 0.6),
#'   V2 = rnorm(n, mean = 2 * h_binary - 1, sd = 0.6),
#'   V3 = discrete_x$V3,
#'   V4 = discrete_x$V4
#' )
#' est_mixed <- learn_sltree(binary_and_gaussian_x, model = "mixed")
#' edges_from_ltree(est_stm)
#' @export
learn_sltree <- function(x, tau1 = 3, tau2 = 5, allow_observed_internal = NULL,
                         epsilon_init = 0.1, step_size = 0.1, relation_tolerance = 0.1,
                         use_near_pd = NULL, maxit = 10000,
                         return_intermediate = FALSE,
                         model = c("stm", "gaussian", "discrete", "mixed",
                                   "semiparametric")) {
  model <- match.arg(model)
  if (model == "semiparametric") {
    model <- "stm"
  }
  use_near_pd <- .resolve_near_pd_setting(model, use_near_pd)
  if (!is.logical(return_intermediate) || length(return_intermediate) != 1L ||
      is.na(return_intermediate)) {
    stop("`return_intermediate` must be one non-missing logical value.")
  }

  model_allows_observed_internal <- model != "mixed"
  if (is.null(allow_observed_internal)) {
    allow_observed_internal <- model_allows_observed_internal
  } else {
    if (!is.logical(allow_observed_internal) || length(allow_observed_internal) != 1L ||
        is.na(allow_observed_internal)) {
      stop("`allow_observed_internal` must be NULL or one non-missing logical value.")
    }
    if (!identical(allow_observed_internal, model_allows_observed_internal)) {
      warning("`allow_observed_internal = ", allow_observed_internal,
              "` overrides the structural restriction for model = '", model,
              "'; the result no longer follows that model's default tree class.",
              call. = FALSE)
    }
  }

  t0 <- proc.time()
  .validate_numeric_data(
    x, min_rows = 2L, min_cols = 3L, reject_constant = TRUE
  )
  variable_types <- .infer_variable_types_validated(x)

  if (model == "gaussian" && any(variable_types != "continuous")) {
    bad <- which(variable_types != "continuous")
    stop("model='gaussian' requires every column to be continuous; ",
         "non-continuous columns: ", paste(bad, collapse = ", "), ".")
  }
  if (model == "discrete" && any(variable_types == "continuous")) {
    bad <- which(variable_types == "continuous")
    stop("model='discrete' requires every column to be binary or ordinal-coded; ",
         "continuous/nonconforming columns: ", paste(bad, collapse = ", "), ".")
  }
  if (model == "mixed" && any(variable_types == "ordinal")) {
    bad <- which(variable_types == "ordinal")
    stop("model='mixed' is the parametric binary/conditional-Gaussian model and ",
         "does not support ordinal columns: ", paste(bad, collapse = ", "),
         ". Use model='stm' for ordinal data.")
  }

  distance_fit <- .estimate_information_distance(
    data = x,
    type = model,
    use_near_pd = use_near_pd,
    maxit = maxit,
    variable_types = if (model == "stm") variable_types else NULL
  )
  correlation <- distance_fit$correlation
  distance <- distance_fit$distance
  t_distance <- unname((proc.time() - t0)[3])

  t1 <- proc.time()
  tree <- sls(distance, allow_observed_internal = allow_observed_internal,
               epsilon_init = epsilon_init, epsilon_max = max(distance),
               tau1 = tau1, tau2 = tau2, step_size = step_size,
               relation_tolerance = relation_tolerance, node_types = variable_types,
               node_names = colnames(x))
  t_sls <- unname((proc.time() - t1)[3])

  if (!return_intermediate) return(tree)

  run_time <- stats::setNames(
    c(t_distance, t_sls, t_distance + t_sls),
    c("distance", "SLS", "all")
  )
  list(tree = tree, correlation = correlation, distance = distance,
       variable_types = variable_types,
       model = model, allow_observed_internal = allow_observed_internal,
       use_near_pd = use_near_pd,
       run_time = run_time)
}


# Default use_near_pd for each model.
.resolve_near_pd_setting <- function(model, use_near_pd = NULL) {
  if (is.null(use_near_pd)) {
    return(identical(model, "stm"))
  }
  if (!is.logical(use_near_pd) || length(use_near_pd) != 1L ||
      is.na(use_near_pd)) {
    stop("`use_near_pd` must be NULL or one non-missing logical value.")
  }
  use_near_pd
}
