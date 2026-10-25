#' sltree: Structure Learning of Semi-Parametric Latent Trees
#'
#' Learns latent tree structures from mixed-type data (continuous, binary,
#' ordinal) via the SLS algorithm. See \code{\link{sls}}, \code{\link{latent_correlation}},
#' \code{\link{learn_sltree}}.
#'
#' @references
#' Zhou, C., Li, N., Wang, S., Wang, X. and Guo, J. (2026).
#' Learning semi-parametric tree models from mixed data.
#' \emph{Artificial Intelligence}, \strong{353}, 104499.
#' \doi{10.1016/j.artint.2026.104499}
#'
#' Zhou, C., Wang, X. and Guo, J. (2020).
#' Learning Mixed Latent Tree Models.
#' \emph{Journal of Machine Learning Research}, \strong{21}(249), 1--35.
#' \url{https://jmlr.org/papers/v21/20-365.html}
#'
#' @keywords internal
#' @useDynLib sltree, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
