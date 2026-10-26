#' sltree: Structure Learning of Semiparametric/Latent Tree Models
#'
#' Learns the structure of Gaussian latent tree models (LTMs),
#' discrete LTMs, parametric mixed LTMs, and semiparametric tree models
#' using information distances and a bottom-up tree reconstruction
#' algorithm. Gaussian models and parametric mixed models with binary
#' and conditionally Gaussian variables use Pearson-correlation-based
#' distances. Discrete models use determinant-based distances and
#' require the same number of observed states across variables.
#' Semiparametric models accommodate continuous, binary, and ordinal
#' variables, with latent correlations estimated under a Gaussian
#' copula model using Kendall's tau and bridge functions. Utilities
#' are provided for data simulation, tree visualization and export,
#' and topology comparison using the Robinson-Foulds distance.
#' The semiparametric framework allowing observed internal nodes is
#' described in Zhou et al. (2026). Earlier structure-learning methods
#' assuming that all observed variables are leaf nodes are described
#' for discrete models in Wang et al. (2017) and for parametric mixed
#' models in Zhou et al. (2020).
#'
#' See \code{\link{sls}}, \code{\link{latent_correlation}},
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
#' Wang, X., Guo, J., Hao, L. and Zhang, N. L. (2017).
#' Spectral methods for learning discrete latent tree models.
#' \emph{Statistics and Its Interface}, \strong{10}(4), 677--698.
#' \doi{10.4310/SII.2017.v10.n4.a11}
#'
#' @keywords internal
#' @useDynLib sltree, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
