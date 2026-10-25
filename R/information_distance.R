#' Computing information distance.
#'
#' \code{information_distance} returns a numeric information-distance matrix.
#' The input of \code{information_distance} is a \eqn{n} by \eqn{p} numeric
#'   matrix or data frame, which is a sample of latent tree models (LTM).
#'   Where \eqn{n} is the sample size and \eqn{p} is the number of observed
#'   variables. The column names of the input will be assigned to column and
#'   row names of the result (information distance matrix).
#'
#' The information distance is additive along path on LTM.
#'   Let \eqn{x,y} be two nodes (random variables) in LTM. For Gaussian and
#'   parametric mixed LTM (\code{type = "gaussian"} or \code{type = "mixed"}),
#'   the information distance \eqn{d_{xy}} between \eqn{x} and \eqn{y}
#'   is \deqn{d_{xy} = -log(|\rho_{xy}|),} where \eqn{\rho_{xy}} is the
#'   Pearson correlation coefficient of \eqn{x} and \eqn{y}.
#'   For discrete LTM (\code{type = "discrete"}), the information distance
#'   \eqn{d_{xy}} is \deqn{d_{xy} = -log(\frac{|\det(P_{xy})|}
#'   {\sqrt{\det(P_{xx})\det(P_{yy})}}),}
#'   where \eqn{P_{xy} = (P(x=a,y=b))} is the joint probability matrix and
#'   \eqn{P_{xx}} and \eqn{P_{yy}} are diagonal marginal probability matrices.
#'   All observed variables must have the same number of observed states
#'   (at least two), so that every joint probability matrix is square.
#'   State counts are determined from distinct values present in each column.
#'   Singular joint matrices yield a distance capped at
#'   \code{-log(.Machine$double.xmin)}.
#'   For semi-parametric tree models (\code{type = "semiparametric"},
#'   the default; alias \code{"stm"}), the information distance \eqn{d_{xy}} is
#'   \deqn{d_{xy} = -log(|\rho_{xy}|),} where \eqn{\rho_{xy}} is the latent Pearson
#'   correlation coefficient of \eqn{x} and \eqn{y}. The coefficient
#'   \eqn{\rho_{xy}} is estimated with a bridge function \eqn{F(\cdot)} that
#'   connects the latent correlation to the population Kendall's \eqn{\tau}
#'   (see \code{\link{latent_correlation}}).
#'
#' @param data a \eqn{n} by \eqn{p} numeric matrix or data frame, with one
#'   numeric vector per data-frame column and no missing or infinite values.
#' @param type A character scalar selecting one of four model classes (matched
#'   by \code{\link[base]{match.arg}}): \code{"semiparametric"} (default;
#'   alias \code{"stm"}) estimates latent Gaussian-copula correlation for
#'   continuous/binary/ordinal data via \code{\link{latent_correlation}};
#'   \code{"gaussian"} and \code{"mixed"} use Pearson correlation for Gaussian
#'   and parametric binary/conditional-Gaussian models, respectively;
#'   \code{"discrete"} uses the absolute-determinant distance and requires the
#'   same number of observed states in every column.
#' @param use_near_pd Logical; for correlation-based distances, whether to
#'   project the estimated correlation matrix onto the nearest positive-definite
#'   correlation matrix before computing the distance. Requires \pkg{Matrix}.
#' @param maxit Maximum number of iterations passed to
#'   \code{\link[Matrix]{nearPD}} when \code{use_near_pd = TRUE}.
#' @return A numeric \eqn{p} by \eqn{p} information-distance matrix; its row and
#'   column names are taken from the columns of \code{data}.
#' @export
#' @examples
#' x <- data.frame(continuous = c(-1.8, -1.2, -0.4, 0.2, 0.8, 1.5),
#'                 binary = c(0, 0, 0, 1, 1, 1),
#'                 ordinal = c(0, 0, 1, 1, 2, 2))
#' information_distance(x, type = "semiparametric")
#'
#' categorical <- data.frame(a = c(0, 0, 1, 1, 2, 2),
#'                           b = c(2, 2, 1, 1, 0, 0))
#' information_distance(categorical, type = "discrete")
information_distance <- function(data,
                                 type = c("semiparametric", "gaussian", "discrete", "mixed", "stm"),
                                 use_near_pd = FALSE, maxit = 10000) {
  type <- match.arg(type)
  if (type == "stm") type <- "semiparametric"
  .validate_numeric_data(
    data,
    arg = "data",
    min_rows = 2L,
    min_cols = if (type == "semiparametric") 2L else 1L,
    reject_constant = TRUE
  )
  variable_types <- if (type == "semiparametric") .infer_variable_types_validated(data) else NULL
  .estimate_information_distance(
    data = data,
    type = type,
    use_near_pd = use_near_pd,
    maxit = maxit,
    variable_types = variable_types
  )$distance
}


# Correlation and information distance.
.estimate_information_distance <- function(
    data, type = c("semiparametric", "gaussian", "discrete", "mixed", "stm"),
    use_near_pd = FALSE, maxit = 10000, variable_types = NULL) {
  type <- match.arg(type)
  if (type == "stm") type <- "semiparametric"

  if (type == "discrete") {
    return(list(
      correlation = NULL,
      distance = .discrete_information_distance(data)
    ))
  }

  if (type == "semiparametric") {
    if (is.null(variable_types)) variable_types <- infer_variable_types(data)
    correlation <- latent_correlation(data, variable_types)
  } else {
    # Positive column scaling preserves Pearson correlation while avoiding
    # overflow and underflow in its squared deviations.
    x <- as.matrix(data)
    scales <- apply(abs(x), 2L, max)
    correlation <- stats::cor(sweep(x, 2L, scales, "/"))
  }

  if (any(!is.finite(correlation)) ||
      any(abs(correlation) > 1 + sqrt(.Machine$double.eps))) {
    stop("Could not estimate finite correlations between -1 and 1 from `data`.")
  }
  correlation <- .project_correlation(
    correlation, use_near_pd = use_near_pd, maxit = maxit
  )
  correlation <- pmin(pmax(correlation, -1), 1)
  distance <- -log(abs(correlation))
  # Cap infinite distances (rho == 0) so that sls() can use the matrix.
  distance[is.infinite(distance)] <- -log(.Machine$double.xmin)
  diag(distance) <- 0

  list(correlation = correlation, distance = distance)
}


.project_correlation <- function(correlation, use_near_pd = FALSE, maxit = 10000) {
  if (!is.logical(use_near_pd) || length(use_near_pd) != 1L || is.na(use_near_pd)) {
    stop("`use_near_pd` must be one non-missing logical value.")
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) ||
      !is.finite(maxit) || maxit < 1 || maxit > .Machine$integer.max ||
      maxit != floor(maxit)) {
    stop("`maxit` must be one positive integer.")
  }
  if (!use_near_pd) return(correlation)

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    warning("use_near_pd=TRUE requires the 'Matrix' package (in Suggests); skipping nearPD.")
    return(correlation)
  }

  as.matrix(Matrix::nearPD(correlation, corr = TRUE, maxit = maxit)$mat)
}


.discrete_information_distance <- function(data) {
  if (!is.matrix(data)) data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  dist_matrix <- matrix(0, nrow = p, ncol = p)
  colnames(dist_matrix) <- rownames(dist_matrix) <- colnames(data)

  state_counts <- vapply(seq_len(p), function(i) {
    length(unique(data[, i]))
  }, integer(1))
  if (any(state_counts < 2L)) {
    stop("Every discrete variable must have at least two observed states.")
  }
  if (length(unique(state_counts)) != 1L) {
    stop("All discrete variables must have the same number of observed states; ",
         "observed state counts: ", paste(state_counts, collapse = ", "), ".")
  }

  # Log scale avoids underflow.
  marginal_log_det <- vapply(seq_len(p), function(i) {
    sum(log(as.numeric(table(data[, i])) / n))
  }, numeric(1))

  for (i in seq_len(p - 1L)) {
    for (j in seq.int(i + 1L, p)) {
      pxy <- unclass(table(data[, i], data[, j])) / n
      log_det <- as.numeric(determinant(pxy, logarithm = TRUE)$modulus)
      log_ratio <- log_det - (marginal_log_det[i] + marginal_log_det[j]) / 2
      if (is.na(log_ratio) || log_ratio == Inf) {
        stop("Could not compute a finite discrete information distance.")
      }
      if (log_ratio > log1p(sqrt(.Machine$double.eps))) {
        stop("The normalized absolute determinant exceeds 1.")
      }
      d <- -min(max(log_ratio, log(.Machine$double.xmin)), 0)
      dist_matrix[i, j] <- dist_matrix[j, i] <- d
    }
  }

  dist_matrix
}


.validate_numeric_data <- function(x, arg = "x", min_rows = 1L,
                                   min_cols = 1L, reject_constant = FALSE) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop(sprintf("`%s` must be a matrix or data frame.", arg))
  }
  if (nrow(x) < min_rows) {
    stop(sprintf("`%s` must have at least %d rows.", arg, min_rows))
  }
  if (ncol(x) < min_cols) {
    stop(sprintf("`%s` must have at least %d columns.", arg, min_cols))
  }
  if (is.data.frame(x) && any(vapply(x, function(column) {
    !is.null(dim(column))
  }, logical(1)))) {
    stop(sprintf("`%s` must have one numeric vector per column (no matrix-valued columns).", arg))
  }

  x_matrix <- as.matrix(x)
  if (!is.numeric(x_matrix)) {
    stop(sprintf("`%s` must contain only numeric values.", arg))
  }
  if (anyNA(x_matrix) || !all(is.finite(x_matrix))) {
    stop(sprintf("`%s` must contain only finite values (no NA or Inf).", arg))
  }

  if (reject_constant) {
    constant <- vapply(seq_len(ncol(x_matrix)), function(i) {
      length(unique(x_matrix[, i])) < 2L
    }, logical(1))
    if (any(constant)) {
      stop(sprintf(
        "`%s` contains constant columns: %s.", arg,
        paste(which(constant), collapse = ", ")
      ))
    }
  }

  invisible(x_matrix)
}


.validate_semiparametric_inputs <- function(x, variable_types, min_cols = 1L) {
  x_matrix <- .validate_numeric_data(
    x, min_rows = 2L, min_cols = min_cols, reject_constant = TRUE
  )

  if (!is.character(variable_types) ||
      length(variable_types) != ncol(x_matrix) || anyNA(variable_types)) {
    stop("`variable_types` must be a non-missing character vector with one value per column.")
  }
  allowed <- c("continuous", "binary", "ordinal")
  if (!all(variable_types %in% allowed)) {
    stop("All variable types must be 'continuous', 'binary', or 'ordinal'.")
  }

  for (i in seq_len(ncol(x_matrix))) {
    values <- sort(unique(x_matrix[, i]))
    if (variable_types[i] == "binary" &&
        !identical(as.numeric(values), c(0, 1))) {
      stop(sprintf("Binary variable in column %d must contain exactly 0 and 1.", i))
    }
    if (variable_types[i] == "ordinal") {
      valid <- length(values) >= 3L && all(values >= 0) &&
        all(values == floor(values)) && values[1] == 0 &&
        identical(as.numeric(values), as.numeric(seq.int(0, max(values))))
      if (!valid) {
        stop(sprintf(
          "Ordinal variable in column %d must have at least three levels coded as consecutive integers starting at 0 (0, 1, ..., k) with no missing levels.",
          i
        ))
      }
    }
  }

  invisible(x_matrix)
}


#' Infer Variable Types in a Dataset
#'
#' Classifies each column as binary, continuous or ordinal from its values.
#'
#' @param x A data frame or matrix where each column represents a variable.
#'   Must not contain NA values or non-numeric data. Binary variables should
#'   only contain 0 or 1. Ordinal variables should only contain consecutive
#'   integers starting from 0 (i.e., 0, 1, 2, ...).
#'
#' @return A character vector of types for all variables with possible values:
#'   "binary", "continuous", or "ordinal".
#'
#' @details
#' \itemize{
#'   \item Binary: exactly 2 unique values that are 0 and 1
#'   \item Ordinal: 3+ unique values that are consecutive integers starting from 0
#'   \item Continuous: all other numeric variables
#' }
#'
#' @examples
#' data <- data.frame(
#'   binary_var = c(0, 1, 0, 1),
#'   ordinal_var = c(0, 1, 2, 1),
#'   continuous_var = c(1.5, 2.3, 4.1, 3.7)
#' )
#' infer_variable_types(data)
#'
#' @export
infer_variable_types <- function(x) {
  .validate_numeric_data(x, min_rows = 1L, min_cols = 1L)
  .infer_variable_types_validated(x)
}


.infer_variable_types_validated <- function(x) {

  variable_types <- rep("", ncol(x))

  for (i in seq_len(ncol(x))) {
    unique_values <- unique(x[, i])
    n_unique <- length(unique_values)

    if (n_unique == 0) {
      stop("Data contains a variable with no values")
    }

    if (n_unique == 1) {
      stop("Data contains a variable with a constant value")
    }

    if (n_unique == 2) {
      variable_types[i] <- if (all(unique_values %in% c(0, 1))) {
        "binary"
      } else {
        "continuous"
      }
    } else {
      sorted_levels <- sort(unique_values)
      max_level <- max(sorted_levels)
      is_nonnegative_integer <-
        all(is.finite(sorted_levels)) &&
        all(sorted_levels >= 0) &&
        all(sorted_levels == floor(sorted_levels))
      is_ordinal <-
        is_nonnegative_integer &&
        min(sorted_levels) == 0 &&
        length(sorted_levels) == max_level + 1 &&
        all(sorted_levels == seq.int(0, max_level))

      if (is_nonnegative_integer && !is_ordinal) {
        warning(sprintf(
          "Column %d looks ordinal but its levels are not consecutive integers starting at 0; treating it as continuous. Recode ordinal variables as 0, 1, ..., k with no missing levels.",
          i))
      }

      variable_types[i] <- if (is_ordinal) "ordinal" else "continuous"
    }
  }

  variable_types
}


#' Latent correlation matrix for mixed data
#'
#' Estimates the latent Pearson correlations of continuous, binary and ordinal
#' variables under a Gaussian copula model.
#'
#' @param x A data frame or matrix where each column represents a variable.
#'   Must not contain NA values. Binary variables should contain only 0 and 1.
#'   Ordinal variables should contain consecutive integers starting from 0.
#' @param variable_types A character vector specifying the type of each variable in x.
#'   Must contain values "continuous", "binary", or "ordinal".
#' @return A symmetric latent Pearson-correlation matrix with one row and column
#'   per variable in \code{x}.
#'
#' @details
#' Each ordinal variable is split into binary indicators, one per cut point.
#' Kendall's tau-a is computed for all pairs and mapped to the latent
#' correlation by inverting the bridge functions of Fan et al. (2017). For
#' ordinal variables the estimates are averaged over cut points.
#'
#' @references
#' Fan, J., Liu, H., Ning, Y. and Zou, H. (2017).
#' High dimensional semiparametric latent graphical model for mixed data.
#' \emph{Journal of the Royal Statistical Society: Series B}, \strong{79}(2), 405--421.
#' \doi{10.1111/rssb.12168}
#'
#' Huang, M., M\"uller, C. L. and Gaynanova, I. (2021).
#' latentcor: An R package for estimating latent correlations from mixed data types.
#' \emph{Journal of Open Source Software}, \strong{6}(65), 3634.
#' \doi{10.21105/joss.03634}
#'
#' @examples
#' set.seed(1)
#' data <- data.frame(
#'   continuous_var = rnorm(100),
#'   binary_var = rep(0:1, length.out = 100),
#'   ordinal_var = rep(0:3, length.out = 100)
#' )
#' variable_types <- c("continuous", "binary", "ordinal")
#' latent_correlation(data, variable_types)
#'
#' @export
latent_correlation <- function(x, variable_types) {
  .validate_semiparametric_inputs(x, variable_types, min_cols = 2L)

  prepared_data <- .prepare_copula_data_validated(x, variable_types)
  tau <- .kendall_tau_a_validated(prepared_data$data)
  correlation <- estimate_latent_correlation_cpp(
    tau,
    prepared_data$source_index,
    variable_types,
    prepared_data$thresholds
  )
  if (!is.null(colnames(x))) {
    dimnames(correlation) <- list(colnames(x), colnames(x))
  }
  correlation
}


#' Expand ordinal variables and compute thresholds
#'
#' An ordinal variable with levels 0, ..., k becomes k binary indicators
#' (>= 1, ..., >= k).
#'
#' @param x Data frame or matrix with binary (0/1), ordinal (0, 1, ...) and
#'   continuous columns.
#' @param variable_types Type of each column: "continuous", "binary" or
#'   "ordinal".
#'
#' @return A list with \code{data} (expanded data), \code{source_index}
#'   (column of \code{x} for each column of \code{data}) and \code{thresholds}
#'   (mean for continuous variables, \code{qnorm(1 - mean)} for each binary
#'   indicator).
#'
#' @examples
#' set.seed(1)
#' data <- data.frame(
#'   cont_var = rnorm(100),
#'   bin_var = rep(0:1, length.out = 100),
#'   ord_var = rep(0:3, length.out = 100)
#' )
#' variable_types <- c("continuous", "binary", "ordinal")
#' result <- prepare_copula_data(data, variable_types)
#'
#' @importFrom stats qnorm
#' @noRd
prepare_copula_data <- function(x, variable_types) {
  .validate_semiparametric_inputs(x, variable_types, min_cols = 1L)
  .prepare_copula_data_validated(x, variable_types)
}


.prepare_copula_data_validated <- function(x, variable_types) {
  n_vars <- ncol(x)

  widths <- ifelse(variable_types == "ordinal",
                   vapply(seq_len(n_vars), function(i) max(x[, i]), numeric(1)),
                   1)
  total_cols <- sum(widths)
  y_list <- vector("list", total_cols)
  source_index <- integer(total_cols)
  thresholds <- vector("list", n_vars)
  out_col <- 0L

  for (i in seq_len(n_vars)) {
    variable_type <- variable_types[i]
    current_var <- x[, i]

    if (variable_type == "continuous") {
      out_col <- out_col + 1L
      source_index[out_col] <- i
      y_list[[out_col]] <- current_var
      thresholds[[i]] <- mean(current_var)

    } else if (variable_type == "binary") {
      out_col <- out_col + 1L
      source_index[out_col] <- i
      y_list[[out_col]] <- current_var
      thresholds[[i]] <- qnorm(1 - mean(current_var))

    } else if (variable_type == "ordinal") {
      max_level <- widths[i]
      threshold_deltas <- numeric(max_level)

      for (j in seq_len(max_level)) {
        binary_cut <- as.numeric(current_var >= j)
        out_col <- out_col + 1L
        y_list[[out_col]] <- binary_cut
        source_index[out_col] <- i
        threshold_deltas[j] <- qnorm(1 - mean(binary_cut))
      }

      thresholds[[i]] <- threshold_deltas
    }
  }

  y_final <- do.call(cbind, y_list)

  list(
    data = y_final,
    source_index = source_index,
    thresholds = thresholds
  )
}


#' Kendall's tau-a matrix
#'
#' Tau-b from \code{pcaPP::cor.fk} is converted to tau-a by
#' tau-a = tau-b * sqrt((n0 - n_x) * (n0 - n_y)) / n0, where n0 = n(n-1)/2
#' and n_x, n_y are the numbers of tied pairs.
#'
#' @param x Numeric data frame or matrix with at least two rows and columns.
#'
#' @return A symmetric matrix of Kendall's tau-a.
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(300), ncol = 3)
#' kendall_tau_a(data)
#'
#' mixed_data <- cbind(
#'   rnorm(100),
#'   sample(1:5, 100, replace = TRUE),
#'   rbinom(100, 1, 0.5)
#' )
#' kendall_tau_a(mixed_data)
#'
#' @importFrom pcaPP cor.fk
#' @noRd
kendall_tau_a <- function(x) {
  .validate_numeric_data(
    x, min_rows = 2L, min_cols = 2L, reject_constant = TRUE
  )
  .kendall_tau_a_validated(x)
}


# Adapted from Kendalltau() in latentcor. See inst/COPYRIGHTS.
.kendall_tau_a_validated <- function(x) {

  compute_tie_correction <- function(x, n) {
    values <- unique(x)
    n_unique <- length(values)

    if (n_unique == n) {
      return(0)
    }
    if (n_unique == 2L && all(values %in% c(0, 1))) {
      n0 <- sum(x == 0)
      n1 <- n - n0
      return(n0 * (n0 - 1) / 2 + n1 * (n1 - 1) / 2)
    }

    run_info <- rle(sort(x))
    tied_lengths <- run_info$lengths[run_info$lengths > 1]

    tie_correction <- sum(tied_lengths * (tied_lengths - 1) / 2)

    return(tie_correction)
  }

  n <- nrow(x)
  n0 <- n * (n - 1) / 2

  tie_corrections <- apply(x, 2, function(x) {
    compute_tie_correction(x, n)
  })

  correction_factors <- sqrt(n0 - tie_corrections)

  kendall_b_matrix <- pcaPP::cor.fk(x)

  kendall_b_lower <- kendall_b_matrix[lower.tri(kendall_b_matrix)]

  lower_tri_indices <- which(lower.tri(kendall_b_matrix), arr.ind = TRUE)
  row_indices <- lower_tri_indices[, 1]
  col_indices <- lower_tri_indices[, 2]

  # tau-a = tau-b * sqrt((n0-n_x)(n0-n_y)) / n0
  correction_products <- (correction_factors[row_indices] *
                            correction_factors[col_indices]) / n0

  kendall_a_lower <- kendall_b_lower * correction_products

  n_vars <- ncol(x)
  result_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)
  result_matrix[lower.tri(result_matrix)] <- kendall_a_lower
  result_matrix <- result_matrix + t(result_matrix)
  diag(result_matrix) <- 1

  if (!is.null(colnames(x))) {
    rownames(result_matrix) <- colnames(x)
    colnames(result_matrix) <- colnames(x)
  }

  result_matrix
}
