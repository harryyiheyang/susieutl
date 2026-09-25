#' Construct a fixed theta-coordinate covariance prior.
#'
#' The theta-to-beta map must use the first ancestry as the reference:
#' `beta = M theta`, with an identity reference row and the named transfer
#' coefficients in the first theta column of the remaining rows. Covariance
#' components are retained at their supplied scale.
#' For automatically generated reference/departure pattern libraries, see
#' [utl_theta_pattern_prior()]. Manual `M`/`Utheta` construction remains fully
#' supported.
#'
#' @param delta Named finite transfer coefficients for non-reference ancestries.
#' @param M Named theta-to-beta coordinate map.
#' @param Utheta Named list of finite, symmetric positive-semidefinite theta
#'   covariance matrices.
#' @param weights Named, fixed component weights that sum to one.
#' @param validation_tol Relative validation tolerance.
#' @param allow_zero_rank Allow a zero-rank covariance component for the exact
#'   fallback route.
#' @return An object of class `utl_theta_prior` containing the supplied inputs,
#'   a full transfer vector, the inverse coordinate map, forward-mapped
#'   observed-scale covariances, and fixed spectral factors `F` with their
#'   ranks.
#' @export
utl_theta_prior <- function(delta, M, Utheta, weights, validation_tol = 1e-10,
                            allow_zero_rank = FALSE) {
  if (!is.numeric(delta) || !is.null(dim(delta)) || length(delta) < 1L ||
      any(!is.finite(delta))) {
    stop("delta must be a finite numeric vector with at least one entry")
  }
  if (is.null(names(delta)) || anyNA(names(delta)) || any(names(delta) == "") ||
      anyDuplicated(names(delta))) {
    stop("delta must have unique, non-empty names for non-reference ancestries")
  }
  if (!is.numeric(validation_tol) || length(validation_tol) != 1L ||
      !is.finite(validation_tol) || validation_tol <= 0) {
    stop("validation_tol must be one finite positive number")
  }
  if (!is.logical(allow_zero_rank) || length(allow_zero_rank) != 1L ||
      is.na(allow_zero_rank)) stop("allow_zero_rank must be TRUE or FALSE")

  d <- length(delta) + 1L
  if (!is.matrix(M) || !is.numeric(M) || nrow(M) != d || ncol(M) != d ||
      any(!is.finite(M))) {
    stop("M must be a finite numeric square matrix with length(delta) + 1 rows and columns")
  }
  ancestry_names <- rownames(M)
  theta_names <- colnames(M)
  if (is.null(ancestry_names) || is.null(theta_names) ||
      anyNA(ancestry_names) || anyNA(theta_names) ||
      any(ancestry_names == "") || any(theta_names == "") ||
      anyDuplicated(ancestry_names) || anyDuplicated(theta_names)) {
    stop("M must have unique, non-empty ancestry row names and theta column names")
  }
  if (!identical(ancestry_names[-1L], names(delta))) {
    stop("M non-reference ancestry row names must match delta names in order")
  }
  M_canonical <- diag(d)
  M_canonical[-1L, 1L] <- delta
  M_scale <- max(1, max(abs(M_canonical)))
  if (max(abs(M - M_canonical)) > validation_tol * M_scale) {
    stop("M must equal the canonical theta-to-beta map beta = M theta")
  }
  M_det <- determinant(M, logarithm = TRUE)
  if (M_det$sign == 0 || !is.finite(M_det$modulus)) {
    stop("M must be invertible")
  }

  if (!is.list(Utheta) || !length(Utheta)) {
    stop("Utheta must be a non-empty named list of covariance matrices")
  }
  component_names <- names(Utheta)
  if (is.null(component_names) || anyNA(component_names) ||
      any(component_names == "") || anyDuplicated(component_names)) {
    stop("Utheta must have unique, non-empty component names")
  }
  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      length(weights) != length(Utheta) || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("weights must be a finite, non-negative numeric vector aligned to Utheta")
  }
  if (is.null(names(weights)) || !identical(names(weights), component_names)) {
    stop("weights names must match Utheta component names in order")
  }
  weight_sum <- sum(weights)
  if (abs(weight_sum - 1) > validation_tol * max(1, abs(weight_sum))) {
    stop("weights must sum to one within validation_tol")
  }

  factor_tol <- 1e-12
  F <- vector("list", length(Utheta))
  rank <- integer(length(Utheta))
  names(F) <- names(rank) <- component_names
  for (k in seq_along(Utheta)) {
    U <- Utheta[[k]]
    component <- component_names[[k]]
    if (!is.matrix(U) || !is.numeric(U) || nrow(U) != d || ncol(U) != d ||
        any(!is.finite(U))) {
      stop("Utheta component '", component, "' must be a finite ", d, " by ", d,
           " numeric matrix")
    }
    if (is.null(rownames(U)) || is.null(colnames(U)) ||
        !identical(rownames(U), theta_names) || !identical(colnames(U), theta_names)) {
      stop("Utheta component '", component,
           "' row and column names must match M theta column names")
    }
    U_scale <- max(1, max(abs(U)))
    if (max(abs(U - t(U))) > validation_tol * U_scale) {
      stop("Utheta component '", component, "' must be symmetric")
    }
    E <- eigen(U, symmetric = TRUE)
    psd_tol <- 100 * .Machine$double.eps * max(1, max(abs(E$values)))
    if (min(E$values) < -psd_tol) {
      stop("Utheta component '", component, "' must be positive semidefinite")
    }
    evals <- E$values
    evals[evals < 0] <- 0
    keep <- evals > 0
    rank[[k]] <- sum(keep)
    if (rank[[k]] == 0L) {
      if (!allow_zero_rank) {
        stop("Utheta component '", component, "' must have positive spectral rank")
      }
      F[[k]] <- matrix(0, d, 0L,
                       dimnames = list(theta_names, character()))
      next
    }
    F[[k]] <- E$vectors[, keep, drop = FALSE] %*%
      diag(sqrt(evals[keep]), nrow = rank[[k]])
    rownames(F[[k]]) <- theta_names
  }

  delta_full <- c(1, delta)
  names(delta_full) <- ancestry_names
  M_inv <- solve(M)
  Ubeta <- lapply(Utheta, function(U) {
    B <- M %*% U %*% t(M)
    dimnames(B) <- list(ancestry_names, ancestry_names)
    B
  })
  names(Ubeta) <- component_names
  initial_weights <- weights
  initial_Utheta <- Utheta
  structure(
    list(delta = delta, delta_full = delta_full, M = M, Utheta = Utheta,
         M_inv = M_inv, Ubeta = Ubeta, weights = weights,
         initial_weights = initial_weights, initial_Utheta = initial_Utheta,
         F = F, rank = rank,
         allow_zero_rank = allow_zero_rank,
         validation_tol = validation_tol, factor_tol = factor_tol,
         psd_tolerance = 100 * .Machine$double.eps,
         ancestry_names = ancestry_names, theta_names = theta_names,
         component_names = component_names,
         component_map = setNames(seq_along(component_names), component_names)),
    class = "utl_theta_prior"
  )
}
