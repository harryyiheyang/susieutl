.utl_pattern_basis_projection <- function(B, validation_tol) {
  S <- svd(B)
  keep <- S$d > validation_tol * max(1, S$d)
  if (!any(keep)) return(matrix(0, nrow(B), nrow(B)))
  Q <- S$u[, keep, drop = FALSE]
  tcrossprod(Q)
}

#' Construct the canonical reference/departure theta pattern library.
#'
#' Builds raw, unnormalized basis matrices for consistent, unrestricted,
#' non-reference-specific, and ancestry-null scientific patterns. Exact duplicate
#' covariance spans are represented as aliases rather than repeated mixture
#' components. The exact zero-variance SER state is metadata and never enters
#' the covariance mixture or its weights.
#'
#' @param delta Strictly named finite transfer coefficients for non-reference
#'   traits, in the desired trait order.
#' @param reference One non-empty reference trait name not present in
#'   `names(delta)`.
#' @param include_AN_ref Whether to request the reference-null pattern.
#' @param weights Optional named non-negative weights aligned exactly to the
#'   final unique non-null components and summing to one. `NULL` gives uniform
#'   weights over those components.
#' @param validation_tol One finite positive validation and duplicate-detection
#'   tolerance.
#' @return An object inheriting from [utl_theta_prior()] with raw unique bases,
#'   the requested-pattern taxonomy, alias map, diagnostics, and implicit
#'   zero-SER metadata.
#' @export
utl_theta_pattern_prior <- function(delta, reference, include_AN_ref = FALSE,
                                    weights = NULL,
                                    validation_tol = 1e-10) {
  if (!is.numeric(delta) || !is.null(dim(delta)) || length(delta) < 1L ||
      any(!is.finite(delta)) || is.null(names(delta)) || anyNA(names(delta)) ||
      any(names(delta) == "") || anyDuplicated(names(delta))) {
    stop("delta must be a strictly named finite numeric vector")
  }
  if (!is.character(reference) || length(reference) != 1L ||
      is.na(reference) || reference == "") {
    stop("reference must be one non-empty character string")
  }
  if (reference %in% names(delta)) {
    stop("reference must not occur in names(delta)")
  }
  if (!is.logical(include_AN_ref) || length(include_AN_ref) != 1L ||
      is.na(include_AN_ref)) {
    stop("include_AN_ref must be TRUE or FALSE")
  }
  if (!is.numeric(validation_tol) || length(validation_tol) != 1L ||
      !is.finite(validation_tol) || validation_tol <= 0) {
    stop("validation_tol must be one finite positive number")
  }
  T <- length(delta) + 1L
  if (T < 2L || T > 5L) stop("automatic theta patterns support T from 2 to 5")
  if (all(delta == 0)) stop("automatic theta patterns require a nonzero delta")

  ancestry_order <- c(reference, names(delta))
  theta_order <- paste0("theta_", ancestry_order)
  M <- diag(T)
  M[-1L, 1L] <- delta
  dimnames(M) <- list(ancestry_order, theta_order)
  e <- diag(T)
  rownames(e) <- colnames(e) <- theta_order
  requested_B <- list()
  requested_family <- character()
  add_pattern <- function(pattern, family, B) {
    B <- as.matrix(B)
    rownames(B) <- theta_order
    if (is.null(colnames(B))) colnames(B) <- paste0(pattern, "_", seq_len(ncol(B)))
    requested_B[[pattern]] <<- B
    requested_family[[pattern]] <<- family
  }
  add_pattern("C", "C", matrix(e[, 1L], T, 1L,
                                  dimnames = list(theta_order, "reference")))
  add_pattern("D", "D", e)
  for (a in names(delta)) {
    j <- match(a, ancestry_order)
    add_pattern(paste0("AS_", a), "AS",
                matrix(e[, j], T, 1L,
                       dimnames = list(theta_order, paste0("AS_", a))))
  }
  for (a in names(delta)) {
    j <- match(a, ancestry_order)
    ca <- e[, 1L] - delta[[a]] * e[, j]
    other <- setdiff(seq.int(2L, T), j)
    B <- cbind(ca, e[, other, drop = FALSE])
    colnames(B) <- c(paste0("reference_minus_", a), theta_order[other])
    add_pattern(paste0("AN_", a), "AN", B)
  }
  if (include_AN_ref) {
    B <- e[, -1L, drop = FALSE]
    add_pattern(paste0("AN_", reference), "AN", B)
  }

  requested_U <- lapply(requested_B, tcrossprod)
  requested_rank <- vapply(requested_B, function(B) {
    sum(svd(B, nu = 0, nv = 0)$d > validation_tol *
          max(1, svd(B, nu = 0, nv = 0)$d))
  }, integer(1))
  canonical <- character()
  aliases <- setNames(character(), character())
  B <- list()
  U <- list()
  for (pattern in names(requested_B)) {
    alias <- NULL
    P <- .utl_pattern_basis_projection(requested_B[[pattern]], validation_tol)
    for (keep in canonical) {
      U_scale <- max(1, max(abs(requested_U[[pattern]])), max(abs(U[[keep]])))
      same_U <- max(abs(requested_U[[pattern]] - U[[keep]])) <=
        validation_tol * U_scale
      P_keep <- .utl_pattern_basis_projection(B[[keep]], validation_tol)
      same_span <- requested_rank[[pattern]] == requested_rank[[keep]] &&
        max(abs(P - P_keep)) <= validation_tol * max(1, max(abs(P_keep)))
      if (same_U || same_span) {
        alias <- keep
        break
      }
    }
    if (is.null(alias)) {
      canonical <- c(canonical, pattern)
      B[[pattern]] <- requested_B[[pattern]]
      U[[pattern]] <- requested_U[[pattern]]
    } else {
      aliases[[pattern]] <- alias
    }
  }
  for (pattern in names(U)) dimnames(U[[pattern]]) <- list(theta_order, theta_order)

  default_weights <- is.null(weights)
  if (default_weights) {
    weights <- setNames(rep(1 / length(U), length(U)), names(U))
  } else if (!is.numeric(weights) || !is.null(dim(weights)) ||
             length(weights) != length(U) || any(!is.finite(weights)) ||
             any(weights < 0) || is.null(names(weights)) ||
             !identical(names(weights), names(U))) {
    stop("weights must be named finite non-negative values aligned to unique patterns")
  }
  prior <- utl_theta_prior(delta, M, U, weights, validation_tol)
  prior$F <- B
  prior$rank <- setNames(unname(requested_rank[names(B)]), names(B))

  patterns <- names(requested_B)
  alias_of <- setNames(rep("", length(patterns)), patterns)
  alias_of[names(aliases)] <- unname(aliases)
  included <- patterns %in% names(U)
  initial_weight <- setNames(numeric(length(patterns)), patterns)
  initial_weight[names(U)] <- prior$weights
  pattern_table <- data.frame(
    pattern = c(patterns, "SER_V0"),
    family = c(unname(requested_family[patterns]), "SER_V0"),
    reference = reference,
    nonnull = c(rep(TRUE, length(patterns)), FALSE),
    included_in_U = c(included, FALSE),
    alias_of = c(unname(alias_of), ""),
    rank = c(unname(requested_rank[patterns]), 0L),
    basis_columns = c(vapply(requested_B[patterns], function(x) {
      paste(colnames(x), collapse = ",")
    }, character(1)), ""),
    initial_weight = c(unname(initial_weight), 0),
    stringsAsFactors = FALSE
  )
  diagnostics <- list(
    T = T, reference = reference, ancestry_order = ancestry_order,
    theta_order = theta_order, delta = delta, include_AN_ref = include_AN_ref,
    requested_nonnull = length(patterns), unique_nonnull = length(U),
    taxonomy_count = length(U) + 1L, alias_count = length(aliases),
    default_weight_initialization = if (default_weights) {
      "uniform_over_unique_nonnull"
    } else {
      "user_supplied"
    }
  )
  null_metadata <- list(
    pattern = "SER_V0", implicit = TRUE, entered_U = FALSE,
    entered_weights = FALSE, semantics = "exact SER prior variance zero"
  )
  prior$B <- B
  prior$U <- prior$Utheta
  prior$pattern_table <- pattern_table
  prior$pattern_aliases <- aliases
  prior$pattern_diagnostics <- diagnostics
  prior$null_metadata <- null_metadata
  class(prior) <- c("utl_theta_pattern_prior", "utl_theta_prior")
  prior
}
