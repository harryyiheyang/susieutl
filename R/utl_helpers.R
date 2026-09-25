.utl_matmul <- function(A, B, transA = FALSE, transB = FALSE) {
  if (requireNamespace("CppMatrix", quietly = TRUE)) {
    CppMatrix::matrixMultiply(A, B, transA = transA, transB = transB)
  } else {
    if (transA) A <- t(A)
    if (transB) B <- t(B)
    A %*% B
  }
}

.utl_matvec <- function(A, b, backend_state) {
  if (requireNamespace("CppMatrix", quietly = TRUE)) {
    out <- CppMatrix::matrixVectorMultiply(A, b)
    backend_state$matvec_count <- backend_state$matvec_count + 1L
    backend_state$cpp_matvec_count <- backend_state$cpp_matvec_count + 1L
    backend_state$name <- "CppMatrix"
    return(out)
  }
  out <- as.vector(A %*% b)
  backend_state$matvec_count <- backend_state$matvec_count + 1L
  out
}

.utl_logsumexp <- function(x) {
  z <- max(x)
  z + log(sum(exp(x - z)))
}

.utl_normalize_ld_eigen <- function(x, variant_names, validation_tol,
                                    context = "LD") {
  if (!is.list(x) || is.null(x$values) || is.null(x$vectors)) {
    stop(context, " must contain numeric values and vectors")
  }
  values <- as.numeric(x$values)
  vectors <- x$vectors
  if (!is.numeric(values) || length(values) < 1L ||
      any(!is.finite(values))) {
    stop(context, "$values must be a non-empty finite numeric vector")
  }
  if (!is.matrix(vectors) || !is.numeric(vectors) ||
      nrow(vectors) != length(variant_names) ||
      ncol(vectors) != length(values) || any(!is.finite(vectors))) {
    stop(context, "$vectors must be a finite numeric matrix with p rows and length(values) columns")
  }
  if (length(values) > length(variant_names)) {
    stop(context, "$values cannot have more entries than the number of variants")
  }
  vnames <- rownames(vectors)
  if (!is.null(vnames)) {
    if (anyNA(vnames) || any(vnames == "") || anyDuplicated(vnames) ||
        !setequal(vnames, variant_names)) {
      stop(context, "$vectors row names must exactly match summary SNPs")
    }
    vectors <- vectors[match(variant_names, vnames), , drop = FALSE]
  } else {
    rownames(vectors) <- variant_names
  }
  values_scale <- max(1, max(abs(values)))
  if (min(values) < -validation_tol * values_scale) {
    stop(context, " must be positive semidefinite")
  }
  values <- pmax(values, 0)
  col_norm <- colSums(vectors^2)
  if (any(!is.finite(col_norm)) ||
      any(abs(col_norm - 1) > 10 * sqrt(validation_tol) *
          max(1, max(col_norm)))) {
    stop(context, "$vectors columns must have unit Euclidean norm")
  }
  keep <- values > validation_tol * values_scale
  if (!any(keep)) stop(context, " must have positive rank")
  values <- values[keep]
  vectors <- vectors[, keep, drop = FALSE]
  diag_ld <- rowSums(sweep(vectors^2, 2L, values, "*"))
  if (any(!is.finite(diag_ld)) || any(diag_ld <= 0)) {
    stop(context, " must have positive diagonal entries")
  }
  rownames(vectors) <- variant_names
  list(values = values, vectors = vectors, rank = length(values),
       diag = diag_ld, variant_names = variant_names)
}

.utl_ld_eigen_from_matrix <- function(A, variant_names, validation_tol,
                                      context = "LD") {
  E <- eigen(A, symmetric = TRUE)
  .utl_normalize_ld_eigen(E, variant_names, validation_tol, context)
}

.utl_ld_factor <- function(ld_eigen, row_scale = NULL) {
  F <- sweep(ld_eigen$vectors, 2L, sqrt(ld_eigen$values), "*")
  if (!is.null(row_scale)) F <- F * as.numeric(row_scale)
  F
}

.utl_eigen_projection <- function(values, vectors, q, y) {
  if (!is.numeric(values) || length(values) < 1L ||
      !is.matrix(vectors) || ncol(vectors) != length(values) ||
      nrow(vectors) != length(q) || any(!is.finite(values)) ||
      any(!is.finite(vectors)) || any(values <= 0) ||
      any(!is.finite(q)) || !is.finite(y)) {
    stop("eigen projection inputs are invalid")
  }
  z <- as.vector(crossprod(vectors, q))
  q_projected <- as.vector(vectors %*% z)
  residual <- q - q_projected
  list(range_residual = sqrt(sum(residual^2)) /
         max(1, sqrt(sum(q^2))),
       projected_schur = as.numeric(y - sum(z^2 / values)))
}
