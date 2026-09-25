#' Construct a fixed theta prior from maximum squared theta estimates.
#'
#' The first ancestry is the reference and `beta = M theta`.  If supplied,
#' `trait_variance` is used directly.  Otherwise `beta_hat` is transformed to
#' theta coordinates and each free-coordinate variance is the maximum squared
#' estimate across variants.  The returned components use the canonical D, C,
#' S, and N basis dictionary without rescaling.
#'
#' @param delta Named finite transfer coefficients for non-reference traits.
#' @param reference One non-empty reference trait name.
#' @param trait_variance Optional named non-negative theta-coordinate variances.
#' @param beta_hat Optional SNP-by-ancestry effect matrix used to derive maxima.
#' @param weights Optional named component weights; `NULL` gives uniform weights.
#' @param validation_tol Positive validation tolerance.
#' @return An object of class `utl_theta_prior` with fixed theta components and
#'   max-theta construction diagnostics.
#' @export
utl_theta_trait_fixed_prior <- function(delta, reference, trait_variance = NULL,
                                        beta_hat = NULL, weights = NULL,
                                        validation_tol = 1e-10) {
  if (!is.numeric(delta) || !is.null(dim(delta)) || !length(delta) ||
      any(!is.finite(delta)) || is.null(names(delta)) || anyNA(names(delta)) ||
      any(names(delta) == "") || anyDuplicated(names(delta))) {
    stop("delta must be a strictly named finite numeric vector")
  }
  if (!is.character(reference) || length(reference) != 1L ||
      is.na(reference) || reference == "" || reference %in% names(delta)) {
    stop("reference must be one non-empty trait absent from delta names")
  }
  T <- length(delta) + 1L
  if (T < 2L || T > 5L) stop("theta fixed prior supports T from 2 to 5")
  if (!is.numeric(validation_tol) || length(validation_tol) != 1L ||
      !is.finite(validation_tol) || validation_tol <= 0) {
    stop("validation_tol must be one finite positive number")
  }
  ancestry_order <- c(reference, names(delta))
  theta_order <- paste0("theta_", ancestry_order)
  M <- diag(T)
  M[-1L, 1L] <- delta
  dimnames(M) <- list(ancestry_order, theta_order)

  theta_hat <- NULL
  if (!is.null(beta_hat)) {
    if (!is.matrix(beta_hat) || !is.numeric(beta_hat) ||
        any(!is.finite(beta_hat)) || is.null(colnames(beta_hat)) ||
        is.null(rownames(beta_hat)) || anyNA(rownames(beta_hat)) ||
        any(rownames(beta_hat) == "") || anyDuplicated(rownames(beta_hat)) ||
        anyNA(colnames(beta_hat)) || anyDuplicated(colnames(beta_hat)) ||
        !identical(sort(colnames(beta_hat)), sort(ancestry_order))) {
      stop("beta_hat must be a finite SNP-by-ancestry matrix with matching names")
    }
    beta_hat <- beta_hat[, ancestry_order, drop = FALSE]
    theta_hat <- beta_hat %*% t(solve(M))
    colnames(theta_hat) <- theta_order
  }
  if (!is.null(trait_variance)) {
    if (!is.numeric(trait_variance) || !is.null(dim(trait_variance)) ||
        length(trait_variance) != T || any(!is.finite(trait_variance)) ||
        any(trait_variance < 0) || is.null(names(trait_variance)) ||
        !setequal(names(trait_variance), theta_order)) {
      stop("trait_variance must be named, finite, non-negative, and theta ordered")
    }
    trait_variance <- trait_variance[theta_order]
    Vfree <- as.numeric(trait_variance)
    names(Vfree) <- theta_order
    V_source <- "trait_variance_external"
  } else {
    if (is.null(theta_hat)) stop("supply trait_variance or beta_hat")
    Vfree <- apply(theta_hat^2, 2L, max)
    names(Vfree) <- theta_order
    V_source <- "beta_hat_max_theta_squared"
  }

  e <- diag(T)
  rownames(e) <- colnames(e) <- theta_order
  B <- list()
  B$D <- e
  for (a in names(delta)) {
    j <- match(a, ancestry_order)
    other <- setdiff(seq.int(2L, T), j)
    B[[paste0("N_", a)]] <- cbind(e[, 1L] - delta[[a]] * e[, j],
                                   e[, other, drop = FALSE])
    colnames(B[[paste0("N_", a)]]) <- c(theta_order[[1L]], theta_order[other])
  }
  B$C <- e[, 1L, drop = FALSE]
  colnames(B$C) <- theta_order[[1L]]
  for (a in names(delta)) {
    j <- match(a, ancestry_order)
    B[[paste0("S_", a)]] <- e[, j, drop = FALSE]
    colnames(B[[paste0("S_", a)]]) <- paste0("theta_", a)
  }
  Utheta <- lapply(B, function(x) {
    z <- x %*% diag(Vfree[colnames(x)], nrow = ncol(x)) %*% t(x)
    dimnames(z) <- list(theta_order, theta_order)
    z
  })
  component_names <- names(Utheta)
  if (is.null(weights)) {
    weights <- setNames(rep(1 / length(Utheta), length(Utheta)), component_names)
  } else {
    if (!is.numeric(weights) || !is.null(dim(weights)) ||
        length(weights) != length(Utheta) || any(!is.finite(weights)) ||
        any(weights < 0) || is.null(names(weights)) ||
        !setequal(names(weights), component_names) || sum(weights) <= 0) {
      stop("weights must be named finite non-negative values aligned to components")
    }
    weights <- weights[component_names]
    weights <- weights / sum(weights)
  }
  prior <- utl_theta_prior(delta, M, Utheta, weights, validation_tol,
                           allow_zero_rank = TRUE)
  prior$B <- B
  prior$Vfree <- Vfree
  prior$Vtheta <- diag(Vfree, nrow = length(Vfree))
  dimnames(prior$Vtheta) <- list(theta_order, theta_order)
  prior$Sigma <- Utheta
  prior$Sigma_beta <- prior$Ubeta
  prior$theta_hat <- theta_hat
  prior$trait_variance_source <- V_source
  prior$trait_fixed_prior <- list(
    reference = reference, delta = delta, Vfree = Vfree,
    source = V_source, weights = weights,
    component_names = component_names,
    formula = "Sigma_k = B_k diag(Vfree) t(B_k)"
  )
  prior
}
