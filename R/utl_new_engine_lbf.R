.utl_new_engine_prior_fingerprint <- function(prior) {
  lapply(lapply(prior$F, tcrossprod), serialize, connection = NULL, version = 3L)
}

.utl_new_engine_cache <- function(Pj, prior, upstream) {
  if (!is.list(Pj) || !length(Pj) || !inherits(prior, "utl_theta_prior")) {
    stop("new engine cache requires precision blocks and a theta prior")
  }
  p <- length(Pj)
  d <- length(prior$theta_names)
  K <- length(prior$F)
  vn <- names(Pj)
  if (is.null(vn)) vn <- as.character(seq_len(p))
  if (anyNA(vn) || any(vn == "") || anyDuplicated(vn)) {
    stop("new engine precision names must be unique and non-empty")
  }
  if (!is.list(upstream) || !is.function(upstream$mv_precompute) ||
      !is.function(upstream$mv_loglik)) {
    stop("new engine cache requires the upstream likelihood functions")
  }
  R <- S <- vector("list", p)
  spd <- TRUE
  for (j in seq_len(p)) {
    P <- Pj[[j]]
    if (!is.matrix(P) || !is.numeric(P) || any(dim(P) != c(d, d)) ||
        any(!is.finite(P))) {
      stop("theta SER precision blocks must be finite d by d matrices")
    }
    # A failed unperturbed factorization selects the existing exact engine.
    C <- tryCatch(chol(P), error = function(e) NULL)
    if (is.null(C)) {
      spd <- FALSE
      next
    }
    A <- chol2inv(C)
    B <- if (all(is.finite(A))) {
      tryCatch(chol(A), error = function(e) NULL)
    } else NULL
    if (is.null(B)) {
      spd <- FALSE
      next
    }
    R[[j]] <- C
    S[[j]] <- A
  }
  common <- spd && all(vapply(S, identical, logical(1), S[[1L]]))
  out <- structure(
    list(route = "exact", p = p, d = d, K = K,
         variant_names = vn, theta_names = prior$theta_names,
         component_names = prior$component_names, upstream = upstream,
         P_chol = R, svs = S, is_common_cov = common, mv_cache = NULL,
         Pj = Pj, component_factors = prior$F, component_rank = prior$rank,
         exact_reason = if (spd) NULL else "non_spd_precision"),
    class = "utl_new_engine_cache"
  )
  attr(out, "prior_fingerprint") <- .utl_new_engine_prior_fingerprint(prior)
  attr(out, "precision_fingerprint") <- lapply(
    Pj, serialize, connection = NULL, version = 3L
  )
  if (!spd) return(out)

  U <- lapply(prior$F, tcrossprod)
  out$mv_cache <- upstream$mv_precompute(S, U, common, max_cache_gb = 8)
  if (is.null(out$mv_cache)) {
    out$exact_reason <- "cache_limit"
  } else {
    out$route <- "cached_mv"
  }
  out
}

.utl_new_engine_betahat <- function(score, cache) {
  if (!is.matrix(score) || !is.numeric(score) || any(!is.finite(score)) ||
      !identical(dim(score), c(cache$p, cache$d))) {
    stop("new engine score must be a finite matrix matching the cache")
  }
  if (!identical(cache$route, "cached_mv")) {
    stop("new engine betahat requires a cached_mv cache")
  }
  if ((!is.null(rownames(score)) &&
       !identical(rownames(score), cache$variant_names)) ||
      (!is.null(colnames(score)) &&
       !identical(colnames(score), cache$theta_names))) {
    stop("new engine score names do not match the cache")
  }
  out <- matrix(0, cache$p, cache$d,
                dimnames = list(cache$variant_names, cache$theta_names))
  for (j in seq_len(cache$p)) {
    R <- cache$P_chol[[j]]
    out[j, ] <- backsolve(R, forwardsolve(t(R), as.numeric(score[j, ])))
  }
  out
}

.utl_new_engine_cache_matches <- function(Pj, prior, cache) {
  if (!is.list(cache) || !inherits(prior, "utl_theta_prior") ||
      !is.list(Pj) || length(Pj) != cache$p ||
      !identical(prior$theta_names, cache$theta_names) ||
      !identical(prior$component_names, cache$component_names) ||
      length(prior$F) != cache$K) {
    return(FALSE)
  }
  vn <- names(Pj)
  if (is.null(vn)) vn <- as.character(seq_along(Pj))
  identical(vn, cache$variant_names) &&
    identical(attr(cache, "prior_fingerprint"),
              .utl_new_engine_prior_fingerprint(prior)) &&
    identical(attr(cache, "precision_fingerprint"),
              lapply(Pj, serialize, connection = NULL, version = 3L))
}

.utl_new_engine_lbf <- function(s, V, weights, cache) {
  if (!identical(cache$route, "cached_mv") || !is.numeric(V) ||
      length(V) != 1L || !is.finite(V) || V <= 0) {
    stop("cached likelihood requires cached_mv precision and a positive scalar V")
  }
  if (!is.numeric(weights) || length(weights) != cache$K ||
      any(!is.finite(weights)) || any(weights < 0) ||
      !identical(names(weights), cache$component_names) ||
      abs(sum(weights) - 1) > 1e-10) {
    stop("weights must be a named non-negative vector summing to one")
  }
  b <- .utl_new_engine_betahat(s, cache)
  BQ <- NULL
  if (cache$is_common_cov) {
    BQ <- lapply(cache$mv_cache$components, function(x) b %*% x$Q)
  }
  ll <- cache$upstream$mv_loglik(b, V, cache$mv_cache, BQ_cache = BQ)
  if (!is.matrix(ll) || !identical(dim(ll), c(cache$p, cache$K + 1L)) ||
      any(!is.finite(ll))) {
    stop("upstream cached likelihood returned invalid values or dimensions")
  }
  # Column one is only the likelihood reference, never an extra prior pattern.
  out <- ll[, -1L, drop = FALSE] - ll[, 1L]
  dimnames(out) <- list(cache$variant_names, cache$component_names)
  out
}
