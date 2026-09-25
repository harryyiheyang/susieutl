.utl_new_engine_component_moments <- function(s, V, cache) {
  if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V <= 0) {
    stop("new engine moments require one positive finite scalar V")
  }
  b <- .utl_new_engine_betahat(s, cache)
  p <- cache$p
  d <- cache$d
  K <- cache$K
  if (!is.list(cache$mv_cache$components) ||
      length(cache$mv_cache$components) != K ||
      length(cache$component_rank) != K ||
      length(cache$component_factors) != K || length(cache$Pj) != p ||
      !is.function(cache$upstream$mv_posterior)) {
    stop("new engine moments cache does not contain the pinned posterior ABI")
  }
  mn <- array(0, c(p, K, d),
              dimnames = list(cache$variant_names, cache$component_names,
                              cache$theta_names))
  cv <- array(0, c(p, K, d, d),
              dimnames = list(cache$variant_names, cache$component_names,
                              cache$theta_names, cache$theta_names))
  pi <- cbind(numeric(p), rep(1, p))
  lowrank_mean <- lowrank_covariance <- NULL
  if (any(cache$component_rank < d)) {
    if (!is.null(cache$common_context) && isTRUE(cache$common_context$ok) &&
        !is.null(cache$common_context_V) &&
        identical(as.numeric(V), as.numeric(cache$common_context_V))) {
      zlow <- utl_cpp_common_component_moments(s, cache$common_context)
    } else {
      zlow <- utl_cpp_component_moments(
        s, cache$Pj, cache$component_factors, V
      )
    }
    lowrank_mean <- array(zlow$component_mean, c(p, K, d))
    lowrank_covariance <- array(zlow$component_covariance, c(p, K, d, d))
  }
  for (k in seq_len(K)) {
    if (cache$component_rank[[k]] < d) {
      mn[, k, ] <- lowrank_mean[, k, ]
      cv[, k, , ] <- lowrank_covariance[, k, , ]
      next
    }
    ec <- cache$mv_cache
    ec$components <- ec$components[k]
    z <- cache$upstream$mv_posterior(
      b, V, ec, pi, em_var_wt = NULL, reduce_params = NULL, BQ_cache = NULL
    )
    if (!is.matrix(z$mu) || !identical(dim(z$mu), c(p, d)) ||
        !is.array(z$mu2) || !identical(dim(z$mu2), c(p, d, d)) ||
        any(!is.finite(z$mu)) || any(!is.finite(z$mu2))) {
      stop("new engine component posterior returned invalid full moments")
    }
    C <- z$mu2
    for (q in seq_len(d)) for (r in seq_len(d)) {
      C[, q, r] <- C[, q, r] - z$mu[, q] * z$mu[, r]
    }
    mn[, k, ] <- z$mu
    cv[, k, , ] <- C
  }
  list(component_mean = mn, component_covariance = cv)
}
