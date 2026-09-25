.utl_theta_top_pattern_reweight <- function(rho, lbf, weights, top_patterns_per_snp) {
  p <- nrow(rho)
  K <- ncol(rho)
  if (top_patterns_per_snp <= 0L || top_patterns_per_snp >= K) {
    return(list(rho = rho, correction = 0))
  }
  top <- min(as.integer(top_patterns_per_snp), K)
  out <- matrix(0, p, K, dimnames = dimnames(rho))
  correction <- 0
  for (j in seq_len(p)) {
    eta <- as.numeric(lbf[j, ]) + log(weights)
    ord <- order(-eta, seq_len(K), na.last = TRUE)
    sel <- ord[seq_len(top)]
    log_full <- .utl_logsumexp(eta)
    log_selected <- .utl_logsumexp(eta[sel])
    logr <- log_selected - log_full
    correction <- correction - sum(rho[j, ]) * logr
    out[j, sel] <- sum(rho[j, ]) * exp(eta[sel] - log_selected)
  }
  list(rho = out, correction = correction)
}

.utl_theta_ser_update_cached_mv <- function(s, Pj, prior, V, moments, weights,
                                            cache, top_patterns_per_snp = 0L) {
  if (!.utl_new_engine_cache_matches(Pj, prior, cache)) {
    stop("new engine cache does not match the current precision/prior state")
  }
  common_ctx <- cache$common_context
  if (!is.null(common_ctx) && isTRUE(common_ctx$ok) &&
      !is.null(cache$common_context_V) &&
      identical(as.numeric(V), as.numeric(cache$common_context_V))) {
    if (moments) {
      zm <- utl_cpp_common_ser_lbf_moments(s, common_ctx)
      lbf <- zm$lbf
      dimnames(lbf) <- list(rownames(s), prior$component_names)
      z <- utl_cpp_ser_reduce(s, Pj, lbf, zm$component_mean,
                              zm$component_covariance, weights, TRUE,
                              top_patterns_per_snp)
      z$mean <- array(z$mean, c(nrow(s), ncol(s)))
      z$second <- array(z$second, c(nrow(s), ncol(s), ncol(s)))
      z$component_mean <- array(zm$component_mean,
                                c(nrow(s), length(prior$F), ncol(s)))
      z$component_covariance <- array(zm$component_covariance,
                                      c(nrow(s), length(prior$F), ncol(s), ncol(s)))
      z$lbf <- lbf
      dimnames(z$rho) <- dimnames(lbf)
      dimnames(z$mean) <- list(rownames(s), prior$theta_names)
      dimnames(z$second) <- list(rownames(s), prior$theta_names,
                                 prior$theta_names)
      dimnames(z$component_mean) <- list(rownames(s), prior$component_names,
                                         prior$theta_names)
      dimnames(z$component_covariance) <- list(
        rownames(s), prior$component_names, prior$theta_names,
        prior$theta_names)
      z$pattern_posterior <- setNames(as.numeric(z$pattern_posterior),
                                      prior$component_names)
      z$pip <- setNames(as.numeric(z$pip), rownames(s))
      return(z)
    }
    zl <- utl_cpp_common_ser_lbf(s, common_ctx)
    lbf <- zl$lbf
    dimnames(lbf) <- list(rownames(s), prior$component_names)
    z <- utl_cpp_ser_reduce(s, Pj, lbf, numeric(), numeric(), weights, FALSE, 0L)
    z$lbf <- lbf
    dimnames(z$rho) <- dimnames(lbf)
    z$pattern_posterior <- setNames(as.numeric(z$pattern_posterior),
                                    prior$component_names)
    z$pip <- setNames(as.numeric(z$pip), rownames(s))
    return(z)
  }
  lbf <- .utl_new_engine_lbf(s, V, weights, cache)
  p <- nrow(s)
  d <- ncol(s)
  K <- length(prior$F)
  variant_names <- rownames(s)
  theta_names <- prior$theta_names
  component_names <- prior$component_names
  if (!is.matrix(lbf) || !identical(dim(lbf), c(p, K)) ||
      any(!is.finite(lbf)) || !identical(rownames(lbf), variant_names) ||
      !identical(colnames(lbf), component_names)) {
    stop("new engine log Bayes factors have invalid dimensions or names")
  }
  lp <- sweep(lbf, 2L, log(weights), "+") - log(p)
  logZ <- .utl_logsumexp(as.vector(lp))
  if (!is.finite(logZ)) stop("theta SER joint normalization is non-finite")
  rho <- exp(lp - logZ)
  dimnames(rho) <- dimnames(lbf)
  if (any(!is.finite(rho)) || abs(sum(rho) - 1) > 1e-10) {
    stop("theta SER joint posterior failed to normalize")
  }
  rw <- .utl_theta_top_pattern_reweight(rho, lbf, weights,
                                        if (moments) top_patterns_per_snp else 0L)
  rho <- rw$rho
  out <- list(lbf = lbf, rho = rho, logZ = logZ,
              pattern_posterior = colSums(rho), pip = rowSums(rho),
              weights = weights, correction = rw$correction)
  if (!moments) return(out)

  cm <- .utl_new_engine_component_moments(s, V, cache)
  required_dims <- c(p, K, d)
  if (!is.array(cm$component_mean) ||
      !identical(dim(cm$component_mean), required_dims) ||
      !is.array(cm$component_covariance) ||
      !identical(dim(cm$component_covariance), c(p, K, d, d)) ||
      any(!is.finite(cm$component_mean)) ||
      any(!is.finite(cm$component_covariance))) {
    stop("new engine component moments have invalid dimensions or values")
  }
  dimnames(cm$component_mean) <- list(variant_names, component_names,
                                      theta_names)
  dimnames(cm$component_covariance) <- list(
    variant_names, component_names, theta_names, theta_names
  )
  mean <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
  second <- array(0, c(p, d, d),
                  dimnames = list(variant_names, theta_names, theta_names))
  e <- 0
  for (j in seq_len(p)) for (k in seq_len(K)) {
    if (rho[j, k] == 0) next
    mu <- as.numeric(cm$component_mean[j, k, ])
    C <- matrix(cm$component_covariance[j, k, , ], d, d)
    Tjk <- C + tcrossprod(mu)
    mean[j, ] <- mean[j, ] + rho[j, k] * mu
    second[j, , ] <- second[j, , ] + rho[j, k] * Tjk
    e <- e + rho[j, k] * (sum(s[j, ] * mu) -
                           0.5 * sum(Pj[[j]] * Tjk))
  }
  KL <- e - logZ + rw$correction
  if (!is.finite(e) || !is.finite(KL) || any(!is.finite(mean)) ||
      any(!is.finite(second))) {
    stop("new engine SER posterior moments are non-finite")
  }
  out$mean <- mean
  out$second <- second
  out$e <- e
  out$KL <- KL
  out$component_mean <- cm$component_mean
  out$component_covariance <- cm$component_covariance
  out
}

.utl_theta_component_posterior <- function(s, P, F, V) {
  d <- length(s)
  if (!is.numeric(V) || (length(V) != 1L && length(V) != d) ||
      any(!is.finite(V)) || any(V < 0)) {
    stop("theta SER prior variance must be a finite scalar or theta-length vector")
  }
  if (length(V) == 1L && V == 0) {
    return(list(lbf = 0, mean = numeric(d), covariance = matrix(0, d, d)))
  }
  if (ncol(F) == 0L) {
    return(list(lbf = 0, mean = numeric(d), covariance = matrix(0, d, d)))
  }
  if (length(V) == d) {
    if (all(V == 0)) {
      return(list(lbf = 0, mean = numeric(d), covariance = matrix(0, d, d)))
    }
    G <- diag(sqrt(V), nrow = d) %*% F
  } else {
    G <- sqrt(V) * F
  }
  K <- diag(ncol(G)) + crossprod(G, P %*% G)
  K <- (K + t(K)) / 2
  R <- chol(K)
  b <- as.vector(crossprod(G, s))
  Kinv_b <- backsolve(R, forwardsolve(t(R), b))
  Kinv_Gt <- backsolve(R, forwardsolve(t(R), t(G)))
  mu <- as.vector(G %*% Kinv_b)
  C <- G %*% Kinv_Gt
  C <- (C + t(C)) / 2
  logdet <- 2 * sum(log(diag(R)))
  lbf <- 0.5 * (sum(b * Kinv_b) - logdet)
  if (!is.finite(lbf) || any(!is.finite(mu)) || any(!is.finite(C))) {
    stop("theta SER component calculation produced non-finite values")
  }
  list(lbf = lbf, mean = mu, covariance = C)
}

.utl_theta_ser_update <- function(s, Pj, prior, V, moments = TRUE,
                                  weights = prior$weights, cache = NULL,
                                  top_patterns_per_snp = 0L) {
  if (!is.matrix(s) || !is.numeric(s) || any(!is.finite(s))) {
    stop("theta SER score must be a finite numeric matrix")
  }
  p <- nrow(s)
  d <- ncol(s)
  K <- length(prior$F)
  component_names <- prior$component_names
  if (p < 1L || d != length(prior$theta_names) || length(Pj) != p) {
    stop("theta SER score dimensions do not match the theta prior")
  }
  if (!is.numeric(V) || (length(V) != 1L && length(V) != d) ||
      any(!is.finite(V)) || any(V < 0)) {
    stop("V must be one finite non-negative number or theta-length vector")
  }
  if (!is.numeric(weights) || length(weights) != K ||
      any(!is.finite(weights)) || any(weights < 0) ||
      is.null(names(weights)) || !identical(names(weights), component_names) ||
      abs(sum(weights) - 1) > 1e-10) {
    stop("weights must be a named non-negative vector summing to one")
  }
  variant_names <- rownames(s)
  if (is.null(variant_names)) variant_names <- as.character(seq_len(p))
  theta_names <- prior$theta_names
  lbf <- matrix(0, p, K, dimnames = list(variant_names, component_names))

  if (length(V) == 1L && V > 0 && !is.null(cache) &&
      identical(cache$route, "cached_mv")) {
    if (is.null(rownames(s))) rownames(s) <- variant_names
    if (is.null(colnames(s))) colnames(s) <- theta_names
    return(.utl_theta_ser_update_cached_mv(
      s, Pj, prior, V, moments, weights, cache, top_patterns_per_snp
    ))
  }

  if ((length(V) == 1L && V == 0) ||
      (length(V) == d && all(V == 0))) {
    rho <- outer(rep(1 / p, p), weights)
    dimnames(rho) <- dimnames(lbf)
    out <- list(lbf = lbf, rho = rho, logZ = 0, e = 0, KL = 0,
                pattern_posterior = colSums(rho), pip = rowSums(rho),
                weights = weights, correction = 0)
    if (moments) {
      out$mean <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
      out$second <- array(0, dim = c(p, d, d),
                          dimnames = list(variant_names, theta_names, theta_names))
      out$component_mean <- array(0, dim = c(p, K, d),
                                  dimnames = list(variant_names, component_names,
                                                  theta_names))
      out$component_covariance <- array(0, dim = c(p, K, d, d),
                                        dimnames = list(variant_names,
                                                        component_names,
                                                        theta_names, theta_names))
    }
    return(out)
  }

  for (j in seq_len(p)) {
    if (!is.matrix(Pj[[j]]) || any(dim(Pj[[j]]) != c(d, d)) ||
        any(!is.finite(Pj[[j]]))) {
      stop("theta SER precision blocks must be finite d by d matrices")
    }
  }
  z <- utl_cpp_ser(s, Pj, prior$F, V, weights, moments,
                   if (moments) top_patterns_per_snp else 0L)
  lbf <- z$lbf
  rho <- z$rho
  dimnames(lbf) <- list(variant_names, component_names)
  dimnames(rho) <- dimnames(lbf)
  if (!is.finite(z$logZ) || any(!is.finite(rho)) ||
      abs(sum(rho) - 1) > 1e-10) {
    stop("theta SER joint posterior failed to normalize")
  }
  out <- list(lbf = lbf, rho = rho, logZ = z$logZ,
              pattern_posterior = setNames(as.numeric(z$pattern_posterior), component_names),
              pip = setNames(as.numeric(z$pip), variant_names),
              weights = weights)
  if (!moments) return(out)
  mean <- matrix(z$mean, p, d,
                 dimnames = list(variant_names, theta_names))
  second <- array(z$second, c(p, d, d),
                  dimnames = list(variant_names, theta_names, theta_names))
  component_mean <- array(z$component_mean, c(p, K, d),
                          dimnames = list(variant_names, component_names,
                                           theta_names))
  component_covariance <- array(z$component_covariance, c(p, K, d, d),
                                dimnames = list(variant_names, component_names,
                                                theta_names, theta_names))
  if (!is.finite(z$e) || !is.finite(z$KL) || any(!is.finite(mean)) ||
      any(!is.finite(second)) || any(!is.finite(component_mean)) ||
      any(!is.finite(component_covariance))) {
    stop("theta SER posterior moments are non-finite")
  }
  out$mean <- mean
  out$second <- second
  out$e <- z$e
  out$KL <- z$KL
  out$component_mean <- component_mean
  out$component_covariance <- component_covariance
  out
}

.utl_mom_sa_validate_control <- function(control) {
  if (!is.list(control) || length(control) != 0L) {
    stop("mom_sa_control must be NULL or an empty list()")
  }
  invisible(list())
}

.utl_theta_mom_sa_raw_components <- function(s, Pj, prior, V, weights,
                                             cache = NULL) {
  p <- nrow(s)
  d <- ncol(s)
  K <- length(prior$component_names)
  vn <- rownames(s)
  if (is.null(vn)) vn <- as.character(seq_len(p))
  tn <- prior$theta_names
  cn <- prior$component_names
  common_ctx <- if (is.null(cache)) NULL else cache$common_context
  common_v <- if (is.null(cache)) NULL else cache$common_context_V
  if (!is.null(common_ctx) && isTRUE(common_ctx$ok) &&
      !is.null(common_v) && identical(as.numeric(V), as.numeric(common_v))) {
    z <- utl_cpp_common_ser_lbf_moments(s, common_ctx)
    lbf <- matrix(as.numeric(z$lbf), p, K,
                  dimnames = list(vn, cn))
    cm <- array(as.numeric(z$component_mean), c(p, K, d),
                dimnames = list(vn, cn, tn))
    cv <- array(as.numeric(z$component_covariance), c(p, K, d, d),
                dimnames = list(vn, cn, tn, tn))
    return(list(lbf = lbf, component_mean = cm,
                component_covariance = cv, route = "common_component"))
  }
  if (!is.null(cache) && identical(cache$route, "cached_mv")) {
    lbf <- .utl_new_engine_lbf(s, V, weights, cache)
    cm <- .utl_new_engine_component_moments(s, V, cache)
    mean <- array(as.numeric(cm$component_mean), c(p, K, d),
                  dimnames = list(vn, cn, tn))
    covariance <- array(as.numeric(cm$component_covariance), c(p, K, d, d),
                        dimnames = list(vn, cn, tn, tn))
    dimnames(lbf) <- list(vn, cn)
    return(list(lbf = lbf, component_mean = mean,
                component_covariance = covariance,
                route = "cached_mv_component"))
  }
  raw <- .utl_theta_ser_update(
    s, Pj, prior, V, moments = TRUE, weights = weights, cache = cache,
    top_patterns_per_snp = 0L
  )
  list(
    lbf = matrix(as.numeric(raw$lbf), p, K, dimnames = list(vn, cn)),
    component_mean = array(as.numeric(raw$component_mean), c(p, K, d),
                           dimnames = list(vn, cn, tn)),
    component_covariance = array(as.numeric(raw$component_covariance),
                                 c(p, K, d, d),
                                 dimnames = list(vn, cn, tn, tn)),
    route = "native_combined_fallback"
  )
}

.utl_theta_mom_sa_update <- function(s, Pj, prior, V, weights,
                                     cache = NULL, retain_audit = TRUE) {
  if (!is.numeric(V) || length(V) != 1L || !is.finite(V) || V <= 0) {
    stop("MOM-SA requires one finite positive scalar V")
  }
  components <- prior$component_names
  expected <- c("D", "N_E", "N_A", "C", "S_E", "S_A")
  if (!identical(components, expected)) {
    stop("MOM-SA requires UTL6 components D, N_E, N_A, C, S_E, S_A")
  }
  raw <- .utl_theta_mom_sa_raw_components(s, Pj, prior, V, weights, cache)
  p <- nrow(s)
  d <- ncol(s)
  K <- length(components)
  vn <- rownames(s)
  if (is.null(vn)) vn <- as.character(seq_len(p))
  tn <- prior$theta_names
  A <- crossprod(prior$M, diag(c(1, 1, 0), nrow = nrow(prior$M)) %*% prior$M)
  A <- (A + t(A)) / 2
  if (!is.matrix(A) || !identical(dim(A), c(d, d)) || any(!is.finite(A))) {
    stop("MOM-SA ancestry tilt matrix is invalid")
  }
  target <- which(components %in% c("D", "N_E"))
  lbf <- matrix(raw$lbf, p, K,
                dimnames = list(vn, components))
  base_mean <- array(raw$component_mean, c(p, K, d),
                     dimnames = list(vn, components, tn))
  base_covariance <- array(raw$component_covariance, c(p, K, d, d),
                           dimnames = list(vn, components, tn, tn))
  tilted_mean <- base_mean
  tilted_covariance <- base_covariance
  t_value <- matrix(NA_real_, p, K,
                    dimnames = list(vn, components))
  z0 <- setNames(rep(NA_real_, K), components)
  for (k in target) {
    C0 <- V * prior$Utheta[[k]]
    z0[[k]] <- sum(A * t(C0))
    if (!is.finite(z0[[k]]) || z0[[k]] <= 0) {
      stop("MOM-SA z0 must be finite and positive for component ", components[[k]])
    }
    for (j in seq_len(p)) {
      m <- as.numeric(base_mean[j, k, ])
      C <- matrix(base_covariance[j, k, , ], d, d)
      T <- C + tcrossprod(m)
      tjk <- sum(A * t(C)) + drop(crossprod(m, A %*% m))
      if (!is.finite(tjk) || tjk <= 0) {
        stop("MOM-SA t must be finite and positive at SNP ", vn[[j]],
             " for component ", components[[k]])
      }
      v <- C %*% A %*% m
      mt <- m + 2 * v / tjk
      Tt <- T + 2 / tjk * (tcrossprod(m, v) + tcrossprod(v, m) +
                             C %*% A %*% C)
      Ct <- (Tt - tcrossprod(mt) + t(Tt - tcrossprod(mt))) / 2
      if (any(!is.finite(c(mt, Ct)))) {
        stop("MOM-SA tilted component moments are non-finite at SNP ", vn[[j]])
      }
      t_value[j, k] <- tjk
      tilted_mean[j, k, ] <- mt
      tilted_covariance[j, k, , ] <- Ct
      lbf[j, k] <- lbf[j, k] + log(tjk) - log(z0[[k]])
    }
  }
  if (any(!is.finite(lbf)) || any(!is.finite(tilted_mean)) ||
      any(!is.finite(tilted_covariance))) {
    stop("MOM-SA adjusted SER inputs are non-finite")
  }
  z <- utl_cpp_ser_reduce(
    s, Pj, lbf, as.numeric(tilted_mean), as.numeric(tilted_covariance),
    weights, TRUE, 0L
  )
  z$lbf <- matrix(as.numeric(z$lbf), p, K,
                  dimnames = list(vn, components))
  z$rho <- matrix(as.numeric(z$rho), p, K,
                  dimnames = list(vn, components))
  z$mean <- matrix(as.numeric(z$mean), p, d,
                   dimnames = list(vn, tn))
  z$second <- array(as.numeric(z$second), c(p, d, d),
                    dimnames = list(vn, tn, tn))
  z$component_mean <- tilted_mean
  z$component_covariance <- tilted_covariance
  z$pattern_posterior <- setNames(as.numeric(z$pattern_posterior), components)
  z$pip <- setNames(as.numeric(z$pip), vn)
  if (!is.finite(z$logZ) || !is.finite(z$e) || !is.finite(z$KL) ||
      any(!is.finite(z$rho)) || any(!is.finite(z$mean)) ||
      any(!is.finite(z$second))) {
    stop("MOM-SA adjusted SER posterior is non-finite")
  }
  z$mom_sa <- list(A = A, target_components = components[target],
                   z0 = z0, t = t_value, component_route = raw$route)
  if (isTRUE(retain_audit)) {
    z$mom_sa <- c(z$mom_sa, list(
      base_lbf = matrix(raw$lbf, p, K, dimnames = list(vn, components)),
      tilted_lbf = z$lbf, base_component_mean = base_mean,
      base_component_covariance = base_covariance,
      tilted_component_mean = tilted_mean,
      tilted_component_covariance = tilted_covariance,
      exact_KL = z$KL
    ))
  }
  z
}

.utl_theta_update_v <- function(s, Pj, prior, V, check_null_threshold,
                                weights = prior$weights, cache = NULL) {
  current_logZ <- .utl_theta_ser_update(
    s, Pj, prior, V, moments = FALSE, weights = weights, cache = cache
  )$logZ
  objective <- function(logV) {
    z <- .utl_theta_ser_update(
      s, Pj, prior, exp(logV), moments = FALSE, weights = weights,
      cache = cache
    )$logZ
    if (!is.finite(z)) stop("theta prior-variance objective is non-finite")
    -z
  }
  opt <- stats::optim(0, objective, method = "Brent", lower = -30, upper = 15)
  if (opt$convergence != 0L || !is.finite(opt$par) || !is.finite(opt$value)) {
    stop("theta prior-variance optimization failed")
  }
  candidate_V <- exp(opt$par)
  candidate_logZ <- .utl_theta_ser_update(
    s, Pj, prior, candidate_V, moments = FALSE, weights = weights,
    cache = cache
  )$logZ
  if (!is.finite(candidate_logZ)) stop("theta prior-variance candidate is non-finite")
  if (candidate_logZ < current_logZ) {
    selected_V <- V
    selected_logZ <- current_logZ
  } else {
    selected_V <- candidate_V
    selected_logZ <- candidate_logZ
  }
  if (check_null_threshold >= selected_logZ) {
    selected_V <- 0
    selected_logZ <- 0
  }
  list(V = selected_V, logZ = selected_logZ, current_logZ = current_logZ,
       candidate_V = candidate_V, candidate_logZ = candidate_logZ)
}
