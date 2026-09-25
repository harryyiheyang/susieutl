.utl_theta_active_slots <- function(prior_variance_mode, V = NULL, tau = NULL,
                                    tau2 = NULL, prior_tol = 1e-9, L = NULL) {
  if (identical(prior_variance_mode, "scalar")) return(as.numeric(V) > 0)
  if (identical(prior_variance_mode, "pooled_theta")) {
    if (is.null(L)) stop("L is required for pooled_theta active-state evaluation")
    return(rep(any(as.numeric(tau) > prior_tol), L))
  }
  if (identical(prior_variance_mode, "multi_theta")) {
    return(colSums(as.matrix(tau2) > prior_tol) > 0)
  }
  stop("unknown theta prior-variance mode")
}

.utl_theta_variance_for_slot <- function(mode, l, V, tau, tau2) {
  if (identical(mode, "scalar")) return(V[[l]])
  if (identical(mode, "pooled_theta")) return(tau)
  tau2[, l]
}

.utl_theta_slot_vhat <- function(slot_second) {
  dd <- dim(slot_second)
  if (length(dd) != 4L) {
    stop("slot theta second moments must be a four-dimensional array")
  }
  L <- dd[[1L]]
  p <- dd[[2L]]
  d <- dd[[3L]]
  out <- matrix(0, d, L,
                dimnames = list(dimnames(slot_second)[[3L]],
                                dimnames(slot_second)[[1L]]))
  for (l in seq_len(L)) for (j in seq_len(p)) {
    out[, l] <- out[, l] + diag(matrix(slot_second[l, j, , ], d, d))
  }
  if (any(!is.finite(out)) || any(out < 0)) {
    stop("slot theta posterior second moments are invalid")
  }
  out
}

.utl_theta_candidate_table <- function(Q, Pj, prior, current, full, eligible,
                                       weights) {
  d <- ncol(Q)
  if (length(current) != d || length(full) != d || length(eligible) != d ||
      any(!is.finite(current)) || any(current < 0) ||
      any(!is.finite(full)) || any(full < 0)) {
    stop("theta variance candidate state is invalid")
  }
  current <- setNames(as.numeric(current), prior$theta_names)
  full <- setNames(as.numeric(full), prior$theta_names)
  values <- list(current = current, full = full)
  for (q in seq_len(d)) if (isTRUE(eligible[[q]])) {
    zq <- full
    zq[[q]] <- 0
    values[[paste0("zero_", prior$theta_names[[q]])]] <- zq
  }
  values$all_zero <- setNames(numeric(d), prior$theta_names)
  candidates <- do.call(rbind, values)
  colnames(candidates) <- prior$theta_names
  logZ <- vapply(values, function(x) {
    .utl_theta_ser_update(Q, Pj, prior, x, moments = FALSE,
                          weights = weights)$logZ
  }, numeric(1))
  if (any(!is.finite(logZ))) stop("theta variance candidate logZ is non-finite")
  chosen_name <- names(logZ)[which.max(logZ)]
  list(candidates = candidates, logZ = logZ, chosen_name = chosen_name,
       chosen = values[[chosen_name]])
}

.utl_theta_pool_vhat <- function(vhat, chosen, absorbed, prior_tol) {
  if (!is.matrix(vhat) || !is.matrix(chosen) ||
      !identical(dim(vhat), dim(chosen)) || length(absorbed) != nrow(vhat)) {
    stop("pooled theta survivor state is invalid")
  }
  tau <- numeric(nrow(vhat))
  survivors <- chosen > prior_tol
  killed <- logical(nrow(vhat))
  for (q in seq_len(nrow(vhat))) {
    if (absorbed[[q]]) next
    keep <- survivors[q, ]
    if (any(keep)) {
      tau[[q]] <- mean(vhat[q, keep])
    } else {
      absorbed[[q]] <- TRUE
      killed[[q]] <- TRUE
    }
  }
  tau[absorbed] <- 0
  names(tau) <- rownames(vhat)
  names(absorbed) <- rownames(vhat)
  names(killed) <- rownames(vhat)
  list(tau = tau, absorbed = absorbed, survivors = survivors, killed = killed)
}

.utl_theta_scalar_candidates <- function(Q, Pj, prior, current, weights,
                                         cache = NULL,
                                         check_null_threshold = 0,
                                         estimate_prior_method = "optim") {
  if (!estimate_prior_method %in% c("optim", "simple")) {
    stop("estimate_prior_method must be optim or simple")
  }
  upstream <- if (!is.null(cache) && !is.null(cache$upstream)) {
    cache$upstream
  } else {
    .utl_new_engine_upstream()
  }
  if (!is.function(upstream$susie_optimize_scalar)) {
    stop("upstream scalar prior-variance optimizer is unavailable")
  }
  V_current <- as.numeric(current)
  if (length(V_current) != 1L || !is.finite(V_current) || V_current < 0) {
    stop("theta scalar prior-variance current value is invalid")
  }
  V_init_used <- if (V_current > 0) V_current else exp(0)
  calls <- 0L
  logZ_at <- function(V) {
    calls <<- calls + 1L
    z <- .utl_theta_ser_update(Q, Pj, prior, V, moments = FALSE,
                               weights = weights, cache = cache)$logZ
    if (!is.finite(z)) stop("theta scalar prior-variance objective is non-finite")
    z
  }
  neg_loglik_fn <- function(logV) -logZ_at(exp(logV))
  loglik_fn <- function(V) logZ_at(V)
  V_selected <- upstream$susie_optimize_scalar(
    V_init = V_init_used, estimate_prior_method = estimate_prior_method,
    neg_loglik_fn = neg_loglik_fn, loglik_fn = loglik_fn,
    optim_init = 0, optim_bounds = c(-30, 15), optim_scale = "log",
    check_null_threshold = check_null_threshold
  )
  if (length(V_selected) != 1L || !is.finite(V_selected) || V_selected < 0) {
    stop("upstream scalar prior-variance optimizer returned an invalid value")
  }
  logZ_current <- logZ_at(V_current)
  logZ_selected <- logZ_at(V_selected)
  logZ_null <- logZ_at(0)
  values <- c(current = V_current, optimized = V_selected, all_zero = 0)
  logZ <- c(current = logZ_current, optimized = logZ_selected,
            all_zero = logZ_null)
  if (V_selected == 0 && V_current == 0) chosen_name <- "current"
  else if (V_selected == V_current) chosen_name <- "current"
  else if (V_selected == 0) chosen_name <- "all_zero"
  else chosen_name <- "optimized"
  list(candidates = values, logZ = logZ, chosen_name = chosen_name,
       chosen = unname(V_selected), optimizer = "susieR::optimize_scalar_prior_variance",
       upstream_call_count = 1L, V_current = V_current,
       V_init_used = V_init_used, V_selected = unname(V_selected),
       logZ_selected = unname(logZ_selected), logZ_null = unname(logZ_null),
       check_null_threshold = check_null_threshold,
       null_selected = isTRUE(V_selected == 0))
}

.utl_theta_pip_condition <- function(message, mode, active, V_state, raw,
                                     slot_pip, iter) {
  structure(
    list(message = message, call = NULL, mode = mode, active = active,
         V = V_state, raw = raw,
         raw_range = if (length(raw)) range(raw) else c(NA_real_, NA_real_),
         slot_range = if (length(slot_pip)) range(slot_pip) else c(NA_real_, NA_real_),
         iter = iter, dims = dim(slot_pip)),
    class = c("utl_theta_pip_error", "error", "condition")
  )
}

.utl_theta_checked_pip <- function(slot_pip, active, mode, V_state, iter) {
  if (!is.matrix(slot_pip) || length(active) != nrow(slot_pip)) {
    stop(.utl_theta_pip_condition("slot PIP dimensions are invalid", mode,
                                  active, V_state, numeric(), slot_pip, iter))
  }
  if (any(!is.finite(slot_pip))) {
    stop(.utl_theta_pip_condition("slot PIP is non-finite", mode, active,
                                  V_state, numeric(), slot_pip, iter))
  }
  if (any(active)) {
    raw <- 1 - apply(1 - slot_pip[active, , drop = FALSE], 2L, prod)
  } else {
    raw <- numeric(ncol(slot_pip))
  }
  tol <- 100 * .Machine$double.eps * max(1, nrow(slot_pip), abs(raw))
  if (any(!is.finite(raw)) || any(raw < -tol) || any(raw > 1 + tol) ||
      any(slot_pip < -tol) || any(slot_pip > 1 + tol)) {
    stop(.utl_theta_pip_condition("posterior inclusion probability is outside tolerance",
                                  mode, active, V_state, raw, slot_pip, iter))
  }
  pip <- pmin(1, pmax(0, raw))
  list(pip = pip,
       diagnostic = list(iter = iter, mode = mode, active = active,
                         tolerance = tol, raw_range = range(raw),
                         slot_range = range(slot_pip), clamped = any(pip != raw)))
}

.utl_theta_pattern_summary <- function(pattern_by_slot, active) {
  if (any(active)) {
    list(posterior = colMeans(pattern_by_slot[active, , drop = FALSE]),
         all_null = FALSE)
  } else {
    list(posterior = setNames(numeric(ncol(pattern_by_slot)),
                              colnames(pattern_by_slot)), all_null = TRUE)
  }
}

.utl_theta_update_weights <- function(slot_lbf, slot_pip, active_slots,
                                      weights_by_slot, null_weights,
                                      update_count, component_names,
                                      null_weight_mode = "global") {
  if (!any(active_slots)) {
    return(list(weights_by_slot = weights_by_slot, null_weights = null_weights,
                full_weights_by_slot = NULL,
                active_components = rep(TRUE, length(component_names)),
                update_count = update_count, pruned = character(),
                prune_event = character(), updated = FALSE))
  }
  p <- ncol(slot_pip)
  ls <- which(active_slots)
  fit_mix <- function(LL, ow, old) {
    if (any(!is.finite(ow)) || any(ow < 0) || sum(ow) <= 0) {
      stop("alpha weights for prior-mixture optimization are invalid")
    }
    mix <- mixsqp::mixsqp(LL, w = ow, log = TRUE,
                          control = list(verbose = FALSE))
    if (is.null(mix$x) || length(mix$x) != ncol(LL) ||
        any(!is.finite(mix$x)) || any(mix$x < -1e-10) ||
        abs(sum(mix$x) - 1) > 1e-7) {
      stop("mixsqp prior-mixture optimization returned invalid weights")
    }
    q <- pmax(as.numeric(mix$x), 0)
    q <- q / sum(q)
    if (sum(q[-1L]) <= 0) {
      return(list(null = q[[1L]], weights = old,
                  full = setNames(q, c("NULL", component_names))))
    }
    q_nonnull <- q[-1L] / sum(q[-1L])
    list(null = q[[1L]], weights = setNames(q_nonnull, component_names),
         full = setNames(q, c("NULL", component_names)))
  }
  out_w <- weights_by_slot
  out_n <- null_weights
  out_full <- NULL
  if (null_weight_mode == "none") {
    LL <- matrix(0, length(ls) * p, length(component_names))
    ow <- numeric(length(ls) * p)
    for (i in seq_along(ls)) {
      l <- ls[[i]]
      jj <- ((i - 1L) * p + 1L):(i * p)
      LL[jj, ] <- matrix(slot_lbf[l, , ], p, length(component_names))
      ow[jj] <- slot_pip[l, ]
    }
    mix <- mixsqp::mixsqp(LL, w = ow, log = TRUE,
                          control = list(verbose = FALSE))
    if (is.null(mix$x) || length(mix$x) != length(component_names) ||
        any(!is.finite(mix$x)) || any(mix$x < -1e-10) ||
        abs(sum(mix$x) - 1) > 1e-7) {
      stop("mixsqp prior-mixture optimization returned invalid weights")
    }
    w_new <- pmax(as.numeric(mix$x), 0)
    w_new <- w_new / sum(w_new)
    out_w[,] <- rep(setNames(w_new, component_names), each = nrow(out_w))
    out_full <- out_w
  } else if (null_weight_mode == "per_ser") {
    for (l in ls) {
      LL <- cbind(0, matrix(slot_lbf[l, , ], p, length(component_names)))
      z <- fit_mix(LL, slot_pip[l, ], weights_by_slot[l, ])
      out_w[l, ] <- z$weights
      out_n[[l]] <- z$null
      if (is.null(out_full)) out_full <- matrix(0, nrow(out_w), length(component_names) + 1L)
      out_full[l, ] <- z$full
    }
  } else {
    LL <- matrix(0, length(ls) * p, length(component_names) + 1L)
    ow <- numeric(length(ls) * p)
    for (i in seq_along(ls)) {
      l <- ls[[i]]
      jj <- ((i - 1L) * p + 1L):(i * p)
      LL[jj, ] <- cbind(0, matrix(slot_lbf[l, , ], p, length(component_names)))
      ow[jj] <- slot_pip[l, ]
    }
    z <- fit_mix(LL, ow, weights_by_slot[1L, ])
    out_w[,] <- rep(z$weights, each = nrow(out_w))
    out_n[] <- z$null
    out_full <- matrix(rep(z$full, each = nrow(out_w)), nrow(out_w), length(component_names) + 1L,
                       byrow = FALSE,
                       dimnames = list(rownames(out_w), names(z$full)))
  }
  update_count <- update_count + 1L
  list(weights_by_slot = out_w, null_weights = out_n,
       full_weights_by_slot = out_full,
       active_components = rep(TRUE, length(component_names)),
       update_count = update_count, pruned = character(),
       prune_event = character(), updated = TRUE)
}

.utl_theta_update_sigma2 <- function(op, theta_bar, slot_mean, slot_second,
                                     active_slots) {
  T <- length(op$n)
  L <- dim(slot_mean)[[1L]]
  p <- dim(slot_mean)[[2L]]
  d <- dim(slot_mean)[[3L]]
  beta_bar <- theta_bar %*% t(op$prior$M)
  out <- numeric(T)
  names(out) <- op$prior$ancestry_names
  for (t in seq_len(T)) {
    b <- beta_bar[, t]
    q2 <- 0
    bslot <- vector("list", L)
    for (l in seq_len(L)) {
      Ml <- matrix(slot_mean[l, , ], p, d)
      bl <- as.vector(Ml %*% op$prior$M[t, ])
      bslot[[l]] <- bl
      if (!active_slots[[l]]) next
      for (j in seq_len(p)) {
        Tlj <- matrix(slot_second[l, j, , ], d, d)
        mb <- as.numeric(op$prior$M[t, ] %*% Tlj %*% op$prior$M[t, ])
        q2 <- q2 + op$XtX_list[[t]][j, j] * mb
      }
    }
    ls <- which(active_slots)
    if (length(ls) > 1L) for (u in seq_len(length(ls) - 1L)) {
      for (v in (u + 1L):length(ls)) {
        Am <- .utl_matvec(op$XtX_list[[t]], bslot[[ls[[v]]]], op$backend)
        q2 <- q2 + 2 * sum(bslot[[ls[[u]]]] * Am)
      }
    }
    erss <- op$yty[[t]] - 2 * sum(op$Xty_list[[t]] * b) + q2
    if (!is.finite(erss) || erss <= 0) {
      stop("estimated residual variance is non-positive for ancestry '",
           op$prior$ancestry_names[[t]], "'")
    }
    out[[t]] <- erss / op$n[[t]]
  }
  out
}

.utl_theta_zero_slot <- function(l, slot_rho, slot_lbf, slot_mean, slot_second,
                                 slot_component_mean,
                                 slot_component_covariance, slot_pip,
                                 pattern_by_slot, slot_KL) {
  slot_rho[l, , ] <- 0
  slot_lbf[l, , ] <- 0
  slot_mean[l, , ] <- 0
  slot_second[l, , , ] <- 0
  slot_component_mean[l, , , ] <- 0
  slot_component_covariance[l, , , , ] <- 0
  slot_pip[l, ] <- 0
  pattern_by_slot[l, ] <- 0
  slot_KL[[l]] <- 0
  list(slot_rho = slot_rho, slot_lbf = slot_lbf, slot_mean = slot_mean,
       slot_second = slot_second, slot_component_mean = slot_component_mean,
       slot_component_covariance = slot_component_covariance,
       slot_pip = slot_pip, pattern_by_slot = pattern_by_slot,
       slot_KL = slot_KL)
}

.utl_theta_ibss <- function(op, L, V, estimate_prior_variance, max_iter, tol,
                            check_null_threshold, verbose,
                            residual_variance = rep(1, length(op$n)),
                            estimate_residual_variance = FALSE,
                            estimate_prior_mixture_weights = FALSE,
                            prune_threshold = 1e-8,
                            prior_variance_mode = "scalar", tau_init = NULL,
                            zeroing_warmup = 3L, prior_tol = 1e-9,
                            fixed_v_null_check = TRUE,
                              weight_prune_start = NULL,
                              weight_prune_interval = 1L,
                              null_weight_mode = "global",
                              top_patterns_per_snp = 0L,
                              mom_sa_control = NULL,
                              mom_sa_source = NULL) {
  modes <- c("scalar", "pooled_theta", "multi_theta")
  if (!prior_variance_mode %in% modes) stop("unknown theta prior-variance mode")
  if (length(null_weight_mode) != 1L ||
      !null_weight_mode %in% c("global", "per_ser", "none")) {
    stop("null_weight_mode must be global, per_ser, or none")
  }
  if (!is.logical(fixed_v_null_check) || length(fixed_v_null_check) != 1L ||
      is.na(fixed_v_null_check)) stop("fixed_v_null_check must be TRUE or FALSE")
  if (!is.null(weight_prune_start) &&
      (!is.numeric(weight_prune_start) || length(weight_prune_start) != 1L ||
       !is.finite(weight_prune_start) || weight_prune_start < 1L ||
       weight_prune_start != as.integer(weight_prune_start))) {
    stop("weight_prune_start must be NULL or one positive integer")
  }
  if (!is.numeric(weight_prune_interval) || length(weight_prune_interval) != 1L ||
      !is.finite(weight_prune_interval) || weight_prune_interval < 1L ||
      weight_prune_interval != as.integer(weight_prune_interval)) {
    stop("weight_prune_interval must be one positive integer")
  }
  if (!is.null(weight_prune_start)) weight_prune_start <- as.integer(weight_prune_start)
  weight_prune_interval <- as.integer(weight_prune_interval)
  p <- op$p
  d <- ncol(op$h)
  K <- length(op$prior$weights)
  T <- length(op$n)
  if (!is.numeric(residual_variance) || length(residual_variance) != T ||
      is.null(names(residual_variance)) ||
      !identical(names(residual_variance), op$prior$ancestry_names) ||
      any(!is.finite(residual_variance)) || any(residual_variance <= 0)) {
    stop("residual_variance must be named, finite, positive, and ancestry ordered")
  }
  slot_names <- paste0("L", seq_len(L))
  variant_names <- op$variant_names
  theta_names <- op$prior$theta_names
  component_names <- op$prior$component_names
  V <- setNames(as.numeric(V), slot_names)
  V_nominal <- V
  tau <- setNames(rep(V[[1L]], d), theta_names)
  tau2 <- matrix(V[[1L]], d, L, dimnames = list(theta_names, slot_names))
  global_absorbed <- setNames(rep(FALSE, d), theta_names)
  local_absorbed <- matrix(FALSE, d, L, dimnames = list(theta_names, slot_names))
  variance_trace <- list(initial = if (prior_variance_mode == "scalar") V else
    if (prior_variance_mode == "pooled_theta") tau else tau2)
  vhat_trace <- list()
  candidate_trace <- list()
  absorbed_trace <- list(initial = if (prior_variance_mode == "pooled_theta")
    global_absorbed else local_absorbed)
  sweep_state_trace <- list()
  pip_diagnostic_trace <- list()

  sigma2 <- residual_variance
  sigma2_trace <- matrix(sigma2, nrow = 1L,
                         dimnames = list("initial", op$prior$ancestry_names))
  weights <- op$prior$initial_weights
  weights <- weights / sum(weights)
  weights_by_slot <- matrix(rep(weights, each = L), L, K,
                            dimnames = list(slot_names, component_names))
  mom_enabled <- !is.null(mom_sa_control)
  if (mom_enabled && is.null(mom_sa_source)) mom_sa_source <- "pilot_alias"
  mom_sa_final <- vector("list", L)
  names(mom_sa_final) <- slot_names
  null_weights <- setNames(rep(0, L), slot_names)
  initial_weights_by_slot <- weights_by_slot
  active_components <- rep(TRUE, K)
  weights_trace <- matrix(weights, nrow = 1L,
                          dimnames = list("initial", component_names))
  weights_by_slot_trace <- list(initial = weights_by_slot)
  null_weight_trace <- list(initial = null_weights)
  full_weights_by_slot <- cbind(null_weights, weights_by_slot)
  colnames(full_weights_by_slot) <- c("NULL", component_names)
  full_weights_by_slot_trace <- list(initial = full_weights_by_slot)
  pruned_components <- character()
  prune_events <- character()
  weight_update_count <- 0L

  slot_rho <- array(0, c(L, p, K),
                    dimnames = list(slot_names, variant_names, component_names))
  slot_lbf <- array(0, c(L, p, K), dimnames = dimnames(slot_rho))
  slot_mean <- array(0, c(L, p, d),
                     dimnames = list(slot_names, variant_names, theta_names))
  slot_second <- array(0, c(L, p, d, d),
                       dimnames = list(slot_names, variant_names,
                                       theta_names, theta_names))
  slot_component_mean <- array(0, c(L, p, K, d),
                               dimnames = list(slot_names, variant_names,
                                               component_names, theta_names))
  slot_component_covariance <- array(0, c(L, p, K, d, d),
                                     dimnames = list(slot_names, variant_names,
                                                     component_names,
                                                     theta_names, theta_names))
  slot_pip <- matrix(0, L, p, dimnames = list(slot_names, variant_names))
  pattern_by_slot <- matrix(0, L, K,
                            dimnames = list(slot_names, component_names))
  slot_KL <- setNames(numeric(L), slot_names)
  slot_top_pattern_kl <- setNames(numeric(L), slot_names)
  theta_bar <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
  trace <- data.frame(iter = integer(), elbo = numeric(), delta = numeric(),
                      d_theta = numeric(), d_pip = numeric(),
                      d_pattern = numeric(), d_weights = numeric(),
                      d_sigma2 = numeric())
  last_elbo <- NA_real_
  last_theta <- NULL
  last_pip <- NULL
  last_pattern <- NULL
  converged <- FALSE
  global <- NULL
  upstream <- if (prior_variance_mode == "scalar") {
    .utl_new_engine_upstream()
  } else NULL
  engine_cache <- NULL
  engine_cache_sigma2 <- NULL
  common_context_tried <- FALSE
  fixed_sigma <- !estimate_residual_variance
  fixed_sc <- if (fixed_sigma) op$scaled(sigma2) else NULL
  fixed_elbo <- !estimate_prior_variance && !estimate_residual_variance &&
    !estimate_prior_mixture_weights
  fixed_elbo_state <- if (fixed_elbo) {
    .utl_theta_fixed_elbo_state(op$n, op$yty, sigma2, L)
  } else NULL

  for (iter in seq_len(max_iter)) {
    sigma_for_t <- sigma2
    weights_for_t <- weights_by_slot
    V_for_t <- V
    tau_for_t <- tau
    tau2_for_t <- tau2
    sc <- if (fixed_sigma) fixed_sc else op$scaled(sigma_for_t)
    if (prior_variance_mode == "scalar" &&
        (is.null(engine_cache_sigma2) ||
         !identical(sigma_for_t, engine_cache_sigma2))) {
      engine_cache <- .utl_new_engine_cache(sc$Pj, op$prior, upstream)
      engine_cache_sigma2 <- sigma_for_t
    }
    if (!common_context_tried && fixed_sigma && prior_variance_mode == "scalar" &&
        !estimate_prior_variance && isTRUE(op$rss_common) &&
        !is.null(engine_cache)) {
      ctx <- utl_cpp_common_ser_context(
        sc$Pj[[1L]], op$prior$F, V_for_t[[1L]], moments = TRUE
      )
      if (isTRUE(ctx$ok)) {
        engine_cache$common_context <- ctx
        engine_cache$common_context_V <- V_for_t[[1L]]
      }
      common_context_tried <- TRUE
    }
    Q_slot <- vector("list", L)
    scalar_candidates <- setNames(vector("list", L), slot_names)
    for (l in seq_len(L)) {
      old <- matrix(slot_mean[l, , ], p, d,
                    dimnames = list(variant_names, theta_names))
      Q <- sc$h - sc$H(theta_bar - old)
      Q_slot[[l]] <- Q
      if (prior_variance_mode == "scalar") {
        if (estimate_prior_variance) {
          use_method <- "optim"
          use_threshold <- if (iter <= zeroing_warmup) -Inf else check_null_threshold
          use_current <- if (estimate_prior_variance) V_for_t[[l]] else V_nominal[[l]]
          cand <- .utl_theta_scalar_candidates(
            Q, sc$Pj, op$prior, use_current, weights_for_t[l, ],
            cache = engine_cache, check_null_threshold = use_threshold,
            estimate_prior_method = use_method
          )
          cand$phase <- use_method
          cand$V_nominal <- V_nominal[[l]]
          V_for_t[[l]] <- cand$chosen
          scalar_candidates[[l]] <- cand
        } else {
          V_for_t[[l]] <- V_nominal[[l]]
          scalar_candidates[[l]] <- list(
            phase = if (fixed_v_null_check && iter <= zeroing_warmup)
              "fixed_warmup" else "fixed_nokill",
            V_nominal = V_nominal[[l]],
            V_current = V_nominal[[l]], V_selected = V_nominal[[l]],
            null_selected = FALSE, check_null_threshold = check_null_threshold,
            optimizer = "not_called", upstream_call_count = 0L
          )
        }
      }
      active_l <- .utl_theta_active_slots(
        prior_variance_mode, V_for_t, tau_for_t, tau2_for_t, prior_tol, L
      )[[l]]
      if (!active_l) {
        if (fixed_elbo) {
          fixed_elbo_state <- .utl_theta_fixed_elbo_replace(
            fixed_elbo_state, l, Q, old, matrix(0, p, d), 0, 0
          )
        }
        z <- .utl_theta_zero_slot(
          l, slot_rho, slot_lbf, slot_mean, slot_second,
          slot_component_mean, slot_component_covariance, slot_pip,
          pattern_by_slot, slot_KL
        )
        list2env(z, environment())
        slot_top_pattern_kl[[l]] <- 0
        theta_bar <- theta_bar - old
        next
      }
      Vl <- .utl_theta_variance_for_slot(prior_variance_mode, l, V_for_t,
                                         tau_for_t, tau2_for_t)
      ser <- if (mom_enabled) {
        .utl_theta_mom_sa_update(
          Q, sc$Pj, op$prior, Vl, weights = weights_for_t[l, ],
          cache = if (prior_variance_mode == "scalar") engine_cache else NULL,
          retain_audit = FALSE
        )
      } else {
        .utl_theta_ser_update(
          Q, sc$Pj, op$prior, Vl, moments = TRUE,
          weights = weights_for_t[l, ],
          cache = if (prior_variance_mode == "scalar") engine_cache else NULL,
          top_patterns_per_snp = top_patterns_per_snp
        )
      }
      if (fixed_elbo) {
        fixed_elbo_state <- .utl_theta_fixed_elbo_replace(
          fixed_elbo_state, l, Q, old, ser$mean, ser$e, ser$KL
        )
      }
      slot_rho[l, , ] <- ser$rho
      slot_lbf[l, , ] <- ser$lbf
      slot_mean[l, , ] <- ser$mean
      slot_second[l, , , ] <- ser$second
      slot_component_mean[l, , , ] <- ser$component_mean
      slot_component_covariance[l, , , , ] <- ser$component_covariance
      slot_pip[l, ] <- ser$pip
      pattern_by_slot[l, ] <- ser$pattern_posterior
      slot_KL[[l]] <- ser$KL
      slot_top_pattern_kl[[l]] <- ser$correction
      theta_bar <- theta_bar - old + ser$mean
    }

    active_for_t <- .utl_theta_active_slots(
      prior_variance_mode, V_for_t, tau_for_t, tau2_for_t, prior_tol, L
    )
    pip_now <- .utl_theta_checked_pip(
      slot_pip, active_for_t, prior_variance_mode,
      list(V = V_for_t, tau = tau_for_t, tau2 = tau2_for_t), iter
    )
    pip <- setNames(pip_now$pip, variant_names)
    pip_diagnostic_trace[[sprintf("sweep_%03d", iter)]] <- pip_now$diagnostic
    pat <- .utl_theta_pattern_summary(pattern_by_slot, active_for_t)
    pattern_posterior <- pat$posterior
    all_null <- pat$all_null
    global <- if (fixed_elbo) {
      list(C0 = fixed_elbo_state$C0, Delta = numeric(L),
           Eloglik = fixed_elbo_state$Eloglik,
           ELBO = fixed_elbo_state$ELBO)
    } else {
      .utl_theta_global_elbo(
        sc$h, sc$H, sc$Pj, theta_bar, slot_mean, slot_second, slot_KL,
        op$n, op$yty, sigma_for_t
      )
    }

    vh <- .utl_theta_slot_vhat(slot_second)
    V_next <- V_for_t
    tau_next <- tau_for_t
    tau2_next <- tau2_for_t
    candidates <- if (prior_variance_mode == "scalar") scalar_candidates else
      setNames(vector("list", L), slot_names)
    if (estimate_prior_variance && prior_variance_mode == "pooled_theta" &&
        any(active_for_t)) {
      if (iter == 1L) {
        for (l in seq_len(L)) candidates[[l]] <- list(phase = "fixed")
      } else if (iter <= zeroing_warmup) {
        tau_next <- rowMeans(vh)
        tau_next[global_absorbed] <- 0
        names(tau_next) <- theta_names
        for (l in seq_len(L)) candidates[[l]] <- list(phase = "direct")
      } else {
        chosen <- matrix(0, d, L, dimnames = list(theta_names, slot_names))
        for (l in seq_len(L)) {
          full <- vh[, l]
          full[global_absorbed] <- 0
          current <- tau_for_t
          current[global_absorbed] <- 0
          eligible <- !global_absorbed & current > prior_tol
          candidates[[l]] <- .utl_theta_candidate_table(
            Q_slot[[l]], sc$Pj, op$prior, current, full, eligible, weights_for_t[l, ]
          )
          candidates[[l]]$phase <- "zero"
          chosen[, l] <- candidates[[l]]$chosen
        }
        pool <- .utl_theta_pool_vhat(vh, chosen, global_absorbed, prior_tol)
        tau_next <- pool$tau
        global_absorbed <- pool$absorbed
      }
    }
    if (estimate_prior_variance && prior_variance_mode == "multi_theta" &&
        any(active_for_t)) {
      if (iter == 1L) {
        for (l in seq_len(L)) candidates[[l]] <- list(phase = "fixed")
      } else if (iter <= zeroing_warmup) {
        tau2_next <- vh
        tau2_next[local_absorbed] <- 0
        for (l in seq_len(L)) candidates[[l]] <- list(phase = "direct")
      } else {
        for (l in seq_len(L)) {
          if (!active_for_t[[l]]) {
            candidates[[l]] <- list(phase = "inactive")
            next
          }
          full <- vh[, l]
          full[local_absorbed[, l]] <- 0
          current <- tau2_for_t[, l]
          current[local_absorbed[, l]] <- 0
          eligible <- !local_absorbed[, l] & current > prior_tol
          candidates[[l]] <- .utl_theta_candidate_table(
            Q_slot[[l]], sc$Pj, op$prior, current, full, eligible, weights_for_t[l, ]
          )
          candidates[[l]]$phase <- "zero"
          tau2_next[, l] <- candidates[[l]]$chosen
          local_absorbed[, l] <- local_absorbed[, l] |
            candidates[[l]]$chosen <= prior_tol
          tau2_next[local_absorbed[, l], l] <- 0
        }
      }
    }

    weights_next <- weights_for_t
    null_weights_next <- null_weights
    full_weights_next <- full_weights_by_slot
    active_components_next <- active_components
    d_weights <- NA_real_
    if (estimate_prior_mixture_weights) {
      wu <- .utl_theta_update_weights(
        slot_lbf, slot_pip, active_for_t, weights_for_t, null_weights,
        weight_update_count, component_names, null_weight_mode
      )
      weights_next <- wu$weights_by_slot
      null_weights_next <- wu$null_weights
      if (!is.null(wu$full_weights_by_slot)) full_weights_next <- wu$full_weights_by_slot
      active_components_next <- rep(TRUE, K)
      weight_update_count <- wu$update_count
      if (wu$updated) {
        weights_by_slot_trace[[length(weights_by_slot_trace) + 1L]] <- weights_next
        null_weight_trace[[length(null_weight_trace) + 1L]] <- null_weights_next
        full_weights_by_slot_trace[[length(full_weights_by_slot_trace) + 1L]] <- full_weights_next
        weights_trace <- rbind(weights_trace, colMeans(weights_next))
        rownames(weights_trace)[nrow(weights_trace)] <-
          sprintf("update_%03d", weight_update_count)
        d_weights <- max(abs(weights_next - weights_for_t))
      }
    }
    sigma_next <- sigma_for_t
    d_sigma2 <- NA_real_
    if (estimate_residual_variance) {
      sigma_next <- .utl_theta_update_sigma2(
        op, theta_bar, slot_mean, slot_second, active_for_t
      )
      d_sigma2 <- max(abs(sigma_next - sigma_for_t))
      sigma2_trace <- rbind(sigma2_trace, sigma_next)
      rownames(sigma2_trace)[nrow(sigma2_trace)] <- sprintf("next_%03d", iter)
    }

    if (iter == 1L) {
      delta <- d_theta <- d_pip <- d_pattern <- NA_real_
    } else {
      delta <- global$ELBO - last_elbo
      d_theta <- max(abs(theta_bar - last_theta))
      d_pip <- max(abs(pip - last_pip))
      d_pattern <- max(abs(pattern_posterior - last_pattern))
      if (any(!is.finite(c(delta, d_theta, d_pip, d_pattern)))) {
        stop("theta IBSS convergence trace is non-finite")
      }
    }
    trace_row <- data.frame(
      iter = iter, elbo = global$ELBO, delta = delta, d_theta = d_theta,
      d_pip = d_pip, d_pattern = d_pattern, d_weights = d_weights,
      d_sigma2 = d_sigma2,
      stringsAsFactors = FALSE
    )
    trace <- rbind(trace, trace_row)
    nm <- sprintf("sweep_%03d", iter)
    vhat_trace[[nm]] <- vh
    candidate_trace[[nm]] <- candidates
    V <- V_next
    tau <- tau_next
    tau2 <- tau2_next
    weights_by_slot <- weights_next
    null_weights <- null_weights_next
    full_weights_by_slot <- full_weights_next
    active_components <- active_components_next
    sigma2 <- sigma_next
    variance_trace[[nm]] <- if (prior_variance_mode == "scalar") V else
      if (prior_variance_mode == "pooled_theta") tau else tau2
    absorbed_trace[[nm]] <- if (prior_variance_mode == "pooled_theta")
      global_absorbed else local_absorbed
    sweep_state_trace[[nm]] <- list(
      V_for_t = V_for_t, tau_for_t = tau_for_t, tau2_for_t = tau2_for_t,
      weights_for_t = weights_for_t, null_weights_for_t = null_weights,
      sigma2_for_t = sigma_for_t,
      V_for_next = V, tau_for_next = tau, tau2_for_next = tau2,
      weights_for_next = weights_by_slot, null_weights_for_next = null_weights,
      full_weights_for_next = full_weights_by_slot,
      sigma2_for_next = sigma2,
      active_for_t = active_for_t
    )
    if (verbose) cat("theta IBSS iteration", iter, "ELBO",
                     format(global$ELBO, digits = 8), "\n")
    if (iter > 1L && delta < 0) {
      warning("theta IBSS global ELBO decreased; continuing because stopping uses ELBO only")
    }
    min_iter <- if (prior_variance_mode == "scalar") zeroing_warmup + 1L else zeroing_warmup + 2L
    legacy_done <- !is.na(delta) && iter >= min_iter &&
      delta >= 0 && delta < tol
    if (legacy_done) {
      converged <- TRUE
      break
    }
    last_elbo <- global$ELBO
    last_theta <- theta_bar
    last_pip <- pip
    last_pattern <- pattern_posterior
  }

  active_slots <- .utl_theta_active_slots(prior_variance_mode, V, tau, tau2,
                                          prior_tol, L)
  sc <- if (fixed_sigma) fixed_sc else op$scaled(sigma2)
  if (prior_variance_mode == "scalar" &&
      (is.null(engine_cache_sigma2) ||
       !identical(sigma2, engine_cache_sigma2))) {
    engine_cache <- .utl_new_engine_cache(sc$Pj, op$prior, upstream)
    engine_cache_sigma2 <- sigma2
  }
  theta_bar <- apply(slot_mean, c(2L, 3L), sum)
  dimnames(theta_bar) <- list(variant_names, theta_names)
  for (l in seq_len(L)) {
    old <- matrix(slot_mean[l, , ], p, d,
                  dimnames = list(variant_names, theta_names))
    Q <- sc$h - sc$H(theta_bar - old)
    if (!active_slots[[l]]) {
      if (fixed_elbo) {
        fixed_elbo_state <- .utl_theta_fixed_elbo_replace(
          fixed_elbo_state, l, Q, old, matrix(0, p, d), 0, 0
        )
      }
      z <- .utl_theta_zero_slot(
        l, slot_rho, slot_lbf, slot_mean, slot_second,
        slot_component_mean, slot_component_covariance, slot_pip,
        pattern_by_slot, slot_KL
      )
      list2env(z, environment())
      slot_top_pattern_kl[[l]] <- 0
      theta_bar <- theta_bar - old
      next
    }
    Vl <- .utl_theta_variance_for_slot(prior_variance_mode, l, V, tau, tau2)
    ser <- if (mom_enabled) {
      .utl_theta_mom_sa_update(
        Q, sc$Pj, op$prior, Vl,
        weights = weights_by_slot[l, ],
        cache = if (prior_variance_mode == "scalar") engine_cache else NULL,
        retain_audit = TRUE
      )
    } else {
      .utl_theta_ser_update(
        Q, sc$Pj, op$prior, Vl,
        weights = weights_by_slot[l, ],
        cache = if (prior_variance_mode == "scalar") engine_cache else NULL,
        top_patterns_per_snp = top_patterns_per_snp
      )
    }
    if (mom_enabled) mom_sa_final[[l]] <- ser$mom_sa
    if (fixed_elbo) {
      fixed_elbo_state <- .utl_theta_fixed_elbo_replace(
        fixed_elbo_state, l, Q, old, ser$mean, ser$e, ser$KL
      )
    }
    slot_rho[l, , ] <- ser$rho
    slot_lbf[l, , ] <- ser$lbf
    slot_mean[l, , ] <- ser$mean
    slot_second[l, , , ] <- ser$second
    slot_component_mean[l, , , ] <- ser$component_mean
    slot_component_covariance[l, , , , ] <- ser$component_covariance
    slot_pip[l, ] <- ser$pip
    pattern_by_slot[l, ] <- ser$pattern_posterior
    slot_KL[[l]] <- ser$KL
    slot_top_pattern_kl[[l]] <- ser$correction
    theta_bar <- theta_bar - old + ser$mean
  }
  global <- .utl_theta_global_elbo(
    sc$h, sc$H, sc$Pj, theta_bar, slot_mean, slot_second, slot_KL,
    op$n, op$yty, sigma2
  )
  pip_final <- .utl_theta_checked_pip(
    slot_pip, active_slots, prior_variance_mode,
    list(V = V, tau = tau, tau2 = tau2), "final_alignment"
  )
  pip <- setNames(pip_final$pip, variant_names)
  pip_diagnostic_trace$final_alignment <- pip_final$diagnostic
  pat <- .utl_theta_pattern_summary(pattern_by_slot, active_slots)
  pattern_posterior <- pat$posterior
  all_null <- pat$all_null
  vhat_trace$final_alignment <- .utl_theta_slot_vhat(slot_second)

  theta_covariance <- array(0, c(p, d, d),
                             dimnames = list(variant_names, theta_names, theta_names))
  for (l in which(active_slots)) for (j in seq_len(p)) {
    Tlj <- matrix(slot_second[l, j, , ], d, d)
    Mlj <- as.vector(slot_mean[l, j, ])
    theta_covariance[j, , ] <- theta_covariance[j, , ] + Tlj - tcrossprod(Mlj)
  }
  theta_variance <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
  for (j in seq_len(p)) theta_variance[j, ] <- diag(theta_covariance[j, , ])
  covariance_tol <- 100 * .Machine$double.eps * max(1, abs(theta_covariance))
  if (any(theta_variance < -covariance_tol)) {
    stop("theta posterior covariance has materially negative diagonal entries")
  }
  fit_converged <- converged
  fit_status <- if (fit_converged) "converged" else "max_iter"
  V_out <- if (prior_variance_mode == "scalar") V else {
    x <- vector("list", L)
    for (l in seq_len(L)) x[[l]] <- if (prior_variance_mode == "pooled_theta")
      tau else tau2[, l]
    names(x) <- slot_names
    x
  }
  weights_out <- if (null_weight_mode == "per_ser") weights_by_slot else
    setNames(weights_by_slot[1L, ], component_names)
  nw <- length(weights_by_slot_trace)
  names(weights_by_slot_trace) <- c("initial",
    if (nw > 1L) paste0("update_", seq_len(nw - 1L)))
  names(null_weight_trace) <- names(weights_by_slot_trace)
  names(full_weights_by_slot_trace) <- names(weights_by_slot_trace)
  out <- list(
    theta_mean = theta_bar, theta_covariance = theta_covariance,
    theta_variance = theta_variance, slot_rho = slot_rho, slot_lbf = slot_lbf,
    slot_mean = slot_mean, slot_second = slot_second, slot_KL = slot_KL,
    slot_component_mean = slot_component_mean,
    slot_component_covariance = slot_component_covariance,
    pip = pip, pip_by_slot = slot_pip,
    pattern_posterior = pattern_posterior,
    pattern_posterior_by_slot = pattern_by_slot,
    all_null = all_null, active_slots = active_slots, V = V_out,
    V_nominal = V_nominal, fixed_v_null_check = fixed_v_null_check,
    weight_prune_start = weight_prune_start, weight_prune_interval = weight_prune_interval,
    variance_trace = variance_trace, slot_vhat = vhat_trace,
    variance_candidates = candidate_trace,
    absorbed_trace = absorbed_trace, sweep_state_trace = sweep_state_trace,
    pip_diagnostic = pip_final$diagnostic,
    pip_diagnostic_trace = pip_diagnostic_trace,
    KL = sum(slot_KL[active_slots]), KL_by_slot = slot_KL,
    initial_weights = op$prior$initial_weights,
    initial_weights_by_slot = initial_weights_by_slot,
    weights = weights_out, weights_by_slot = weights_by_slot,
    weights_trace = weights_trace, weights_by_slot_trace = weights_by_slot_trace,
    null_weight = if (null_weight_mode == "per_ser") null_weights else null_weights[[1L]],
    null_weight_trace = null_weight_trace,
    full_mixture_weights = if (null_weight_mode == "none") weights_by_slot else
      if (null_weight_mode == "per_ser") full_weights_by_slot else full_weights_by_slot[1L, ],
    full_mixture_weights_by_slot_trace = full_weights_by_slot_trace,
    null_weight_mode = null_weight_mode, active_weights = weights_out,
    active_mask = rep(TRUE, K), weight_update_count = weight_update_count,
    prune_events = prune_events, pruned_components = pruned_components,
    pattern_pruning = "disabled_by_user_contract",
    residual_variance = sigma2, residual_variance_trace = sigma2_trace,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_mixture_weights = estimate_prior_mixture_weights,
    elbo = global$ELBO, eloglik = global$Eloglik, Delta = global$Delta,
    C0 = global$C0, trace = trace, delta = trace$delta[[nrow(trace)]],
     converged = fit_converged, status = fit_status,
    L = L, backend = op$backend$name,
    backend_calls = list(matvec_count = op$backend$matvec_count,
                         cpp_matvec_count = op$backend$cpp_matvec_count),
    covariance_negative_tolerance = covariance_tol
  )
  if (top_patterns_per_snp > 0L) {
    out$top_patterns_per_snp <- top_patterns_per_snp
    out$slot_top_pattern_kl <- slot_top_pattern_kl
  }
  if (mom_enabled) {
    out$weight_model <- "mom_sa"
    out$estimate_pattern_weight <- TRUE
    out$pattern_weight_method <- "mom_sa_q1_NE_D"
    out$prior_weight_estimation <- FALSE
    out$mom_sa_control_source <- mom_sa_source
    out$mom_sa_slot <- mom_sa_final
    out$mom_sa_objective_disclosure <-
      "MOM-SA effects and exact KL use the specified local moment tilt; posterior quantities are conditional on plug-in pattern weights and do not propagate weight uncertainty."
  }
  if (prior_variance_mode == "pooled_theta") {
    out$tau <- tau
    out$global_absorbed <- global_absorbed
  }
  if (prior_variance_mode == "multi_theta") {
    out$tau2 <- tau2
    out$local_absorbed <- local_absorbed
  }
  out
}
