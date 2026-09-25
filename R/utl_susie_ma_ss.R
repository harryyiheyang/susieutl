#' Fit a theta-coordinate multi-ancestry SuSiE model from sufficient statistics.
#'
#' The implementation fits finite multi-ancestry summary statistics with
#' native theta covariance components and exact joint SNP-by-component SER
#' updates.  Residual variances and, when requested, component mixture weights
#' are fit locally without mutating the supplied prior.
#'
#' @param XtX_list Named ancestry-specific LD crossproduct matrices.
#' @param Xty_list Named ancestry-specific genotype-response crossproducts.
#' @param yty Named ancestry-specific response sums of squares.
#' @param n Named ancestry-specific sample sizes.
#' @param prior A prior returned by [utl_theta_pattern_prior()] or constructed
#'   manually by [utl_theta_prior()].
#' @param L Number of IBSS single-effect slots.
#' @param V One positive initial prior-variance scale, broadcast according to
#'   `prior_variance_mode`.
#' @param estimate_prior_variance Whether to optimize each slot scale.
#' @param max_iter Maximum number of IBSS sweeps.
#' @param tol Non-negative global-ELBO increment used for convergence.
#' @param check_null_threshold Threshold used when comparing the selected scale
#'   with the exact zero-scale model.
#' @param prior_tol Final activity tolerance. After final sequential alignment,
#'   a slot with no retained prior-variance coordinate greater than prior_tol is
#'   inactive and contributes zero posterior inclusion mass, credible sets,
#'   theta moments, and pattern posterior mass.
#' @param verbose Print one line per IBSS sweep.
#' @param residual_variance Optional named starting residual variance vector.
#' @param estimate_residual_variance Whether to update residual variances from
#'   expected residual sums of squares.
#' @param estimate_prior_mixture_weights Whether to update component weights
#'   using `mixsqp`.
#' @param prune_threshold Deprecated compatibility argument; pattern removal is
#'   disabled by the user contract and this value is ignored after validation.
#' @param prior_variance_mode Per-slot scalar, shared theta-coordinate, or
#'   per-slot theta-coordinate variance mode.
#' @param tau_init Withdrawn argument. Non-`NULL` values are an error.
#' @param zeroing_warmup Number of completed sweeps before theta-coordinate
#'   modes may compare coordinate-zero and all-zero candidates.
#' @param input_mode Exact sufficient-statistics validation or RSS pseudo-SS
#'   validation. RSS pseudo-SS records range and Schur compatibility diagnostics
#'   in the result without stopping or warning on augmented-Gram mismatch.
#' @param fixed_v_null_check For fixed scalar V, compare the nominal V with
#'   the exact null each sweep; set FALSE to retain the historical fixed-V path.
#' @param weight_prune_start Deprecated compatibility argument retained for
#'   validation; permanent pattern removal is disabled and this value is ignored.
#' @param weight_prune_interval Deprecated compatibility argument retained for
#'   validation; permanent pattern removal is disabled and this value is ignored.
#' @param null_weight_mode NULL-likelihood weight scheme: `global` (the
#'   default), `per_ser`, or legacy-compatible `none`. NULL is auxiliary and
#'   is not included in theta posterior `rho`.
#' @param top_patterns_per_snp Optional positive number of patterns retained
#'   within each SNP for posterior moment reweighting; `NULL` keeps all.
#' @param mom_sa_control Optional empty list alias enabling the MOM-SA
#'   moment-tilt update. `NULL` is the default categorical SER path unless
#'   `estimate_pattern_weight = TRUE`.
#' @param estimate_pattern_weight Whether to enable the formal MOM-SA
#'   q1 pattern-weight estimator. This uses a fixed prior and no weight
#'   uncertainty propagation; `mom_sa_control = list()` remains an alias.
#' @return An object of class `utl_theta_ma_ss_fit`.
#' @export
utl_susie_ma_ss <- function(XtX_list, Xty_list, yty, n, prior, L = 5L, V = 1,
                            estimate_prior_variance = TRUE, max_iter = 100L,
                            tol = 1e-4, check_null_threshold = 0,
                            prior_tol = 1e-9, verbose = TRUE,
                            residual_variance = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_mixture_weights = FALSE,
                            prune_threshold = 1e-8,
                            prior_variance_mode = c("scalar", "pooled_theta",
                                                    "multi_theta"),
                            tau_init = NULL, zeroing_warmup = 3L,
                            input_mode = c("exact_ss", "rss_pseudo_ss"),
                            fixed_v_null_check = TRUE,
                            weight_prune_start = NULL,
                             weight_prune_interval = 1L,
                             null_weight_mode = c("global", "per_ser", "none"),
                             input_certificate = NULL,
                             top_patterns_per_snp = NULL,
                             mom_sa_control = NULL,
                             estimate_pattern_weight = FALSE) {
  if (!is.numeric(L) || length(L) != 1L || !is.finite(L) ||
      L != as.integer(L) || L < 1L) {
    stop("L must be one positive integer")
  }
  L <- as.integer(L)
  if (length(prior_variance_mode) == 1L &&
      identical(prior_variance_mode, "slot_theta")) {
    stop("slot_theta has been renamed; use multi_theta")
  }
  prior_variance_mode <- match.arg(prior_variance_mode,
                                    c("scalar", "pooled_theta", "multi_theta"))
  if (!is.null(tau_init)) {
    stop("tau_init has been withdrawn; use the scalar V initializer")
  }
  if (!is.numeric(V) || !is.null(dim(V)) || length(V) != 1L ||
      !is.finite(V) || V <= 0) {
    stop("V must be one finite positive initial value")
  }
  V <- rep(as.numeric(V), L)
  names(V) <- paste0("L", seq_len(L))
  if (!is.logical(estimate_prior_variance) || length(estimate_prior_variance) != 1L ||
      is.na(estimate_prior_variance)) {
    stop("estimate_prior_variance must be TRUE or FALSE")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter != as.integer(max_iter) || max_iter < 1L) {
    stop("max_iter must be one positive integer")
  }
  max_iter <- as.integer(max_iter)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("tol must be one finite non-negative number")
  }
  if (!is.numeric(check_null_threshold) || length(check_null_threshold) != 1L ||
      is.na(check_null_threshold) ||
      (is.infinite(check_null_threshold) && check_null_threshold > 0)) {
    stop("check_null_threshold must be finite or -Inf")
  }
  if (!is.numeric(prior_tol) || length(prior_tol) != 1L || !is.finite(prior_tol) ||
      prior_tol < 0) {
    stop("prior_tol must be one finite non-negative number")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  if (!is.logical(estimate_prior_mixture_weights) ||
      length(estimate_prior_mixture_weights) != 1L ||
      is.na(estimate_prior_mixture_weights)) {
    stop("estimate_prior_mixture_weights must be TRUE or FALSE")
  }
  if (!is.numeric(prune_threshold) || length(prune_threshold) != 1L ||
      !is.finite(prune_threshold) || prune_threshold < 0) {
    stop("prune_threshold must be one finite non-negative number")
  }
  if (!is.numeric(zeroing_warmup) || length(zeroing_warmup) != 1L ||
      !is.finite(zeroing_warmup) || zeroing_warmup != as.integer(zeroing_warmup) ||
      zeroing_warmup < 0) {
    stop("zeroing_warmup must be one non-negative integer")
  }
  zeroing_warmup <- as.integer(zeroing_warmup)
  if (!is.logical(fixed_v_null_check) || length(fixed_v_null_check) != 1L ||
      is.na(fixed_v_null_check)) {
    stop("fixed_v_null_check must be TRUE or FALSE")
  }
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
  input_mode <- match.arg(input_mode)
  null_weight_mode <- match.arg(null_weight_mode,
                                 c("global", "per_ser", "none"))
  mom_control <- if (is.null(mom_sa_control)) NULL else
    .utl_mom_sa_validate_control(mom_sa_control)
  if (!is.logical(estimate_pattern_weight) ||
      length(estimate_pattern_weight) != 1L ||
      is.na(estimate_pattern_weight)) {
    stop("estimate_pattern_weight must be TRUE or FALSE")
  }
  mom_enabled <- isTRUE(estimate_pattern_weight) || !is.null(mom_control)
  mom_source <- if (isTRUE(estimate_pattern_weight) && !is.null(mom_control)) {
    "both"
  } else if (isTRUE(estimate_pattern_weight)) {
    "formal"
  } else if (!is.null(mom_control)) {
    "pilot_alias"
  } else {
    NULL
  }
  top_patterns_per_snp <- if (is.null(top_patterns_per_snp)) 0L else
    as.integer(top_patterns_per_snp)

  op <- utl_ma_ss_operator(XtX_list, Xty_list, yty, n, prior,
                           input_mode = input_mode,
                           input_certificate = input_certificate)
  legacy_defaults <- missing(residual_variance) && missing(estimate_residual_variance)
  if (is.null(residual_variance)) {
    residual_variance <- setNames(rep(1, length(n)), prior$ancestry_names)
  } else if (!is.numeric(residual_variance) || length(residual_variance) != length(n) ||
             is.null(names(residual_variance)) ||
             !identical(names(residual_variance), prior$ancestry_names) ||
             any(!is.finite(residual_variance)) || any(residual_variance <= 0)) {
    stop("residual_variance must be a named, finite, positive vector ordered by ancestry")
  }
  if (legacy_defaults) estimate_residual_variance <- FALSE
  if (!is.logical(estimate_residual_variance) ||
      length(estimate_residual_variance) != 1L ||
      is.na(estimate_residual_variance)) {
    stop("estimate_residual_variance must be TRUE or FALSE")
  }
  if (mom_enabled) {
    expected_components <- c("D", "N_E", "N_A", "C", "S_E", "S_A")
    if (!identical(prior$component_names, expected_components) ||
        !identical(names(prior$Utheta), expected_components)) {
      stop("mom_sa_control requires the six UTL6 components D, N_E, N_A, C, S_E, S_A")
    }
    if (isTRUE(estimate_prior_variance) ||
        !identical(prior_variance_mode, "scalar") ||
        isTRUE(estimate_prior_mixture_weights) ||
        isTRUE(estimate_residual_variance) ||
        !identical(null_weight_mode, "none") ||
        !identical(input_mode, "exact_ss") ||
        top_patterns_per_snp > 0L || any(abs(residual_variance - 1) > 0)) {
      stop("MOM-SA pattern weights require fixed positive scalar V, direct SS, unit residual variance, all mixture updates disabled, and all patterns")
    }
  }
  if (!mom_enabled) {
  out <- .utl_theta_ibss(
    op, L, V, estimate_prior_variance, max_iter, tol, check_null_threshold,
    verbose, residual_variance = residual_variance,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_mixture_weights = estimate_prior_mixture_weights,
    prune_threshold = prune_threshold,
    prior_variance_mode = prior_variance_mode, tau_init = tau_init,
    zeroing_warmup = zeroing_warmup, prior_tol = prior_tol,
    fixed_v_null_check = fixed_v_null_check,
     weight_prune_start = weight_prune_start,
      weight_prune_interval = weight_prune_interval,
      null_weight_mode = null_weight_mode,
      top_patterns_per_snp = top_patterns_per_snp
   )
  } else {
    out <- .utl_theta_ibss(
      op, L, V, estimate_prior_variance, max_iter, tol, check_null_threshold,
      verbose, residual_variance = residual_variance,
      estimate_residual_variance = estimate_residual_variance,
      estimate_prior_mixture_weights = estimate_prior_mixture_weights,
      prune_threshold = prune_threshold,
      prior_variance_mode = prior_variance_mode, tau_init = tau_init,
      zeroing_warmup = zeroing_warmup, prior_tol = prior_tol,
      fixed_v_null_check = fixed_v_null_check,
      weight_prune_start = weight_prune_start,
      weight_prune_interval = weight_prune_interval,
      null_weight_mode = null_weight_mode,
      top_patterns_per_snp = top_patterns_per_snp,
      mom_sa_control = list(),
      mom_sa_source = mom_source
    )
  }
  out$M <- prior$M
  out$prior <- prior
  out$XtX_list <- op$XtX_list
  out$Xty_list <- op$Xty_list
  out$n <- op$n
  out$yty <- op$yty
  out$sigma2 <- out$residual_variance
  out$estimate_prior_variance <- estimate_prior_variance
  out$max_iter <- max_iter
  out$tol <- tol
  out$check_null_threshold <- check_null_threshold
  out$prior_tol <- prior_tol
  out$prior_variance_mode <- prior_variance_mode
  out$zeroing_warmup <- zeroing_warmup
  out$fixed_v_null_check <- fixed_v_null_check
  out$weight_prune_start <- weight_prune_start
  out$weight_prune_interval <- weight_prune_interval
  out$null_weight_mode <- null_weight_mode
  if (!is.null(mom_control)) out$mom_sa_control <- mom_control
  if (mom_enabled) out$estimate_pattern_weight <- TRUE
  if (top_patterns_per_snp > 0L) {
    out$top_patterns_per_snp <- top_patterns_per_snp
  }
  out$input_mode <- input_mode
  out$input_qc <- op$input_qc
  class(out) <- "utl_theta_ma_ss_fit"
  if (!mom_enabled) {
    out$ancestry_lbf <- utl_ancestry_lbf(out, bf_gate = 0)
    out$lfsr <- utl_lfsr(out)
  }
  out
}

#' @export
print.utl_theta_ma_ss_fit <- function(x, ...) {
  cat("UTL theta multi-ancestry sufficient-statistics fit\n")
  cat("  p:", nrow(x$theta_mean), "slots:", x$L,
      "status:", x$status, "top SNP:", names(x$pip)[which.max(x$pip)], "\n")
  invisible(x)
}
