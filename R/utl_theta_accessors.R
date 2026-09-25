.utl_check_theta_fit <- function(fit) {
  if (!inherits(fit, "utl_theta_ma_ss_fit")) {
    stop("fit must be a utl_theta_ma_ss_fit")
  }
  invisible(TRUE)
}

#' Extract theta-coordinate SNP posterior inclusion probabilities.
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @param by_slot Return the slot-by-SNP probabilities when `TRUE`.
#' @return A named numeric vector or slot-by-SNP matrix.
#' @export
utl_theta_pip <- function(fit, by_slot = FALSE) {
  .utl_check_theta_fit(fit)
  if (!is.logical(by_slot) || length(by_slot) != 1L || is.na(by_slot)) {
    stop("by_slot must be TRUE or FALSE")
  }
  if (by_slot) fit$pip_by_slot else fit$pip
}

#' Extract fixed theta-prior component posterior probabilities.
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @param by_slot Return the slot-by-component probabilities when `TRUE`.
#' @return A named numeric vector or slot-by-component matrix.
#' @export
utl_theta_pattern_posterior <- function(fit, by_slot = FALSE) {
  .utl_check_theta_fit(fit)
  if (!is.logical(by_slot) || length(by_slot) != 1L || is.na(by_slot)) {
    stop("by_slot must be TRUE or FALSE")
  }
  if (by_slot) fit$pattern_posterior_by_slot else fit$pattern_posterior
}

#' Extract theta-coordinate posterior effects.
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @return A list containing theta posterior means, full covariances, variances,
#'   and slot-level moments.
#' @export
utl_theta_effects <- function(fit) {
  .utl_check_theta_fit(fit)
  list(mean = fit$theta_mean, covariance = fit$theta_covariance,
       variance = fit$theta_variance, slot_mean = fit$slot_mean,
       slot_second = fit$slot_second)
}

#' Extract theta-coordinate fitting status and diagnostics.
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @return A list with convergence, objective, and backend diagnostics.
#' @export
utl_theta_status <- function(fit) {
  .utl_check_theta_fit(fit)
  out <- list(
    status = fit$status, converged = fit$converged, elbo = fit$elbo,
    eloglik = fit$eloglik, delta = fit$delta, trace = fit$trace, V = fit$V,
    KL = fit$KL, KL_by_slot = fit$KL_by_slot, Delta = fit$Delta,
    backend = fit$backend, backend_calls = fit$backend_calls,
    sigma2 = fit$sigma2, residual_variance_trace = fit$residual_variance_trace,
    prior_variance_mode = fit$prior_variance_mode, active_slots = fit$active_slots,
    all_null = fit$all_null, variance_trace = fit$variance_trace,
    slot_vhat = fit$slot_vhat, variance_candidates = fit$variance_candidates,
    absorbed_trace = fit$absorbed_trace,
    sweep_state_trace = fit$sweep_state_trace,
    pip_diagnostic = fit$pip_diagnostic,
    pip_diagnostic_trace = fit$pip_diagnostic_trace,
    initial_weights = fit$initial_weights, weights = fit$weights,
    weights_trace = fit$weights_trace, active_weights = fit$active_weights,
    active_mask = fit$active_mask,
    weight_update_count = fit$weight_update_count,
    prune_events = fit$prune_events,
    pruned_components = fit$pruned_components
  )
  if (identical(fit$weight_model, "mom_sa") ||
      identical(fit$weight_model, "mom_sa_pilot")) {
    out$weight_model <- fit$weight_model
    if (!is.null(fit$estimate_pattern_weight)) {
      out$estimate_pattern_weight <- fit$estimate_pattern_weight
    }
    if (!is.null(fit$pattern_weight_method)) {
      out$pattern_weight_method <- fit$pattern_weight_method
    }
    if (!is.null(fit$prior_weight_estimation)) {
      out$prior_weight_estimation <- fit$prior_weight_estimation
    }
    if (!is.null(fit$mom_sa_control_source)) {
      out$mom_sa_control_source <- fit$mom_sa_control_source
    }
    if (!is.null(fit$mom_sa_slot)) out$mom_sa_slot <- fit$mom_sa_slot
    if (!is.null(fit$mom_sa_objective_disclosure)) {
      out$mom_sa_objective_disclosure <- fit$mom_sa_objective_disclosure
    }
  }
  if (identical(fit$prior_variance_mode, "pooled_theta")) {
    out$tau <- fit$tau
    out$global_absorbed <- fit$global_absorbed
  }
  if (identical(fit$prior_variance_mode, "multi_theta")) {
    out$tau2 <- fit$tau2
    out$local_absorbed <- fit$local_absorbed
  }
  out
}
