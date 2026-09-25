.utl_new_engine_upstream <- function() {
  if (!requireNamespace("mvsusieR", quietly = TRUE)) {
    stop("mvsusieR (>= 0.3) is required for the new engine")
  }
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("susieR (>= 0.16.5) is required for the new engine")
  }

  f <- list(
    mv_precompute = utils::getFromNamespace("precompute_eigen_cache", "mvsusieR"),
    mv_loglik = utils::getFromNamespace("loglik_precomputed", "mvsusieR"),
    mv_posterior = utils::getFromNamespace("posterior_precomputed", "mvsusieR"),
    susie_optimize_scalar = utils::getFromNamespace(
      "optimize_scalar_prior_variance", "susieR"
    )
  )
  expected <- list(
    mv_precompute = c("svs", "V_structure", "is_common_cov", "max_cache_gb"),
    mv_loglik = c("betahat", "V_scalar", "eigen_cache", "BQ_cache"),
    mv_posterior = c("betahat", "V_scalar", "eigen_cache", "pi_V_post",
                     "em_var_wt", "reduce_params", "BQ_cache"),
    susie_optimize_scalar = c("V_init", "estimate_prior_method",
                              "neg_loglik_fn", "loglik_fn", "optim_init",
                              "optim_bounds", "optim_scale",
                              "check_null_threshold")
  )
  for (nm in names(f)) {
    if (!is.function(f[[nm]]) ||
        !identical(names(formals(f[[nm]])), expected[[nm]])) {
      stop("upstream function formals do not match the pinned adapter ABI: ", nm)
    }
  }
  f
}
