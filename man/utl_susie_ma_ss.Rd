\name{utl_susie_ma_ss}
\alias{utl_susie_ma_ss}
\alias{print.utl_theta_ma_ss_fit}
\title{Fit theta-coordinate multi-ancestry SuSiE from sufficient statistics}
\usage{
utl_susie_ma_ss(XtX_list, Xty_list, yty, n, prior, L = 5L, V = 1,
  estimate_prior_variance = TRUE, max_iter = 100L, tol = 1e-4,
  check_null_threshold = 0, prior_tol = 1e-9, verbose = TRUE,
  residual_variance = NULL, estimate_residual_variance = TRUE,
  estimate_prior_mixture_weights = FALSE, prune_threshold = 1e-8,
  prior_variance_mode = c("scalar", "pooled_theta", "multi_theta"),
  tau_init = NULL,
  zeroing_warmup = 3L,
  input_mode = c("exact_ss", "rss_pseudo_ss"),
  fixed_v_null_check = TRUE, weight_prune_start = NULL,
  weight_prune_interval = 1L,
  null_weight_mode = c("global", "per_ser", "none"),
  mom_sa_control = NULL, estimate_pattern_weight = FALSE)

\method{print}{utl_theta_ma_ss_fit}(x, ...)
}
\arguments{
  \item{XtX_list}{Named ancestry-specific LD crossproduct matrices.}
  \item{Xty_list}{Named ancestry-specific genotype-response crossproducts.}
  \item{yty}{Named ancestry-specific response sums of squares.}
  \item{n}{Named ancestry-specific sample sizes.}
  \item{prior}{A prior returned by \code{utl_theta_pattern_prior()} or constructed manually by \code{utl_theta_prior()}.}
  \item{L}{Number of IBSS single-effect slots.}
  \item{V}{One finite positive initial prior-variance scale, broadcast according to the selected mode.}
  \item{estimate_prior_variance}{Whether to optimize slot scales.}
  \item{max_iter}{Maximum number of IBSS sweeps.}
  \item{tol}{Global-ELBO convergence tolerance.}
  \item{check_null_threshold}{Zero-scale comparison threshold.}
  \item{prior_tol}{Final activity tolerance. After final sequential alignment, a slot with no retained prior-variance coordinate greater than prior_tol is inactive and contributes zero posterior inclusion mass, credible sets, theta moments, and pattern posterior mass.}
  \item{verbose}{Whether to print sweep diagnostics.}
  \item{residual_variance}{Optional named positive starting residual variance vector.}
  \item{estimate_residual_variance}{Whether to update residual variances from expected residual sums of squares.}
  \item{estimate_prior_mixture_weights}{Whether to optimize component weights with \code{mixsqp}.}
  \item{prune_threshold}{Deprecated compatibility argument; pattern removal is disabled by the user contract and this value is ignored after validation.}
  \item{prior_variance_mode}{Per-slot scalar, shared theta-coordinate, or per-slot theta-coordinate prior variance.}
  \item{tau_init}{Withdrawn argument. Non-\code{NULL} values are an error.}
  \item{zeroing_warmup}{Number of completed sweeps before theta-coordinate modes may compare coordinate-zero and all-zero candidates.}
  \item{input_mode}{Input validation mode. \code{exact_ss} validates exact sufficient-statistics inputs; \code{rss_pseudo_ss} records range and Schur compatibility diagnostics without stopping or warning on augmented-Gram mismatch.}
  \item{fixed_v_null_check}{For fixed scalar V, compare nominal V with the exact null each sweep; FALSE retains the historical fixed-V path.}
  \item{weight_prune_start}{Deprecated compatibility argument retained for validation; permanent pattern removal is disabled and this value is ignored.}
  \item{weight_prune_interval}{Deprecated compatibility argument retained for validation; permanent pattern removal is disabled and this value is ignored.}
  \item{null_weight_mode}{NULL-likelihood weight scheme: \code{global} (the default), \code{per_ser}, or legacy-compatible \code{none}. NULL is auxiliary and is not included in theta posterior \code{rho}.}
  \item{mom_sa_control}{Optional empty-list alias for the MOM-SA q1 pattern-weight estimator. \code{NULL} retains the default categorical SER path unless \code{estimate_pattern_weight = TRUE}.}
  \item{estimate_pattern_weight}{Whether to enable formal MOM-SA q1 pattern weights. The estimator uses all six UTL6 patterns, a fixed positive scalar \code{V}, direct sufficient statistics, and no prior-variance, residual-variance, or mixture-weight updates. Fitted posterior quantities are conditional on the plug-in weights and do not propagate their uncertainty.}
  \item{x}{A fitted object.}
  \item{...}{Unused.}
}
\value{\code{utl_susie_ma_ss} returns a \code{utl_theta_ma_ss_fit} object.}
\description{Fits theta-coordinate multi-ancestry summary-statistic SuSiE with fixed covariance components, optional residual-variance estimation, and optional live mixture-weight optimization. The formal MOM-SA opt-in returns \code{weight_model = "mom_sa"}; \code{fit$weights} retains the fixed prior and \code{fit$mom_sa_slot} records the exact NE/D moment-tilt diagnostics.}
