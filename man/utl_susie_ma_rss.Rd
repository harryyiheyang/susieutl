\name{utl_susie_ma_rss}
\alias{utl_susie_ma_rss}
\title{Fit theta-coordinate SuSiE from per-SNP RSS summaries}
\usage{
  utl_susie_ma_rss(summary_stats, R_list = NULL, prior, L = 5L, V = 1,
  estimate_prior_variance = TRUE, max_iter = 100L, tol = 1e-4,
  check_null_threshold = 0, prior_tol = 1e-9, verbose = TRUE,
  residual_variance = NULL, estimate_residual_variance = TRUE,
  estimate_prior_mixture_weights = FALSE, prune_threshold = 1e-8,
  prior_variance_mode = c("scalar", "pooled_theta", "multi_theta"),
  tau_init = NULL, zeroing_warmup = 3L,
  fixed_v_null_check = TRUE, weight_prune_start = NULL,
  weight_prune_interval = 1L, null_weight_mode = c("global", "per_ser", "none"),
  fiteigen = NULL, top_patterns_per_snp = NULL,
  estimate_pattern_weight = FALSE)
}
\arguments{
  \item{summary_stats}{An ancestry-named list of data frames with \code{SNP}, \code{Zscore}, and \code{N}.}
  \item{R_list}{An ancestry-named list of LD matrices or eigen-factor lists. Each eigen-factor list contains \code{values} and \code{vectors}.}
  \item{prior}{A \code{utl_theta_prior} object.}
  \item{L}{Number of IBSS single-effect slots.}
  \item{V}{Initial prior-variance scale.}
  \item{estimate_prior_variance}{Whether to optimize prior variance.}
  \item{max_iter}{Maximum number of IBSS sweeps.}
  \item{tol}{Global-ELBO convergence tolerance.}
  \item{check_null_threshold}{Zero-scale comparison threshold.}
  \item{prior_tol}{Final activity tolerance.}
  \item{verbose}{Whether to print sweep diagnostics.}
  \item{residual_variance}{Optional named residual variance vector.}
  \item{estimate_residual_variance}{Whether to estimate residual variances.}
  \item{estimate_prior_mixture_weights}{Whether to estimate mixture weights.}
  \item{prune_threshold}{Deprecated compatibility argument; pattern removal is disabled by the user contract and this value is ignored after validation.}
  \item{prior_variance_mode}{Prior-variance mode passed to the sufficient-statistics fit.}
  \item{tau_init}{Withdrawn argument passed to the sufficient-statistics API.}
  \item{zeroing_warmup}{Completed sweeps before zeroing candidates.}
  \item{fixed_v_null_check}{For fixed scalar V, compare nominal V with the exact null each sweep; FALSE retains the historical fixed-V path.}
  \item{weight_prune_start}{Deprecated compatibility argument retained for validation; permanent pattern removal is disabled and this value is ignored.}
  \item{weight_prune_interval}{Deprecated compatibility argument retained for validation; permanent pattern removal is disabled and this value is ignored.}
  \item{null_weight_mode}{NULL-likelihood weight scheme: \code{global}, \code{per_ser}, or \code{none}.}
  \item{fiteigen}{Optional final-position alias for an ancestry-named list of LD eigen-factor objects. \code{R_list} and \code{fiteigen} cannot both be supplied.}
  \item{top_patterns_per_snp}{Optional positive number of patterns retained within each SNP for posterior moment reweighting; \code{NULL} keeps all.}
  \item{estimate_pattern_weight}{Whether to enable formal MOM-SA pattern weights. RSS input is rejected for this opt-in; use \code{utl_susie_ma_ss()} with direct sufficient statistics.}
}
\value{A \code{utl_theta_ma_ss_fit} object with an \code{rss_bridge} contract, ancestry-specific rounded-median \code{n_t}/\code{yty_t}, and per-SNP RSS conversion vectors.}
\description{Converts exact per-SNP RSS quantities without rounding, computes \code{n_t = as.integer(round(median(N_t)))} and \code{yty_t = n_t - 1} separately for each ancestry, and fits using \code{utl_susie_ma_ss(input_mode = "rss_pseudo_ss")}.}
