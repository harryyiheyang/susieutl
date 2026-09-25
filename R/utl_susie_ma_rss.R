#' Fit theta-coordinate multi-ancestry SuSiE from per-SNP RSS summaries.
#'
#' This bridge converts ancestry-specific per-SNP Z scores and sample sizes to
#' sufficient statistics without rounding the per-SNP sample sizes.  It then
#' calls [utl_susie_ma_ss()] in RSS pseudo-sufficient-statistics mode, using
#' one rounded median sample size and response total per ancestry.
#'
#' @param summary_stats Named list of data frames with `SNP`, `Zscore`, and `N`.
#' @param R_list Named ancestry-specific LD matrices or eigen-factor lists.
#'   Each eigen-factor list contains `values` and `vectors`.
#' @param prior A [utl_theta_pattern_prior()] or [utl_theta_prior()] object.
#' @param L Number of IBSS single-effect slots.
#' @param V Initial prior-variance scale.
#' @param estimate_prior_variance Whether to optimize prior variance.
#' @param max_iter Maximum number of IBSS sweeps.
#' @param tol Global-ELBO convergence tolerance.
#' @param check_null_threshold Zero-scale comparison threshold.
#' @param prior_tol Final activity tolerance.
#' @param verbose Whether to print sweep diagnostics.
#' @param residual_variance Optional named residual variance vector.
#' @param estimate_residual_variance Whether to estimate residual variances.
#' @param estimate_prior_mixture_weights Whether to estimate mixture weights.
#' @param prune_threshold Deprecated compatibility argument; pattern removal is
#'   disabled by the user contract and this value is ignored after validation.
#' @param prior_variance_mode Prior-variance mode.
#' @param tau_init Withdrawn argument forwarded to the sufficient-statistics API.
#' @param zeroing_warmup Completed sweeps before zeroing candidates.
#' @param fixed_v_null_check For fixed scalar V, compare nominal V with the
#'   exact null each sweep; FALSE retains the historical fixed-V path.
#' @param weight_prune_start Deprecated compatibility argument retained for
#'   validation; permanent pattern removal is disabled and this value is ignored.
#' @param weight_prune_interval Deprecated compatibility argument retained for
#'   validation; permanent pattern removal is disabled and this value is ignored.
#' @param null_weight_mode NULL-likelihood weight scheme: `global` (the
#'   default), `per_ser`, or legacy-compatible `none`.
#' @param top_patterns_per_snp Optional positive number of patterns retained
#'   within each SNP for posterior moment reweighting; `NULL` keeps all.
#' @param estimate_pattern_weight Whether to enable formal MOM-SA pattern
#'   weights. RSS input is rejected for this opt-in; use [utl_susie_ma_ss()].
#' @param fiteigen Optional final-position alias for an ancestry-named list of
#'   LD eigen-factor objects. `R_list` and `fiteigen` cannot both be supplied.
#' @return A `utl_theta_ma_ss_fit` object with an `rss_bridge` contract.
#' @export
utl_susie_ma_rss <- function(summary_stats, R_list = NULL, prior, L = 5L, V = 1,
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
                              fixed_v_null_check = TRUE,
                              weight_prune_start = NULL,
                              weight_prune_interval = 1L,
                              null_weight_mode = c("global", "per_ser", "none"),
                              fiteigen = NULL,
                              top_patterns_per_snp = NULL,
                              estimate_pattern_weight = FALSE) {
  if (!inherits(prior, "utl_theta_prior")) stop("prior must be a utl_theta_prior")
  if (!is.logical(estimate_pattern_weight) ||
      length(estimate_pattern_weight) != 1L ||
      is.na(estimate_pattern_weight)) {
    stop("estimate_pattern_weight must be TRUE or FALSE")
  }
  if (isTRUE(estimate_pattern_weight)) {
    stop("estimate_pattern_weight requires direct sufficient statistics via utl_susie_ma_ss")
  }
  an <- prior$ancestry_names
  if (!is.list(summary_stats) || is.null(names(summary_stats)) ||
      !identical(names(summary_stats), an)) {
    stop("summary_stats must be an ancestry-named list ordered by prior")
  }
  if (!is.null(R_list) && !is.null(fiteigen)) {
    stop("R_list and fiteigen cannot both be supplied")
  }
  if (is.null(R_list) && is.null(fiteigen)) {
    stop("one of R_list or fiteigen must be supplied")
  }
  R_input <- if (is.null(R_list)) fiteigen else R_list
  if (!is.list(R_input) || is.null(names(R_input)) || !identical(names(R_input), an)) {
    stop("R_list/fiteigen must be an ancestry-named list ordered by prior")
  }
  required <- c("SNP", "Zscore", "N")
  canonical <- NULL
  ss <- vector("list", length(an))
  names(ss) <- an
  RR <- vector("list", length(an))
  names(RR) <- an
  for (a in an) {
    dat <- summary_stats[[a]]
    if (!is.data.frame(dat) || !all(required %in% names(dat))) {
      stop("summary_stats for ancestry '", a, "' must contain SNP, Zscore, and N")
    }
    snp <- as.character(dat$SNP)
    if (length(snp) < 1L || anyNA(snp) || any(snp == "") || anyDuplicated(snp)) {
      stop("summary_stats SNP identifiers are invalid for ancestry '", a, "'")
    }
    if (is.null(canonical)) canonical <- snp
    if (!setequal(snp, canonical) || length(snp) != length(canonical)) {
      stop("summary_stats SNP sets do not exactly overlap")
    }
    ii <- match(canonical, snp)
    dat <- dat[ii, , drop = FALSE]
    z <- as.numeric(dat$Zscore)
    N <- as.numeric(dat$N)
    if (any(!is.finite(z)) || any(!is.finite(N)) || any(N <= 2)) {
      stop("Zscore and N must be finite, with N greater than two")
    }
    dat$SNP <- canonical
    dat$Zscore <- z
    dat$N <- N
    ss[[a]] <- dat
    R0 <- R_input[[a]]
    if (is.matrix(R0)) {
      R <- R0
      if (!is.numeric(R) || any(dim(R) != c(length(canonical), length(canonical))) ||
          is.null(rownames(R)) || is.null(colnames(R)) || anyDuplicated(rownames(R)) ||
          anyDuplicated(colnames(R)) || !setequal(rownames(R), canonical) ||
          !setequal(colnames(R), canonical)) {
        stop("R_list matrix variant names must exactly match summary SNPs for ancestry '", a, "'")
      }
      ir <- match(canonical, rownames(R))
      ic <- match(canonical, colnames(R))
      R <- R[ir, ic, drop = FALSE]
      if (any(!is.finite(R))) stop("R_list contains non-finite values for ancestry '", a, "'")
      R_scale <- max(1, max(abs(R)))
      if (max(abs(R - t(R))) > prior$validation_tol * R_scale) {
        stop("R_list must be symmetric for ancestry '", a, "'")
      }
      if (any(diag(R) <= 0)) {
        stop("R_list diagonal must be positive for ancestry '", a, "'")
      }
    } else if (is.list(R0) && !is.null(R0$values) && !is.null(R0$vectors)) {
      values <- as.numeric(R0$values)
      vectors <- R0$vectors
      if (!is.numeric(values) || length(values) < 1L ||
          any(!is.finite(values)) || !is.matrix(vectors) ||
          !is.numeric(vectors) || nrow(vectors) != length(canonical) ||
          ncol(vectors) != length(values) || any(!is.finite(vectors))) {
        stop("fiteigen entry for ancestry '", a,
             "' must contain finite values and a matching finite vectors matrix")
      }
      vnames <- rownames(vectors)
      if (!is.null(vnames)) {
        if (anyNA(vnames) || any(vnames == "") || anyDuplicated(vnames) ||
            !setequal(vnames, canonical)) {
          stop("fiteigen vectors row names must exactly match summary SNPs for ancestry '", a, "'")
        }
        vectors <- vectors[match(canonical, vnames), , drop = FALSE]
      } else {
        rownames(vectors) <- canonical
      }
      R <- vectors %*% sweep(t(vectors), 1L, values, "*")
    } else {
      stop("R_list/fiteigen entry for ancestry '", a,
           "' must be a matrix or a list with values and vectors")
    }
    dimnames(R) <- list(canonical, canonical)
    RR[[a]] <- R
  }
  N_list <- lapply(ss, function(x) as.numeric(x$N))
  n_by_ancestry <- setNames(vapply(N_list, stats::median, numeric(1)), an)
  c_by_ancestry <- n_by_ancestry - 1
  z_adjusted <- vector("list", length(an))
  names(z_adjusted) <- an
  N_diagnostic <- do.call(rbind, lapply(an, function(a) {
    N_a <- N_list[[a]]
    nt <- n_by_ancestry[[a]]
    ratio <- N_a / nt
    outlier <- ratio < 0.9 | ratio > 1.1
    summary_row <- data.frame(
      row_type = "summary", SNP = "__summary__", trait = a,
      rawN = nt, nt = nt, ratio = 1, direction = "summary",
      stringsAsFactors = FALSE
    )
    outlier_rows <- data.frame(
      row_type = rep("outlier", sum(outlier)), SNP = ss[[a]]$SNP[outlier],
      trait = rep(a, sum(outlier)), rawN = N_a[outlier],
      nt = rep(nt, sum(outlier)), ratio = ratio[outlier],
      direction = ifelse(ratio[outlier] < 0.9, "low", "high"),
      stringsAsFactors = FALSE
    )
    rbind(summary_row, outlier_rows)
  }))
  b_list <- list()
  XtX_list <- list()
  Xty_list <- list()
  yty <- c_by_ancestry
  R2_list <- sigma2_list <- d_list <- se_list <- list()
  for (a in an) {
    z <- ss[[a]]$Zscore
    N_a <- ss[[a]]$N
    nt <- n_by_ancestry[[a]]
    c_a <- c_by_ancestry[[a]]
    den <- z^2 + nt - 2
    if (any(!is.finite(den)) || any(den <= 0)) stop("RSS bridge denominator is invalid for ancestry '", a, "'")
    z_adj <- z * sqrt(c_a / den)
    b <- z_adj / sqrt(c_a)
    se <- rep(1 / sqrt(c_a), length(z))
    R2 <- z^2 / den
    sigma2 <- (nt - 1) / den
    d <- rep(c_a, length(z))
    z_adjusted[[a]] <- z_adj
    XtX <- c_a * RR[[a]]
    Xty <- sqrt(c_a) * z_adj
    dimnames(XtX) <- list(canonical, canonical)
    names(Xty) <- canonical
    names(b) <- canonical
    names(se) <- canonical
    names(R2) <- canonical
    names(sigma2) <- canonical
    names(d) <- canonical
    b_list[[a]] <- b
    se_list[[a]] <- se
    R2_list[[a]] <- R2
    sigma2_list[[a]] <- sigma2
    d_list[[a]] <- d
    XtX_list[[a]] <- XtX
    Xty_list[[a]] <- Xty
  }
  fit <- utl_susie_ma_ss(
    XtX_list = XtX_list, Xty_list = Xty_list, yty = yty, n = n_by_ancestry,
    prior = prior, L = L, V = V, estimate_prior_variance = estimate_prior_variance,
    max_iter = max_iter, tol = tol, check_null_threshold = check_null_threshold,
    prior_tol = prior_tol, verbose = verbose, residual_variance = residual_variance,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_mixture_weights = estimate_prior_mixture_weights,
    prune_threshold = prune_threshold, prior_variance_mode = prior_variance_mode,
    tau_init = tau_init, zeroing_warmup = zeroing_warmup,
    fixed_v_null_check = fixed_v_null_check,
    weight_prune_start = weight_prune_start,
    weight_prune_interval = weight_prune_interval,
    null_weight_mode = null_weight_mode,
    top_patterns_per_snp = top_patterns_per_snp,
    input_certificate = NULL,
    input_mode = "rss_pseudo_ss")
  fit$rss_bridge <- list(
    contract = "full_overlap_per_snp_n_rss_pseudo_ss",
    canonical = canonical, N_policy = "n_t=median(raw_N_t); c_t=n_t-1; yty_t=c_t",
    n_scalar_by_ancestry = n_by_ancestry,
    precision_by_ancestry = c_by_ancestry,
    response_scale = "constant-N RSS pseudo-SS; b=z_adjusted/sqrt(c), se=1/sqrt(c), d=c",
    N = N_list, b = b_list, se = se_list, R2 = R2_list,
    sigma2 = sigma2_list, d = d_list, z_adjusted = z_adjusted,
    N_diagnostic = N_diagnostic,
    alignment = "exact SNP-name alignment; no intersect/drop",
    R_input = R_input
  )
  fit
}
