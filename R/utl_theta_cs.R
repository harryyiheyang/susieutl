#' Build multi-LD credible sets from slot-level posterior inclusion masses.
#'
#' The set construction is performed independently for each slot.  For every
#' pair of variants, purity uses the maximum absolute LD over ancestries; the
#' reported set purity is then the minimum of those pairwise maxima. Duplicate
#' sets and slots with zero prior variance are omitted.
#'
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @param R_list Named list of ancestry-specific LD matrices.
#' @param coverage Desired cumulative slot mass.
#' @param min_abs_corr Minimum multi-LD purity.
#' @return A list with a set table, slot masses, and filtering diagnostics.
#' @export
utl_credible_sets <- function(fit, R_list, coverage = 0.95,
                              min_abs_corr = 0.5) {
  .utl_check_theta_fit(fit)
  if (!is.numeric(coverage) || length(coverage) != 1L || !is.finite(coverage) ||
      coverage <= 0 || coverage > 1) stop("coverage must be in (0, 1]")
  if (!is.numeric(min_abs_corr) || length(min_abs_corr) != 1L ||
      !is.finite(min_abs_corr) || min_abs_corr < 0 || min_abs_corr > 1) {
    stop("min_abs_corr must be in [0, 1]")
  }
  an <- fit$prior$ancestry_names
  vn <- names(fit$pip)
  p <- length(vn)
  if (!is.list(R_list) || is.null(names(R_list)) ||
      !identical(names(R_list), an)) {
    stop("R_list must be a named list ordered by fit ancestries")
  }
  for (t in seq_along(an)) {
    R <- R_list[[t]]
    if (!is.matrix(R) || !is.numeric(R) || any(dim(R) != c(p, p)) ||
        any(!is.finite(R)) || is.null(rownames(R)) ||
        !identical(rownames(R), vn) || !identical(colnames(R), vn)) {
      stop("each R_list matrix must be finite p by p LD with fit variant names")
    }
  }
  alpha <- fit$pip_by_slot
  if (!is.matrix(alpha) || nrow(alpha) != fit$L || ncol(alpha) != p) {
    stop("fit pip_by_slot has incompatible dimensions")
  }
  rows <- list()
  set_key <- character()
  dropped <- data.frame(slot = character(), reason = character(),
                        stringsAsFactors = FALSE)
  nr <- 0L
  active_slots <- .utl_theta_active_slots(
    fit$prior_variance_mode,
    V = if (fit$prior_variance_mode == "scalar") fit$V else NULL,
    tau = if (fit$prior_variance_mode == "pooled_theta") fit$tau else NULL,
    tau2 = if (fit$prior_variance_mode == "multi_theta") fit$tau2 else NULL,
    prior_tol = fit$prior_tol, L = fit$L
  )
  if (!is.logical(active_slots) || length(active_slots) != fit$L ||
      anyNA(active_slots)) {
    stop("fit active-slot metadata is invalid")
  }
  for (l in seq_len(fit$L)) {
    if (!active_slots[[l]] || max(alpha[l, ]) <= 0) {
      dropped <- rbind(dropped, data.frame(slot = paste0("L", l),
                                           reason = "inactive_prior_variance",
                                           stringsAsFactors = FALSE))
      next
    }
    ii <- order(-alpha[l, ], vn)
    cs <- ii[cumsum(alpha[l, ii]) <= coverage]
    first <- ii[[1L]]
    if (!first %in% cs) cs <- c(first, cs)
    if (sum(alpha[l, cs]) < coverage && length(cs) < p) {
      nxt <- setdiff(ii, cs)[[1L]]
      cs <- c(cs, nxt)
      while (sum(alpha[l, cs]) < coverage && length(cs) < p) {
        cs <- c(cs, setdiff(ii, cs)[[1L]])
      }
    }
    pur <- 1
    if (length(cs) > 1L) {
      for (u in seq_len(length(cs) - 1L)) {
        for (v in (u + 1L):length(cs)) {
          pair_purity <- max(vapply(seq_along(an), function(t) {
            abs(R_list[[t]][cs[[u]], cs[[v]]])
          }, numeric(1)))
          pur <- min(pur, pair_purity)
        }
      }
    }
    key <- paste(vn[sort(cs)], collapse = ",")
    if (pur < min_abs_corr) {
      dropped <- rbind(dropped, data.frame(slot = paste0("L", l),
                                           reason = "purity_below_threshold",
                                           stringsAsFactors = FALSE))
      next
    }
    if (key %in% set_key) {
      dropped <- rbind(dropped, data.frame(slot = paste0("L", l),
                                           reason = "duplicate_set",
                                           stringsAsFactors = FALSE))
      next
    }
    nr <- nr + 1L
    set_key[[nr]] <- key
    rows[[nr]] <- data.frame(cs = paste(vn[cs], collapse = ","),
                             slot = paste0("L", l),
                             size = length(cs),
                             coverage = sum(alpha[l, cs]), purity = pur,
                             stringsAsFactors = FALSE)
  }
  sets <- if (length(rows)) do.call(rbind, rows) else
    data.frame(cs = character(), slot = character(), size = integer(),
               coverage = numeric(), purity = numeric(), stringsAsFactors = FALSE)
  list(sets = sets, alpha = alpha, coverage = coverage,
       min_abs_corr = min_abs_corr, dropped = dropped)
}
