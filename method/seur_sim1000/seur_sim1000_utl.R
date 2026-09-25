ar <- commandArgs(trailingOnly = TRUE)
get_arg <- function(x) {
  i <- match(x, ar)
  if (is.na(i) || i == length(ar)) stop("missing command-line value for ", x)
  ar[[i + 1L]]
}

rep_root <- normalizePath(get_arg("--rep-root"), winslash = "/", mustWork = TRUE)
r_lib <- normalizePath(get_arg("--r-lib"), winslash = "/", mustWork = TRUE)
arm <- get_arg("--arm")
if (arm != "UTL6") stop("arm must be UTL6; the reference-specific arm has been discontinued")
.libPaths(c(r_lib, .libPaths()))
if (!requireNamespace("SuSiEUTL", quietly = TRUE)) stop("SuSiEUTL is unavailable")
expected <- normalizePath(file.path(r_lib, "SuSiEUTL"), winslash = "/", mustWork = TRUE)
actual <- normalizePath(find.package("SuSiEUTL"), winslash = "/", mustWork = TRUE)
if (actual != expected) stop("SuSiEUTL was not loaded from the campaign runtime")

ip <- file.path(rep_root, "input", "input.rds")
dat <- readRDS(ip)
if (!identical(dat$p, 500L) || !identical(dat$L, 5L) ||
    !identical(dat$causal, c(209L, 280L, 386L)) ||
    !identical(dat$ancestry, c("EUR", "EAS", "AFR"))) {
  stop("UTL simulation input contract failed")
}
od <- file.path(rep_root, "methods", tolower(arm))
if (dir.exists(od)) stop("refusing to overwrite UTL output: ", od)
dir.create(od, recursive = TRUE)
lib <- dat$libraries[[arm]]
if (is.null(lib) || length(lib$Utheta) != 6L ||
    !identical(names(lib$Utheta), c("D", "N_E", "N_A", "C", "S_E", "S_A"))) {
  stop("UTL prior library contract failed")
}
prior <- SuSiEUTL::utl_theta_prior(dat$delta, dat$M, lib$Utheta, lib$weights)
sigma2 <- setNames(rep(1, 3L), dat$ancestry)
t0 <- proc.time()[["elapsed"]]
fit <- SuSiEUTL::utl_susie_ma_ss(
  dat$XtX_list, dat$Xty_list, dat$yty, dat$n, prior,
  L = 5L, V = 0.00015, estimate_prior_variance = TRUE,
  max_iter = 100L, tol = 1e-4, check_null_threshold = 0,
  prior_tol = 1e-9, verbose = FALSE, residual_variance = sigma2,
  estimate_residual_variance = FALSE,
  estimate_prior_mixture_weights = TRUE,
  prior_variance_mode = "multi_theta", zeroing_warmup = 5L
)
elapsed <- proc.time()[["elapsed"]] - t0
if (!isTRUE(fit$converged) || !is.finite(fit$elbo) || any(!is.finite(fit$pip)) ||
    any(fit$pip < 0 | fit$pip > 1) || fit$weight_update_count < 1L ||
    !"update_001" %in% rownames(fit$weights_trace)) {
  stop(arm, " returned invalid posterior output")
}
saveRDS(fit, file.path(od, "fit.rds"))

vn <- names(fit$pip)
rk <- rank(-fit$pip, ties.method = "min")
write.table(data.frame(variant = vn, variant_row = seq_len(dat$p),
                       causal = as.integer(seq_len(dat$p) %in% dat$causal),
                       pip = as.numeric(fit$pip), rank = rk),
            file.path(od, "pip.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)

theta_rows <- vector("list", 3L)
beta_rows <- vector("list", 3L)
bm <- fit$theta_mean %*% t(dat$M)
bv <- matrix(0, dat$p, 3L, dimnames = list(vn, dat$ancestry))
for (j in seq_len(dat$p)) {
  S <- matrix(fit$theta_covariance[j, , ], 3L, 3L)
  bv[j, ] <- diag(dat$M %*% S %*% t(dat$M))
}
for (q in seq_len(3L)) {
  theta_rows[[q]] <- data.frame(
    variant = vn, variant_row = seq_len(dat$p), coordinate = dat$theta_names[[q]],
    causal = as.integer(seq_len(dat$p) %in% dat$causal),
    truth = dat$truth_theta[, q], mean = fit$theta_mean[, q],
    variance = fit$theta_variance[, q], bias = fit$theta_mean[, q] - dat$truth_theta[, q],
    se2 = (fit$theta_mean[, q] - dat$truth_theta[, q])^2,
    stringsAsFactors = FALSE)
  beta_rows[[q]] <- data.frame(
    variant = vn, variant_row = seq_len(dat$p), ancestry = dat$ancestry[[q]],
    causal = as.integer(seq_len(dat$p) %in% dat$causal),
    truth = dat$truth_beta[, q], mean = bm[, q], variance = bv[, q],
    bias = bm[, q] - dat$truth_beta[, q], se2 = (bm[, q] - dat$truth_beta[, q])^2,
    stringsAsFactors = FALSE)
}
write.table(do.call(rbind, theta_rows), file.path(od, "theta_effects.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, beta_rows), file.path(od, "beta_effects.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

lf <- SuSiEUTL::utl_lfsr(fit)
ltheta <- vector("list", 3L)
lbeta <- vector("list", 3L)
for (q in seq_len(3L)) {
  z <- as.numeric(lf$variant$theta[, q])
  est <- sign(fit$theta_mean[, q])
  ts <- sign(dat$truth_theta[, q])
  ltheta[[q]] <- data.frame(
    variant = vn, variant_row = seq_len(dat$p), coordinate = dat$theta_names[[q]],
    causal = as.integer(seq_len(dat$p) %in% dat$causal), value = z,
    call_0.10 = as.integer(z <= .10), call_0.05 = as.integer(z <= .05),
    truth_sign = ts, correct = as.integer(z <= .10 & ts != 0 & est == ts),
    false_call = as.integer(z <= .10 & (ts == 0 | est != ts)),
    stringsAsFactors = FALSE)
  z <- as.numeric(lf$variant$beta[, q])
  est <- sign(bm[, q])
  ts <- sign(dat$truth_beta[, q])
  lbeta[[q]] <- data.frame(
    variant = vn, variant_row = seq_len(dat$p), ancestry = dat$ancestry[[q]],
    causal = as.integer(seq_len(dat$p) %in% dat$causal), value = z,
    call_0.10 = as.integer(z <= .10), call_0.05 = as.integer(z <= .05),
    truth_sign = ts, correct = as.integer(z <= .10 & ts != 0 & est == ts),
    false_call = as.integer(z <= .10 & (ts == 0 | est != ts)),
    stringsAsFactors = FALSE)
}
lt <- do.call(rbind, ltheta)
lb <- do.call(rbind, lbeta)
if (nrow(lt) != 1500L || nrow(lb) != 1500L ||
    any(!is.finite(lt$value)) || any(!is.finite(lb$value)) ||
    any(lt$value < 0 | lt$value > 1) || any(lb$value < 0 | lb$value > 1)) {
  stop(arm, " public variant lfsr contract failed")
}
write.table(lt, file.path(od, "lfsr_theta.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(lb, file.path(od, "lfsr_beta.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

cs <- SuSiEUTL:::utl_credible_sets(fit, dat$R_list, coverage = .95,
                                    min_abs_corr = .5)
write.table(cs$sets, file.path(od, "cs.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
pat <- as.data.frame(as.table(fit$pattern_posterior_by_slot),
                     stringsAsFactors = FALSE)
names(pat) <- c("slot", "pattern", "posterior")
pat$scope <- "slot"
po <- data.frame(slot = "all_active", pattern = names(fit$pattern_posterior),
                 posterior = as.numeric(fit$pattern_posterior), scope = "overall")
write.table(rbind(pat, po), file.path(od, "pattern.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
wt <- as.data.frame(as.table(fit$weights_trace), stringsAsFactors = FALSE)
names(wt) <- c("update", "pattern", "weight")
write.table(wt, file.path(od, "weights_trace.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
vr <- list()
i <- 0L
for (s in names(fit$variance_trace)) {
  x <- fit$variance_trace[[s]]
  for (q in seq_len(nrow(x))) for (l in seq_len(ncol(x))) {
    i <- i + 1L
    vr[[i]] <- data.frame(step = s, coordinate = rownames(x)[[q]],
                          slot = colnames(x)[[l]], variance = x[q, l])
  }
}
write.table(do.call(rbind, vr), file.path(od, "variance_trace.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(
  method = arm, status = "success", converged = fit$converged,
  fit_status = fit$status, iterations = nrow(fit$trace), elbo = fit$elbo,
  elapsed_sec = elapsed, input_md5 = unname(tools::md5sum(ip)),
  package_version = as.character(utils::packageVersion("SuSiEUTL")),
  V_policy = "oracle_initial_0.00015_adaptive_multi_theta",
  weight_updates = fit$weight_update_count, active_slots = sum(fit$active_slots),
  stringsAsFactors = FALSE), file.path(od, "status.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE)
