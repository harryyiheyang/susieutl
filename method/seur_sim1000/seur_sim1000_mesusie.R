ar <- commandArgs(trailingOnly = TRUE)
get_arg <- function(x) {
  i <- match(x, ar)
  if (is.na(i) || i == length(ar)) stop("missing command-line value for ", x)
  ar[[i + 1L]]
}

root <- normalizePath(get_arg("--rep-root"), winslash = "/", mustWork = TRUE)
lib <- normalizePath(get_arg("--r-lib"), winslash = "/", mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
if (!requireNamespace("MESuSiE", quietly = TRUE) ||
    as.character(utils::packageVersion("MESuSiE")) != "1.0") {
  stop("MESuSiE 1.0 is unavailable")
}
dat <- readRDS(file.path(root, "input", "input.rds"))
od <- file.path(root, "methods", "mesusie")
if (dir.exists(od)) stop("refusing to overwrite MESuSiE output: ", od)
dir.create(od, recursive = TRUE)
an <- c("EUR", "EAS", "AFR")
vn <- sprintf("v%04d", seq_len(500L))
ss <- vector("list", 3L)
names(ss) <- an
for (a in an) {
  ss[[a]] <- data.frame(
    SNP = vn, Beta = dat$z_list[[a]] / sqrt(dat$n[[a]]),
    Se = rep(1 / sqrt(dat$n[[a]]), 500L), Z = dat$z_list[[a]],
    N = rep(dat$n[[a]], 500L), REF = "A", ALT = "G")
}
t0 <- proc.time()[["elapsed"]]
fit <- MESuSiE::meSuSie_core(dat$R_list, ss, L = 5L)
elapsed <- proc.time()[["elapsed"]] - t0
pip <- as.numeric(fit$pip)
if (length(pip) != 500L || any(!is.finite(pip)) || any(pip < 0 | pip > 1)) {
  stop("MESuSiE native PIP is invalid")
}
rk <- rank(-pip, ties.method = "min")
write.table(data.frame(variant = vn, variant_row = seq_len(500L),
                       causal = as.integer(seq_len(500L) %in% dat$causal),
                       pip = pip, rank = rk), file.path(od, "pip.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
if (!is.null(fit$pip_config)) {
  pc <- as.matrix(fit$pip_config)
  if (nrow(pc) != 500L || any(!is.finite(pc)) || any(pc < 0 | pc > 1)) {
    stop("MESuSiE native configuration posterior is invalid")
  }
  colnames(pc) <- fit$name_config
  write.table(data.frame(variant = vn, pc, check.names = FALSE),
              file.path(od, "configuration_posterior.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  write.table(data.frame(metric = "configuration_posterior",
                         interpretation = "configuration_not_false_sign"),
              file.path(od, "configuration_status.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
}
saveRDS(fit$cs, file.path(od, "cs.rds"))
cs_tab <- data.frame(cs_id = character(), variants = character(), size = integer(),
                     category = character(), stringsAsFactors = FALSE)
if (length(fit$cs$cs)) {
  cs_tab <- data.frame(
    cs_id = names(fit$cs$cs),
    variants = vapply(fit$cs$cs, function(x) paste(vn[x], collapse = ","),
                      character(1)),
    size = vapply(fit$cs$cs, length, integer(1)),
    category = unname(fit$cs$cs_category[names(fit$cs$cs)]))
}
write.table(cs_tab, file.path(od, "cs.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)

moment_ok <- is.list(fit$EB1) && is.list(fit$EB2) && length(fit$EB1) == 5L &&
  length(fit$EB2) == 5L &&
  all(vapply(fit$EB1, function(x) is.matrix(x) && all(dim(x) == c(500L, 3L)),
             logical(1))) &&
  all(vapply(fit$EB2, function(x) is.matrix(x) && all(dim(x) == c(500L, 3L)),
             logical(1)))
metric_status <- "unavailable_native_moment_contract"
if (moment_ok) {
  m1 <- simplify2array(fit$EB1)
  m2 <- simplify2array(fit$EB2)
  vv <- m2 - m1^2
  vt <- 100 * .Machine$double.eps * max(1, max(abs(m2)), max(abs(m1^2)))
  moment_ok <- all(is.finite(m1)) && all(is.finite(m2)) && min(vv) >= -vt
}
if (moment_ok) {
  mu <- apply(m1, c(1L, 2L), sum)
  va <- apply(vv, c(1L, 2L), sum)
  va[va < 0] <- 0
  metric_status <- "available_native_posterior_moments"
  rows <- list()
  lr <- list()
  for (t in seq_len(3L)) {
    truth <- dat$truth_beta[, t]
    sd <- sqrt(va[, t])
    lf <- rep(1, 500L)
    nz <- sd > 0
    lf[nz] <- stats::pnorm(-abs(mu[nz, t] / sd[nz]))
    lf[!nz & mu[, t] != 0] <- 0
    rows[[t]] <- data.frame(variant = vn, variant_row = seq_len(500L),
                            ancestry = an[[t]], causal = as.integer(seq_len(500L) %in% dat$causal),
                            truth = truth, mean = mu[, t], variance = va[, t],
                            bias = mu[, t] - truth, se2 = (mu[, t] - truth)^2,
                            metric_status = metric_status)
    lr[[t]] <- data.frame(variant = vn, variant_row = seq_len(500L), ancestry = an[[t]],
                          causal = as.integer(seq_len(500L) %in% dat$causal), value = lf,
                          call_0.10 = as.integer(lf <= .10), call_0.05 = as.integer(lf <= .05),
                          truth_sign = sign(truth),
                          correct = as.integer(lf <= .10 & truth != 0 & sign(mu[, t]) == sign(truth)),
                          false_call = as.integer(lf <= .10 & (truth == 0 | sign(mu[, t]) != sign(truth))),
                          metric_status = metric_status)
  }
  write.table(do.call(rbind, rows), file.path(od, "effects.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  write.table(do.call(rbind, lr), file.path(od, "lfsr.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
} else {
  write.table(data.frame(metric_status = metric_status), file.path(od, "effects.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(data.frame(metric_status = metric_status), file.path(od, "lfsr.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
saveRDS(fit, file.path(od, "fit.rds"))
write.table(data.frame(method = "mesusie", status = "success", elapsed_sec = elapsed,
                       input_md5 = unname(tools::md5sum(file.path(root, "input", "input.rds"))),
                       version = as.character(utils::packageVersion("MESuSiE")),
                       commit = "5a3dd323e5d07dbd190952937f29b7a0bcf5eee0",
                       V_policy = "native_auto_default", metric_status = metric_status,
                       call = "MESuSiE::meSuSie_core(R_list,summary_stat_list,L=5)"),
            file.path(od, "status.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
