ar <- commandArgs(trailingOnly = TRUE)
get_arg <- function(x) {
  i <- match(x, ar)
  if (is.na(i) || i == length(ar)) stop("missing command-line value for ", x)
  ar[[i + 1L]]
}

root <- normalizePath(get_arg("--rep-root"), winslash = "/", mustWork = TRUE)
lib <- normalizePath(get_arg("--r-lib"), winslash = "/", mustWork = TRUE)
a <- get_arg("--ancestry")
if (!a %in% c("EUR", "EAS", "AFR")) stop("invalid ancestry")
.libPaths(c(lib, .libPaths()))
if (!requireNamespace("susieR", quietly = TRUE) ||
    as.character(utils::packageVersion("susieR")) != "0.16.5") {
  stop("susieR 0.16.5 is unavailable")
}
dat <- readRDS(file.path(root, "input", "input.rds"))
od <- file.path(root, "methods", paste0("susie_", a))
if (dir.exists(od)) stop("refusing to overwrite single-SuSiE output: ", od)
dir.create(od, recursive = TRUE)
t0 <- proc.time()[["elapsed"]]
fit <- susieR::susie_rss(z = dat$z_list[[a]], R = dat$R_list[[a]],
                         n = dat$n[[a]], L = 5L)
elapsed <- proc.time()[["elapsed"]] - t0
pip <- as.numeric(susieR::susie_get_pip(fit))
if (length(pip) != 500L || any(!is.finite(pip)) || any(pip < 0 | pip > 1)) {
  stop("single-SuSiE native PIP is invalid")
}
vn <- sprintf("v%04d", seq_len(500L))
rk <- rank(-pip, ties.method = "min")
mu <- as.numeric(susieR::susie_get_posterior_mean(fit))
sd <- as.numeric(susieR::susie_get_posterior_sd(fit))
if (length(mu) != 500L || length(sd) != 500L || any(!is.finite(mu)) ||
    any(!is.finite(sd)) || any(sd < 0)) stop("single-SuSiE native moments are invalid")
lf <- rep(1, 500L)
nz <- sd > 0
lf[nz] <- stats::pnorm(-abs(mu[nz] / sd[nz]))
lf[!nz & mu != 0] <- 0
truth <- dat$truth_beta[, a]
write.table(data.frame(variant = vn, variant_row = seq_len(500L),
                       causal = as.integer(seq_len(500L) %in% dat$causal),
                       pip = pip, rank = rk), file.path(od, "pip.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(variant = vn, variant_row = seq_len(500L), ancestry = a,
                       causal = as.integer(seq_len(500L) %in% dat$causal),
                       truth = truth, mean = mu, variance = sd^2,
                       bias = mu - truth, se2 = (mu - truth)^2,
                       metric_status = "available_native_posterior_moments"),
            file.path(od, "effects.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
write.table(data.frame(variant = vn, variant_row = seq_len(500L), ancestry = a,
                       causal = as.integer(seq_len(500L) %in% dat$causal), value = lf,
                       call_0.10 = as.integer(lf <= .10),
                       call_0.05 = as.integer(lf <= .05), truth_sign = sign(truth),
                       correct = as.integer(lf <= .10 & truth != 0 & sign(mu) == sign(truth)),
                       false_call = as.integer(lf <= .10 & (truth == 0 | sign(mu) != sign(truth))),
                       metric_status = "available_native_posterior_moments"),
            file.path(od, "lfsr.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
cs <- susieR::susie_get_cs(fit, Xcorr = dat$R_list[[a]])
saveRDS(cs, file.path(od, "cs.rds"))
cs_tab <- data.frame(cs_id = character(), variants = character(), size = integer(),
                     stringsAsFactors = FALSE)
if (length(cs$cs)) {
  cs_tab <- data.frame(cs_id = names(cs$cs),
                       variants = vapply(cs$cs, function(x) paste(vn[x], collapse = ","),
                                         character(1)),
                       size = vapply(cs$cs, length, integer(1)))
}
write.table(cs_tab, file.path(od, "cs.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
saveRDS(fit, file.path(od, "fit.rds"))
write.table(data.frame(method = paste0("susie_", a), status = "success",
                       elapsed_sec = elapsed,
                       input_md5 = unname(tools::md5sum(file.path(root, "input", "input.rds"))),
                       version = as.character(utils::packageVersion("susieR")),
                       commit = "69f6f38d70467f4a57a74200fbe98f71de8a70d1",
                       V_policy = "native_auto_default",
                       call = "susieR::susie_rss(z,R,n,L=5)"),
            file.path(od, "status.tsv"), sep = "\t", row.names = FALSE,
            quote = FALSE)
