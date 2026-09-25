ar <- commandArgs(trailingOnly = TRUE)
preflight <- "--preflight" %in% ar
get_arg <- function(x) {
  i <- match(x, ar)
  if (is.na(i) || i == length(ar)) stop("missing command-line value for ", x)
  ar[[i + 1L]]
}

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
src <- file.path(repo, "tmp", "research_sampling_models", "XMAP", "data",
                 "example_data.RData")
if (!file.exists(src)) stop("XMAP example_data.RData is missing: ", src)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
load(src)
if (!exists("R1") || !exists("R2")) stop("XMAP data must contain R1 and R2")

p <- 500L
L <- 5L
V0 <- 0.00015
b <- 0.8 * sqrt(V0)
causal <- c(209L, 280L, 386L)
an <- c("EUR", "EAS", "AFR")
tn <- c("theta0", "thetaE", "thetaA")
n <- c(EUR = 1000000, EAS = 150000, AFR = 150000)
delta <- c(EAS = 1, AFR = 1)
M <- rbind(EUR = c(1, 0, 0), EAS = c(1, 1, 0), AFR = c(1, 0, 1))
colnames(M) <- tn
cases <- c("C_weak", "C_standard", "near_C", "original_mixed",
           "EUR_specific", "EAS_specific", "AFR_specific", "EUR_null",
           "EAS_null", "AFR_null")

B <- list(
  D = diag(3L),
  N_E = cbind(c(1, -1, 0), c(0, 0, 1)),
  N_A = cbind(c(1, 0, -1), c(0, 1, 0)),
  C = matrix(c(1, 0, 0), 3L, 1L),
  S_E = matrix(c(0, 1, 0), 3L, 1L),
  S_A = matrix(c(0, 0, 1), 3L, 1L)
)
for (k in names(B)) {
  rownames(B[[k]]) <- tn
  colnames(B[[k]]) <- paste0(k, "_b", seq_len(ncol(B[[k]])))
}
U <- lapply(B, tcrossprod)
for (k in names(U)) dimnames(U[[k]]) <- list(tn, tn)
u6 <- c("D", "N_E", "N_A", "C", "S_E", "S_A")
libraries <- list(
  UTL6 = list(B = B[u6], Utheta = U[u6],
              weights = setNames(rep(1 / 6, 6L), u6))
)
for (k in names(U)) {
  ev <- eigen(U[[k]], symmetric = TRUE, only.values = TRUE)$values
  et <- 100 * .Machine$double.eps * max(1, max(abs(ev)))
  if (min(ev) < -et) stop("prior covariance is not machine-PSD: ", k)
}
if (anyDuplicated(vapply(U, digest::digest, character(1), algo = "sha256"))) {
  stop("UTL6 contains a duplicated prior covariance")
}

truth_for_case <- function(case) {
  beta <- matrix(0, p, 3L, dimnames = list(sprintf("v%04d", seq_len(p)), an))
  if (case == "original_mixed") {
    beta[209L, ] <- sqrt(V0) * c(.9, .7, 0)
    beta[280L, ] <- sqrt(V0) * c(0, 0, .8)
    beta[386L, ] <- sqrt(V0) * c(.8, .8, .8)
  } else {
    z <- switch(case,
      C_weak = rep(.0065, 3L),
      C_standard = rep(b, 3L),
      near_C = b * c(1, 1.1, .9),
      EUR_specific = b * c(1, 0, 0),
      EAS_specific = b * c(0, 1, 0),
      AFR_specific = b * c(0, 0, 1),
      EUR_null = b * c(0, 1, 1),
      EAS_null = b * c(1, 0, 1),
      AFR_null = b * c(1, 1, 0),
      stop("unknown case: ", case)
    )
    beta[causal, ] <- matrix(rep(z, each = length(causal)), length(causal), 3L)
  }
  theta <- t(solve(M, t(beta)))
  dimnames(theta) <- list(rownames(beta), tn)
  list(beta = beta, theta = theta)
}

if (!is.matrix(R1) || !is.matrix(R2) || any(dim(R1) != c(p, p)) ||
    any(dim(R2) != c(p, p))) stop("XMAP R1/R2 must be 500 by 500")
R1 <- (R1 + t(R1)) / 2
R2 <- (R2 + t(R2)) / 2
diag(R1) <- diag(R2) <- 1
AR <- 0.9 ^ abs(outer(seq_len(p), seq_len(p), "-"))
R_list <- list(EUR = R2, EAS = R1, AFR = (0.7 * R2 + 0.3 * R1) * AR)
R_list$AFR <- (R_list$AFR + t(R_list$AFR)) / 2
diag(R_list$AFR) <- 1
vn <- sprintf("v%04d", seq_len(p))
F_list <- vector("list", 3L)
names(F_list) <- an
ld_psd <- vector("list", 3L)
for (t in seq_along(an)) {
  dimnames(R_list[[t]]) <- list(vn, vn)
  ee <- eigen(R_list[[t]], symmetric = TRUE)
  et <- 100 * .Machine$double.eps * max(1, max(abs(ee$values)))
  if (min(ee$values) < -et) stop("LD has a substantive negative eigenvalue: ", an[[t]])
  F_list[[t]] <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), p)
  ld_psd[[t]] <- data.frame(
    ancestry = an[[t]], min_eigen_raw = min(ee$values),
    max_eigen_raw = max(ee$values), machine_psd_tol = et,
    clipped_machine_negative = sum(ee$values < 0), stringsAsFactors = FALSE
  )
}
ld_psd <- do.call(rbind, ld_psd)

if (preflight) {
  for (case in cases) {
    tr <- truth_for_case(case)
    if (any(!is.finite(tr$beta)) || any(!is.finite(tr$theta)) ||
        !isTRUE(all.equal(unname(tr$theta %*% t(M)), unname(tr$beta),
                          tolerance = 1e-14))) {
      stop("truth transformation failed: ", case)
    }
  }
  quit(save = "no", status = 0L)
}

rep_root <- normalizePath(get_arg("--rep-root"), winslash = "/", mustWork = FALSE)
case <- get_arg("--case")
case_id <- as.integer(get_arg("--case-id"))
rep_id <- as.integer(get_arg("--rep"))
seed <- as.integer(get_arg("--seed"))
if (dir.exists(rep_root)) stop("replicate root must be absent: ", rep_root)
if (is.na(case_id) || case_id < 1L || case_id > 10L || cases[[case_id]] != case) {
  stop("case/case_id contract failed")
}
if (is.na(rep_id) || rep_id < 1L || rep_id > 100L) stop("rep must be 1..100")
if (seed != 20270000L + 100L * (case_id - 1L) + rep_id) stop("seed formula changed")

tr <- truth_for_case(case)
set.seed(seed)
z_list <- vector("list", 3L)
XtX_list <- vector("list", 3L)
Xty_list <- vector("list", 3L)
names(z_list) <- names(XtX_list) <- names(Xty_list) <- an
for (t in seq_along(an)) {
  mu <- as.vector(sqrt(n[[t]]) * R_list[[t]] %*% tr$beta[, t])
  z_list[[t]] <- setNames(mu + as.vector(F_list[[t]] %*% rnorm(p)), vn)
  XtX_list[[t]] <- n[[t]] * R_list[[t]]
  dimnames(XtX_list[[t]]) <- list(vn, vn)
  Xty_list[[t]] <- setNames(sqrt(n[[t]]) * z_list[[t]], vn)
}
yty <- setNames(as.numeric(n), an)
prior_rows <- list()
i <- 0L
for (arm in names(libraries)) for (k in names(libraries[[arm]]$Utheta)) {
  i <- i + 1L
  prior_rows[[i]] <- data.frame(
    arm = arm, component = k, rank = qr(libraries[[arm]]$B[[k]])$rank,
    B_sha256 = digest::digest(libraries[[arm]]$B[[k]], algo = "sha256"),
    U_sha256 = digest::digest(libraries[[arm]]$Utheta[[k]], algo = "sha256"),
    initial_weight = libraries[[arm]]$weights[[k]], stringsAsFactors = FALSE
  )
}
dat <- list(case = case, case_id = case_id, rep = rep_id,
            pair_id = sprintf("%s-rep%03d", case, rep_id), seed = seed,
            p = p, L = L, V0 = V0, b = b, causal = causal, n = n,
            ancestry = an, theta_names = tn, delta = delta, M = M,
            truth_beta = tr$beta, truth_theta = tr$theta,
            libraries = libraries, prior_library = do.call(rbind, prior_rows),
            R_list = R_list, z_list = z_list, XtX_list = XtX_list,
            Xty_list = Xty_list, yty = yty, ld_psd = ld_psd,
            xmap_path = normalizePath(src, winslash = "/", mustWork = TRUE),
            xmap_md5 = unname(tools::md5sum(src)))
dir.create(file.path(rep_root, "input"), recursive = TRUE)
dir.create(file.path(rep_root, "truth"))
dir.create(file.path(rep_root, "status"))
saveRDS(dat, file.path(rep_root, "input", "input.rds"))
tb <- data.frame(variant = rep(vn, 3L), variant_row = rep(seq_len(p), 3L),
                 ancestry = rep(an, each = p), truth = as.vector(tr$beta),
                 causal = as.integer(rep(seq_len(p), 3L) %in% causal))
tt <- data.frame(variant = rep(vn, 3L), variant_row = rep(seq_len(p), 3L),
                 coordinate = rep(tn, each = p), truth = as.vector(tr$theta),
                 causal = as.integer(rep(seq_len(p), 3L) %in% causal))
write.table(tb, file.path(rep_root, "truth", "truth_beta.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(tt, file.path(rep_root, "truth", "truth_theta.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(dat$prior_library, file.path(rep_root, "truth", "prior_library.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ld_psd, file.path(rep_root, "status", "ld_psd.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
ip <- file.path(rep_root, "input", "input.rds")
ih <- unname(tools::md5sum(ip))
write.table(data.frame(path = "input/input.rds", md5 = ih,
                       size = format(file.info(ip)$size, scientific = FALSE,
                                     trim = TRUE)),
            file.path(rep_root, "status", "input_hashes.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(case = case, case_id = case_id, rep = rep_id, seed = seed,
                       method = "producer", status = "success", input_md5 = ih),
            file.path(rep_root, "status", "producer_status.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
