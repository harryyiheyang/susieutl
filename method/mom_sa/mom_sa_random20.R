args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (length(hit) != 1L) stop("one ", name, "=... argument is required")
  value <- sub(paste0("^", name, "="), "", hit)
  if (!nzchar(value)) stop(name, " must not be empty")
  value
}

mode <- get_arg("--mode")
if (!mode %in% c("fit", "collect")) stop("--mode must be fit or collect")
out_dir <- normalizePath(get_arg("--out"), winslash = "/", mustWork = FALSE)
file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
if (length(file_arg) != 1L) stop("run this script with Rscript so --file is available")
analysis_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
pkg_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(pkg_dir, "..", ".."), mustWork = TRUE)

rep_ids <- c(7L, 8L, 16L, 17L, 21L, 32L, 35L, 41L, 45L, 54L,
             56L, 58L, 63L, 65L, 69L, 72L, 74L, 82L, 84L, 87L)
if (length(rep_ids) != 20L || anyDuplicated(rep_ids) ||
    any(rep_ids < 6L | rep_ids > 100L)) stop("random20 replicate list is invalid")
rep_position <- match(rep_ids, rep_ids)

if (mode == "fit") {
  rep_text <- get_arg("--rep")
  if (!grepl("^[0-9]{3}$", rep_text)) stop("--rep must be a three-digit replicate")
  rep_id <- as.integer(rep_text)
  list_position <- match(rep_id, rep_ids)
  if (is.na(list_position)) stop("replicate is not in the frozen random20 list")
  if (dir.exists(file.path(out_dir, sprintf("rep%03d", rep_id)))) {
    stop("replicate output directory already exists: ",
         file.path(out_dir, sprintf("rep%03d", rep_id)))
  }
}

relative_path <- function(x, root) {
  y <- normalizePath(x, winslash = "/", mustWork = FALSE)
  r <- normalizePath(root, winslash = "/", mustWork = FALSE)
  prefix <- paste0(r, "/")
  if (startsWith(y, prefix)) substring(y, nchar(prefix) + 1L) else y
}

fit_args <- function(dat, prior) {
  residual_variance <- setNames(rep(1, length(dat$n)), names(dat$n))
  list(
    XtX_list = dat$XtX_list, Xty_list = dat$Xty_list, yty = dat$yty,
    n = dat$n, prior = prior, L = 5L, V = as.numeric(dat$V0),
    estimate_prior_variance = FALSE, max_iter = 100L, tol = 1e-4,
    verbose = FALSE, residual_variance = residual_variance,
    estimate_residual_variance = FALSE,
    estimate_prior_mixture_weights = FALSE,
    prior_variance_mode = "scalar", input_mode = "exact_ss",
    null_weight_mode = "none", top_patterns_per_snp = 0L,
    fixed_v_null_check = FALSE
  )
}

if (mode == "fit") {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
  library(SuSiEUTL)
  input_path <- file.path(root_dir, "tmp", "seur_6vs7_sim1000", "original_mixed",
                          sprintf("rep%03d", rep_id), "input", "input.rds")
  if (!file.exists(input_path)) stop("random20 input is missing: ", input_path)
  input_md5 <- unname(tools::md5sum(input_path))
  dat <- readRDS(input_path)
  if (!identical(dat$case, "original_mixed") ||
      as.integer(dat$rep) != rep_id || as.integer(dat$p) != 500L ||
      as.integer(dat$L) != 5L || !is.numeric(dat$V0) || length(dat$V0) != 1L ||
      !is.finite(dat$V0) || dat$V0 <= 0 ||
      !identical(names(dat$n), c("EUR", "EAS", "AFR")) ||
      !isTRUE(all.equal(as.numeric(dat$n), c(1000000, 150000, 150000), tolerance = 0)) ||
      !is.numeric(dat$seed) || length(dat$seed) != 1L ||
      dat$seed != 20270300L + rep_id) {
    stop("random20 input contract failed for rep", sprintf("%03d", rep_id))
  }
  components <- c("D", "N_E", "N_A", "C", "S_E", "S_A")
  lib <- dat$libraries$UTL6
  if (is.null(lib) || !identical(names(lib$Utheta), components) ||
      !identical(names(lib$weights), components) ||
      !isTRUE(all.equal(unname(lib$weights), rep(1 / 6, 6L), tolerance = 0))) {
    stop("random20 UTL6 prior contract failed for rep", sprintf("%03d", rep_id))
  }
  prior <- utl_theta_prior(dat$delta, dat$M, lib$Utheta, lib$weights)
  base <- fit_args(dat, prior)
  parameter_hash <- digest::digest(
    list(L = base$L, V = base$V, n = base$n,
         estimate_prior_variance = base$estimate_prior_variance,
         max_iter = base$max_iter, tol = base$tol,
         residual_variance = base$residual_variance,
         estimate_residual_variance = base$estimate_residual_variance,
         estimate_prior_mixture_weights = base$estimate_prior_mixture_weights,
         prior_variance_mode = base$prior_variance_mode,
         input_mode = base$input_mode, null_weight_mode = base$null_weight_mode,
         top_patterns_per_snp = base$top_patterns_per_snp,
         fixed_v_null_check = base$fixed_v_null_check), algo = "sha256"
  )
  case_dir <- file.path(out_dir, sprintf("rep%03d", rep_id))
  dir.create(case_dir, recursive = TRUE)
  arm_order <- if (list_position %% 2L == 1L) c("gaussian", "mom") else c("mom", "gaussian")
  metadata <- data.frame(
    rep = rep_id, list_position = list_position, case = dat$case,
    input_path = relative_path(input_path, root_dir), input_md5 = input_md5,
    seed = as.integer(dat$seed), p = as.integer(dat$p), L = as.integer(base$L),
    V = base$V, n_EUR = dat$n[["EUR"]], n_EAS = dat$n[["EAS"]],
    n_AFR = dat$n[["AFR"]], tol = base$tol, max_iter = base$max_iter,
    arm_order = paste(arm_order, collapse = ","), parameter_hash = parameter_hash,
    stringsAsFactors = FALSE
  )
  write.table(metadata, file.path(case_dir, "metadata.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  ledger_rows <- list()
  for (arm in arm_order) {
    arm_args <- base
    if (arm == "gaussian") {
      arm_args$estimate_pattern_weight <- FALSE
      arm_args$mom_sa_control <- NULL
    } else {
      arm_args$estimate_pattern_weight <- TRUE
    }
    start <- proc.time()[["elapsed"]]
    fit <- do.call(utl_susie_ma_ss, arm_args)
    elapsed <- unname(proc.time()[["elapsed"]] - start)
    fit_path <- file.path(case_dir, paste0(arm, "_fit.rds"))
    runtime_path <- file.path(case_dir, paste0(arm, "_runtime.rds"))
    saveRDS(fit, fit_path)
    saveRDS(list(runtime_seconds = elapsed, fit_call_only = TRUE,
                arm = arm, rep = rep_id), runtime_path)
    fit_status <- as.character(fit$status)
    iterations <- if (is.data.frame(fit$trace)) nrow(fit$trace) else NA_integer_
    converged <- isTRUE(fit$converged)
    weight_model <- if (!is.null(fit$weight_model)) as.character(fit$weight_model) else "legacy"
    ledger_rows[[length(ledger_rows) + 1L]] <- data.frame(
      rep = rep_id, list_position = list_position, arm = arm,
      arm_order = paste(arm_order, collapse = ","), input_md5 = input_md5,
      seed = as.integer(dat$seed), p = as.integer(dat$p),
      n_EUR = dat$n[["EUR"]], n_EAS = dat$n[["EAS"]], n_AFR = dat$n[["AFR"]],
      V = base$V, L = base$L, max_iter = base$max_iter, tol = base$tol,
      status = "success", fit_status = fit_status,
      error = "", fit_path = relative_path(fit_path, out_dir),
      runtime_path = relative_path(runtime_path, out_dir),
      runtime_seconds = elapsed, iterations = iterations, converged = converged,
      weight_model = weight_model, parameter_hash = parameter_hash,
      stringsAsFactors = FALSE
    )
    ledger <- do.call(rbind, ledger_rows)
    write.table(ledger, file.path(case_dir, "ledger.tsv"), sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  cat("mom_sa_random20 fit rep", sprintf("%03d", rep_id), " completed\n", sep = "")
  quit(save = "no", status = 0L)
}

if (!dir.exists(out_dir)) stop("random20 output directory is missing for collect: ", out_dir)
library(SuSiEUTL)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")

make_missing_ledger <- function(rep_id, list_position, arm, arm_order) {
  data.frame(
    rep = rep_id, list_position = list_position, arm = arm,
    arm_order = paste(arm_order, collapse = ","), input_md5 = NA_character_,
    seed = NA_integer_, p = NA_integer_, n_EUR = NA_real_, n_EAS = NA_real_,
    n_AFR = NA_real_, V = NA_real_, L = NA_integer_, max_iter = NA_integer_,
    tol = NA_real_, status = "missing", fit_status = "", error = "missing ledger",
    fit_path = "", runtime_path = "", runtime_seconds = NA_real_,
    iterations = NA_integer_, converged = NA, weight_model = "legacy",
    parameter_hash = NA_character_, stringsAsFactors = FALSE
  )
}

format_values <- function(x) {
  if (!length(x)) return(NA_character_)
  paste(vapply(x, function(z) if (is.na(z)) "NA" else sprintf("%.17g", z), character(1)),
        collapse = ";")
}

collect_credible_sets <- function(fit, dat) {
  cs <- SuSiEUTL:::utl_credible_sets(fit, dat$R_list, coverage = .99,
                                      min_abs_corr = .5)
  sets <- cs$sets
  members <- if (nrow(sets)) strsplit(sets$cs, ",", fixed = TRUE) else list()
  member_union <- if (length(members)) unique(unlist(members, use.names = FALSE)) else character()
  list(sets = sets, members = members, member_union = member_union)
}

owner_scalar <- function(x) {
  if (!length(x)) return(NA_real_)
  if (length(x) == 1L || length(unique(x)) == 1L) as.numeric(x[[1L]]) else NA_real_
}

reference_rows <- function(fit, dat, rep_id, arm, ledger_row) {
  rho <- fit$slot_rho
  if (!is.array(rho) || length(dim(rho)) != 3L ||
      any(!is.finite(rho))) stop("slot_rho is invalid for rep", sprintf("%03d", rep_id), " ", arm)
  alpha <- apply(rho, c(1L, 2L), sum)
  vn <- dimnames(rho)[[2L]]
  cn <- dimnames(rho)[[3L]]
  sn <- dimnames(rho)[[1L]]
  if (is.null(vn) || is.null(cn) || is.null(sn)) stop("slot_rho has no required dimnames")
  pip <- as.numeric(fit$pip)
  if (length(pip) != length(vn) || is.null(names(fit$pip)) ||
      !identical(names(fit$pip), vn) || any(!is.finite(pip))) {
    stop("PIP is invalid for reference collection")
  }
  cs_info <- collect_credible_sets(fit, dat)
  refs <- data.frame(
    reference = c("v0209", "v0280", "v0386"),
    variant_index = c(209L, 280L, 386L),
    expected_component = c("N_A", "S_A", "C"), stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(refs))
  for (i in seq_len(nrow(refs))) {
    j <- refs$variant_index[[i]]
    k <- match(refs$expected_component[[i]], cn)
    ne <- match("N_E", cn)
    d <- match("D", cn)
    if (is.na(k) || is.na(ne) || is.na(d) || is.na(vn[[j]])) {
      stop("reference dimensions are incompatible")
    }
    a <- alpha[, j]
    owner_value <- max(a)
    owners <- which(a == owner_value)
    joint <- rho[, j, k]
    conditional <- ifelse(a > 0, joint / a, NA_real_)
    nd_joint <- rho[, j, ne] + rho[, j, d]
    nd_conditional <- ifelse(a > 0, nd_joint / a, NA_real_)
    use_nd <- identical(refs$expected_component[[i]], "S_A")
    owner_pattern_posterior <- vapply(owners, function(l) sum(rho[l, , k]), numeric(1))
    nd_pattern_posterior <- vapply(owners, function(l) sum(rho[l, , ne] + rho[l, , d]), numeric(1))
    rows[[i]] <- data.frame(
      rep = rep_id, arm = arm, reference = refs$reference[[i]],
      variant_index = j, variant = vn[[j]], variant_pip = pip[[j]],
      expected_component = refs$expected_component[[i]],
      reference_cs_covered = vn[[j]] %in% cs_info$member_union,
      owner_argmax_slots = paste(sn[owners], collapse = ";"),
      owner_argmax_count = length(owners),
      owner_slot_alpha = owner_scalar(a[owners]),
      owner_slot_alpha_values = format_values(a[owners]),
      owner_slot_conditional_pattern = owner_scalar(conditional[owners]),
      owner_slot_conditional_pattern_values = format_values(conditional[owners]),
      owner_slot_joint_pattern = owner_scalar(joint[owners]),
      owner_slot_joint_pattern_values = format_values(joint[owners]),
      owner_slot_pattern_posterior = owner_scalar(owner_pattern_posterior),
      owner_slot_pattern_posterior_values = format_values(owner_pattern_posterior),
      all_slot_reference_component_joint = sum(joint),
      ne_d_owner_slot_conditional_pattern = if (use_nd) owner_scalar(nd_conditional[owners]) else NA_real_,
      ne_d_owner_slot_conditional_pattern_values = if (use_nd) format_values(nd_conditional[owners]) else NA_character_,
      ne_d_owner_slot_joint_pattern = if (use_nd) owner_scalar(nd_joint[owners]) else NA_real_,
      ne_d_owner_slot_joint_pattern_values = if (use_nd) format_values(nd_joint[owners]) else NA_character_,
      ne_d_owner_slot_pattern_posterior = if (use_nd) owner_scalar(nd_pattern_posterior) else NA_real_,
      ne_d_owner_slot_pattern_posterior_values = if (use_nd) format_values(nd_pattern_posterior) else NA_character_,
      status = ledger_row$status, fit_status = ledger_row$fit_status,
      input_md5 = ledger_row$input_md5, parameter_hash = ledger_row$parameter_hash,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

case_row <- function(fit, dat, rep_id, arm, ledger_row) {
  pip <- as.numeric(fit$pip)
  vn <- names(fit$pip)
  causal <- paste0("v", sprintf("%04d", as.integer(dat$causal)))
  is_causal <- vn %in% causal
  if (length(pip) != length(vn) || any(!is.finite(pip))) stop("PIP is invalid")
  truth <- as.numeric(is_causal)
  cs_info <- collect_credible_sets(fit, dat)
  sets <- cs_info$sets
  members <- cs_info$members
  member_union <- cs_info$member_union
  causal_covered <- sum(causal %in% member_union)
  noncausal_in_set <- sum(setdiff(member_union, causal) %in% vn)
  set_has_causal <- if (nrow(sets)) vapply(members, function(x) any(x %in% causal), logical(1)) else logical()
  set_has_noncausal <- if (nrow(sets)) vapply(members, function(x) any(x %in% setdiff(vn, causal)), logical(1)) else logical()
  fit_weight_model <- if (is.null(fit$weight_model)) "legacy" else as.character(fit$weight_model)
  beta_hat <- fit$theta_mean %*% t(dat$M)
  beta_truth <- dat$truth_beta
  if (!is.matrix(beta_hat) || !is.matrix(beta_truth) ||
      is.null(rownames(beta_hat)) || is.null(colnames(beta_hat)) ||
      is.null(rownames(beta_truth)) || is.null(colnames(beta_truth)) ||
      !setequal(rownames(beta_hat), rownames(beta_truth)) ||
      !setequal(colnames(beta_hat), colnames(beta_truth))) {
    stop("beta estimate and truth dimensions are not alignable")
  }
  beta_hat <- beta_hat[rownames(beta_truth), colnames(beta_truth), drop = FALSE]
  beta_error <- beta_hat - beta_truth
  causal_index <- as.integer(dat$causal)
  if (any(!is.finite(causal_index)) || any(causal_index < 1L) ||
      any(causal_index > nrow(beta_error))) stop("causal indices are invalid for beta RMSE")
  beta_causal_error <- beta_error[causal_index, , drop = FALSE]
  beta_noncausal_error <- beta_error[-causal_index, , drop = FALSE]
  data.frame(
    rep = rep_id, arm = arm, status = ledger_row$status,
    fit_status = ledger_row$fit_status, converged = ledger_row$converged,
    weight_model = fit_weight_model, input_md5 = ledger_row$input_md5,
    parameter_hash = ledger_row$parameter_hash,
    pip_mean_causal = mean(pip[is_causal]), pip_mean_noncausal = mean(pip[!is_causal]),
    pip_max_noncausal = max(pip[!is_causal]), pip_causal_min = min(pip[is_causal]),
    pip_rmse_all = sqrt(mean((pip - truth)^2)),
    pip_rmse_causal = sqrt(mean((pip[is_causal] - 1)^2)),
    pip_rmse_noncausal = sqrt(mean(pip[!is_causal]^2)),
    beta_rmse_all = sqrt(mean(beta_error^2)),
    beta_rmse_causal = sqrt(mean(beta_causal_error^2)),
    beta_rmse_noncausal = sqrt(mean(beta_noncausal_error^2)),
    cs_count_all = nrow(sets), cs_size_mean_all = if (nrow(sets)) mean(sets$size) else NA_real_,
    cs_coverage_mean_all = if (nrow(sets)) mean(sets$coverage) else NA_real_,
    cs_purity_min_all = if (nrow(sets)) min(sets$purity) else NA_real_,
    cs_causal_covered_count = causal_covered, cs_causal_coverage = causal_covered / length(causal),
    cs_noncausal_in_set_count = noncausal_in_set,
    cs_causal_set_count = sum(set_has_causal), cs_noncausal_set_count = sum(set_has_noncausal),
    cs_causal_size_mean = if (any(set_has_causal)) mean(sets$size[set_has_causal]) else NA_real_,
    cs_noncausal_size_mean = if (any(set_has_noncausal)) mean(sets$size[set_has_noncausal]) else NA_real_,
    runtime_seconds = ledger_row$runtime_seconds, iterations = ledger_row$iterations,
    stringsAsFactors = FALSE
  )
}

ledger_all <- list()
reference_all <- list()
case_all <- list()
for (list_position in seq_along(rep_ids)) {
  rep_id <- rep_ids[[list_position]]
  arm_order <- if (list_position %% 2L == 1L) c("gaussian", "mom") else c("mom", "gaussian")
  case_dir <- file.path(out_dir, sprintf("rep%03d", rep_id))
  ledger_path <- file.path(case_dir, "ledger.tsv")
  led <- if (file.exists(ledger_path)) read.delim(ledger_path, stringsAsFactors = FALSE,
                                                  check.names = FALSE) else NULL
  for (arm in c("gaussian", "mom")) {
    row <- if (!is.null(led) && any(led$arm == arm)) led[which(led$arm == arm)[[1L]], , drop = FALSE] else
      make_missing_ledger(rep_id, list_position, arm, arm_order)
    if (nrow(row) != 1L) stop("ledger has invalid row count for rep", sprintf("%03d", rep_id), " ", arm)
    ledger_all[[length(ledger_all) + 1L]] <- row
    fit_path <- if (nzchar(row$fit_path)) file.path(out_dir, row$fit_path) else ""
    if (identical(row$status, "success") && file.exists(fit_path)) {
      input_path <- file.path(root_dir, "tmp", "seur_6vs7_sim1000", "original_mixed",
                              sprintf("rep%03d", rep_id), "input", "input.rds")
      if (!file.exists(input_path)) stop("input missing while collecting rep", sprintf("%03d", rep_id))
      dat <- readRDS(input_path)
      fit <- readRDS(fit_path)
      reference_all[[length(reference_all) + 1L]] <- reference_rows(fit, dat, rep_id, arm, row)
      case_all[[length(case_all) + 1L]] <- case_row(fit, dat, rep_id, arm, row)
    } else {
      reference_all[[length(reference_all) + 1L]] <- data.frame(
        rep = rep_id, arm = arm,
        reference = c("v0209", "v0280", "v0386"),
        variant_index = c(209L, 280L, 386L),
        variant = c("v0209", "v0280", "v0386"),
        variant_pip = NA_real_,
        expected_component = c("N_A", "S_A", "C"),
        reference_cs_covered = NA,
        owner_argmax_slots = NA_character_, owner_argmax_count = NA_integer_,
        owner_slot_alpha = NA_real_, owner_slot_alpha_values = NA_character_,
        owner_slot_conditional_pattern = NA_real_, owner_slot_conditional_pattern_values = NA_character_,
        owner_slot_joint_pattern = NA_real_, owner_slot_joint_pattern_values = NA_character_,
        owner_slot_pattern_posterior = NA_real_, owner_slot_pattern_posterior_values = NA_character_,
        all_slot_reference_component_joint = NA_real_,
        ne_d_owner_slot_conditional_pattern = NA_real_,
        ne_d_owner_slot_conditional_pattern_values = NA_character_,
        ne_d_owner_slot_joint_pattern = NA_real_,
        ne_d_owner_slot_joint_pattern_values = NA_character_,
        ne_d_owner_slot_pattern_posterior = NA_real_,
        ne_d_owner_slot_pattern_posterior_values = NA_character_,
        status = row$status, fit_status = row$fit_status,
        input_md5 = row$input_md5, parameter_hash = row$parameter_hash,
        stringsAsFactors = FALSE
      )
      case_all[[length(case_all) + 1L]] <- data.frame(
        rep = rep_id, arm = arm, status = row$status, fit_status = row$fit_status,
        converged = row$converged, weight_model = row$weight_model,
        input_md5 = row$input_md5, parameter_hash = row$parameter_hash,
         pip_mean_causal = NA_real_, pip_mean_noncausal = NA_real_, pip_max_noncausal = NA_real_,
         pip_causal_min = NA_real_, pip_rmse_all = NA_real_, pip_rmse_causal = NA_real_,
         pip_rmse_noncausal = NA_real_, beta_rmse_all = NA_real_, beta_rmse_causal = NA_real_,
         beta_rmse_noncausal = NA_real_, cs_count_all = NA_real_, cs_size_mean_all = NA_real_,
        cs_coverage_mean_all = NA_real_, cs_purity_min_all = NA_real_,
        cs_causal_covered_count = NA_real_, cs_causal_coverage = NA_real_,
        cs_noncausal_in_set_count = NA_real_, cs_causal_set_count = NA_real_,
        cs_noncausal_set_count = NA_real_, cs_causal_size_mean = NA_real_,
        cs_noncausal_size_mean = NA_real_, runtime_seconds = row$runtime_seconds,
        iterations = row$iterations, stringsAsFactors = FALSE
      )
    }
  }
}

ledger_table <- do.call(rbind, ledger_all)
reference_table <- do.call(rbind, reference_all)
case_table <- do.call(rbind, case_all)
write.table(ledger_table, file.path(out_dir, "ledger.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(reference_table, file.path(out_dir, "reference_metrics.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(case_table, file.path(out_dir, "case_summary.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

metric_names <- c(
  "pip_mean_causal", "pip_mean_noncausal", "pip_max_noncausal",
  "pip_causal_min", "pip_rmse_all", "pip_rmse_causal", "pip_rmse_noncausal",
  "beta_rmse_all", "beta_rmse_causal", "beta_rmse_noncausal",
  "cs_count_all", "cs_size_mean_all", "cs_coverage_mean_all", "cs_purity_min_all",
  "cs_causal_covered_count", "cs_causal_coverage", "cs_noncausal_in_set_count",
  "cs_causal_set_count", "cs_noncausal_set_count", "cs_causal_size_mean",
  "cs_noncausal_size_mean", "runtime_seconds", "iterations"
)
higher_is_better <- setNames(c(
  TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
  TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE,
  FALSE, FALSE, FALSE
), metric_names)
paired_rows <- vector("list", length(metric_names))
for (i in seq_along(metric_names)) {
  metric <- metric_names[[i]]
  g <- case_table[case_table$arm == "gaussian", c("rep", metric), drop = FALSE]
  m <- case_table[case_table$arm == "mom", c("rep", metric), drop = FALSE]
  names(g)[[2L]] <- "gaussian"
  names(m)[[2L]] <- "mom"
  pair <- merge(g, m, by = "rep", all = TRUE, sort = TRUE)
  valid <- is.finite(pair$gaussian) & is.finite(pair$mom)
  diff <- pair$mom[valid] - pair$gaussian[valid]
  n_pair <- length(diff)
  p_value <- if (n_pair >= 2L) stats::wilcox.test(diff, mu = 0, exact = FALSE)$p.value else NA_real_
  win <- if (higher_is_better[[metric]]) diff > 0 else diff < 0
  loss <- if (higher_is_better[[metric]]) diff < 0 else diff > 0
  tie <- diff == 0
  paired_rows[[i]] <- data.frame(
    metric = metric, higher_is_better = higher_is_better[[metric]],
    valid_paired_n = n_pair,
    gaussian_failures = sum(!is.finite(g[["gaussian"]])),
    mom_failures = sum(!is.finite(m[["mom"]])),
    mom_minus_gaussian_median = if (n_pair) median(diff) else NA_real_,
    mom_minus_gaussian_iqr = if (n_pair) IQR(diff) else NA_real_,
    mom_minus_gaussian_range_min = if (n_pair) min(diff) else NA_real_,
    mom_minus_gaussian_range_max = if (n_pair) max(diff) else NA_real_,
    wilcoxon_p = p_value, mom_win_n = sum(win), mom_tie_n = sum(tie),
    mom_loss_n = sum(loss), mom_win_tie_loss = paste(sum(win), sum(tie), sum(loss), sep = "/"),
    stringsAsFactors = FALSE
  )
}
paired_table <- do.call(rbind, paired_rows)
write.table(paired_table, file.path(out_dir, "paired_summary.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("mom_sa_random20 collect completed\n")
