args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  hit <- args[startsWith(args, paste0(name, "="))]
  if (length(hit) != 1L) stop("one ", name, "=... argument is required")
  value <- sub(paste0("^", name, "="), "", hit)
  if (!nzchar(value)) stop(name, " must not be empty")
  value
}

cases_dir <- normalizePath(get_arg("--cases"), winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(get_arg("--out"), winslash = "/", mustWork = FALSE)
file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
if (length(file_arg) != 1L) stop("run this script with Rscript so --file is available")
analysis_dir <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
pkg_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(pkg_dir, "..", ".."), mustWork = TRUE)
if (dir.exists(out_dir)) stop("LFSR output directory already exists: ", out_dir)

rep_ids <- c(7L, 8L, 16L, 17L, 21L, 32L, 35L, 41L, 45L, 54L,
             56L, 58L, 63L, 65L, 69L, 72L, 74L, 82L, 84L, 87L)
if (length(rep_ids) != 20L || anyDuplicated(rep_ids) ||
    any(rep_ids < 6L | rep_ids > 100L)) stop("LFSR random20 replicate list is invalid")
dir.create(out_dir, recursive = TRUE)
library(SuSiEUTL)
if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")

refs <- c("v0209", "v0280", "v0386")
ref_indices <- c(209L, 280L, 386L)
arms <- c("gaussian", "mom")
format_values <- function(x) {
  if (!length(x)) return(NA_character_)
  paste(vapply(x, function(z) if (is.na(z)) "NA" else sprintf("%.17g", z), character(1)),
        collapse = ";")
}
scalar_value <- function(x) {
  if (!length(x)) return(NA_real_)
  if (length(x) == 1L || length(unique(x)) == 1L) as.numeric(x[[1L]]) else NA_real_
}
append_missing_owner <- function(dat, rep_id, arm, status, fit_status, input_md5) {
  vn <- rownames(dat$truth_beta)
  bn <- colnames(dat$truth_beta)
  rows <- vector("list", length(refs) * length(bn))
  z <- 0L
  for (i in seq_along(refs)) {
    j <- ref_indices[[i]]
    for (a in seq_along(bn)) {
      z <- z + 1L
      rows[[z]] <- data.frame(
        rep = rep_id, arm = arm, reference = refs[[i]],
        variant = vn[[j]], ancestry = bn[[a]], truth_beta = dat$truth_beta[j, a],
        truth_zero = isTRUE(dat$truth_beta[j, a] == 0),
        owner_slots = NA_character_, owner_tie_count = NA_integer_,
        owner_alpha_max = NA_real_, owner_conditional_lfsr = NA_real_,
        owner_conditional_lfsr_values = NA_character_, status = status,
        fit_status = fit_status, input_md5 = input_md5, fit_sha256 = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
append_missing_variant <- function(dat, rep_id, arm, status, fit_status, input_md5) {
  vn <- rownames(dat$truth_beta)
  coords <- list(theta = colnames(dat$M), beta = colnames(dat$truth_beta))
  rows <- vector("list", length(vn) * (length(coords$theta) + length(coords$beta)))
  z <- 0L
  for (sp in names(coords)) {
    for (a in seq_along(coords[[sp]])) {
      for (j in seq_along(vn)) {
        z <- z + 1L
        truth_zero <- if (sp == "beta") isTRUE(dat$truth_beta[j, coords[[sp]][[a]]] == 0) else NA
        rows[[z]] <- data.frame(
          rep = rep_id, arm = arm, space = sp, coordinate = coords[[sp]][[a]],
          variant = vn[[j]], variant_index = j, causal = j %in% as.integer(dat$causal),
          truth_zero = truth_zero, variant_lfsr = NA_real_, status = status,
          fit_status = fit_status, input_md5 = input_md5, fit_sha256 = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

owner_rows <- list()
variant_rows <- list()
metadata_rows <- list()
gate_rows <- list()
for (rep_id in rep_ids) {
  input_path <- file.path(root_dir, "tmp", "seur_6vs7_sim1000", "original_mixed",
                          sprintf("rep%03d", rep_id), "input", "input.rds")
  if (!file.exists(input_path)) stop("random20 input is missing: ", input_path)
  input_md5 <- unname(tools::md5sum(input_path))
  dat <- readRDS(input_path)
  if (!is.matrix(dat$truth_beta) || is.null(rownames(dat$truth_beta)) ||
      is.null(colnames(dat$truth_beta)) || !is.matrix(dat$M) ||
      !identical(colnames(dat$truth_beta), rownames(dat$M)) ||
      is.null(colnames(dat$M)) ||
      !identical(as.integer(dat$causal), ref_indices)) {
    stop("random20 LFSR input contract failed for rep", sprintf("%03d", rep_id))
  }
  case_dir <- file.path(cases_dir, sprintf("rep%03d", rep_id))
  case_out <- file.path(out_dir, sprintf("rep%03d", rep_id))
  dir.create(case_out, recursive = TRUE)
  ledger_path <- file.path(case_dir, "ledger.tsv")
  ledger <- if (file.exists(ledger_path)) read.delim(
    ledger_path, stringsAsFactors = FALSE, check.names = FALSE
  ) else NULL
  for (arm in arms) {
    fit_path <- file.path(case_dir, paste0(arm, "_fit.rds"))
    ledger_row <- if (!is.null(ledger) && any(ledger$arm == arm)) {
      ledger[which(ledger$arm == arm)[[1L]], , drop = FALSE]
    } else NULL
    fit_ok <- !is.null(ledger_row) && identical(as.character(ledger_row$status), "success") &&
      file.exists(fit_path)
    status <- if (fit_ok) "success" else "missing"
    fit_status <- if (fit_ok) as.character(ledger_row$fit_status) else
      if (is.null(ledger_row)) "" else as.character(ledger_row$fit_status)
    fit_sha256 <- NA_character_
    lfsr <- NULL
    if (fit_ok) {
      fit_sha256 <- digest::digest(fit_path, algo = "sha256", file = TRUE)
      fit <- readRDS(fit_path)
      lfsr <- utl_lfsr(fit, space = "both", what = "all")
      gate <- utl_ancestry_lbf(fit, bf_gate = 0)
      rho <- fit$slot_rho
      if (!is.array(rho) || length(dim(rho)) != 3L || any(!is.finite(rho))) {
        stop("slot_rho is invalid for LFSR rep", sprintf("%03d", rep_id), " ", arm)
      }
      alpha <- apply(rho, c(1L, 2L), sum)
      sn <- dimnames(rho)[[1L]]
      vn <- dimnames(rho)[[2L]]
      cn <- dimnames(rho)[[3L]]
      if (is.null(sn) || is.null(vn) || is.null(cn) ||
          !identical(vn, rownames(dat$truth_beta))) {
        stop("slot_rho names are incompatible for LFSR rep", sprintf("%03d", rep_id))
      }
      saveRDS(list(rep = rep_id, arm = arm, input_md5 = input_md5,
                   fit_sha256 = fit_sha256, lfsr = lfsr,
                   lfsr_output_semantics =
                     "supersedes_initial_lfsr_prior_support_bug"),
              file.path(case_out, paste0(arm, "_lfsr.rds")))
      owner <- lfsr$conditional$beta
      if (!is.array(owner) || !identical(dim(owner)[[1L]], dim(rho)[[1L]]) ||
          !identical(dim(owner)[[2L]], dim(rho)[[2L]]) ||
          !identical(dimnames(owner)[[2L]], vn) ||
          !identical(dimnames(owner)[[3L]], colnames(dat$truth_beta))) {
        stop("conditional beta LFSR dimensions are incompatible")
      }
      gate_method <- if (is.null(gate$gate_method)) {
        "legacy_residualized_wakefield"
      } else gate$gate_method
      raw_single <- matrix(1, nrow(alpha), ncol(gate$lbf_outcome),
                           dimnames = dimnames(gate$lbf_outcome))
      for (a in seq_len(ncol(gate$lbf_outcome))) {
        raw_single[, a] <- rowSums(alpha * lfsr$conditional$beta[, , a])
      }
      for (l in seq_len(nrow(gate$lbf_outcome))) for (a in seq_len(ncol(gate$lbf_outcome))) {
        gate_rows[[length(gate_rows) + 1L]] <- data.frame(
          rep = rep_id, arm = arm, slot = rownames(gate$lbf_outcome)[[l]],
          ancestry = colnames(gate$lbf_outcome)[[a]],
          raw_single = raw_single[l, a],
          gated_single = lfsr$single$beta[l, a],
          lbf_outcome = gate$lbf_outcome[l, a], pass = gate$gate[l, a],
          threshold = gate$bf_gate, gate_method = gate_method,
          status = status, fit_status = fit_status, input_md5 = input_md5,
          fit_sha256 = fit_sha256, stringsAsFactors = FALSE
        )
      }
      for (i in seq_along(refs)) {
        j <- ref_indices[[i]]
        amax <- max(alpha[, j])
        owners <- which(alpha[, j] >= amax - 1e-12)
        for (a in seq_along(cn <- colnames(dat$truth_beta))) {
          vals <- as.numeric(owner[owners, j, a])
          owner_rows[[length(owner_rows) + 1L]] <- data.frame(
            rep = rep_id, arm = arm, reference = refs[[i]], variant = vn[[j]],
            ancestry = cn[[a]], truth_beta = dat$truth_beta[j, a],
            truth_zero = isTRUE(dat$truth_beta[j, a] == 0),
            owner_slots = paste(sn[owners], collapse = ";"),
            owner_tie_count = length(owners), owner_alpha_max = amax,
            owner_conditional_lfsr = scalar_value(vals),
            owner_conditional_lfsr_values = format_values(vals),
            status = status, fit_status = fit_status, input_md5 = input_md5,
            fit_sha256 = fit_sha256, stringsAsFactors = FALSE
          )
        }
      }
      for (sp in c("theta", "beta")) {
        arr <- lfsr$variant[[sp]]
        coords <- dimnames(arr)[[2L]]
        for (a in seq_along(coords)) {
          truth_zero <- if (sp == "beta") dat$truth_beta[, coords[[a]]] == 0 else rep(NA, length(vn))
          for (j in seq_along(vn)) {
            variant_rows[[length(variant_rows) + 1L]] <- data.frame(
              rep = rep_id, arm = arm, space = sp, coordinate = coords[[a]],
              variant = vn[[j]], variant_index = j, causal = j %in% ref_indices,
              truth_zero = truth_zero[[j]], variant_lfsr = arr[j, a],
              status = status, fit_status = fit_status, input_md5 = input_md5,
              fit_sha256 = fit_sha256, stringsAsFactors = FALSE
            )
          }
        }
      }
      metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
        rep = rep_id, arm = arm, status = status, fit_status = fit_status,
        weight_model = if (is.null(fit$weight_model)) "legacy" else as.character(fit$weight_model),
        input_md5 = input_md5,
        fit_sha256 = fit_sha256, lfsr_path = file.path(sprintf("rep%03d", rep_id),
                                                       paste0(arm, "_lfsr.rds")),
        beta_single_supported = if (!is.null(lfsr$metadata)) lfsr$metadata$beta_single_supported else NA,
        beta_single_unsupported_reason = if (!is.null(lfsr$metadata) &&
                                             !is.null(lfsr$metadata$beta_single_unsupported_reason)) {
          lfsr$metadata$beta_single_unsupported_reason
        } else "",
        lfsr_output_semantics = "supersedes_initial_lfsr_prior_support_bug",
        stringsAsFactors = FALSE
      )
    } else {
      owner_rows[[length(owner_rows) + 1L]] <- append_missing_owner(
        dat, rep_id, arm, status, fit_status, input_md5
      )
      variant_rows[[length(variant_rows) + 1L]] <- append_missing_variant(
        dat, rep_id, arm, status, fit_status, input_md5
      )
      metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
        rep = rep_id, arm = arm, status = status, fit_status = fit_status,
        weight_model = if (arm == "mom") "mom_sa" else "legacy",
        input_md5 = input_md5, fit_sha256 = NA_character_, lfsr_path = "",
        beta_single_supported = NA, beta_single_unsupported_reason = "",
        lfsr_output_semantics = "supersedes_initial_lfsr_prior_support_bug",
        stringsAsFactors = FALSE
      )
      for (slot in paste0("L", seq_len(5L))) for (ancestry in colnames(dat$truth_beta)) {
        gate_rows[[length(gate_rows) + 1L]] <- data.frame(
          rep = rep_id, arm = arm, slot = slot, ancestry = ancestry,
          raw_single = NA_real_, gated_single = NA_real_,
          lbf_outcome = NA_real_, pass = NA, threshold = 0,
          gate_method = NA_character_, status = status, fit_status = fit_status,
          input_md5 = input_md5, fit_sha256 = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
owner_table <- do.call(rbind, owner_rows)
variant_table <- do.call(rbind, variant_rows)
metadata_table <- do.call(rbind, metadata_rows)
gate_table <- do.call(rbind, gate_rows)
write.table(owner_table, file.path(out_dir, "conditional_owner_lfsr.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(variant_table, file.path(out_dir, "variant_lfsr.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(metadata_table, file.path(out_dir, "metadata.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gate_table, file.path(out_dir, "beta_single_gate.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

paired_rows <- list()
pair_summary <- function(g, m, metric, space, coordinate, key_names) {
  names(g)[names(g) == "value"] <- "gaussian"
  names(m)[names(m) == "value"] <- "mom"
  pair <- merge(g, m, by = key_names, all = TRUE, sort = TRUE)
  valid <- is.finite(pair$gaussian) & is.finite(pair$mom)
  diff <- pair$mom[valid] - pair$gaussian[valid]
  n_pair <- length(diff)
  p_value <- if (n_pair >= 2L) stats::wilcox.test(diff, mu = 0, exact = FALSE)$p.value else NA_real_
  win <- diff < 0
  loss <- diff > 0
  tie <- diff == 0
  data.frame(
    metric = metric, space = space, coordinate = coordinate,
    expected_pair_n = nrow(pair), valid_paired_n = n_pair,
    gaussian_failures = sum(!is.finite(pair$gaussian)),
    mom_failures = sum(!is.finite(pair$mom)),
    mom_minus_gaussian_median = if (n_pair) median(diff) else NA_real_,
    mom_minus_gaussian_iqr = if (n_pair) IQR(diff) else NA_real_,
    mom_minus_gaussian_range_min = if (n_pair) min(diff) else NA_real_,
    mom_minus_gaussian_range_max = if (n_pair) max(diff) else NA_real_,
    lower_is_better = TRUE, wilcoxon_p = p_value,
    mom_win_n = sum(win), mom_tie_n = sum(tie), mom_loss_n = sum(loss),
    mom_win_tie_loss = paste(sum(win), sum(tie), sum(loss), sep = "/"),
    stringsAsFactors = FALSE
  )
}
for (sp in c("theta", "beta")) {
  x <- variant_table[variant_table$space == sp, c("rep", "arm", "coordinate", "variant", "variant_lfsr")]
  names(x)[names(x) == "variant_lfsr"] <- "value"
  for (coordinate in unique(x$coordinate)) {
    y <- x[x$coordinate == coordinate, c("rep", "arm", "variant", "value")]
    g <- y[y$arm == "gaussian", c("rep", "variant", "value")]
    m <- y[y$arm == "mom", c("rep", "variant", "value")]
    paired_rows[[length(paired_rows) + 1L]] <- pair_summary(
      g, m, "variant_lfsr", sp, coordinate, c("rep", "variant")
    )
  }
}
for (coordinate in unique(owner_table$ancestry)) {
  y <- owner_table[owner_table$ancestry == coordinate,
                   c("rep", "arm", "reference", "owner_conditional_lfsr")]
  names(y)[names(y) == "owner_conditional_lfsr"] <- "value"
  g <- y[y$arm == "gaussian", c("rep", "reference", "value")]
  m <- y[y$arm == "mom", c("rep", "reference", "value")]
  paired_rows[[length(paired_rows) + 1L]] <- pair_summary(
    g, m, "conditional_owner_lfsr", "beta", coordinate, c("rep", "reference")
  )
}
paired_table <- do.call(rbind, paired_rows)
write.table(paired_table, file.path(out_dir, "lfsr_paired_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("mom_sa_random20 post hoc LFSR completed\n")
