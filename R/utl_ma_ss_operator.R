.utl_validate_ma_ss_inputs <- function(XtX_list, Xty_list, yty, n, prior,
                                       input_mode = c("exact_ss", "rss_pseudo_ss"),
                                       input_certificate = NULL) {
  input_mode <- match.arg(input_mode)
  if (!inherits(prior, "utl_theta_prior")) {
    stop("prior must be a utl_theta_prior")
  }
  ancestry_names <- prior$ancestry_names
  T <- length(ancestry_names)
  if (!is.numeric(n) || !is.null(dim(n)) || length(n) != T ||
      is.null(names(n)) || !identical(names(n), ancestry_names)) {
    stop("n must be a named vector ordered as prior ancestry_names")
  }
  if (any(!is.finite(n)) || any(n <= 1) ||
      (!identical(input_mode, "rss_pseudo_ss") && any(n != round(n)))) {
    stop("n must contain finite integer-like values greater than one")
  }
  if (!is.numeric(yty) || !is.null(dim(yty)) || length(yty) != T ||
      is.null(names(yty)) || !identical(names(yty), ancestry_names)) {
    stop("yty must be a named vector ordered as prior ancestry_names")
  }
  if (any(!is.finite(yty))) stop("yty must be finite")
  if (any(yty <= 0)) stop("yty must be positive")
  if (!is.list(XtX_list) || is.null(names(XtX_list)) ||
      !identical(names(XtX_list), ancestry_names)) {
    stop("XtX_list must be a named list ordered by prior ancestry_names")
  }
  if (!is.list(Xty_list) || is.null(names(Xty_list)) ||
      !identical(names(Xty_list), ancestry_names)) {
    stop("Xty_list must be a named list ordered by prior ancestry_names")
  }
  if (!is.null(input_certificate)) {
    if (!is.list(input_certificate) ||
        !identical(input_certificate$type, "rss") ||
        !is.list(input_certificate$values) ||
        !identical(names(input_certificate$values), ancestry_names) ||
        !is.list(input_certificate$vectors) ||
        !identical(names(input_certificate$vectors), ancestry_names) ||
        !is.list(input_certificate$Dhalf) ||
        !identical(names(input_certificate$Dhalf), ancestry_names) ||
        is.null(input_certificate$source) ||
        !identical(names(input_certificate$source), ancestry_names) ||
        is.null(input_certificate$rank) ||
        !identical(names(input_certificate$rank), ancestry_names) ||
        is.null(input_certificate$ld_diag) ||
        !identical(names(input_certificate$ld_diag), ancestry_names)) {
      stop("input_certificate is invalid")
    }
  }

  p <- NULL
  variant_names <- NULL
  input_qc <- vector("list", T)
  for (t in seq_len(T)) {
    ancestry <- ancestry_names[[t]]
    A <- XtX_list[[t]]
    b <- Xty_list[[t]]
    if (!is.matrix(A) || !is.numeric(A) || nrow(A) != ncol(A) || nrow(A) < 1L) {
      stop("XtX for ancestry '", ancestry, "' must be a non-empty square numeric matrix")
    }
    if (any(!is.finite(A))) {
      stop("XtX for ancestry '", ancestry,
           "' contains missing or non-finite summary statistics")
    }
    A_scale <- max(1, max(abs(A)))
    if (max(abs(A - t(A))) > prior$validation_tol * A_scale) {
      stop("XtX for ancestry '", ancestry, "' must be symmetric")
    }
    if (any(diag(A) <= 0)) {
      stop("XtX for ancestry '", ancestry,
           "' must be positive semidefinite with positive diagonal entries; monomorphic variants are not supported")
    }
    A_names <- rownames(A)
    if (is.null(A_names) || is.null(colnames(A)) || !identical(A_names, colnames(A)) ||
        anyNA(A_names) || any(A_names == "") || anyDuplicated(A_names)) {
      stop("XtX for ancestry '", ancestry,
           "' must have matching, unique, non-empty variant row and column names")
    }
    if (!is.numeric(b) || !is.null(dim(b)) || length(b) != nrow(A) ||
        any(!is.finite(b))) {
      stop("Xty for ancestry '", ancestry,
           "' must be a finite numeric vector aligned to XtX")
    }
    if (is.null(names(b)) || !identical(names(b), A_names)) {
      stop("Xty for ancestry '", ancestry,
           "' names must match XtX variant names in order")
    }
    if (is.null(p)) {
      p <- nrow(A)
      variant_names <- A_names
    } else if (nrow(A) != p || !identical(A_names, variant_names)) {
      stop("XtX and Xty must use the same variant names and order across ancestries")
    }

    if (identical(input_mode, "rss_pseudo_ss") && is.null(input_certificate)) {
      input_qc[[t]] <- data.frame(
        ancestry = ancestry, input_mode = input_mode,
        qc_method = "not_computed",
        computed = FALSE, enforced = FALSE,
        source = NA_character_, rank = NA_integer_,
        ld_diag_min = NA_real_, ld_diag_max = NA_real_,
        diag_min = NA_real_, diag_max = NA_real_,
        range_ok = NA, range_residual = NA_real_,
        range_residual_scaled = NA_real_, projected_schur = NA_real_,
        schur_ok = NA, joint_compatible = NA,
        schur_complement = NA_real_, schur_complement_scaled = NA_real_,
        augmented_min_eig_raw = NA_real_, augmented_min_eig_scaled = NA_real_,
        rss_ld_mismatch = NA, warning_message = "",
        stringsAsFactors = FALSE
      )
      next
    }

    cert_source <- "xtx_matrix"
    ld_diag <- NULL
    if (!is.null(input_certificate)) {
      values <- input_certificate$values[[ancestry]]
      vectors <- input_certificate$vectors[[ancestry]]
      Dhalf <- input_certificate$Dhalf[[ancestry]]
      cert_source <- as.character(input_certificate$source[[ancestry]])
      ld_diag <- input_certificate$ld_diag[[ancestry]]
      if (!is.numeric(values) || length(values) < 1L ||
          any(!is.finite(values)) || any(values <= 0) ||
          !is.matrix(vectors) || !is.numeric(vectors) ||
          nrow(vectors) != nrow(A) || ncol(vectors) != length(values) ||
          any(!is.finite(vectors))) {
        stop("input_certificate eigen fields for ancestry '", ancestry, "' are invalid")
      }
      if (!is.numeric(Dhalf) || length(Dhalf) != nrow(A) ||
          any(!is.finite(Dhalf)) || any(Dhalf <= 0)) {
        stop("input_certificate Dhalf for ancestry '", ancestry, "' is invalid")
      }
      if (!is.numeric(ld_diag) || length(ld_diag) != nrow(A) ||
          any(!is.finite(ld_diag)) || any(ld_diag <= 0)) {
        stop("input_certificate LD diagonal for ancestry '", ancestry,
             "' must be a finite positive vector")
      }
      if (!is.null(rownames(vectors)) && !identical(rownames(vectors), A_names)) {
        stop("input_certificate eigen row names must match XtX for ancestry '",
             ancestry, "'")
      }
      A_rank <- as.integer(input_certificate$rank[[ancestry]])
      if (!is.finite(A_rank) || A_rank < 1L || A_rank != length(values)) {
        stop("input_certificate rank does not match its eigen fields for ancestry '",
             ancestry, "'")
      }
      R_diag <- rowSums(sweep(vectors^2, 2L, values, "*"))
      if (max(abs(R_diag - ld_diag)) > 100 * prior$validation_tol *
          max(1, max(abs(ld_diag)))) {
        stop("input_certificate LD diagonal does not match eigen fields for ancestry '",
             ancestry, "'")
      }
      A_diag_from_cert <- Dhalf^2 * R_diag
      if (max(abs(A_diag_from_cert - diag(A))) >
          100 * prior$validation_tol * max(1, max(abs(diag(A))))) {
        stop("input_certificate eigen fields do not match XtX for ancestry '",
             ancestry, "'")
      }
      q <- b / Dhalf
    } else {
      E <- .utl_ld_eigen_from_matrix(A, A_names, prior$validation_tol,
                                     paste0("XtX for ancestry '", ancestry, "'"))
      values <- E$values
      vectors <- E$vectors
      A_rank <- E$rank
      q <- b
    }
    proj <- .utl_eigen_projection(values, vectors, q, yty[[t]])
    range_ok <- proj$range_residual <= prior$validation_tol
    projected_schur <- proj$projected_schur
    schur_scale <- max(1, abs(yty[[t]]),
                       abs(yty[[t]] - projected_schur))
    schur_ok <- isTRUE(range_ok) &&
      projected_schur >= -prior$validation_tol * schur_scale
    joint_compatible <- isTRUE(range_ok) && isTRUE(schur_ok)
    schur_raw <- if (isTRUE(range_ok)) projected_schur else NA_real_
    schur_scaled <- if (isTRUE(range_ok)) projected_schur / yty[[t]] else NA_real_
    input_qc[[t]] <- data.frame(
      ancestry = ancestry, input_mode = input_mode,
      qc_method = "eigen_projection",
      computed = TRUE, enforced = FALSE,
      source = cert_source, rank = A_rank,
      ld_diag_min = if (is.null(ld_diag)) NA_real_ else min(ld_diag),
      ld_diag_max = if (is.null(ld_diag)) NA_real_ else max(ld_diag),
      diag_min = min(diag(A)), diag_max = max(diag(A)),
      range_ok = range_ok,
      range_residual = proj$range_residual,
      range_residual_scaled = proj$range_residual,
      projected_schur = projected_schur,
      schur_ok = schur_ok,
      joint_compatible = joint_compatible,
      schur_complement = schur_raw,
      schur_complement_scaled = schur_scaled,
      augmented_min_eig_raw = NA_real_,
      augmented_min_eig_scaled = NA_real_,
      rss_ld_mismatch = !joint_compatible,
      warning_message = "",
      stringsAsFactors = FALSE
    )
  }
  list(p = p, variant_names = variant_names,
       input_qc = do.call(rbind, input_qc))
}

#' Construct a multi-ancestry sufficient-statistics operator in theta space.
#'
#' This constructor validates centered, standardized ancestry-specific summary
#' statistics and returns `h`, the per-variant theta precision blocks, and a
#' matrix-free Hessian action. It does not fit a single-effect model.
#'
#' @param XtX_list Named ancestry-specific LD crossproduct matrices.
#' @param Xty_list Named ancestry-specific genotype-response crossproducts.
#' @param yty Named ancestry-specific response sums of squares.
#' @param n Named ancestry-specific sample sizes.
#' @param prior A [utl_theta_prior()] object.
#' @return An object of class `utl_ma_ss_operator` with `h`, `Pj`, and `H`.
utl_ma_ss_operator <- function(XtX_list, Xty_list, yty, n, prior,
                               input_mode = c("exact_ss", "rss_pseudo_ss"),
                               input_certificate = NULL) {
  input_mode <- match.arg(input_mode)
  dat <- .utl_validate_ma_ss_inputs(XtX_list, Xty_list, yty, n, prior,
                                     input_mode = input_mode,
                                     input_certificate = input_certificate)
  p <- dat$p
  variant_names <- dat$variant_names
  ancestry_names <- prior$ancestry_names
  theta_names <- prior$theta_names
  T <- length(ancestry_names)
  rss_common <- identical(input_mode, "rss_pseudo_ss") && is.null(input_certificate)
  d <- length(theta_names)
  M <- prior$M

  h0 <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
  Pj0 <- vector("list", p)
  names(Pj0) <- variant_names
  for (t in seq_len(T)) {
    m <- M[t, ]
    h0 <- h0 + tcrossprod(Xty_list[[t]], m)
  }
  if (rss_common) {
    Pcommon0 <- matrix(0, d, d,
                       dimnames = list(theta_names, theta_names))
    for (t in seq_len(T)) {
      m <- M[t, ]
      Pcommon0 <- Pcommon0 + yty[[t]] * tcrossprod(m)
    }
    Pj0 <- rep(list(Pcommon0), p)
    names(Pj0) <- variant_names
  } else {
    for (j in seq_len(p)) {
      P <- matrix(0, d, d, dimnames = list(theta_names, theta_names))
      for (t in seq_len(T)) {
        m <- M[t, ]
        P <- P + diag(XtX_list[[t]])[[j]] * tcrossprod(m)
      }
      Pj0[[j]] <- P
    }
  }

  backend <- new.env(parent = emptyenv())
  backend$name <- "base"
  backend$matvec_count <- 0L
  backend$cpp_matvec_count <- 0L
  H_scaled <- function(X, sigma2 = rep(1, T)) {
    if (!is.matrix(X) || !is.numeric(X) || nrow(X) != p || ncol(X) != d ||
        any(!is.finite(X))) {
      stop("X must be a finite numeric p by n_theta matrix")
    }
    if (!is.numeric(sigma2) || length(sigma2) != T || any(!is.finite(sigma2)) ||
        any(sigma2 <= 0)) {
      stop("sigma2 must be a finite positive vector with one value per ancestry")
    }
    out <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
    for (t in seq_len(T)) {
      m <- M[t, ]
      Axm <- .utl_matvec(XtX_list[[t]], as.vector(X %*% m), backend)
      out <- out + tcrossprod(Axm, m) / sigma2[[t]]
    }
    out
  }

  scaled <- function(sigma2 = rep(1, T)) {
    if (!is.numeric(sigma2) || length(sigma2) != T ||
        is.null(names(sigma2)) || !identical(names(sigma2), ancestry_names) ||
        any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
      stop("sigma2 must be a named, finite, positive vector ordered by ancestry")
    }
    h <- h0 / sigma2[[1L]]
    Pj <- lapply(Pj0, function(P) P / sigma2[[1L]])
    if (rss_common) {
      Pcommon <- matrix(0, d, d,
                        dimnames = list(theta_names, theta_names))
      for (t in seq_len(T)) {
        m <- M[t, ]
        Pcommon <- Pcommon + yty[[t]] * tcrossprod(m) / sigma2[[t]]
      }
      Pj <- rep(list(Pcommon), p)
      names(Pj) <- variant_names
    }
    if (T > 1L) {
      h <- matrix(0, p, d, dimnames = list(variant_names, theta_names))
      if (!rss_common) {
        Pj <- lapply(seq_len(p), function(j) {
          matrix(0, d, d, dimnames = list(theta_names, theta_names))
        })
      }
      for (t in seq_len(T)) {
        m <- M[t, ]
        h <- h + tcrossprod(Xty_list[[t]], m) / sigma2[[t]]
        if (!rss_common) {
          for (j in seq_len(p)) {
            Pj[[j]] <- Pj[[j]] + diag(XtX_list[[t]])[[j]] *
              tcrossprod(m) / sigma2[[t]]
          }
        }
      }
    }
    names(Pj) <- variant_names
    list(h = h, Pj = Pj, H = function(X) H_scaled(X, sigma2), sigma2 = sigma2)
  }
  H <- function(X) H_scaled(X, rep(1, T))

  structure(
    list(h = h0, Pj = Pj0, H = H, scaled = scaled, h0 = h0, Pj0 = Pj0,
         XtX_list = XtX_list, Xty_list = Xty_list, n = n, yty = yty,
         prior = prior, backend = backend, p = p, variant_names = variant_names,
         rss_common = rss_common,
         ancestry_names = ancestry_names, input_mode = input_mode,
         input_qc = dat$input_qc),
    class = "utl_ma_ss_operator"
  )
}
