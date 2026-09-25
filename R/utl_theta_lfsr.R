.utl_lfsr_sign_probability <- function(mu, variance) {
  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) ||
      !is.numeric(variance) || length(variance) != 1L ||
      !is.finite(variance) || variance < 0) {
    stop("lfsr requires finite means and non-negative variances")
  }
  if (variance == 0) {
    if (mu > 0) return(c(minus = 0, zero = 0, plus = 1))
    if (mu < 0) return(c(minus = 1, zero = 0, plus = 0))
    return(c(minus = 0, zero = 1, plus = 0))
  }
  pminus <- stats::pnorm(0, mean = mu, sd = sqrt(variance))
  c(minus = pminus, zero = 0, plus = 1 - pminus)
}

.utl_lfsr_declared_support <- function(fit, contrast) {
  if (!is.matrix(contrast) || !is.numeric(contrast) ||
      any(!is.finite(contrast))) {
    stop("lfsr contrast is invalid")
  }
  L <- fit$L
  K <- dim(fit$slot_rho)[[3L]]
  d <- ncol(fit$M)
  q <- nrow(contrast)
  Ulist <- fit$prior$Utheta
  if (!is.list(Ulist) || length(Ulist) != K) {
    stop("lfsr requires declared prior covariance components")
  }
  if (is.null(fit$V)) stop("lfsr requires declared prior slot variances")
  factor_tol <- fit$prior$factor_tol
  if (is.null(factor_tol)) factor_tol <- 0
  if (!is.numeric(factor_tol) || length(factor_tol) != 1L ||
      !is.finite(factor_tol) || factor_tol < 0) {
    stop("lfsr prior factor tolerance is invalid")
  }
  r <- max(factor_tol, 100 * .Machine$double.eps)
  xmin <- .Machine$double.xmin
  zero <- array(FALSE, c(L, K, q))
  vscale <- array(0, c(L, K, q))
  for (l in seq_len(L)) {
    Vl <- if (is.list(fit$V)) fit$V[[l]] else fit$V[[l]]
    if (!is.numeric(Vl) || !(length(Vl) %in% c(1L, d)) ||
        any(!is.finite(Vl)) || any(Vl < 0)) {
      stop("lfsr prior slot variance is invalid for slot ", l)
    }
    for (k in seq_len(K)) {
      U <- Ulist[[k]]
      if (!is.matrix(U) || !identical(dim(U), c(d, d)) ||
          any(!is.finite(U))) {
        stop("lfsr declared prior covariance is invalid for component ", k)
      }
      C0 <- if (length(Vl) == 1L) {
        Vl[[1L]] * U
      } else {
        D <- diag(sqrt(Vl), nrow = d, ncol = d)
        D %*% U %*% D
      }
      C0 <- (C0 + t(C0)) / 2
      cscale <- norm(C0, type = "2")
      if (!is.finite(cscale)) stop("lfsr declared prior covariance is invalid")
      for (a in seq_len(q)) {
        z <- as.numeric(contrast[a, ])
        q0 <- drop(crossprod(z, C0 %*% z))
        vs <- cscale * sum(z^2)
        if (!is.finite(q0) || !is.finite(vs) || vs < 0) {
          stop("lfsr declared prior contrast variance is invalid")
        }
        if (q0 < 0) {
          if (abs(q0) <= 100 * .Machine$double.eps * vs) {
            q0 <- 0
          } else {
            stop("lfsr declared prior contrast covariance is not PSD")
          }
        }
        vscale[l, k, a] <- vs
        zero[l, k, a] <- q0 == 0
      }
    }
  }
  list(zero = zero, vscale = vscale, r = r, xmin = xmin)
}

.utl_lfsr_check_declared_zero <- function(mu, variance, vscale, r, xmin) {
  if (any(!is.finite(mu)) || any(!is.finite(variance))) {
    stop("lfsr declared structural zero has non-finite posterior moments")
  }
  variance_tol <- 10 * r * max(vscale, xmin)
  mean_tol <- 10 * r * max(sqrt(vscale), xmin)
  if (any(abs(variance) > variance_tol) || any(abs(mu) > mean_tol)) {
    stop("lfsr declared structural zero disagrees with posterior moments")
  }
  list(mu = rep(0, length(mu)), variance = rep(0, length(variance)))
}

.utl_lfsr_space <- function(fit, contrast, coordinate_names) {
  L <- fit$L
  p <- length(fit$pip)
  K <- dim(fit$slot_rho)[[3L]]
  q <- nrow(contrast)
  vn <- names(fit$pip)
  sn <- dimnames(fit$slot_rho)[[1L]]
  alpha <- apply(fit$slot_rho, c(1L, 2L), sum)
  if (any(!is.finite(alpha)) || any(alpha < 0) ||
      any(alpha > 1 + 100 * .Machine$double.eps)) {
    stop("lfsr requires finite slot inclusion probabilities in [0,1]")
  }
  pminus <- pzero <- array(
    0, c(L, p, q), dimnames = list(sn, vn, coordinate_names)
  )
  rho <- fit$slot_rho
  alpha3 <- array(rep(alpha, K), c(L, p, K))
  gamma <- sweep(rho, c(1L, 2L), alpha, "/")
  positive_alpha <- alpha3 > 0
  if (any(!is.finite(gamma[positive_alpha])) || any(gamma[positive_alpha] < 0) ||
      any(abs(apply(gamma, c(1L, 2L), sum)[alpha > 0] - 1) >
          100 * .Machine$double.eps)) {
    stop("conditional pattern probabilities must be finite and sum to one")
  }
  gamma[!positive_alpha] <- 0
  gamma4 <- array(gamma, c(L, p, K, q))

  cm <- matrix(fit$slot_component_mean, nrow = L * p * K,
               ncol = ncol(fit$M))
  mu4 <- array(cm %*% t(contrast), c(L, p, K, q))
  cv <- matrix(fit$slot_component_covariance, nrow = L * p * K,
               ncol = ncol(fit$M) * ncol(fit$M))
  variance4 <- array(0, c(L, p, K, q))
  for (a in seq_len(q)) {
    ma <- contrast[a, ]
    variance4[, , , a] <- array(
      cv %*% as.vector(outer(ma, ma)), c(L, p, K)
    )
  }
  mu <- as.vector(mu4)
  variance <- as.vector(variance4)
  if (any(!is.finite(mu)) || any(!is.finite(variance))) {
    stop("lfsr requires finite means and non-negative variances")
  }
  pminus4 <- pzero4 <- numeric(length(mu))
  positive_variance <- variance > 0
  pminus4[positive_variance] <- stats::pnorm(
    0, mean = mu[positive_variance], sd = sqrt(variance[positive_variance])
  )
  zero_variance <- !positive_variance
  pminus4[zero_variance & mu > 0] <- 0
  pminus4[zero_variance & mu < 0] <- 1
  pzero4[zero_variance & mu == 0] <- 1
  pminus4 <- array(pminus4, c(L, p, K, q))
  pzero4 <- array(pzero4, c(L, p, K, q))
  support <- .utl_lfsr_declared_support(fit, contrast)
  for (l in seq_len(L)) for (k in seq_len(K)) for (a in seq_len(q)) {
    if (support$zero[l, k, a]) {
      checked <- .utl_lfsr_check_declared_zero(
        mu4[l, , k, a], variance4[l, , k, a],
        support$vscale[l, k, a], support$r, support$xmin
      )
      pminus4[l, , k, a] <- checked$mu
      pzero4[l, , k, a] <- 1
    } else if (any(variance4[l, , k, a] < 0)) {
      stop("lfsr requires finite means and non-negative variances")
    }
  }
  gm <- aperm(gamma4 * pminus4, c(1L, 2L, 4L, 3L))
  gz <- aperm(gamma4 * pzero4, c(1L, 2L, 4L, 3L))
  pminus <- array(rowSums(matrix(gm, nrow = L * p * q, ncol = K)),
                  c(L, p, q), dimnames = list(sn, vn, coordinate_names))
  pzero <- array(rowSums(matrix(gz, nrow = L * p * q, ncol = K)),
                 c(L, p, q), dimnames = list(sn, vn, coordinate_names))
  for (a in seq_len(q)) {
    za <- pzero[, , a]
    za[alpha == 0] <- 1
    pzero[, , a] <- za
  }
  probability_tol <- 100 * .Machine$double.eps
  probability_sum <- pminus + pzero
  if (any(!is.finite(pminus)) || any(!is.finite(pzero)) ||
      any(pminus < 0) || any(pzero < 0) ||
      any(probability_sum > 1 + probability_tol)) {
    stop("conditional lfsr probabilities are invalid")
  }
  pplus <- array(pmax(0, 1 - pminus - pzero), dim(pminus),
                 dimnames = dimnames(pminus))
  conditional <- ashr::compute_lfsr(pminus, pzero)
  probability_sum <- pminus + pzero + pplus
  if (any(!is.finite(probability_sum)) ||
      any(abs(probability_sum - 1) > 1e-10) ||
      any(!is.finite(conditional)) || any(conditional < 0 | conditional > 1)) {
    stop("conditional lfsr probabilities are invalid")
  }
  single <- matrix(1, L, q, dimnames = list(sn, coordinate_names))
  active_alpha <- rowSums(alpha) > 0
  for (a in seq_len(q)) {
    sa <- rowSums(alpha * conditional[, , a])
    single[active_alpha, a] <- sa[active_alpha]
  }
  variant <- matrix(1, p, q, dimnames = list(vn, coordinate_names))
  for (a in seq_len(q)) {
    variant[, a] <- pmax(1e-20,
                          1 - apply(alpha * (1 - conditional[, , a]), 2L, max))
  }
  list(conditional = conditional, single = single, variant = variant,
       probability = list(minus = pminus, zero = pzero, plus = pplus),
       alpha = alpha)
}

.utl_mom_lfsr_space <- function(fit, contrast, coordinate_names) {
  L <- fit$L
  p <- length(fit$pip)
  K <- dim(fit$slot_rho)[[3L]]
  d <- ncol(fit$M)
  q <- nrow(contrast)
  vn <- names(fit$pip)
  sn <- dimnames(fit$slot_rho)[[1L]]
  cn <- dimnames(fit$slot_rho)[[3L]]
  if (!is.array(fit$slot_rho) || length(dim(fit$slot_rho)) != 3L ||
      !identical(dim(fit$slot_rho), c(L, p, K)) || is.null(sn) ||
      is.null(vn) || is.null(cn)) {
    stop("MOM lfsr requires named slot pattern probabilities")
  }
  alpha <- apply(fit$slot_rho, c(1L, 2L), sum)
  if (any(!is.finite(alpha)) || any(alpha < 0) ||
      any(alpha > 1 + 100 * .Machine$double.eps)) {
    stop("MOM lfsr requires finite slot inclusion probabilities in [0,1]")
  }
  if (!is.matrix(contrast) || !identical(dim(contrast), c(q, d)) ||
      any(!is.finite(contrast))) stop("MOM lfsr contrast is invalid")
  if (!is.list(fit$mom_sa_slot) || length(fit$mom_sa_slot) != L) {
    stop("MOM lfsr requires per-slot MOM audit data")
  }
  if (is.null(names(fit$mom_sa_slot))) names(fit$mom_sa_slot) <- sn
  pminus4 <- pzero4 <- array(0, c(L, p, K, q))
  support <- .utl_lfsr_declared_support(fit, contrast)
  target_names <- c("D", "N_E")
  probability_tol <- 1e-12
  for (l in seq_len(L)) {
    audit <- fit$mom_sa_slot[[l]]
    required <- c("A", "target_components", "z0", "t",
                  "base_component_mean", "base_component_covariance")
    if (!is.list(audit) || any(vapply(required, function(n) is.null(audit[[n]]), logical(1)))) {
      stop("MOM lfsr audit is incomplete for slot ", sn[[l]])
    }
    A <- audit$A
    if (!is.matrix(A) || !identical(dim(A), c(d, d)) || any(!is.finite(A))) {
      stop("MOM lfsr audit A is invalid for slot ", sn[[l]])
    }
    A_tol <- 1e-10 * max(1, max(abs(A)))
    if (max(abs(A - t(A))) > A_tol) stop("MOM lfsr audit A is not symmetric")
    target <- as.character(audit$target_components)
    if (!setequal(target, target_names)) {
      stop("MOM lfsr audit targets must be D and N_E")
    }
    z0 <- audit$z0
    tmat <- audit$t
    cm <- audit$base_component_mean
    cv <- audit$base_component_covariance
    if (!is.numeric(z0) || is.null(names(z0)) ||
        !is.matrix(tmat) || !identical(dim(tmat), c(p, K)) ||
        !is.array(cm) || !identical(dim(cm), c(p, K, d)) ||
        !is.array(cv) || !identical(dim(cv), c(p, K, d, d)) ||
        any(!is.finite(cm)) || any(!is.finite(cv))) {
      stop("MOM lfsr audit moments have invalid dimensions for slot ", sn[[l]])
    }
    if (!identical(dimnames(tmat)[[1L]], vn) ||
        !identical(dimnames(tmat)[[2L]], cn) ||
        !identical(dimnames(cm)[[1L]], vn) ||
        !identical(dimnames(cm)[[2L]], cn) ||
        !identical(dimnames(cv)[[1L]], vn) ||
        !identical(dimnames(cv)[[2L]], cn)) {
      stop("MOM lfsr audit moment names are incompatible")
    }
    for (k in seq_len(K)) {
      is_target <- cn[[k]] %in% target_names
      if (is_target && (!is.finite(z0[[cn[[k]]]]) || z0[[cn[[k]]]] <= 0)) {
        stop("MOM lfsr audit z0 is invalid for ", cn[[k]], " in slot ", sn[[l]])
      }
      for (j in seq_len(p)) {
        m <- as.numeric(cm[j, k, ])
        C <- matrix(cv[j, k, , ], d, d)
        C_tol <- 1e-10 * max(1, max(abs(C)))
        if (max(abs(C - t(C))) > C_tol) {
          stop("MOM lfsr component covariance is not symmetric")
        }
        tval <- if (is_target) tmat[j, k] else NA_real_
        if (is_target) {
          if (!is.finite(tval) || tval <= 0) {
            stop("MOM lfsr saved tilt moment is invalid for ", vn[[j]], " ", cn[[k]])
          }
          tcheck <- sum(diag(A %*% C)) + drop(crossprod(m, A %*% m))
          t_tol <- 1e-8 * max(1e-12, abs(tval), abs(tcheck))
          if (!is.finite(tcheck) || abs(tcheck - tval) > t_tol) {
            stop("MOM lfsr saved tilt moment does not match audit moments")
          }
        }
        for (a in seq_len(q)) {
          z <- as.numeric(contrast[a, ])
          mu <- drop(crossprod(z, m))
          dvec <- drop(C %*% z)
          s2 <- drop(crossprod(z, dvec))
          tau_v <- 100 * .Machine$double.eps *
            max(1, sum(z^2) * norm(C, type = "I"))
          if (!is.finite(mu) || !is.finite(s2) || s2 < -tau_v) {
            stop("MOM lfsr contrast variance is invalid")
          }
          if (s2 < 0) s2 <- 0
          if (support$zero[l, k, a]) {
            checked <- .utl_lfsr_check_declared_zero(
              mu, s2, support$vscale[l, k, a], support$r, support$xmin
            )
            pminus4[l, j, k, a] <- 0
            pzero4[l, j, k, a] <- 1
          } else if (s2 == 0) {
            prob <- .utl_lfsr_sign_probability(mu, 0)
            pminus4[l, j, k, a] <- prob[["minus"]]
            pzero4[l, j, k, a] <- prob[["zero"]]
          } else if (!is_target) {
            prob <- .utl_lfsr_sign_probability(mu, s2)
            pminus4[l, j, k, a] <- prob[["minus"]]
            pzero4[l, j, k, a] <- prob[["zero"]]
          } else {
            sd <- sqrt(s2)
            zz <- -mu / sd
            Phi <- stats::pnorm(zz)
            phi <- stats::dnorm(zz)
            B0 <- sum(diag(A %*% (C - tcrossprod(dvec) / s2))) +
              drop(crossprod(m, A %*% m))
            B1 <- 2 * drop(crossprod(m, A %*% dvec)) / s2
            B2 <- drop(crossprod(dvec, A %*% dvec)) / s2^2
            pminus <- (B0 * Phi - B1 * sd * phi +
                       B2 * s2 * (Phi - zz * phi)) / tval
            if (!is.finite(pminus) || pminus < -probability_tol ||
                pminus > 1 + probability_tol) {
              stop("MOM lfsr signed probability is out of bounds")
            }
            if (pminus < 0) pminus <- 0
            if (pminus > 1) pminus <- 1
            pminus4[l, j, k, a] <- pminus
            pzero4[l, j, k, a] <- 0
          }
        }
      }
    }
  }
  alpha3 <- array(rep(alpha, K), c(L, p, K))
  gamma <- sweep(fit$slot_rho, c(1L, 2L), alpha, "/")
  positive_alpha <- alpha3 > 0
  if (any(!is.finite(gamma[positive_alpha])) || any(gamma[positive_alpha] < 0) ||
      any(abs(apply(gamma, c(1L, 2L), sum)[alpha > 0] - 1) >
          100 * .Machine$double.eps)) {
    stop("conditional pattern probabilities must be finite and sum to one")
  }
  gamma[!positive_alpha] <- 0
  gamma4 <- array(gamma, c(L, p, K, q))
  pminus <- array(rowSums(matrix(aperm(gamma4 * pminus4, c(1L, 2L, 4L, 3L)),
                                 nrow = L * p * q, ncol = K)),
                  c(L, p, q), dimnames = list(sn, vn, coordinate_names))
  pzero <- array(rowSums(matrix(aperm(gamma4 * pzero4, c(1L, 2L, 4L, 3L)),
                                nrow = L * p * q, ncol = K)),
                 c(L, p, q), dimnames = list(sn, vn, coordinate_names))
  for (a in seq_len(q)) {
    za <- pzero[, , a]
    za[alpha == 0] <- 1
    pzero[, , a] <- za
  }
  if (any(!is.finite(pminus)) || any(!is.finite(pzero)) ||
      any(pminus < 0) || any(pzero < 0) ||
      any(pminus + pzero > 1 + probability_tol)) {
    stop("MOM lfsr probabilities are invalid")
  }
  pplus <- array(pmax(0, 1 - pminus - pzero), dim(pminus),
                 dimnames = dimnames(pminus))
  probability_sum <- pminus + pzero + pplus
  if (any(!is.finite(probability_sum)) ||
      any(abs(probability_sum - 1) > 1e-10)) {
    stop("MOM lfsr probabilities do not sum to one")
  }
  conditional <- ashr::compute_lfsr(pminus, pzero)
  if (any(!is.finite(conditional)) || any(conditional < 0 | conditional > 1)) {
    stop("MOM lfsr conditional rates are invalid")
  }
  single <- matrix(1, L, q, dimnames = list(sn, coordinate_names))
  active_alpha <- rowSums(alpha) > 0
  for (a in seq_len(q)) {
    sa <- rowSums(alpha * conditional[, , a])
    single[active_alpha, a] <- sa[active_alpha]
  }
  variant <- matrix(1, p, q, dimnames = list(vn, coordinate_names))
  for (a in seq_len(q)) {
    variant[, a] <- pmax(1e-20,
                          1 - apply(alpha * (1 - conditional[, , a]), 2L, max))
  }
  list(conditional = conditional, single = single, variant = variant,
       probability = list(minus = pminus, zero = pzero, plus = pplus),
       alpha = alpha)
}

.utl_mom_ancestry_lbf <- function(fit, bf_gate = 0) {
  method <- "mom_sa_prior_second_moment_matched_normal_approx"
  L <- fit$L
  p <- length(fit$pip)
  T <- nrow(fit$M)
  d <- ncol(fit$M)
  K <- length(fit$prior$Utheta)
  expected_components <- c("D", "N_E", "N_A", "C", "S_E", "S_A")
  if (!is.matrix(fit$M) || !identical(dim(fit$M), c(T, d)) ||
      any(!is.finite(fit$M)) || is.null(rownames(fit$M))) {
    stop("MOM beta gate requires a finite named ancestry map")
  }
  if (!is.array(fit$slot_rho) || !identical(dim(fit$slot_rho), c(L, p, K)) ||
      is.null(dimnames(fit$slot_rho)[[1L]]) ||
      is.null(dimnames(fit$slot_rho)[[2L]]) ||
      is.null(dimnames(fit$slot_rho)[[3L]])) {
    stop("MOM beta gate requires named slot pattern probabilities")
  }
  component_names <- dimnames(fit$slot_rho)[[3L]]
  if (!identical(component_names, expected_components) ||
      !is.list(fit$prior$Utheta) ||
      !identical(names(fit$prior$Utheta), expected_components) ||
      length(fit$prior$Utheta) != K) {
    stop("MOM beta gate requires the UTL6 prior component order")
  }
  alpha <- apply(fit$slot_rho, c(1L, 2L), sum)
  if (any(!is.finite(alpha)) || any(alpha < 0) ||
      any(alpha > 1 + 100 * .Machine$double.eps)) {
    stop("MOM beta gate requires finite slot inclusion probabilities in [0,1]")
  }
  active <- fit$active_slots
  if (is.null(active)) active <- rowSums(alpha) > 0
  if (!is.logical(active) || length(active) != L || anyNA(active)) {
    stop("MOM beta gate active slot state is invalid")
  }
  if (!is.numeric(fit$weights) || length(fit$weights) != K ||
      !identical(names(fit$weights), expected_components) ||
      any(!is.finite(fit$weights)) || any(fit$weights < 0) ||
      abs(sum(fit$weights) - 1) > 1e-10 ||
      any(abs(fit$weights - rep(1 / K, K)) > 1e-10)) {
    stop("MOM beta gate requires uniform fixed prior weights")
  }
  weights_by_slot <- fit$weights_by_slot
  if (!is.matrix(weights_by_slot) ||
      !identical(dim(weights_by_slot), c(L, K)) ||
      is.null(rownames(weights_by_slot)) ||
      !identical(colnames(weights_by_slot), expected_components) ||
      any(!is.finite(weights_by_slot)) || any(weights_by_slot < 0) ||
      any(abs(rowSums(weights_by_slot) - 1) > 1e-10) ||
      any(abs(sweep(weights_by_slot, 2L, fit$weights, "-") ) > 1e-10)) {
    stop("MOM beta gate requires final weights_by_slot matching fixed prior weights")
  }
  if (is.null(fit$V) || is.null(fit$mom_sa_slot) ||
      !is.list(fit$mom_sa_slot) || length(fit$mom_sa_slot) != L) {
    stop("MOM beta gate requires per-slot prior variance and MOM audit")
  }
  if (is.null(fit$theta_mean) || is.null(fit$slot_mean) ||
      is.null(fit$XtX_list) || is.null(fit$Xty_list) ||
      is.null(fit$residual_variance)) {
    stop("MOM beta gate requires retained residualized likelihood state")
  }
  sn <- dimnames(fit$slot_rho)[[1L]]
  vn <- dimnames(fit$slot_rho)[[2L]]
  an <- rownames(fit$M)
  out <- array(0, c(L, p, T), dimnames = list(sn, vn, an))
  beta_mean <- fit$theta_mean %*% t(fit$M)
  if (!is.matrix(beta_mean) || !identical(dim(beta_mean), c(p, T)) ||
      any(!is.finite(beta_mean))) {
    stop("MOM beta gate theta means have incompatible dimensions")
  }
  for (l in seq_len(L)) {
    if (!active[[l]]) next
    Vl <- if (is.list(fit$V)) fit$V[[l]] else fit$V[[l]]
    if (!is.numeric(Vl) || length(Vl) != 1L || !is.finite(Vl) || Vl <= 0) {
      stop("MOM beta gate requires one finite positive scalar V for slot ", sn[[l]])
    }
    audit <- fit$mom_sa_slot[[l]]
    required <- c("A", "target_components", "z0")
    if (!is.list(audit) || any(vapply(required, function(n) is.null(audit[[n]]), logical(1)))) {
      stop("MOM beta gate audit is incomplete for slot ", sn[[l]])
    }
    A <- audit$A
    if (!is.matrix(A) || !identical(dim(A), c(d, d)) || any(!is.finite(A))) {
      stop("MOM beta gate audit A is invalid for slot ", sn[[l]])
    }
    A <- (A + t(A)) / 2
    if (max(abs(A - t(A))) > 1e-10 * max(1, max(abs(A)))) {
      stop("MOM beta gate audit A is not symmetric for slot ", sn[[l]])
    }
    if (!setequal(as.character(audit$target_components), c("D", "N_E"))) {
      stop("MOM beta gate audit targets must be D and N_E for slot ", sn[[l]])
    }
    z0 <- audit$z0
    if (!is.numeric(z0) || is.null(names(z0))) {
      stop("MOM beta gate audit z0 is invalid for slot ", sn[[l]])
    }
    Cstar <- vector("list", K)
    for (k in seq_len(K)) {
      U <- fit$prior$Utheta[[k]]
      if (!is.matrix(U) || !identical(dim(U), c(d, d)) || any(!is.finite(U))) {
        stop("MOM beta gate prior covariance is invalid for component ", component_names[[k]])
      }
      C0 <- Vl[[1L]] * (U + t(U)) / 2
      if (component_names[[k]] %in% c("D", "N_E")) {
        z0k <- z0[[component_names[[k]]]]
        denom <- sum(diag(A %*% C0))
        if (!is.finite(z0k) || z0k <= 0 || !is.finite(denom) || denom <= 0 ||
            abs(z0k - denom) > 1e-8 * max(1e-12, abs(z0k), abs(denom))) {
          stop("MOM beta gate audit z0 does not match declared prior for ",
               component_names[[k]], " in slot ", sn[[l]])
        }
        Cstar[[k]] <- C0 + 2 * C0 %*% A %*% C0 / denom
      } else {
        Cstar[[k]] <- C0
      }
      if (any(!is.finite(Cstar[[k]]))) {
        stop("MOM beta gate matched prior covariance is invalid")
      }
      Cstar[[k]] <- (Cstar[[k]] + t(Cstar[[k]])) / 2
    }
    beta_slot <- matrix(fit$slot_mean[l, , ], p, d) %*% t(fit$M)
    if (any(!is.finite(beta_slot))) stop("MOM beta gate slot means are invalid")
    beta_without <- beta_mean - beta_slot
    for (t in seq_len(T)) {
      a <- fit$M[t, ]
      W <- sum(vapply(seq_len(K), function(k) {
        fit$weights_by_slot[l, k] * drop(crossprod(a, Cstar[[k]] %*% a))
      }, numeric(1)))
      if (!is.finite(W) || W < 0) stop("MOM beta gate prior variance is invalid")
      XtX <- fit$XtX_list[[t]]
      Xty <- fit$Xty_list[[t]]
      Dg <- diag(XtX)
      sigma2 <- fit$residual_variance[[t]]
      if (!is.matrix(XtX) || !identical(dim(XtX), c(p, p)) ||
          !is.numeric(Xty) || length(Xty) != p ||
          any(!is.finite(Dg)) || any(Dg <= 0) ||
          !is.finite(sigma2) || sigma2 <= 0) {
        stop("MOM beta gate likelihood state is invalid for ancestry ", an[[t]])
      }
      residual <- Xty - as.numeric(XtX %*% beta_without[, t])
      bhat <- residual / Dg
      se2 <- sigma2 / Dg
      if (any(!is.finite(bhat)) || any(!is.finite(se2)) || any(se2 <= 0)) {
        stop("MOM beta gate residualized likelihood is invalid")
      }
      if (W > 0) {
        z2 <- bhat^2 / se2
        tau <- W / se2
        out[l, , t] <- .5 * z2 * tau / (1 + tau) - .5 * log1p(tau)
      }
    }
  }
  outcome <- matrix(0, L, T, dimnames = list(sn, an))
  for (l in seq_len(L)) for (t in seq_len(T)) {
    outcome[l, t] <- sum(alpha[l, ] * out[l, , t])
  }
  gate <- outcome >= bf_gate
  list(lbf = out, lbf_outcome = outcome, bf_gate = bf_gate, gate = gate,
       gate_status = ifelse(active, ifelse(gate, "pass", "fail"), "inactive"),
       gate_method = method, gate_approximate = TRUE,
       gate_disclosure = "MOM prior second-moment-matched normal approximation")
}

#' Extract theta and observed-beta local false sign rates.
#'
#' Conditional rates integrate all prior-pattern components within each SER
#' and SNP. Exact-zero component mass is retained as false-sign uncertainty.
#' Single-effect rates average conditional rates over SNP posterior inclusion
#' probabilities. Variant rates use the most confident SER without combining
#' independent SERs into a cross-slot sign distribution.
#'
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @param space One of `"theta"`, `"beta"`, or `"both"`.
#' @param what One of `"conditional"`, `"single"`, `"variant"`, or `"all"`.
#' @return Requested theta-coordinate and observed-beta lfsr arrays. Beta
#'   single-effect rates failing the fit's ancestry gate are one; MOM gate
#'   metadata identifies its second-moment-matched normal approximation.
#' @export
utl_lfsr <- function(fit, space = c("both", "theta", "beta"),
                     what = c("all", "conditional", "single", "variant")) {
  .utl_check_theta_fit(fit)
  space <- match.arg(space)
  what <- match.arg(what)
  is_mom <- identical(fit$weight_model, "mom_sa") ||
    identical(fit$weight_model, "mom_sa_pilot")
  if (is_mom) {
    theta <- .utl_mom_lfsr_space(
      fit, diag(ncol(fit$M)), colnames(fit$M)
    )
    beta <- .utl_mom_lfsr_space(fit, fit$M, rownames(fit$M))
    ancestry_lbf <- utl_ancestry_lbf(fit, bf_gate = 0)
    beta$single[!ancestry_lbf$gate] <- 1
    metadata <- list(
      model = as.character(fit$weight_model),
      beta_single_supported = TRUE,
      beta_single_gate_method = ancestry_lbf$gate_method,
      beta_single_gate_threshold = ancestry_lbf$bf_gate,
      beta_single_gate_approximate = ancestry_lbf$gate_approximate,
      beta_single_gate_disclosure = ancestry_lbf$gate_disclosure,
      beta_single_gate_status = ancestry_lbf$gate_status
    )
    ans <- list(
      conditional = list(theta = theta$conditional, beta = beta$conditional),
      single = list(theta = theta$single, beta = beta$single),
      variant = list(theta = theta$variant, beta = beta$variant)
    )
    if (space == "theta") ans <- lapply(ans, `[[`, "theta")
    if (space == "beta") ans <- lapply(ans, `[[`, "beta")
    if (what == "all") {
      if (space == "both") ans$metadata <- metadata else attr(ans, "metadata") <- metadata
      return(ans)
    }
    selected <- ans[[what]]
    attr(selected, "metadata") <- metadata
    return(selected)
  }
  if (is.null(fit$slot_component_mean) ||
      is.null(fit$slot_component_covariance)) {
    stop("fit does not contain component posterior moments required for lfsr")
  }
  d <- ncol(fit$M)
  theta_contrast <- diag(d)
  dimnames(theta_contrast) <- list(colnames(fit$M), colnames(fit$M))
  beta_contrast <- fit$M
  theta <- .utl_lfsr_space(fit, theta_contrast, colnames(fit$M))
  beta <- .utl_lfsr_space(fit, beta_contrast, rownames(fit$M))
  ancestry_lbf <- fit$ancestry_lbf
  if (is.null(ancestry_lbf)) ancestry_lbf <- utl_ancestry_lbf(fit, bf_gate = 0)
  if (is.null(ancestry_lbf$gate) ||
      !identical(dim(ancestry_lbf$gate), c(fit$L, nrow(fit$M)))) {
    stop("fit ancestry-LBF gate has incompatible dimensions")
  }
  beta$single[!ancestry_lbf$gate] <- 1
  ans <- list(
    conditional = list(theta = theta$conditional, beta = beta$conditional),
    single = list(theta = theta$single, beta = beta$single),
    variant = list(theta = theta$variant, beta = beta$variant)
  )
  if (space == "theta") ans <- lapply(ans, `[[`, "theta")
  if (space == "beta") ans <- lapply(ans, `[[`, "beta")
  if (what == "all") return(ans)
  ans[[what]]
}

#' Extract posterior effects in observed beta coordinates.
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @return A list containing beta means, full covariance, variances, and slot means.
#' @export
utl_beta_effects <- function(fit) {
  .utl_check_theta_fit(fit)
  M <- fit$M
  p <- nrow(fit$theta_mean)
  T <- nrow(M)
  bm <- fit$theta_mean %*% t(M)
  dimnames(bm) <- list(rownames(fit$theta_mean), rownames(M))
  bc <- array(0, c(p, T, T), dimnames = list(rownames(fit$theta_mean),
                                               rownames(M), rownames(M)))
  bv <- matrix(0, p, T, dimnames = list(rownames(fit$theta_mean), rownames(M)))
  for (j in seq_len(p)) {
    C <- M %*% matrix(fit$theta_covariance[j, , ], ncol(M), ncol(M)) %*% t(M)
    bc[j, , ] <- C
    bv[j, ] <- diag(C)
  }
  sm <- lapply(seq_len(fit$L), function(l) {
    x <- matrix(fit$slot_mean[l, , ], p, ncol(M)) %*% t(M)
    dimnames(x) <- list(rownames(fit$theta_mean), rownames(M))
    x
  })
  names(sm) <- dimnames(fit$slot_mean)[[1L]]
  list(mean = bm, covariance = bc, variance = bv, slot_mean = sm)
}

#' Compute ancestry-specific observed-scale log Bayes factors.
#'
#' For each SER and ancestry, the final pattern-weighted theta prior covariance
#' is projected to one observed-beta prior variance before evaluating the
#' univariate log Bayes factor at every SNP. The outcome score is the
#' SNP-posterior-weighted mean log Bayes factor.
#'
#' @param fit A `utl_theta_ma_ss_fit` object.
#' @param bf_gate Finite log-Bayes-factor threshold for the ancestry gate.
#' @return A list with `lbf` (`slot` by `variant` by `ancestry`),
#'   `lbf_outcome` (`slot` by `ancestry`), `bf_gate`, logical `gate`, and
#'   character `gate_status`.
#' @export
utl_ancestry_lbf <- function(fit, bf_gate = 0) {
  .utl_check_theta_fit(fit)
  if (identical(fit$weight_model, "mom_sa") ||
      identical(fit$weight_model, "mom_sa_pilot")) {
    if (!is.numeric(bf_gate) || length(bf_gate) != 1L || !is.finite(bf_gate)) {
      stop("bf_gate must be one finite number")
    }
    cached <- fit$ancestry_lbf
    if (is.list(cached) && identical(cached$gate_method,
                                     "mom_sa_prior_second_moment_matched_normal_approx") &&
        identical(cached$bf_gate, bf_gate) &&
        is.array(cached$lbf) &&
        identical(dim(cached$lbf), c(fit$L, length(fit$pip), nrow(fit$M))) &&
        is.matrix(cached$lbf_outcome) &&
        identical(dim(cached$lbf_outcome), c(fit$L, nrow(fit$M)))) {
      return(cached)
    }
    return(.utl_mom_ancestry_lbf(fit, bf_gate = bf_gate))
  }
  if (!is.numeric(bf_gate) || length(bf_gate) != 1L || !is.finite(bf_gate)) {
    stop("bf_gate must be one finite number")
  }
  if (is.null(fit$XtX_list) || is.null(fit$Xty_list) ||
      is.null(fit$slot_mean) || is.null(fit$residual_variance)) {
    stop("fit does not retain state required for ancestry LBF")
  }
  L <- fit$L
  p <- length(fit$pip)
  T <- nrow(fit$M)
  d <- ncol(fit$M)
  K <- length(fit$prior$Utheta)
  if (!is.list(fit$prior$F) || length(fit$prior$F) != K) {
    stop("fit prior factors are unavailable for ancestry LBF")
  }
  weights <- fit$weights
  if (!is.numeric(weights) || length(weights) != K ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-10) {
    stop("final prior mixture weights must be finite, non-negative, and sum to one")
  }
  alpha <- apply(fit$slot_rho, c(1L, 2L), sum)
  beta_mean <- fit$theta_mean %*% t(fit$M)
  out <- array(0, c(L, p, T),
               dimnames = list(dimnames(fit$slot_rho)[[1L]],
                               names(fit$pip), rownames(fit$M)))
  for (l in seq_len(L)) {
    Vl <- fit$V[[l]]
    if (!is.numeric(Vl) || any(!is.finite(Vl)) || any(Vl < 0) ||
        !length(Vl) %in% c(1L, d)) {
      stop("final prior variance has incompatible dimensions")
    }
    beta_slot <- matrix(fit$slot_mean[l, , ], p, d) %*% t(fit$M)
    beta_without <- beta_mean - beta_slot
    for (t in seq_len(T)) {
      m <- fit$M[t, ]
      W <- 0
      for (k in seq_len(K)) {
        if (length(Vl) == 1L) {
          G <- sqrt(Vl) * fit$prior$F[[k]]
        } else {
          D <- diag(sqrt(Vl), d, d)
          G <- D %*% fit$prior$F[[k]]
        }
        projected_factor <- as.numeric(crossprod(G, m))
        W <- W + weights[[k]] * sum(projected_factor^2)
      }
      if (!is.finite(W) || W < 0) {
        stop("ancestry prior variance must be finite and non-negative")
      }
      Dg <- diag(fit$XtX_list[[t]])
      if (length(Dg) != p || any(!is.finite(Dg)) || any(Dg <= 0)) {
        stop("ancestry XtX diagonal must be finite and positive")
      }
      sigma2 <- fit$residual_variance[[t]]
      if (!is.finite(sigma2) || sigma2 <= 0) {
        stop("ancestry residual variance must be finite and positive")
      }
      residual <- fit$Xty_list[[t]] -
        as.numeric(fit$XtX_list[[t]] %*% beta_without[, t])
      bhat <- residual / Dg
      se2 <- sigma2 / Dg
      if (W > 0) {
        z2 <- bhat^2 / se2
        tau <- W / se2
        out[l, , t] <- 0.5 * z2 * tau / (1 + tau) - 0.5 * log1p(tau)
      }
    }
  }
  outcome <- matrix(0, L, T,
                    dimnames = list(dimnames(out)[[1L]], dimnames(out)[[3L]]))
  for (l in seq_len(L)) for (t in seq_len(T)) {
    outcome[l, t] <- sum(alpha[l, ] * out[l, , t])
  }
  gate <- outcome >= bf_gate
  list(lbf = out, lbf_outcome = outcome, bf_gate = bf_gate, gate = gate,
       gate_status = ifelse(gate, "pass", "fail"))
}
