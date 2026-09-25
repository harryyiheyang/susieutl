.utl_theta_global_elbo <- function(h, H, Pj, theta_bar, slot_mean, slot_second,
                                   slot_KL, n, yty, sigma2 = rep(1, length(n))) {
  p <- nrow(h)
  d <- ncol(h)
  L <- dim(slot_mean)[[1L]]
  if (!is.numeric(sigma2) || length(sigma2) != length(n) ||
      any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
    stop("theta ELBO sigma2 must be finite and positive")
  }
  if (any(!is.finite(theta_bar)) || any(!is.finite(slot_mean)) ||
      any(!is.finite(slot_second)) || any(!is.finite(slot_KL))) {
    stop("theta ELBO inputs must be finite")
  }
  Delta <- numeric(L)
  for (l in seq_len(L)) {
    Ml <- matrix(slot_mean[l, , ], p, d)
    H_Ml <- H(Ml)
    if (any(!is.finite(H_Ml))) stop("theta Hessian action is non-finite")
    tr_second <- 0
    for (j in seq_len(p)) {
      Tlj <- matrix(slot_second[l, j, , ], d, d)
      tr_second <- tr_second + sum(Pj[[j]] * Tlj)
    }
    Delta[[l]] <- tr_second - sum(Ml * H_Ml)
  }
  H_theta <- H(theta_bar)
  if (any(!is.finite(H_theta)) || any(!is.finite(Delta))) {
    stop("theta ELBO quadratic terms are non-finite")
  }
  C0 <- -0.5 * sum(n * log(2 * pi * sigma2) + yty / sigma2)
  Eloglik <- C0 + sum(h * theta_bar) -
    0.5 * (sum(theta_bar * H_theta) + sum(Delta))
  ELBO <- Eloglik - sum(slot_KL)
  if (!is.finite(C0) || !is.finite(Eloglik) || !is.finite(ELBO)) {
    stop("theta global ELBO is non-finite")
  }
  list(C0 = C0, Delta = Delta, Eloglik = Eloglik, ELBO = ELBO)
}

.utl_theta_fixed_elbo_state <- function(n, yty, sigma2, L) {
  C0 <- -0.5 * sum(n * log(2 * pi * sigma2) + yty / sigma2)
  list(C0 = C0, Eloglik = C0, ELBO = C0,
       A = numeric(L), KL = numeric(L))
}

.utl_theta_fixed_elbo_replace <- function(state, l, Q, Mold, Mnew,
                                          ser_e, KLnew) {
  Aold <- state$A[[l]]
  Anew <- 2 * (sum(Q * Mnew) - ser_e)
  dELL <- sum(Q * (Mnew - Mold)) - 0.5 * (Anew - Aold)
  state$Eloglik <- state$Eloglik + dELL
  state$ELBO <- state$ELBO + dELL - (KLnew - state$KL[[l]])
  state$A[[l]] <- Anew
  state$KL[[l]] <- KLnew
  state
}
