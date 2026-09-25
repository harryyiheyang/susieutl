\name{utl_lfsr}
\alias{utl_lfsr}
\title{Extract Theta and Observed-Beta Local False Sign Rates}
\usage{
utl_lfsr(
  fit,
  space = c("both", "theta", "beta"),
  what = c("all", "conditional", "single", "variant")
)
}
\arguments{
\item{fit}{A \code{utl_theta_ma_ss_fit} object.}

\item{space}{One of \code{"theta"}, \code{"beta"}, or \code{"both"}.}

\item{what}{One of \code{"conditional"}, \code{"single"},
\code{"variant"}, or \code{"all"}.}
}
\value{
Requested theta-coordinate and observed-beta lfsr arrays. Gaussian beta
single-effect rates failing the fit's ancestry-LBF gate are one. MOM beta
single-effect rates use the MOM second-moment-matched normal ancestry gate and
retain explicit gate metadata.
}
\description{
Conditional rates integrate every prior-pattern component within each SER
and SNP. Exact-zero component mass remains false-sign uncertainty.
For both Gaussian and MOM fits, structural support is determined from the
declared prior covariance \code{Utheta} and slot variance. Only an exactly
zero prior contrast variance is an atom; a positive variance remains
continuous, including when it is very small. Declared zero support is checked
against the posterior moments before its zero atom is assigned.
The conditional calculation uses the stable posterior negative-probability
and point-zero-mass definition implemented by \pkg{ashr}.
Single-effect rates average conditional rates over SNP posterior inclusion
probabilities. Variant rates use the most confident SER and do not construct a
cross-SER sign distribution.
For MOM-SA fits, this accessor is a post hoc calculation from the saved MOM
component audit. Its theta conditional, single-effect, and variant rates are
available. MOM observed-beta single-effect rates use the MOM second-moment-
matched normal approximation and the existing residualized Wakefield gate.
MOM observed-beta conditional and variant rates remain available. The gate
is approximate and is reported by \code{utl_ancestry_lbf}; it uses fixed prior
weights rather than posterior pattern probabilities.
}
