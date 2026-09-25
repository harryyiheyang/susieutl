\name{utl_credible_sets}
\alias{utl_credible_sets}
\title{Build multi-LD credible sets from slot posterior masses}
\usage{
utl_credible_sets(fit, R_list, coverage = 0.95, min_abs_corr = 0.5)
}
\arguments{
  \item{fit}{A \code{utl_theta_ma_ss_fit} object.}
  \item{R_list}{Named ancestry-specific LD matrices in fit ancestry order.}
  \item{coverage}{Desired cumulative slot mass.}
  \item{min_abs_corr}{Minimum multi-LD purity.}
}
\value{A list containing the credible-set table, slot masses, and filtering diagnostics.}
\description{Constructs one 95-percent-style credible set per active slot. For each variant pair, purity is the maximum absolute LD across ancestries, followed by the minimum over pairs; duplicate sets and inactive prior-variance slots are filtered.}
