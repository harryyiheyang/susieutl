\name{utl_theta_pattern_posterior}
\alias{utl_theta_pattern_posterior}
\title{Extract theta prior-component posterior probabilities}
\usage{
utl_theta_pattern_posterior(fit, by_slot = FALSE)
}
\arguments{
  \item{fit}{A \code{utl_theta_ma_ss_fit} object.}
  \item{by_slot}{Whether to return slot-by-component probabilities.}
}
\value{A named numeric vector or slot-by-component matrix.}
\description{Returns posterior probabilities for fixed theta covariance components.}
