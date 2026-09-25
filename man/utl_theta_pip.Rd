\name{utl_theta_pip}
\alias{utl_theta_pip}
\title{Extract theta-coordinate SNP posterior inclusion probabilities}
\usage{
utl_theta_pip(fit, by_slot = FALSE)
}
\arguments{
  \item{fit}{A \code{utl_theta_ma_ss_fit} object.}
  \item{by_slot}{Whether to return slot-by-SNP probabilities.}
}
\value{A named numeric vector or slot-by-SNP matrix.}
\description{Returns posterior inclusion probabilities from a theta-coordinate fit.}
