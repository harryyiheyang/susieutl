\name{utl_theta_prior}
\alias{utl_theta_prior}
\title{Construct a fixed theta-coordinate covariance prior}
\usage{
utl_theta_prior(delta, M, Utheta, weights, validation_tol = 1e-10,
  allow_zero_rank = FALSE)
}
\arguments{
  \item{delta}{Named finite transfer coefficients for non-reference ancestries.}
  \item{M}{Named theta-coordinate map.}
  \item{Utheta}{Named list of theta covariance matrices.}
  \item{weights}{Named fixed component weights summing to one.}
  \item{validation_tol}{Positive relative validation tolerance.}
  \item{allow_zero_rank}{Logical compatibility option allowing a zero-rank
  covariance component for the exact fallback route. Defaults to code{FALSE}.}
}
\value{An object of class \code{utl_theta_prior}.}
\description{Validates a canonical theta-coordinate map and fixed positive-semidefinite covariance components. For automatically generated reference/departure pattern libraries, see \code{utl_theta_pattern_prior()}. Manual \code{M}/\code{Utheta} construction remains fully supported.}
