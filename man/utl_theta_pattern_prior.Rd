\name{utl_theta_pattern_prior}
\alias{utl_theta_pattern_prior}
\title{Construct the canonical reference/departure theta pattern library}
\usage{
utl_theta_pattern_prior(delta, reference, include_AN_ref = FALSE,
  weights = NULL, validation_tol = 1e-10)
}
\arguments{
  \item{delta}{Strictly named finite transfer coefficients for non-reference traits, in the desired trait order.}
  \item{reference}{One non-empty reference trait name not present in \code{names(delta)}.}
  \item{include_AN_ref}{Whether to request the reference-null pattern.}
  \item{weights}{Optional named non-negative weights aligned exactly to the final unique non-null components and summing to one. \code{NULL} gives uniform weights.}
  \item{validation_tol}{One finite positive validation and duplicate-detection tolerance.}
}
\value{An object of classes \code{utl_theta_pattern_prior} and \code{utl_theta_prior}, including raw unique bases, covariance components, taxonomy, aliases, diagnostics, and implicit zero-SER metadata.}
\description{Builds raw, unnormalized basis matrices for consistent, unrestricted, non-reference-specific, and ancestry-null scientific patterns. Duplicate covariance spans become aliases; \code{SER_V0} remains metadata and does not enter \code{Utheta} or its weights.}
