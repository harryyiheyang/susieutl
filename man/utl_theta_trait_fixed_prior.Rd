\name{utl_theta_trait_fixed_prior}
\alias{utl_theta_trait_fixed_prior}
\title{Construct a fixed theta prior from maximum squared theta estimates}
\usage{
utl_theta_trait_fixed_prior(delta, reference, trait_variance = NULL,
  beta_hat = NULL, weights = NULL, validation_tol = 1e-10)
}
\arguments{
  \item{delta}{Named finite transfer coefficients for non-reference traits.}
  \item{reference}{One non-empty reference trait name.}
  \item{trait_variance}{Optional named non-negative theta-coordinate variances;
  when supplied they take precedence over \code{beta_hat}.}
  \item{beta_hat}{Optional finite SNP-by-ancestry effect matrix. Its columns
  are transformed to theta coordinates and the maximum squared estimate is
  used for each free coordinate.}
  \item{weights}{Optional named non-negative component weights. \code{NULL}
  gives uniform weights; supplied weights are normalized.}
  \item{validation_tol}{Positive validation tolerance.}
}
\value{An object of class \code{utl_theta_prior}, with \code{B},
\code{Vfree}, \code{Vtheta}, \code{Sigma}, and construction diagnostics.}
\description{
Constructs fixed theta-coordinate covariance components using the canonical
\code{D}, \code{C}, \code{S_*}, and \code{N_*} dictionary. For each component,
\code{Sigma_k = B_k diag(Vfree) t(B_k)} without additional rescaling. If
\code{trait_variance} is absent, \code{Vfree} is the per-coordinate maximum
squared transformed effect across variants.
}
