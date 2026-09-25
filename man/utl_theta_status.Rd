\name{utl_theta_status}
\alias{utl_theta_status}
\title{Extract theta-coordinate fit diagnostics}
\usage{
utl_theta_status(fit)
}
\arguments{
  \item{fit}{A \code{utl_theta_ma_ss_fit} object.}
}
\value{A list of convergence, objective, scale, backend, mixture-weight, and
optional MOM-SA diagnostics.}
\description{Returns status information for a theta-coordinate sufficient-statistics fit.
Categorical \code{weights*} fields retain their fixed-prior meaning unless
prior mixture-weight estimation is enabled. For formal MOM-SA fits,
\code{weight_model},
\code{estimate_pattern_weight}, \code{pattern_weight_method},
\code{prior_weight_estimation}, \code{mom_sa_control_source}, and
\code{mom_sa_slot} describe the fixed-prior q1 estimator.}
