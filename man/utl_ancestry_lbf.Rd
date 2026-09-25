\name{utl_ancestry_lbf}
\alias{utl_ancestry_lbf}
\title{Compute Ancestry-Specific Observed-Scale Log Bayes Factors}
\usage{
utl_ancestry_lbf(fit, bf_gate = 0)
}
\arguments{
\item{fit}{A \code{utl_theta_ma_ss_fit} object.}

\item{bf_gate}{Finite log-Bayes-factor threshold for the ancestry gate.}
}
\value{
A list containing \code{lbf}, an array indexed by slot, variant, and
ancestry; \code{lbf_outcome}, a slot-by-ancestry matrix; \code{bf_gate};
logical \code{gate}; and character \code{gate_status}.
}
\description{
For each SER and ancestry, final pattern weights are first used to form one
projected observed-beta prior variance. The function then evaluates the
univariate log Bayes factor at each SNP using the residual that excludes the
current SER. The outcome score is the SNP-posterior-weighted mean log Bayes
factor.
For MOM-SA fits, it uses the saved prior and MOM audit to form a second-moment-matched
normal approximation for the ancestry gate, followed by the residualized
Wakefield log Bayes factor and alpha-weighted outcome. The approximation uses
the fixed final prior weights and is reported with \code{gate_method}.
}
