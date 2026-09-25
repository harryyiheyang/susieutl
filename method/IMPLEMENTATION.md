# Implementation path

Pattern and trait priors feed the SS/RSS entry points, which call the shared
operator. The operator runs IBSS with SER updates and the optional MOM-SA
moment-based pattern-weight update. Accessors expose theta effects, PIPs,
credible sets, LFSRs, and ancestry log Bayes factors.

Native bindings connect `R/RcppExports.R` through `src/RcppExports.cpp` to the
active kernels in `src/snp_kernels.cpp`.

The default method is the fixed-weight categorical Gaussian comparator;
MOM-SA is opt-in and is the current research line. The single-Gaussian
simplex-QP and GEE routes are retired and are not published.
