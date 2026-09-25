#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

using arma::mat;
using arma::vec;

inline int worker_count() {
#ifdef _OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

std::vector<mat> as_matrices(const Rcpp::List& x, const char* what) {
  std::vector<mat> out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    Rcpp::NumericMatrix z = Rcpp::as<Rcpp::NumericMatrix>(x[i]);
    out[i] = mat(z.begin(), z.nrow(), z.ncol(), false, true);
    if (!out[i].is_finite()) {
      Rcpp::stop(std::string(what) + " contains non-finite values");
    }
  }
  return out;
}

mat copy_numeric_matrix(const Rcpp::NumericMatrix& x) {
  mat out(x.nrow(), x.ncol());
  std::copy(x.begin(), x.end(), out.begin());
  return out;
}

void check_same_square(const std::vector<mat>& x, arma::uword d,
                       const char* what) {
  for (const mat& z : x) {
    if (z.n_rows != d || z.n_cols != d) {
      Rcpp::stop(std::string(what) + " has incompatible dimensions");
    }
  }
}

bool chol_upper(const mat& x, mat& r) {
  return arma::chol(r, (x + x.t()) * 0.5);
}

struct ComponentPost {
  double lbf;
  vec mean;
  mat cov;
  bool ok;
  ComponentPost() : lbf(0.0), ok(true) {}
};

ComponentPost one_exact(const vec& score, const mat& P, const mat& F,
                        const vec& v, bool want_moments) {
  ComponentPost out;
  const arma::uword d = score.n_elem;
  const arma::uword r = F.n_cols;
  out.mean.zeros(d);
  out.cov.zeros(d, d);
  if (r == 0 || arma::accu(v) == 0.0) return out;
  mat G = F;
  for (arma::uword i = 0; i < d; ++i) G.row(i) *= std::sqrt(v(i));
  mat K = arma::eye<mat>(r, r) + G.t() * P * G;
  mat R;
  if (!chol_upper(K, R)) {
    out.ok = false;
    return out;
  }
  vec b = G.t() * score;
  vec kb = arma::solve(arma::trimatu(R),
                       arma::solve(arma::trimatl(R.t()), b));
  out.lbf = 0.5 * (arma::dot(b, kb) -
                   2.0 * arma::accu(arma::log(R.diag())));
  if (!std::isfinite(out.lbf)) {
    out.ok = false;
    return out;
  }
  if (!want_moments) return out;
  vec mu = G * kb;
  mat Rt = arma::solve(arma::trimatu(R),
                       arma::solve(arma::trimatl(R.t()), G.t()));
  mat C = G * Rt;
  out.mean = mu;
  out.cov = (C + C.t()) * 0.5;
  if (!out.mean.is_finite() || !out.cov.is_finite()) out.ok = false;
  return out;
}

inline double logsumexp_ordered(const std::vector<double>& x) {
  double z = -std::numeric_limits<double>::infinity();
  for (double a : x) z = std::max(z, a);
  double s = 0.0;
  for (double a : x) s += std::exp(a - z);
  return z + std::log(s);
}

Rcpp::List make_ser_result(const std::vector<double>& lbf,
                           const std::vector<double>& cm,
                           const std::vector<double>& cc,
                           const Rcpp::NumericMatrix& score,
                           const std::vector<mat>& P,
                           const Rcpp::NumericVector& weights,
                           bool moments,
                           int topcount) {
  const arma::uword p = score.nrow();
  const arma::uword d = score.ncol();
  const arma::uword K = weights.size();
  Rcpp::NumericMatrix out_lbf(p, K);
  Rcpp::NumericMatrix out_rho(p, K);
  Rcpp::NumericVector out_pip(p), out_pat(K);
  std::vector<double> lp(p * K);
  std::vector<double> full_rho(p * K, 0.0), full_pip(p, 0.0);
  for (arma::uword k = 0; k < K; ++k) {
    if (!(weights[k] >= 0.0) || !std::isfinite(weights[k]))
      Rcpp::stop("weights must be finite and non-negative");
  }
  double wsum = 0.0;
  for (double w : weights) wsum += w;
  if (!std::isfinite(wsum) || wsum <= 0.0)
    Rcpp::stop("weights must have positive sum");
  for (arma::uword k = 0; k < K; ++k) {
    for (arma::uword j = 0; j < p; ++j) {
      const double z = lbf[j + p * k];
      out_lbf(j, k) = z;
      lp[j + p * k] = weights[k] > 0.0 ?
        z + std::log(weights[k]) - std::log((double)p) :
        -std::numeric_limits<double>::infinity();
    }
  }
  const double logz = logsumexp_ordered(lp);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < K; ++k) {
      const double r = std::exp(lp[j + p * k] - logz);
      full_rho[j + p * k] = r;
      full_pip[j] += r;
    }
  }
  double correction = 0.0;
  if (topcount > 0 && topcount < (int)K) {
    const int top = std::min(topcount, (int)K);
    std::vector<int> order(K);
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> eta(K), selected(top);
    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword k = 0; k < K; ++k)
        eta[k] = lbf[j + p * k] + std::log(weights[k]);
      std::stable_sort(order.begin(), order.end(),
                       [&eta](int a, int b) { return eta[a] > eta[b]; });
      const double log_full = logsumexp_ordered(eta);
      for (int q = 0; q < top; ++q) selected[q] = eta[order[q]];
      const double log_selected = logsumexp_ordered(selected);
      const double logr = log_selected - log_full;
      correction -= full_pip[j] * logr;
      for (arma::uword k = 0; k < K; ++k) out_rho(j, k) = 0.0;
      for (int q = 0; q < top; ++q) {
        const int k = order[q];
        out_rho(j, k) = full_pip[j] *
          std::exp(eta[k] - log_selected);
        out_pat[k] += out_rho(j, k);
      }
      out_pip[j] = full_pip[j];
    }
  } else {
    for (arma::uword j = 0; j < p; ++j) {
      out_pip[j] = full_pip[j];
      for (arma::uword k = 0; k < K; ++k) {
        out_rho(j, k) = full_rho[j + p * k];
        out_pat[k] += out_rho(j, k);
      }
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("lbf") = out_lbf,
    Rcpp::Named("rho") = out_rho,
    Rcpp::Named("logZ") = logz,
    Rcpp::Named("pattern_posterior") = out_pat,
    Rcpp::Named("pip") = out_pip,
    Rcpp::Named("weights") = weights,
    Rcpp::Named("correction") = correction);
  if (!moments) return out;
  Rcpp::NumericVector mean(p * d), second(p * d * d);
  Rcpp::NumericVector cmean(cm.size()), ccov(cc.size());
  std::copy(cm.begin(), cm.end(), cmean.begin());
  std::copy(cc.begin(), cc.end(), ccov.begin());
  double e = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < K; ++k) {
      const double r = out_rho(j, k);
      for (arma::uword q = 0; q < d; ++q) {
        const double mu = cmean[j + p * (k + K * q)];
        mean[j + p * q] += r * mu;
        for (arma::uword h = 0; h < d; ++h) {
          const double cv = ccov[j + p * (k + K * (q + d * h))];
          second[j + p * (q + d * h)] +=
            r * (cv + mu * cmean[j + p * (k + K * h)]);
        }
      }
      mat T(d, d);
      vec mu(d);
      for (arma::uword q = 0; q < d; ++q)
        mu(q) = cmean[j + p * (k + K * q)];
      for (arma::uword q = 0; q < d; ++q) {
        for (arma::uword h = 0; h < d; ++h)
          T(q, h) = ccov[j + p * (k + K * (q + d * h))] + mu(q) * mu(h);
      }
      vec sj(d);
      for (arma::uword q = 0; q < d; ++q) sj(q) = score(j, q);
      e += r * (arma::dot(sj, mu) - 0.5 * arma::accu(P[j] % T));
    }
  }
  out["mean"] = mean;
  out["second"] = second;
  out["e"] = e;
  out["KL"] = e - logz + correction;
  out["component_mean"] = cmean;
  out["component_covariance"] = ccov;
  return out;
}

} // namespace
// [[Rcpp::export]]
Rcpp::List utl_cpp_ser(const Rcpp::NumericMatrix& score,
                       const Rcpp::List& Pj,
                       const Rcpp::List& Flist,
                       const Rcpp::NumericVector& V,
                       const Rcpp::NumericVector& weights,
                       bool moments = true,
                       int topcount = 0) {
  const arma::uword p = score.nrow();
  const arma::uword d = score.ncol();
  const arma::uword K = Flist.size();
  if (p < 1 || d < 1 || Pj.size() != (R_xlen_t)p || K < 1 ||
      (V.size() != 1 && V.size() != (R_xlen_t)d) ||
      weights.size() != (R_xlen_t)K) {
    Rcpp::stop("C++ SER dimensions are incompatible");
  }
  std::vector<mat> P = as_matrices(Pj, "precision blocks");
  std::vector<mat> F = as_matrices(Flist, "component factors");
  check_same_square(P, d, "precision blocks");
  for (const mat& z : F) if (z.n_rows != d)
    Rcpp::stop("component factors have incompatible dimensions");
  vec vv(V.begin(), V.size());
  if (!vv.is_finite() || arma::any(vv < 0.0))
    Rcpp::stop("V must be finite and non-negative");
  if (vv.n_elem == 1) vv = arma::vec(d, arma::fill::value(vv(0)));
  if (!score.begin() || !std::all_of(score.begin(), score.end(),
                                      [](double x){ return std::isfinite(x); }))
    Rcpp::stop("score must be finite");
  std::vector<double> lbf(p * K), cm(p * K * d, 0.0),
    cc(p * K * d * d, 0.0);
  std::vector<int> bad(p, 0);
  mat sc(p, d);
  std::copy(score.begin(), score.end(), sc.begin());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int jj = 0; jj < (int)p; ++jj) {
    try {
      vec sj(d);
      for (arma::uword q = 0; q < d; ++q) sj(q) = sc(jj, q);
      for (arma::uword k = 0; k < K; ++k) {
        ComponentPost z = one_exact(sj, P[jj], F[k], vv, moments);
        if (!z.ok) { bad[jj] = 1; continue; }
        lbf[jj + p * k] = z.lbf;
        if (moments) {
          for (arma::uword q = 0; q < d; ++q) {
            cm[jj + p * (k + K * q)] = z.mean(q);
            for (arma::uword h = 0; h < d; ++h)
              cc[jj + p * (k + K * (q + d * h))] = z.cov(q, h);
          }
        }
      }
    } catch (...) {
      bad[jj] = 1;
    }
  }
  for (int x : bad) if (x) Rcpp::stop("C++ SER GLS factorization failed");
  Rcpp::List out = make_ser_result(lbf, cm, cc, score, P, weights, moments,
                                   topcount);
  out["omp_threads"] = worker_count();
  return out;
}

// Component posterior moments without recomputing lbf or the joint mixture.
// The SNP loop is owned by OpenMP and the R layer only reshapes the result.
// [[Rcpp::export]]
Rcpp::List utl_cpp_component_moments(const Rcpp::NumericMatrix& score,
                                     const Rcpp::List& Pj,
                                     const Rcpp::List& Flist,
                                     const Rcpp::NumericVector& V) {
  const int p = score.nrow(), d = score.ncol(), K = Flist.size();
  if (p < 1 || d < 1 || Pj.size() != p || K < 1 ||
      (V.size() != 1 && V.size() != d))
    Rcpp::stop("C++ component moment dimensions are incompatible");
  std::vector<mat> P = as_matrices(Pj, "precision blocks");
  std::vector<mat> F = as_matrices(Flist, "component factors");
  check_same_square(P, d, "precision blocks");
  for (const mat& z : F) if (z.n_rows != (arma::uword)d)
    Rcpp::stop("component factors have incompatible dimensions");
  vec vv(V.begin(), V.size());
  if (!vv.is_finite() || arma::any(vv < 0.0))
    Rcpp::stop("V must be finite and non-negative");
  if (vv.n_elem == 1) vv = arma::vec(d, arma::fill::value(vv(0)));
  std::vector<double> cm((std::size_t)p*K*d, 0.0),
    cc((std::size_t)p*K*d*d, 0.0);
  std::vector<int> bad(p, 0);
  mat sc(p, d);
  std::copy(score.begin(), score.end(), sc.begin());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int jj = 0; jj < p; ++jj) {
    try {
      vec sj(d);
      for (int q = 0; q < d; ++q) sj(q) = sc(jj, q);
      for (int k = 0; k < K; ++k) {
        if (F[k].n_cols == 0 || arma::accu(vv) == 0.0) continue;
        mat G = F[k];
        for (int q = 0; q < d; ++q) G.row(q) *= std::sqrt(vv(q));
        mat Kmat = arma::eye<mat>(G.n_cols, G.n_cols) + G.t() * P[jj] * G;
        mat R;
        if (!chol_upper(Kmat, R)) { bad[jj] = 1; continue; }
        vec b = G.t() * sj;
        vec kb = arma::solve(arma::trimatu(R),
                             arma::solve(arma::trimatl(R.t()), b));
        mat KinvGt = arma::solve(arma::trimatu(R),
                                 arma::solve(arma::trimatl(R.t()), G.t()));
        vec mu = G * kb;
        mat C = (G * KinvGt + (G * KinvGt).t()) * 0.5;
        if (!mu.is_finite() || !C.is_finite()) { bad[jj] = 1; continue; }
        for (int q = 0; q < d; ++q) {
          cm[jj + p * (k + K*q)] = mu(q);
          for (int r = 0; r < d; ++r)
            cc[jj + p * (k + K*(q + d*r))] = C(q,r);
        }
      }
    } catch (...) { bad[jj] = 1; }
  }
  for (int x : bad) if (x) Rcpp::stop("C++ component posterior moments failed");
  return Rcpp::List::create(Rcpp::Named("component_mean") =
                              Rcpp::wrap(cm),
                            Rcpp::Named("component_covariance") =
                              Rcpp::wrap(cc),
                            Rcpp::Named("omp_threads") = worker_count());
}

// Fixed-order joint SER reducer used by the native cached route.
// [[Rcpp::export]]
Rcpp::List utl_cpp_ser_reduce(const Rcpp::NumericMatrix& score,
                              const Rcpp::List& Pj,
                              const Rcpp::NumericMatrix& lbf,
                              const Rcpp::NumericVector& component_mean,
                              const Rcpp::NumericVector& component_covariance,
                              const Rcpp::NumericVector& weights,
                              bool moments = true,
                              int topcount = 0) {
  const int p = score.nrow(), d = score.ncol(), K = weights.size();
  if (lbf.nrow() != p || lbf.ncol() != K ||
      (moments && (component_mean.size() != (R_xlen_t)p*K*d ||
                   component_covariance.size() != (R_xlen_t)p*K*d*d)))
    Rcpp::stop("C++ SER reducer dimensions are incompatible");
  std::vector<mat> P = as_matrices(Pj, "precision blocks");
  check_same_square(P, d, "precision blocks");
  std::vector<double> ll((std::size_t)p*K), cm(component_mean.size()),
    cc(component_covariance.size());
  for (int k = 0; k < K; ++k) for (int j = 0; j < p; ++j)
    ll[j + p*k] = lbf(j,k);
  if (moments) {
    std::copy(component_mean.begin(), component_mean.end(), cm.begin());
    std::copy(component_covariance.begin(), component_covariance.end(), cc.begin());
  }
  return make_ser_result(ll, cm, cc, score, P, weights, moments, topcount);
}

// Precompute the common precision posterior algebra for a fixed RSS fit.
// The score-dependent work remains in the OpenMP kernels below.
// [[Rcpp::export]]
Rcpp::List utl_cpp_common_ser_context(const Rcpp::NumericMatrix& P0,
                                      const Rcpp::List& Flist,
                                      const Rcpp::NumericVector& V,
                                      bool moments = true) {
  mat P = copy_numeric_matrix(P0);
  std::vector<mat> F = as_matrices(Flist, "component factors");
  const int d = P.n_rows;
  const int K = F.size();
  if (P.n_cols != (arma::uword)d || d < 1 || K < 1 ||
      (V.size() != 1 && V.size() != d)) {
    Rcpp::stop("common SER context dimensions are incompatible");
  }
  for (const mat& z : F) if (z.n_rows != (arma::uword)d)
    Rcpp::stop("component factors have incompatible dimensions");
  vec vv(V.begin(), V.size());
  if (!vv.is_finite() || arma::any(vv < 0.0))
    Rcpp::stop("V must be finite and non-negative");
  if (vv.n_elem == 1) vv = arma::vec(d, arma::fill::value(vv(0)));

  Rcpp::List components(K);
  bool ok = true;
  for (int k = 0; k < K; ++k) {
    mat G = F[k];
    for (int q = 0; q < d; ++q) G.row(q) *= std::sqrt(vv(q));
    Rcpp::List one = Rcpp::List::create(
      Rcpp::Named("G") = Rcpp::wrap(G),
      Rcpp::Named("rank") = (int)G.n_cols);
    if (G.n_cols > 0 && arma::accu(vv) != 0.0) {
      mat Kmat = arma::eye<mat>(G.n_cols, G.n_cols) + G.t() * P * G;
      mat R;
      if (!chol_upper(Kmat, R)) {
        ok = false;
        components[k] = one;
        continue;
      }
      one["R"] = Rcpp::wrap(R);
      one["logdet"] = 2.0 * arma::accu(arma::log(R.diag()));
      if (moments) {
        mat KinvGt = arma::solve(arma::trimatu(R),
                                 arma::solve(arma::trimatl(R.t()), G.t()));
        one["C"] = Rcpp::wrap((G * KinvGt + (G * KinvGt).t()) * 0.5);
      }
    } else {
      one["R"] = Rcpp::wrap(arma::eye<mat>(G.n_cols, G.n_cols));
      one["logdet"] = 0.0;
      if (moments) one["C"] = Rcpp::wrap(arma::zeros<mat>(d, d));
    }
    components[k] = one;
  }
  return Rcpp::List::create(Rcpp::Named("ok") = ok,
                            Rcpp::Named("moments") = moments,
                            Rcpp::Named("components") = components,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("K") = K,
                            Rcpp::Named("omp_threads") = worker_count());
}

// Evaluate common-P component log Bayes factors with fixed factorizations.
// [[Rcpp::export]]
Rcpp::List utl_cpp_common_ser_lbf(const Rcpp::NumericMatrix& score,
                                  const Rcpp::List& context) {
  const int p = score.nrow();
  const int d = score.ncol();
  Rcpp::List components = context["components"];
  const int K = components.size();
  std::vector<mat> G(K), R(K);
  std::vector<double> logdet(K, 0.0);
  for (int k = 0; k < K; ++k) {
    Rcpp::List one = components[k];
    G[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["G"]));
    R[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["R"]));
    logdet[k] = Rcpp::as<double>(one["logdet"]);
  }
  std::vector<double> sc(score.begin(), score.end());
  std::vector<double> out((std::size_t)p * K, 0.0);
  std::vector<int> bad(p, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    try {
      vec sj(d);
      for (int q = 0; q < d; ++q) sj(q) = sc[j + p * q];
      for (int k = 0; k < K; ++k) {
        if (G[k].n_cols == 0) continue;
        vec b = G[k].t() * sj;
        vec kb = arma::solve(arma::trimatu(R[k]),
                             arma::solve(arma::trimatl(R[k].t()), b));
        out[j + p * k] = 0.5 * (arma::dot(b, kb) - logdet[k]);
      }
    } catch (...) { bad[j] = 1; }
  }
  for (int x : bad) if (x) Rcpp::stop("common SER log Bayes factor failed");
  Rcpp::NumericMatrix lbf(p, K);
  std::copy(out.begin(), out.end(), lbf.begin());
  return Rcpp::List::create(Rcpp::Named("lbf") = lbf,
                            Rcpp::Named("omp_threads") = worker_count());
}

// Evaluate common-P component posterior moments with fixed factorizations.
// [[Rcpp::export]]
Rcpp::List utl_cpp_common_component_moments(const Rcpp::NumericMatrix& score,
                                            const Rcpp::List& context) {
  const int p = score.nrow();
  const int d = score.ncol();
  Rcpp::List components = context["components"];
  const int K = components.size();
  std::vector<mat> G(K), R(K), C(K);
  for (int k = 0; k < K; ++k) {
    Rcpp::List one = components[k];
    G[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["G"]));
    R[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["R"]));
    C[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["C"]));
  }
  std::vector<double> sc(score.begin(), score.end());
  std::vector<double> cm((std::size_t)p * K * d, 0.0);
  std::vector<double> cc((std::size_t)p * K * d * d, 0.0);
  std::vector<int> bad(p, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    try {
      vec sj(d);
      for (int q = 0; q < d; ++q) sj(q) = sc[j + p * q];
      for (int k = 0; k < K; ++k) {
        if (G[k].n_cols == 0) continue;
        vec b = G[k].t() * sj;
        vec kb = arma::solve(arma::trimatu(R[k]),
                             arma::solve(arma::trimatl(R[k].t()), b));
        vec mu = G[k] * kb;
        if (!mu.is_finite() || !C[k].is_finite()) { bad[j] = 1; continue; }
        for (int q = 0; q < d; ++q) {
          cm[j + p * (k + K * q)] = mu(q);
          for (int h = 0; h < d; ++h)
            cc[j + p * (k + K * (q + d * h))] = C[k](q, h);
        }
      }
    } catch (...) { bad[j] = 1; }
  }
  for (int x : bad) if (x) Rcpp::stop("common SER component moments failed");
  Rcpp::NumericVector rmean = Rcpp::wrap(cm);
  Rcpp::NumericVector rcov = Rcpp::wrap(cc);
  rmean.attr("dim") = Rcpp::IntegerVector::create(p, K, d);
  rcov.attr("dim") = Rcpp::IntegerVector::create(p, K, d, d);
  return Rcpp::List::create(Rcpp::Named("component_mean") = rmean,
                            Rcpp::Named("component_covariance") = rcov,
                            Rcpp::Named("omp_threads") = worker_count());
}

// Evaluate common-P log Bayes factors and component posterior moments together.
// [[Rcpp::export]]
Rcpp::List utl_cpp_common_ser_lbf_moments(const Rcpp::NumericMatrix& score,
                                          const Rcpp::List& context) {
  const int p = score.nrow();
  const int d = score.ncol();
  Rcpp::List components = context["components"];
  const int K = components.size();
  std::vector<mat> G(K), R(K), C(K);
  std::vector<double> logdet(K, 0.0);
  for (int k = 0; k < K; ++k) {
    Rcpp::List one = components[k];
    G[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["G"]));
    R[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["R"]));
    C[k] = copy_numeric_matrix(Rcpp::as<Rcpp::NumericMatrix>(one["C"]));
    logdet[k] = Rcpp::as<double>(one["logdet"]);
  }
  std::vector<double> sc(score.begin(), score.end());
  std::vector<double> out((std::size_t)p * K, 0.0);
  std::vector<double> cm((std::size_t)p * K * d, 0.0);
  std::vector<double> cc((std::size_t)p * K * d * d, 0.0);
  std::vector<int> bad(p, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int j = 0; j < p; ++j) {
    try {
      vec sj(d);
      for (int q = 0; q < d; ++q) sj(q) = sc[j + p * q];
      for (int k = 0; k < K; ++k) {
        if (G[k].n_cols == 0) continue;
        vec b = G[k].t() * sj;
        vec kb = arma::solve(arma::trimatu(R[k]),
                             arma::solve(arma::trimatl(R[k].t()), b));
        out[j + p * k] = 0.5 * (arma::dot(b, kb) - logdet[k]);
        vec mu = G[k] * kb;
        if (!mu.is_finite() || !C[k].is_finite()) { bad[j] = 1; continue; }
        for (int q = 0; q < d; ++q) {
          cm[j + p * (k + K*q)] = mu(q);
          for (int h = 0; h < d; ++h)
            cc[j + p * (k + K * (q + d*h))] = C[k](q, h);
        }
      }
    } catch (...) { bad[j] = 1; }
  }
  for (int x : bad) if (x) Rcpp::stop("common SER lbf and moments failed");
  Rcpp::NumericMatrix lbf(p, K);
  std::copy(out.begin(), out.end(), lbf.begin());
  Rcpp::NumericVector rmean = Rcpp::wrap(cm);
  Rcpp::NumericVector rcov = Rcpp::wrap(cc);
  rmean.attr("dim") = Rcpp::IntegerVector::create(p, K, d);
  rcov.attr("dim") = Rcpp::IntegerVector::create(p, K, d, d);
  return Rcpp::List::create(Rcpp::Named("lbf") = lbf,
                            Rcpp::Named("component_mean") = rmean,
                            Rcpp::Named("component_covariance") = rcov,
                            Rcpp::Named("omp_threads") = worker_count());
}

// [[Rcpp::export]]
Rcpp::List utl_cpp_v_cache(const Rcpp::NumericMatrix& M,
                           const Rcpp::NumericMatrix& Bk,
                           const Rcpp::NumericMatrix& bhat,
                           const Rcpp::NumericMatrix& N) {
  const int p=bhat.nrow(), T=bhat.ncol(), d=M.ncol(), r=Bk.ncol();
  if (M.nrow()!=T || Bk.nrow()!=d || N.nrow()!=p || N.ncol()!=T)
    Rcpp::stop("V cache dimensions are incompatible");
  mat mm = copy_numeric_matrix(M);
  mat bb = copy_numeric_matrix(Bk);
  mat yy = copy_numeric_matrix(bhat);
  mat nn = copy_numeric_matrix(N);
  if (!mm.is_finite() || !bb.is_finite() || !yy.is_finite() || !nn.is_finite() ||
      arma::any(arma::vectorise(nn) <= 0.0))
    Rcpp::stop("V cache inputs must be finite with positive sample sizes");
  mat Ak=mm*bb;
  bool constant=true;
  for(int j=1;j<p && constant;++j) for(int t=0;t<T;++t) if(nn(j,t)!=nn(0,t)) {constant=false;break;}
  std::vector<double> sv((std::size_t)p*r, 0.0), uv((std::size_t)p*r, 0.0);
  if(constant) {
    mat H=Ak.t()*arma::diagmat(nn.row(0))*Ak;
    mat S;
    if(!arma::inv_sympd(S,H)) Rcpp::stop("V cache GLS system is singular");
    vec ev; mat P;
    if(!arma::eig_sym(ev,P,S)) Rcpp::stop("V cache eigendecomposition failed");
    ev=arma::reverse(ev); P=arma::fliplr(P);
    mat q=S*Ak.t()*arma::diagmat(nn.row(0))*yy.t();
    mat uu=q.t()*P;
    for(int j=0;j<p;++j) for(int k=0;k<r;++k){sv[j+p*k]=ev(k);uv[j+p*k]=uu(j,k);}
  } else {
    std::vector<int> bad(p, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int j=0;j<p;++j){
      try {
        mat H=Ak.t()*arma::diagmat(nn.row(j))*Ak, S;
        if(!arma::inv_sympd(S,H)) { bad[j] = 1; continue; }
        vec ev; mat P;
        if(!arma::eig_sym(ev,P,S)) { bad[j] = 1; continue; }
        ev=arma::reverse(ev); P=arma::fliplr(P);
        vec g=S*Ak.t()*arma::diagmat(nn.row(j))*yy.row(j).t();
        vec uu=P.t()*g;
        if (!ev.is_finite() || !uu.is_finite()) { bad[j] = 1; continue; }
        for(int k=0;k<r;++k){sv[j+p*k]=ev(k);uv[j+p*k]=uu(k);}
      } catch (...) { bad[j] = 1; }
    }
    for (int x : bad) if (x)
      Rcpp::stop("V cache GLS system is singular or eigendecomposition failed");
  }
  Rcpp::NumericMatrix s(p,r), u(p,r);
  std::copy(sv.begin(), sv.end(), s.begin());
  std::copy(uv.begin(), uv.end(), u.begin());
  return Rcpp::List::create(Rcpp::Named("N_constant")=constant,
                            Rcpp::Named("s")=s,Rcpp::Named("u")=u,
                            Rcpp::Named("rank")=r,
                            Rcpp::Named("omp_threads")=worker_count());
}

// [[Rcpp::export]]
Rcpp::NumericVector utl_cpp_v_loglik(const Rcpp::NumericMatrix& s,
                                     const Rcpp::NumericMatrix& u,
                                     double V) {
  if (!std::isfinite(V) || V < 0.0 || s.nrow()!=u.nrow() || s.ncol()!=u.ncol())
    Rcpp::stop("V log likelihood inputs are invalid");
  const int p=s.nrow(), r=s.ncol();
  std::vector<double> sv(s.begin(), s.end()), uv(u.begin(), u.end());
  std::vector<double> lbf(p), dlbf(p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int j=0;j<p;++j){
    double z=0.0,dz=0.0;
    for(int k=0;k<r;++k){
      double sk=sv[j + p*k], uk=uv[j + p*k], sp=sk+V;
      if(!(sk>0.0) || !(sp>0.0) || !std::isfinite(sk) || !std::isfinite(uk)) continue;
      z += std::log(sk)-std::log(sp)+uk*uk*(1.0/sk-1.0/sp);
      dz += -1.0/sp+uk*uk/(sp*sp);
    }
    lbf[j]=0.5*z; dlbf[j]=0.5*dz;
  }
  double mx=-std::numeric_limits<double>::infinity();
  for(double z:lbf) mx=std::max(mx,z-std::log((double)p));
  double den=0.0, num=0.0;
  for(int j=0;j<p;++j){double a=std::exp(lbf[j]-std::log((double)p)-mx);den+=a;num+=a*dlbf[j];}
  return Rcpp::NumericVector::create(mx+std::log(den),num/den);
}

// Per-trait initial V estimates, with the SNP reduction performed outside R.
// [[Rcpp::export]]
Rcpp::NumericVector utl_cpp_vraw(const Rcpp::NumericMatrix& bhat,
                                  const Rcpp::NumericMatrix& N) {
  if (bhat.nrow() < 1 || bhat.ncol() < 1 ||
      bhat.nrow() != N.nrow() || bhat.ncol() != N.ncol())
    Rcpp::stop("V initial estimate dimensions are incompatible");
  const int p = bhat.nrow(), T = bhat.ncol();
  std::vector<double> bv(bhat.begin(), bhat.end()), nv(N.begin(), N.end());
  std::vector<double> out(T, -std::numeric_limits<double>::infinity());
  std::vector<int> bad(T, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int t = 0; t < T; ++t) {
    double z = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < p; ++j) {
      const double n = nv[j + p*t], b = bv[j + p*t];
      if (!std::isfinite(n) || n <= 0.0 || !std::isfinite(b)) {
        bad[t] = 1; continue;
      }
      z = std::max(z, b * b - 1.0 / n);
    }
    out[t] = z;
  }
  for (int x : bad) if (x) Rcpp::stop("V initial estimate inputs are invalid");
  return Rcpp::wrap(out);
}
