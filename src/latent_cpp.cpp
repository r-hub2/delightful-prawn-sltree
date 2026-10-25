#include <Rcpp.h>
#include <cmath>
#include <float.h>
#include <string>
#include <vector>
using namespace Rcpp;

static inline double Phi(double x) { return R::pnorm(x, 0.0, 1.0, 1, 0); }

// Gauss-Legendre nodes and weights for bvnd().
static const double GL3_X[3] = {
  0.9324695142031522, 0.6612093864662647, 0.2386191860831970
};
static const double GL3_W[3] = {
  0.1713244923791705, 0.3607615730481384, 0.4679139345726904
};
static const double GL6_X[6] = {
  0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
  0.5873179542866171, 0.3678314989981802, 0.1252334085114692
};
static const double GL6_W[6] = {
  0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
  0.2031674267230659, 0.2334925365383547, 0.2491470458134029
};
static const double GL10_X[10] = {
  0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
  0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
  0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
  0.07652652113349733
};
static const double GL10_W[10] = {
  0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
  0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
  0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
  0.1527533871307259
};

// P(X > dh, Y > dk) for a standard bivariate normal, translated from BVND
// in TVPACK by Alan Genz (Statistics and Computing 14, 251-260, 2004).
// Copyright (C) 2013, Alan Genz. See inst/COPYRIGHTS.
static double bvnd(double dh, double dk, double r) {
  const double TWOPI = 6.283185307179586;
  const double* x;
  const double* w;
  int lg;
  if (std::fabs(r) < 0.3) {
    lg = 3; x = GL3_X; w = GL3_W;
  } else if (std::fabs(r) < 0.75) {
    lg = 6; x = GL6_X; w = GL6_W;
  } else {
    lg = 10; x = GL10_X; w = GL10_W;
  }

  double h = dh, k = dk, hk = h * k, bvn = 0.0;
  if (std::fabs(r) < 0.925) {
    double hs = (h * h + k * k) / 2.0, asr = std::asin(r);
    for (int i = 0; i < lg; i++) {
      for (int is = -1; is <= 1; is += 2) {
        double sn = std::sin(asr * (is * x[i] + 1.0) / 2.0);
        bvn += w[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
      }
    }
    bvn = bvn * asr / (2.0 * TWOPI) + Phi(-h) * Phi(-k);
  } else {
    if (r < 0) { k = -k; hk = -hk; }
    if (std::fabs(r) < 1.0) {
      double as = (1.0 - r) * (1.0 + r), a = std::sqrt(as);
      double bs = (h - k) * (h - k);
      double c = (4.0 - hk) / 8.0, d = (12.0 - hk) / 16.0;
      double asr = -(bs / as + hk) / 2.0;
      if (asr > -100.0)
        bvn = a * std::exp(asr) *
              (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);
      if (-hk < 100.0) {
        double b = std::sqrt(bs);
        bvn -= std::exp(-hk / 2.0) * std::sqrt(TWOPI) * Phi(-b / a) * b *
               (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }
      a = a / 2.0;
      for (int i = 0; i < lg; i++) {
        for (int is = -1; is <= 1; is += 2) {
          double xs = a * (is * x[i] + 1.0); xs = xs * xs;
          double rs = std::sqrt(1.0 - xs);
          asr = -(bs / xs + hk) / 2.0;
          if (asr > -100.0) {
            bvn += a * w[i] * std::exp(asr) *
                   (std::exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
                    (1.0 + c * xs * (1.0 + d * xs)));
          }
        }
      }
      bvn = -bvn / TWOPI;
    }
    if (r > 0) {
      bvn += Phi(-std::fmax(h, k));
    } else {
      bvn = -bvn;
      if (k > h) bvn += Phi(k) - Phi(h);
    }
  }
  return bvn;
}

// Lower-orthant bivariate normal CDF P(X <= a, Y <= b).
static inline double Phi2(double a, double b, double r) {
  return bvnd(-a, -b, r);
}

// Ported from Brent_fmin() in R's optimize.c. See inst/COPYRIGHTS.
template <typename F>
static double Brent_fmin(double ax, double bx, F f, double tol) {
  const double cc = (3.0 - std::sqrt(5.0)) * 0.5;
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  eps = DBL_EPSILON;
  tol1 = eps + 1.0;
  eps = std::sqrt(eps);

  a = ax; b = bx; v = a + cc * (b - a); w = v; x = v;
  d = 0.0; e = 0.0;
  fx = f(x); fv = fx; fw = fx;
  tol3 = tol / 3.0;

  for (;;) {
    xm = (a + b) * 0.5;
    tol1 = eps * std::fabs(x) + tol3;
    t2 = tol1 * 2.0;
    if (std::fabs(x - xm) <= t2 - (b - a) * 0.5) break;
    p = 0.0; q = 0.0; r = 0.0;
    if (std::fabs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.0;
      if (q > 0.0) p = -p; else q = -q;
      r = e; e = d;
    }
    if (std::fabs(p) >= std::fabs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      if (x < xm) e = b - x; else e = a - x;
      d = cc * e;
    } else {
      d = p / q; u = x + d;
      if (u - a < t2 || b - u < t2) { d = tol1; if (x >= xm) d = -d; }
    }
    if (std::fabs(d) >= tol1) u = x + d;
    else if (d > 0.0) u = x + tol1;
    else u = x - tol1;
    fu = f(u);
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  return x;
}

static const double OPT_TOL = 1.2207031250000000e-04; // .Machine$double.eps^0.25
static const double LO = -0.9999, HI = 0.9999;

// Bridge functions of Fan, Liu, Ning and Zou (2017, JRSS-B).
// Binary-binary bridge objective.
static inline double bb_obj(double x, double tau, double d1, double d2) {
  double v = 2.0 * Phi2(d1, d2, x) - 2.0 * Phi(d1) * Phi(d2) - tau;
  return v * v;
}
// Binary-continuous bridge objective.
static inline double bc_obj(double x, double tau, double d) {
  double v = 4.0 * Phi2(d, 0.0, x / std::sqrt(2.0)) - 2.0 * Phi(d) - tau;
  return v * v;
}

static inline double bb_solve(double tau, double d1, double d2) {
  return Brent_fmin(LO, HI, [&](double x){ return bb_obj(x, tau, d1, d2); }, OPT_TOL);
}
static inline double bc_solve(double tau, double d) {
  return Brent_fmin(LO, HI, [&](double x){ return bc_obj(x, tau, d); }, OPT_TOL);
}

// [[Rcpp::export]]
NumericMatrix estimate_latent_correlation_cpp(
    NumericMatrix kendall_tau, IntegerVector source_index,
    CharacterVector variable_types, List thresholds) {
  int p = variable_types.size();
  std::vector< std::vector<int> > idv(p);          // 0-based expanded cols per variable
  for (int c = 0; c < source_index.size(); ++c) idv[source_index[c] - 1].push_back(c);

  enum VarType { CONTINUOUS, BINARY, ORDINAL, UNKNOWN };
  std::vector<VarType> type(p, UNKNOWN);
  for (int i = 0; i < p; ++i) {
    std::string value = as<std::string>(variable_types[i]);
    if (value == "continuous") type[i] = CONTINUOUS;
    else if (value == "binary") type[i] = BINARY;
    else if (value == "ordinal") type[i] = ORDINAL;
  }

  NumericMatrix R(p, p);
  for (int i = 0; i < p; ++i) R(i, i) = 1.0;
  const double HALFPI = M_PI / 2.0;

  for (int i = 0; i < p - 1; ++i) {
    NumericVector di = thresholds[i];
    const std::vector<int>& ci = idv[i];
    for (int j = i + 1; j < p; ++j) {
      NumericVector dj = thresholds[j];
      const std::vector<int>& cj = idv[j];
      double val = 0.0;

      if (type[i] == CONTINUOUS && type[j] == CONTINUOUS) {
        val = std::sin(HALFPI * kendall_tau(ci[0], cj[0]));
      } else if (type[i] == BINARY && type[j] == CONTINUOUS) {
        val = bc_solve(kendall_tau(ci[0], cj[0]), di[0]);
      } else if (type[i] == CONTINUOUS && type[j] == BINARY) {
        val = bc_solve(kendall_tau(ci[0], cj[0]), dj[0]);
      } else if (type[i] == BINARY && type[j] == BINARY) {
        val = bb_solve(kendall_tau(ci[0], cj[0]), di[0], dj[0]);
      } else if (type[i] == CONTINUOUS && type[j] == ORDINAL) {
        int L = (int)cj.size(); double s = 0.0;
        for (int h = 0; h < L; ++h) s += bc_solve(kendall_tau(ci[0], cj[h]), dj[h]);
        val = s / L;
      } else if (type[i] == ORDINAL && type[j] == CONTINUOUS) {
        int L = (int)ci.size(); double s = 0.0;
        for (int h = 0; h < L; ++h) s += bc_solve(kendall_tau(ci[h], cj[0]), di[h]);
        val = s / L;
      } else if (type[i] == BINARY && type[j] == ORDINAL) {
        int L = (int)cj.size(); double s = 0.0;
        for (int h = 0; h < L; ++h) s += bb_solve(kendall_tau(ci[0], cj[h]), di[0], dj[h]);
        val = s / L;
      } else if (type[i] == ORDINAL && type[j] == BINARY) {
        int L = (int)ci.size(); double s = 0.0;
        for (int h = 0; h < L; ++h) s += bb_solve(kendall_tau(ci[h], cj[0]), di[h], dj[0]);
        val = s / L;
      } else if (type[i] == ORDINAL && type[j] == ORDINAL) {
        int L1 = (int)ci.size(), L2 = (int)cj.size(); double s = 0.0;
        for (int h = 0; h < L1; ++h)
          for (int k = 0; k < L2; ++k)
            s += bb_solve(kendall_tau(ci[h], cj[k]), di[h], dj[k]);
        val = s / (L1 * L2);
      }
      R(i, j) = val; R(j, i) = val;
    }
  }
  return R;
}
