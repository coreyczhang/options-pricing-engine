#pragma once

#include <vector>

namespace ope::svi {

/// Raw SVI (Gatheral 2004) total-variance smile for one expiry:
///
///   w(k) = a + b * (rho * (k - m) + sqrt((k - m)^2 + sigma^2)),   k = ln(K / F)
///
/// Five parameters: a sets the level, b the overall slope of the wings, rho
/// their asymmetry, m shifts the smile in k, sigma rounds off the vertex. The
/// wings are asymptotically linear in k, with slopes b (1 - rho) on the left
/// and b (1 + rho) on the right; Roger Lee's moment formula caps both at 2.
struct Params {
    double a;
    double b;
    double rho;
    double m;
    double sigma;

    double w(double k) const;    ///< total variance
    double dw(double k) const;   ///< dw/dk
    double d2w(double k) const;  ///< d2w/dk2
};

/// Durrleman's condition: a smile is free of butterfly arbitrage iff
///
///   g(k) = (1 - k w'/(2 w))^2 - (w'^2 / 4) (1/w + 1/4) + w''/2  >=  0
///
/// for all k (g is proportional to the risk-neutral density of ln S_T).
/// Returns the minimum of g over [k_lo, k_hi] on an n-point grid and where it
/// occurs; g = -inf if w is not positive somewhere on the grid.
struct DurrlemanMin {
    double g;
    double k;
};
DurrlemanMin min_durrleman_g(const Params& p, double k_lo, double k_hi, int n = 2001);

struct Fit {
    Params params;
    double expiry;
    double rmse_vol;       ///< RMS of (fitted vol - quoted vol) over the points
    double max_abs_vol;    ///< worst |fitted vol - quoted vol|
    int points;
};

/// Least-squares fit of raw SVI to one expiry's implied total variances w_i
/// at log-moneyness k_i.
///
/// Method (Zeliade / De Marco & Martini quasi-explicit): for fixed (m, sigma),
/// substituting y = (k - m) / sigma makes w = a + d y + c sqrt(y^2 + 1) linear
/// in (a, d, c) = (a, rho b sigma, b sigma), so the inner problem is a
/// 3-variable weighted least squares over a polytope:
///
///   0 <= a <= max w_i            (level non-negative, not above the data)
///   |d| <= c                     (|rho| <= 1, b >= 0)
///   c + |d| <= 2 sigma           (Lee: wing slopes b (1 +/- rho) <= 2)
///
/// solved exactly by enumerating active constraint sets. The outer problem in
/// (m, ln sigma) is 2-dimensional: a coarse grid, then Nelder-Mead from the
/// best grid point. a >= 0 is sufficient (not necessary) for w > 0; it rules
/// out some valid smiles in exchange for a convex inner problem.
///
/// Weights 1 / w_i make the objective approximately the squared error in
/// *vol*: d(vol) = dw / (2 vol T), so (d vol)^2 ~ dw^2 / (4 w T).
///
/// The constraints keep the wings sane but do NOT guarantee no butterfly
/// arbitrage; check min_durrleman_g on the result. Throws on fewer than 5
/// points, mismatched sizes, non-finite input, or non-positive w or expiry.
Fit fit(const std::vector<double>& k, const std::vector<double>& w, double expiry);

}  // namespace ope::svi
