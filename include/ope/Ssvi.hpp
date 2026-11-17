#pragma once

#include <vector>

#include "ope/Svi.hpp"

namespace ope::ssvi {

/// SSVI (Gatheral & Jacquier 2014): a whole surface from three global
/// parameters and one number per expiry, the ATM total variance theta(T):
///
///   w(k, theta) = theta/2 (1 + rho phi k + sqrt((phi k + rho)^2 + 1 - rho^2))
///   phi(theta)  = eta / (theta^gamma (1 + theta)^(1 - gamma))      (power law)
///
/// rho is the skew (shared by every expiry), eta the overall curvature, gamma
/// how fast the smile flattens as theta grows. The surface is free of static
/// arbitrage if
///   - theta(T) is non-decreasing                        (no calendar arbitrage)
///   - eta (1 + |rho|) <= 2 and 0 <= gamma <= 1/2        (no butterfly arbitrage)
/// (their Section 4). Under these bounds theta phi(theta) = eta (theta /
/// (1 + theta))^(1 - gamma) is increasing and below eta, and theta phi^2 <=
/// eta^2, so their sufficient butterfly conditions theta phi (1 + |rho|) < 4
/// and theta phi^2 (1 + |rho|) <= 4 hold, and their calendar condition reduces
/// to theta non-decreasing. Each slice is a raw SVI slice (Lemma 3.2): raw().
struct Params {
    double rho;
    double eta;
    double gamma;

    double phi(double theta) const;
    double w(double k, double theta) const;
    /// The same slice as raw SVI: a = theta (1 - rho^2)/2, b = theta phi/2,
    /// m = -rho/phi, sigma = sqrt(1 - rho^2)/phi.
    svi::Params raw(double theta) const;
};

/// One expiry's quotes: implied total variances w at log-moneyness k.
struct SliceData {
    double expiry;
    std::vector<double> k;
    std::vector<double> w;
};

struct Fit {
    Params params;
    std::vector<double> expiries;  ///< increasing
    std::vector<double> thetas;    ///< fitted theta at each expiry, non-decreasing
    double rmse_vol;               ///< RMS of (fitted vol - quoted vol) over all points
    double max_abs_vol;            ///< worst |fitted vol - quoted vol|
    int points;

    /// theta at any T > 0: linear in T between fitted expiries, theta_1 T / T_1
    /// before the first, theta_n T / T_n after the last (ATM vol held
    /// constant). Every T then lies on the same SSVI surface, so the whole
    /// surface keeps the no-arbitrage guarantees, not only the fitted slices.
    double theta(double T) const;
    double w(double k, double T) const { return params.w(k, theta(T)); }
};

/// Joint least-squares fit of SSVI to several expiries at once.
///
/// Minimizes sum over all points of (w_model - w_quote)^2 / (w_quote T):
/// approximately squared error in vol, since d vol = dw / (2 vol T). (svi::fit
/// weights by 1/w, which differs only by a constant within one expiry; across
/// expiries 1/w would weight each slice by T.) Subject to the no-arbitrage
/// constraints above, which every fit satisfies by construction.
///
/// Method: the 3 globals by Nelder-Mead in unconstrained coordinates
/// (rho = sin u, gamma = (1 + sin v)/4, eta = (1 + sin z)/(1 + |rho|), which map
/// R^3 onto exactly the arbitrage-free box), started from the best few of a
/// 7 x 3 x 3 grid. For fixed globals the thetas are profiled out: each slice's
/// theta by a 1-D scan plus golden-section search, then pool-adjacent-violators
/// to make them non-decreasing (merging neighbouring slices that want
/// decreasing theta and refitting one common theta for the pool). PAV is exact
/// when each slice's objective is convex in theta, which is not proven here.
///
/// Throws on no slices, a slice with fewer than 3 points, mismatched sizes,
/// non-finite input, non-positive w or expiry, or repeated expiries.
Fit fit(std::vector<SliceData> slices);

}  // namespace ope::ssvi
