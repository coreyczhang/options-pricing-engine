#pragma once

#include "ope/Option.hpp"

namespace ope::iv {

/// Why a solve stopped. Anything other than Converged or AtLowerBound means no
/// volatility reproduces the quote, and Result::vol is NaN.
enum class Status {
    Converged,        ///< |vol step| or the bracket width fell below Settings::vol_tol
    AtLowerBound,     ///< price == no-arbitrage lower bound; the implied vol is 0
    BelowLowerBound,  ///< price < max(+/-(S e^{-qT} - K e^{-rT}), 0): arbitrage
    AboveUpperBound,  ///< price >= S e^{-qT} (call) or K e^{-rT} (put): arbitrage
    VolAboveMax,      ///< a solution, if any, needs vol > Settings::max_vol
    NotConverged,     ///< hit Settings::max_iter (not expected; see tests)
};

struct Settings {
    double vol_tol = 1e-12;  ///< absolute tolerance on sigma
    double max_vol = 10.0;   ///< 1000% vol; above this a quote is treated as junk
    int max_iter = 100;
};

struct Result {
    double vol;      ///< implied sigma, or NaN if status says there is none
    Status status;
    int iterations;  ///< Newton/bisection steps after bracketing
};

/// Black-Scholes implied volatility: the sigma at which bs::price(opt, mkt)
/// equals `price`.
///
/// The BS price is strictly increasing in sigma (vega > 0), running from the
/// no-arbitrage lower bound at sigma = 0 to the upper bound as sigma -> inf, so
/// a quote strictly inside the bounds has exactly one implied vol. Quotes on or
/// outside the bounds are reported through Status rather than thrown, because
/// they are routine in real option chains (stale or crossed quotes) and a
/// surface builder needs to skip them, not abort.
///
/// Method: Newton's method on sigma using the analytic vega, safeguarded by a
/// bracket [lo, hi] that always contains the root. A Newton step that would
/// leave the bracket, or a vanishing vega, falls back to bisection, so the
/// solve always converges; near the root Newton is quadratic.
///
/// `mkt.vol` is the starting guess (e.g. a neighbouring strike's vol when
/// building a surface). Pass 0 for an automatic start: the vol at which vega
/// peaks, sigma* = sqrt(2 |ln(S e^{-qT} / K e^{-rT})| / T) (Manaster-Koehler),
/// or for an at-the-forward option the Brenner-Subrahmanyam approximation
/// sigma ~ sqrt(2 pi / T) * price / S e^{-qT}. The guess changes only the
/// iteration count, never the answer.
///
/// Accuracy limit: the vol error is roughly (price rounding error) / vega, so
/// when vega is tiny — deep in or out of the money, or very short-dated — the
/// quote pins down sigma only loosely. The solver then still returns a sigma
/// that reproduces the price to rounding, but many nearby sigmas do too.
///
/// Throws std::invalid_argument on invalid contract/market inputs (see
/// ope::validate), a non-finite or negative price, or T == 0 (at expiry the
/// price no longer depends on vol, so implied vol is undefined).
Result solve(const EuropeanOption& opt, const Market& mkt, double price,
             const Settings& settings = {});

}  // namespace ope::iv
