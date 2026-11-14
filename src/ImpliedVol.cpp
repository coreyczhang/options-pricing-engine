#include "ope/ImpliedVol.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "ope/BlackScholes.hpp"

namespace ope::iv {
namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kTwoPi = 6.28318530717958647693;
constexpr double kEps = std::numeric_limits<double>::epsilon();

Result fail(Status s) { return {kNaN, s, 0}; }

// Starting point when the caller gives none. sigma* is where vega peaks (the
// price's inflection point in sigma); at the forward that is sigma -> 0, so use
// the ATM approximation C ~ S e^{-qT} sigma sqrt(T / 2 pi) instead.
double automatic_guess(double fwd_spot, double pv_strike, double T, double price) {
    const double log_moneyness = std::log(fwd_spot / pv_strike);
    if (log_moneyness != 0.0) return std::sqrt(2.0 * std::abs(log_moneyness) / T);
    return std::sqrt(kTwoPi / T) * price / fwd_spot;
}

}  // namespace

Result solve(const EuropeanOption& opt, const Market& mkt, double price,
             const Settings& settings) {
    validate(opt, mkt);
    if (!std::isfinite(price) || price < 0.0)
        throw std::invalid_argument("price must be finite and >= 0");
    if (opt.expiry == 0.0)
        throw std::invalid_argument("implied vol is undefined at expiry (T == 0)");
    if (!(settings.vol_tol > 0.0) || !(settings.max_vol > 0.0) || settings.max_iter < 1)
        throw std::invalid_argument("invalid implied-vol settings");

    const double T = opt.expiry;
    const bool call = opt.type == OptionType::Call;
    const double fwd_spot = mkt.spot * std::exp(-mkt.div_yield * T);
    const double pv_strike = opt.strike * std::exp(-mkt.rate * T);
    const double lower = std::max(call ? fwd_spot - pv_strike : pv_strike - fwd_spot, 0.0);
    const double upper = call ? fwd_spot : pv_strike;

    // In the money, lower is a difference of two numbers of size
    // ~max(S e^{-qT}, K e^{-rT}), so it (and any BS price near it) carries
    // rounding error of a few ulps of that size. Quotes within that band are
    // treated as on the bound: otherwise a deep-ITM price computed by bs::price
    // itself could round to just under intrinsic and be rejected as an
    // arbitrage. Out of the money the bound is exactly 0 and erfc keeps tiny
    // prices accurate, so no band: a 1e-13 OTM quote still has an implied vol.
    const double bound_tol = lower > 0.0 ? 4.0 * kEps * std::max(fwd_spot, pv_strike) : 0.0;
    if (price < lower - bound_tol) return fail(Status::BelowLowerBound);
    if (price >= upper) return fail(Status::AboveUpperBound);
    if (price <= lower + bound_tol) return {0.0, Status::AtLowerBound, 0};

    Market m = mkt;
    auto excess = [&](double sigma) {  // BS price minus quote: increasing in sigma
        m.vol = sigma;
        return bs::price(opt, m) - price;
    };

    // Bracket the root: excess(0) < 0 is guaranteed by price > lower; grow hi
    // until excess(hi) >= 0. Three doublings from 1.0 already reach 800% vol.
    double lo = 0.0;
    double hi = std::min(1.0, settings.max_vol);
    while (excess(hi) < 0.0) {
        if (hi >= settings.max_vol) return fail(Status::VolAboveMax);
        lo = hi;
        hi = std::min(2.0 * hi, settings.max_vol);
    }

    double x = mkt.vol > 0.0 ? mkt.vol : automatic_guess(fwd_spot, pv_strike, T, price);
    if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);

    // Stop once the residual is indistinguishable from rounding noise in
    // bs::price: further steps would only chase that noise. The price is a
    // difference of two terms, S e^{-qT} N(+/-d1) and K e^{-rT} N(+/-d2), so its
    // absolute rounding error is ~eps * (sum of their magnitudes) — large next
    // to the price itself whenever the terms nearly cancel (deep ITM, or
    // near the forward at low vol). That sum is |V - 2 S delta| for calls and
    // puts alike. The denorm_min term covers subnormal prices (< 2.2e-308),
    // whose precision is absolute, not relative.
    //
    // (Two earlier versions scaled the noise by the price, then by price +
    // 2 * intrinsic. Both underestimated it in some regime, so those solves
    // bisected on noise down to vol_tol, taking a libm-dependent number of
    // iterations: 44 on Linux/glibc vs 30 on macOS.)
    double prev_step = hi - lo;
    for (int iter = 1; iter <= settings.max_iter; ++iter) {
        m.vol = x;
        const double model = bs::price(opt, m);
        const Greeks g = bs::greeks(opt, m);
        const double f = model - price;
        const double noise = 8.0 * kEps * std::abs(model - 2.0 * mkt.spot * g.delta) +
                             8.0 * std::numeric_limits<double>::denorm_min();
        if (std::abs(f) <= noise) return {x, Status::Converged, iter};
        (f < 0.0 ? lo : hi) = x;

        // Newton unless it would leave the bracket or isn't at least halving
        // the step each iteration (|f / vega| > |prev_step| / 2); then bisect.
        // The second test is what guarantees bisection-speed progress when
        // Newton crawls, e.g. on the flat tails of the price curve.
        const double vega = g.vega;
        double next = x - f / vega;
        if (!(vega > 0.0) || !(next > lo && next < hi) ||
            std::abs(2.0 * f) > std::abs(prev_step * vega))
            next = 0.5 * (lo + hi);

        prev_step = next - x;
        if (std::abs(prev_step) <= settings.vol_tol || hi - lo <= settings.vol_tol)
            return {next, Status::Converged, iter};
        x = next;
    }
    return {x, Status::NotConverged, settings.max_iter};
}

}  // namespace ope::iv
