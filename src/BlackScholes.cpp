#include "ope/BlackScholes.hpp"

#include <algorithm>
#include <cmath>

namespace ope::bs {
namespace {

constexpr double kInvSqrt2Pi = 0.39894228040143267794;  // 1 / sqrt(2 pi)
constexpr double kInvSqrt2 = 0.70710678118654752440;    // 1 / sqrt(2)

// N(x) via erfc rather than 0.5 * (1 + erf(x)): for large negative x the erf
// form cancels catastrophically (1 + (-1 + tiny)), erfc keeps full relative
// precision deep in the tails, which matters for far-OTM prices.
double norm_cdf(double x) { return 0.5 * std::erfc(-x * kInvSqrt2); }

double norm_pdf(double x) { return kInvSqrt2Pi * std::exp(-0.5 * x * x); }

// Quantities shared by price() and greeks().
struct Terms {
    double fwd_spot;    // S e^{-qT}: PV of receiving the stock at T
    double pv_strike;   // K e^{-rT}: PV of paying the strike at T
    double sqrt_t;
    bool degenerate;    // T == 0 or sigma == 0: no diffusion
    double d1;          // only meaningful when !degenerate
    double d2;
};

Terms terms(const EuropeanOption& opt, const Market& mkt) {
    validate(opt, mkt);
    Terms t{};
    const double T = opt.expiry;
    t.fwd_spot = mkt.spot * std::exp(-mkt.div_yield * T);
    t.pv_strike = opt.strike * std::exp(-mkt.rate * T);
    t.sqrt_t = std::sqrt(T);
    const double sig_sqrt_t = mkt.vol * t.sqrt_t;
    t.degenerate = !(sig_sqrt_t > 0.0);
    if (!t.degenerate) {
        // ln(S e^{-qT} / K e^{-rT}) == ln(S/K) + (r - q) T.
        t.d1 = std::log(t.fwd_spot / t.pv_strike) / sig_sqrt_t + 0.5 * sig_sqrt_t;
        t.d2 = t.d1 - sig_sqrt_t;
    }
    return t;
}

// Weight on the "call finishes in the money" branch of the deterministic
// limit: 1 if the forward is above the strike, 0 below, 1/2 exactly at it.
double call_itm_weight(const Terms& t) {
    if (t.fwd_spot > t.pv_strike) return 1.0;
    if (t.fwd_spot < t.pv_strike) return 0.0;
    return 0.5;
}

}  // namespace

double price(const EuropeanOption& opt, const Market& mkt) {
    const Terms t = terms(opt, mkt);
    const bool call = opt.type == OptionType::Call;

    if (t.degenerate) {
        const double fwd_intrinsic = t.fwd_spot - t.pv_strike;
        return std::max(call ? fwd_intrinsic : -fwd_intrinsic, 0.0);
    }
    if (call) return t.fwd_spot * norm_cdf(t.d1) - t.pv_strike * norm_cdf(t.d2);
    return t.pv_strike * norm_cdf(-t.d2) - t.fwd_spot * norm_cdf(-t.d1);
}

Greeks greeks(const EuropeanOption& opt, const Market& mkt) {
    const Terms t = terms(opt, mkt);
    const bool call = opt.type == OptionType::Call;
    const double T = opt.expiry;
    const double r = mkt.rate;
    const double q = mkt.div_yield;
    const double disc_q = std::exp(-q * T);
    Greeks g{};

    if (t.degenerate) {
        // Derivatives of max(+/-(S e^{-qT} - K e^{-rT}), 0) on its active branch.
        const double w = call ? call_itm_weight(t) : 1.0 - call_itm_weight(t);
        const double sign = call ? 1.0 : -1.0;
        g.delta = sign * w * disc_q;
        g.gamma = 0.0;
        g.vega = 0.0;
        g.theta = sign * w * (q * t.fwd_spot - r * t.pv_strike);
        g.rho = sign * w * T * t.pv_strike;
        return g;
    }

    const double pdf_d1 = norm_pdf(t.d1);
    // Shared by call and put: the time-decay of the diffusion (optionality).
    const double theta_vol = -t.fwd_spot * pdf_d1 * mkt.vol / (2.0 * t.sqrt_t);

    g.gamma = disc_q * pdf_d1 / (mkt.spot * mkt.vol * t.sqrt_t);
    g.vega = t.fwd_spot * pdf_d1 * t.sqrt_t;

    if (call) {
        const double n_d1 = norm_cdf(t.d1);
        const double n_d2 = norm_cdf(t.d2);
        g.delta = disc_q * n_d1;
        g.theta = theta_vol - r * t.pv_strike * n_d2 + q * t.fwd_spot * n_d1;
        g.rho = T * t.pv_strike * n_d2;
    } else {
        const double n_md1 = norm_cdf(-t.d1);
        const double n_md2 = norm_cdf(-t.d2);
        g.delta = -disc_q * n_md1;
        g.theta = theta_vol + r * t.pv_strike * n_md2 - q * t.fwd_spot * n_md1;
        g.rho = -T * t.pv_strike * n_md2;
    }
    return g;
}

}  // namespace ope::bs
