#pragma once

#include <algorithm>

namespace ope {

enum class OptionType { Call, Put };

/// A European option contract: exercisable only at expiry.
struct EuropeanOption {
    OptionType type;
    double strike;  ///< K > 0
    double expiry;  ///< T >= 0, in years
};

/// Market state under Black-Scholes assumptions: constant rate, constant
/// volatility, continuous dividend yield (all annualized, continuously
/// compounded).
struct Market {
    double spot;              ///< S > 0
    double rate;              ///< r
    double vol;               ///< sigma >= 0
    double div_yield = 0.0;   ///< q
};

/// First- and second-order sensitivities of the option price V.
///
/// Units are the "raw derivative" convention — no per-1% or per-day scaling:
///   delta = dV/dS
///   gamma = d2V/dS2
///   vega  = dV/dsigma   (per 1.00 of vol; divide by 100 for per vol-point)
///   theta = dV/dt       (per year of calendar time passing, i.e. -dV/dT;
///                        divide by 365 for per-day)
///   rho   = dV/dr       (per 1.00 of rate; divide by 100 for per 1%)
struct Greeks {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
};

/// Payoff at expiry for terminal spot S_T.
inline double payoff(OptionType type, double strike, double spot_at_expiry) {
    return type == OptionType::Call ? std::max(spot_at_expiry - strike, 0.0)
                                    : std::max(strike - spot_at_expiry, 0.0);
}

/// Throws std::invalid_argument on S <= 0, K <= 0, T < 0, sigma < 0, or any
/// non-finite input. Shared by every pricer so they reject the same inputs.
void validate(const EuropeanOption& opt, const Market& mkt);

}  // namespace ope
