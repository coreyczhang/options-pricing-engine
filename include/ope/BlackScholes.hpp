#pragma once

#include "ope/Option.hpp"

namespace ope::bs {

/// Closed-form Black-Scholes(-Merton) price of a European option, with a
/// continuous dividend yield q:
///
///   d1 = [ln(S/K) + (r - q + sigma^2/2) T] / (sigma sqrt(T)),  d2 = d1 - sigma sqrt(T)
///   C  = S e^{-qT} N(d1) - K e^{-rT} N(d2)
///   P  = K e^{-rT} N(-d2) - S e^{-qT} N(-d1)
///
/// Degenerate inputs (T == 0 or sigma == 0) have no diffusion, so the formula's
/// d1/d2 are undefined; we return the exact deterministic limit instead:
/// the discounted intrinsic value of the forward, max(+/-(S e^{-qT} - K e^{-rT}), 0).
///
/// Throws std::invalid_argument on invalid inputs (see ope::validate).
double price(const EuropeanOption& opt, const Market& mkt);

/// Analytic Greeks (see Greeks for units). In the degenerate case the Greeks
/// are the derivatives of the deterministic limit price: gamma and vega are 0,
/// and at the exact kink (S e^{-qT} == K e^{-rT}) delta/theta/rho take the
/// midpoint of the one-sided values — matching N(d1) -> 1/2 as sigma -> 0.
Greeks greeks(const EuropeanOption& opt, const Market& mkt);

}  // namespace ope::bs
