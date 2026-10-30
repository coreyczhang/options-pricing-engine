#pragma once

#include <cstddef>
#include <cstdint>

#include "ope/Option.hpp"

namespace ope::mc {

struct Config {
    /// Number of simulated terminal prices (payoff evaluations), >= 2. With
    /// antithetic sampling this must be even and >= 4: paths/2 draws Z, each
    /// used as both Z and -Z (and >= 2 pairs, so a sample variance exists).
    /// Holding `paths` fixed keeps plain-vs-antithetic comparisons at equal
    /// payoff work.
    std::size_t paths;
    std::uint64_t seed = 42;
    /// Antithetic variates: pair every draw Z with -Z.
    bool antithetic = false;
};

/// A Monte Carlo estimate with its own error bar.
struct Result {
    double price;
    /// Standard error of `price`, computed over the *independent* sampling
    /// units: single payoffs for plain MC, pair averages (f(Z) + f(-Z)) / 2
    /// for antithetic — the two halves of a pair are correlated, so treating
    /// them as 2 independent samples would misstate the error.
    /// price +/- 1.96 * std_error is an approximate 95% confidence interval
    /// for the true (model) price.
    double std_error;
    std::size_t paths;
};

/// Monte Carlo price of a European option under the same GBM model as
/// ope::bs, so the two must agree up to sampling error.
///
/// A European payoff depends only on S_T, and GBM's S_T is exactly lognormal:
///   S_T = S exp((r - q - sigma^2/2) T + sigma sqrt(T) Z),  Z ~ N(0, 1)
/// so each path is one normal draw — no time-stepping and no discretization
/// bias. The only error is sampling error, which shrinks as 1/sqrt(paths).
///
/// Antithetic variates: since Z and -Z are equally likely, the pair average
/// is still unbiased, and its variance is Var(f(Z)) (1 + rho) / 2 where
/// rho = corr(f(Z), f(-Z)). At equal payoff work the variance reduction factor
/// is therefore 1 / (1 + rho): large when the payoff is close to linear in Z
/// (rho -> -1), ~1 when it is not. For a monotone payoff such as a vanilla
/// call or put, rho <= 0, so antithetic never hurts.
///
/// Deterministic for a given seed. Throws std::invalid_argument on invalid
/// market/option inputs (see ope::validate), paths < 2, or (antithetic) odd
/// paths or paths < 4.
Result price(const EuropeanOption& opt, const Market& mkt, const Config& cfg);

/// Bump sizes for finite-difference Greeks.
///
/// Bias vs noise: a central difference is biased by O(h^2), so smaller is
/// more accurate — but noise depends on the Greek. With common random numbers
/// the per-path difference for delta/vega/rho/theta stays bounded as h -> 0,
/// so small bumps are cheap. Gamma's second difference is nonzero only for
/// paths ending within ~h of the strike, where it is ~1/h, so its variance
/// grows like 1/h: gamma needs the larger spot bump. Defaults are chosen so
/// the O(h^2) bias is well under the sampling error at ~1e6 paths.
struct Bumps {
    double spot_rel = 0.01;       ///< spot bump as a fraction of S, in (0, 1)
    double vol = 1e-3;            ///< absolute vol bump
    double time = 1.0 / 365.0;    ///< absolute expiry bump (one calendar day)
    double rate = 1e-4;           ///< absolute rate bump (1 bp)
};

/// Finite-difference Greeks with a standard error on each.
struct GreeksResult {
    Greeks value;
    Greeks std_error;   ///< over independent units, as for Result::std_error
    std::size_t paths;
};

/// Bump-and-reprice Greeks on the Monte Carlo engine, in the same units as
/// ope::bs::greeks.
///
/// Common random numbers: every bumped payoff is evaluated on the *same*
/// draw Z as the base payoff, and each Greek is the sample mean of the
/// per-path finite differences. Pricing the bumps with independent draws
/// instead would leave both prices' full noise in a difference divided by h:
/// variance ~ 1/h^2, which explodes as the bump shrinks.
///
/// Central differences throughout, except a forward difference where the
/// down-bump would leave the valid domain (vol < bump, expiry < bump).
///
/// What this estimates is the finite difference of the model price, which
/// differs from the true derivative by the O(h^2) bump bias; see the tests
/// for how the two error sources are checked separately.
///
/// Throws std::invalid_argument on the same inputs as price(), or on
/// non-positive bumps / spot_rel >= 1.
GreeksResult greeks(const EuropeanOption& opt, const Market& mkt, const Config& cfg,
                    const Bumps& bumps = {});

}  // namespace ope::mc
