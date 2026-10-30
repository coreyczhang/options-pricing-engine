#include "ope/MonteCarlo.hpp"

#include <cmath>
#include <stdexcept>

#include "ope/NormalRng.hpp"

namespace ope::mc {
namespace {

// Welford's online mean/variance: one pass, and no catastrophic
// cancellation from sum(x^2) - n*mean^2 when the payoff's spread is small
// relative to its level (deep ITM, low vol).
class RunningStats {
public:
    void add(double x) {
        ++n_;
        const double delta = x - mean_;
        mean_ += delta / static_cast<double>(n_);
        m2_ += delta * (x - mean_);
    }
    double mean() const { return mean_; }
    // Standard error of the mean. Callers guarantee n >= 2.
    double std_error() const {
        const double n = static_cast<double>(n_);
        return std::sqrt(m2_ / (n - 1.0) / n);
    }

private:
    std::size_t n_ = 0;
    double mean_ = 0.0;
    double m2_ = 0.0;
};

void validate_config(const Config& cfg) {
    if (cfg.paths < 2) throw std::invalid_argument("mc: need at least 2 paths");
    if (cfg.antithetic && (cfg.paths % 2 != 0 || cfg.paths < 4))
        throw std::invalid_argument("mc: antithetic sampling needs an even path count >= 4");
}

// One independent sampling unit per iteration: a single draw (plain) or a
// pair average over Z and -Z (antithetic). Statistics are taken over these
// units, which is what makes std_error correct in both modes.
std::size_t units(const Config& cfg) { return cfg.antithetic ? cfg.paths / 2 : cfg.paths; }

}  // namespace

Result price(const EuropeanOption& opt, const Market& mkt, const Config& cfg) {
    validate(opt, mkt);
    validate_config(cfg);

    const double T = opt.expiry;
    const double drift = (mkt.rate - mkt.div_yield - 0.5 * mkt.vol * mkt.vol) * T;
    const double vol_sqrt_t = mkt.vol * std::sqrt(T);
    const double discount = std::exp(-mkt.rate * T);

    auto payoff_at = [&](double z) {
        return payoff(opt.type, opt.strike, mkt.spot * std::exp(drift + vol_sqrt_t * z));
    };

    NormalRng rng(cfg.seed);
    RunningStats stats;
    for (std::size_t i = 0, n = units(cfg); i < n; ++i) {
        const double z = rng.next();
        stats.add(cfg.antithetic ? 0.5 * (payoff_at(z) + payoff_at(-z)) : payoff_at(z));
    }
    // Discount once at the end: it is a constant factor on every path.
    return Result{discount * stats.mean(), discount * stats.std_error(), cfg.paths};
}

GreeksResult greeks(const EuropeanOption& opt, const Market& mkt, const Config& cfg,
                    const Bumps& bumps) {
    validate(opt, mkt);
    validate_config(cfg);
    if (!(bumps.spot_rel > 0.0 && bumps.spot_rel < 1.0))
        throw std::invalid_argument("mc: spot_rel bump must be in (0, 1)");
    if (!(bumps.vol > 0.0 && bumps.time > 0.0 && bumps.rate > 0.0))
        throw std::invalid_argument("mc: bumps must be > 0");

    const double S = mkt.spot;
    const double K = opt.strike;
    const double T = opt.expiry;
    const double r = mkt.rate;
    const double q = mkt.div_yield;
    const double sigma = mkt.vol;

    // Bumped parameter values. Central difference unless the down-bump would
    // leave the domain, in which case the "down" point is the base point.
    const double h_s = bumps.spot_rel * S;
    const double vol_up = sigma + bumps.vol;
    const double vol_dn = sigma >= bumps.vol ? sigma - bumps.vol : sigma;
    const double t_up = T + bumps.time;
    const double t_dn = T >= bumps.time ? T - bumps.time : T;
    const double r_up = r + bumps.rate;
    const double r_dn = r - bumps.rate;  // rates may go negative

    // A bumped (vol, expiry, rate) scenario, with everything that does not
    // depend on the draw z hoisted out of the path loop. Spot needs no
    // scenario: S_T is linear in S, so spot bumps just rescale the base S_T.
    struct Scenario {
        double drift, vol_sqrt_t, discount;
        // Discounted payoff on draw z.
        double value(const EuropeanOption& o, double spot, double z) const {
            return discount * payoff(o.type, o.strike, spot * std::exp(drift + vol_sqrt_t * z));
        }
    };
    auto scenario = [&](double vol, double t, double rate) {
        return Scenario{(rate - q - 0.5 * vol * vol) * t, vol * std::sqrt(t),
                        std::exp(-rate * t)};
    };
    const Scenario base = scenario(sigma, T, r);
    const Scenario s_vol_up = scenario(vol_up, T, r), s_vol_dn = scenario(vol_dn, T, r);
    const Scenario s_t_up = scenario(sigma, t_up, r), s_t_dn = scenario(sigma, t_dn, r);
    const Scenario s_r_up = scenario(sigma, T, r_up), s_r_dn = scenario(sigma, T, r_dn);

    // Per-draw finite differences, all on the same z (common random numbers).
    struct Contrib {
        double delta, gamma, vega, theta, rho;
    };
    auto contrib = [&](double z) {
        const double g = std::exp(base.drift + base.vol_sqrt_t * z);  // S_T / S
        const double mid = base.discount * payoff(opt.type, K, S * g);
        const double up = base.discount * payoff(opt.type, K, (S + h_s) * g);
        const double dn = base.discount * payoff(opt.type, K, (S - h_s) * g);

        Contrib c{};
        c.delta = (up - dn) / (2.0 * h_s);
        c.gamma = (up - 2.0 * mid + dn) / (h_s * h_s);
        c.vega = (s_vol_up.value(opt, S, z) - s_vol_dn.value(opt, S, z)) / (vol_up - vol_dn);
        // theta = dV/dt = -dV/dT
        c.theta = -(s_t_up.value(opt, S, z) - s_t_dn.value(opt, S, z)) / (t_up - t_dn);
        c.rho = (s_r_up.value(opt, S, z) - s_r_dn.value(opt, S, z)) / (r_up - r_dn);
        return c;
    };

    NormalRng rng(cfg.seed);
    RunningStats s_delta, s_gamma, s_vega, s_theta, s_rho;
    for (std::size_t i = 0, n = units(cfg); i < n; ++i) {
        const double z = rng.next();
        Contrib c = contrib(z);
        if (cfg.antithetic) {
            const Contrib m = contrib(-z);
            c = {0.5 * (c.delta + m.delta), 0.5 * (c.gamma + m.gamma),
                 0.5 * (c.vega + m.vega), 0.5 * (c.theta + m.theta), 0.5 * (c.rho + m.rho)};
        }
        s_delta.add(c.delta);
        s_gamma.add(c.gamma);
        s_vega.add(c.vega);
        s_theta.add(c.theta);
        s_rho.add(c.rho);
    }

    GreeksResult out{};
    out.value = {s_delta.mean(), s_gamma.mean(), s_vega.mean(), s_theta.mean(), s_rho.mean()};
    out.std_error = {s_delta.std_error(), s_gamma.std_error(), s_vega.std_error(),
                     s_theta.std_error(), s_rho.std_error()};
    out.paths = cfg.paths;
    return out;
}

}  // namespace ope::mc
