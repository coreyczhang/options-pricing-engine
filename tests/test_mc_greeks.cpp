#include "ope/MonteCarlo.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "ope/BlackScholes.hpp"

namespace ope {
namespace {

// Bump-and-reprice MC Greeks carry two separate errors:
//   1. bump bias: a finite difference is not the derivative (O(h^2) central);
//   2. sampling noise: the MC estimate of that finite difference.
// Testing MC directly against the analytic Greek with one loose tolerance
// would blur the two. Instead:
//   * MC vs the *same finite difference applied to the BS price* — any gap is
//     pure sampling noise, so a tight 4-sigma test applies;
//   * that BS finite difference vs the analytic Greek — deterministic bump
//     bias, checked against a small tolerance.
// Thresholds are ~4 standard errors, as in test_monte_carlo.cpp.
constexpr double kSigmas = 4.0;

struct Case {
    OptionType type;
    double S, K, T, r, vol, q;
};

std::vector<Case> cases() {
    return {
        {OptionType::Call, 100.0, 100.0, 1.0, 0.05, 0.20, 0.00},   // ATM
        {OptionType::Put, 100.0, 100.0, 1.0, 0.05, 0.20, 0.00},
        {OptionType::Call, 100.0, 120.0, 0.5, 0.03, 0.25, 0.02},   // OTM call
        {OptionType::Put, 100.0, 120.0, 0.5, 0.03, 0.25, 0.02},    // ITM put
        {OptionType::Call, 100.0, 80.0, 2.0, 0.01, 0.40, 0.03},    // ITM, long, high vol
        {OptionType::Call, 49.0, 50.0, 20.0 / 52.0, 0.05, 0.20, 0.0},  // Hull Greeks example
    };
}

EuropeanOption option(const Case& c) { return {c.type, c.K, c.T}; }
Market market(const Case& c) { return {c.S, c.r, c.vol, c.q}; }

// The finite differences mc::greeks computes, applied to the exact BS price.
Greeks bs_fd(const Case& c, const mc::Bumps& b) {
    auto px = [&](double S, double T, double r, double vol) {
        return bs::price(EuropeanOption{c.type, c.K, T}, Market{S, r, vol, c.q});
    };
    const double h = b.spot_rel * c.S;
    const double vol_dn = c.vol >= b.vol ? c.vol - b.vol : c.vol;
    const double t_dn = c.T >= b.time ? c.T - b.time : c.T;
    const double up = px(c.S + h, c.T, c.r, c.vol);
    const double mid = px(c.S, c.T, c.r, c.vol);
    const double dn = px(c.S - h, c.T, c.r, c.vol);
    Greeks g{};
    g.delta = (up - dn) / (2 * h);
    g.gamma = (up - 2 * mid + dn) / (h * h);
    g.vega = (px(c.S, c.T, c.r, c.vol + b.vol) - px(c.S, c.T, c.r, vol_dn)) /
             (c.vol + b.vol - vol_dn);
    g.theta = -(px(c.S, c.T + b.time, c.r, c.vol) - px(c.S, t_dn, c.r, c.vol)) /
              (c.T + b.time - t_dn);
    g.rho = (px(c.S, c.T, c.r + b.rate, c.vol) - px(c.S, c.T, c.r - b.rate, c.vol)) /
            (2 * b.rate);
    return g;
}

std::string label(const Case& c) {
    return std::string(c.type == OptionType::Call ? "call" : "put") +
           " S=" + std::to_string(c.S) + " K=" + std::to_string(c.K) +
           " T=" + std::to_string(c.T);
}

// ---- 1. Sampling noise: MC vs the same finite difference of BS -----------

class McGreeks : public ::testing::TestWithParam<bool> {};

TEST_P(McGreeks, MatchBlackScholesFiniteDifferenceWithinStdErrors) {
    const bool antithetic = GetParam();
    const mc::Bumps bumps{};
    for (const Case& c : cases()) {
        SCOPED_TRACE(label(c));
        const mc::GreeksResult mc =
            mc::greeks(option(c), market(c), mc::Config{400'000, 17, antithetic}, bumps);
        const Greeks target = bs_fd(c, bumps);
        EXPECT_NEAR(mc.value.delta, target.delta, kSigmas * mc.std_error.delta);
        EXPECT_NEAR(mc.value.gamma, target.gamma, kSigmas * mc.std_error.gamma);
        EXPECT_NEAR(mc.value.vega, target.vega, kSigmas * mc.std_error.vega);
        EXPECT_NEAR(mc.value.theta, target.theta, kSigmas * mc.std_error.theta);
        EXPECT_NEAR(mc.value.rho, target.rho, kSigmas * mc.std_error.rho);
    }
}

INSTANTIATE_TEST_SUITE_P(PlainAndAntithetic, McGreeks, ::testing::Bool(),
                         [](const auto& info) { return info.param ? "Antithetic" : "Plain"; });

// ---- 2. Bump bias: BS finite difference vs analytic Greek ------------------

TEST(McGreeksBumpBias, DefaultBumpsAreCloseToAnalyticGreeks) {
    // With the default bumps the deterministic bias is well under 0.5% of
    // each Greek on these contracts — far below the MC noise at practical
    // path counts, so MC Greeks are effectively unbiased estimates of the
    // analytic ones.
    for (const Case& c : cases()) {
        SCOPED_TRACE(label(c));
        const Greeks fd = bs_fd(c, mc::Bumps{});
        const Greeks exact = bs::greeks(option(c), market(c));
        auto tol = [](double x) { return 1e-6 + 5e-3 * std::abs(x); };
        EXPECT_NEAR(fd.delta, exact.delta, tol(exact.delta));
        EXPECT_NEAR(fd.gamma, exact.gamma, tol(exact.gamma));
        EXPECT_NEAR(fd.vega, exact.vega, tol(exact.vega));
        EXPECT_NEAR(fd.theta, exact.theta, tol(exact.theta));
        EXPECT_NEAR(fd.rho, exact.rho, tol(exact.rho));
    }
}

TEST(McGreeksBumpBias, CombinedErrorAgainstAnalyticGreeks) {
    // End-to-end: MC Greek vs the analytic Greek, allowing 4 se of noise plus
    // the (known, deterministic) bump bias.
    const mc::Bumps bumps{};
    for (const Case& c : cases()) {
        SCOPED_TRACE(label(c));
        const mc::GreeksResult mc =
            mc::greeks(option(c), market(c), mc::Config{400'000, 23}, bumps);
        const Greeks exact = bs::greeks(option(c), market(c));
        const Greeks fd = bs_fd(c, bumps);
        auto band = [](double se, double fd_v, double exact_v) {
            return kSigmas * se + std::abs(fd_v - exact_v);
        };
        EXPECT_NEAR(mc.value.delta, exact.delta, band(mc.std_error.delta, fd.delta, exact.delta));
        EXPECT_NEAR(mc.value.gamma, exact.gamma, band(mc.std_error.gamma, fd.gamma, exact.gamma));
        EXPECT_NEAR(mc.value.vega, exact.vega, band(mc.std_error.vega, fd.vega, exact.vega));
        EXPECT_NEAR(mc.value.theta, exact.theta, band(mc.std_error.theta, fd.theta, exact.theta));
        EXPECT_NEAR(mc.value.rho, exact.rho, band(mc.std_error.rho, fd.rho, exact.rho));
    }
}

// ---- 3. Why common random numbers -----------------------------------------

TEST(McGreeksCrn, DeltaNoiseStaysBoundedAsBumpShrinksWithCrn) {
    // With CRN, delta's per-path difference is bounded (the payoff is
    // Lipschitz in S), so its std error barely moves as h shrinks 100x.
    // Pricing up/down with independent seeds instead leaves both prices' full
    // noise in (V_up - V_dn) / 2h: std error ~ 1/h, i.e. 100x worse here.
    const Case c = cases()[0];
    constexpr std::size_t kPaths = 100'000;
    auto crn_se = [&](double rel) {
        return mc::greeks(option(c), market(c), mc::Config{kPaths, 5}, mc::Bumps{rel})
            .std_error.delta;
    };
    auto independent_se = [&](double rel) {
        const double h = rel * c.S;
        Market up = market(c), dn = market(c);
        up.spot += h;
        dn.spot -= h;
        const double se_up = mc::price(option(c), up, mc::Config{kPaths, 1001}).std_error;
        const double se_dn = mc::price(option(c), dn, mc::Config{kPaths, 2002}).std_error;
        return std::sqrt(se_up * se_up + se_dn * se_dn) / (2 * h);
    };

    const double crn_big = crn_se(1e-2), crn_small = crn_se(1e-4);
    const double ind_big = independent_se(1e-2), ind_small = independent_se(1e-4);
    EXPECT_NEAR(crn_small / crn_big, 1.0, 0.1);
    EXPECT_NEAR(ind_small / ind_big, 100.0, 10.0);
    // Even at the larger 1% bump CRN is far ahead. Expected ratio:
    //   independent: sqrt(2) sd(payoff) / (2h)  = sqrt(2) * 14.5 / 2  ~ 10.3
    //   CRN:         sd(per-path delta)         ~ sd(e^{-rT} S_T/S 1{S_T>K}) ~ 0.55
    // so ~18x — and it grows 10x for every 10x smaller bump.
    EXPECT_NEAR(ind_big / crn_big, 18.5, 2.0);
}

TEST(McGreeksCrn, GammaNoiseGrowsLikeOneOverSqrtBump) {
    // Gamma's second difference is ~1/h on the ~h-wide band of paths that
    // end near the strike, and 0 elsewhere: variance ~ 1/h, se ~ 1/sqrt(h).
    // Shrinking the bump 10x should cost ~sqrt(10) in std error even with CRN.
    const Case c = cases()[0];
    auto gamma_se = [&](double rel) {
        return mc::greeks(option(c), market(c), mc::Config{400'000, 9}, mc::Bumps{rel})
            .std_error.gamma;
    };
    EXPECT_NEAR(gamma_se(1e-3) / gamma_se(1e-2), std::sqrt(10.0), 0.15 * std::sqrt(10.0));
}

// ---- Degenerate inputs & validation ---------------------------------------

TEST(McGreeksDegenerate, ZeroVolUsesForwardVolBumpAndIsDeterministic) {
    // sigma = 0: every path is the forward, so delta is exact with zero error;
    // the vol down-bump would go negative, so vega falls back to a forward
    // difference instead of throwing.
    const EuropeanOption opt{OptionType::Call, 90.0, 1.0};
    const Market mkt{100.0, 0.05, 0.0, 0.01};
    const mc::GreeksResult g = mc::greeks(opt, mkt, mc::Config{1'000});
    EXPECT_NEAR(g.value.delta, std::exp(-0.01), 1e-12);
    EXPECT_EQ(g.std_error.delta, 0.0);
    EXPECT_NEAR(g.value.gamma, 0.0, 1e-10);  // (up - 2 mid + dn) / h^2: roundoff only
    EXPECT_TRUE(std::isfinite(g.value.vega));
}

TEST(McGreeksDegenerate, ExpiryShorterThanBumpUsesForwardDifference) {
    const EuropeanOption opt{OptionType::Call, 100.0, 0.5 / 365.0};
    const Market mkt{100.0, 0.05, 0.2};
    const mc::GreeksResult g = mc::greeks(opt, mkt, mc::Config{10'000});
    EXPECT_TRUE(std::isfinite(g.value.theta));
    EXPECT_LT(g.value.theta, 0.0);  // a long ATM call loses value as time passes
}

TEST(McGreeksValidation, RejectsBadBumps) {
    const EuropeanOption opt{OptionType::Call, 100.0, 1.0};
    const Market mkt{100.0, 0.05, 0.2};
    const mc::Config cfg{1'000};
    EXPECT_THROW(mc::greeks(opt, mkt, cfg, mc::Bumps{0.0}), std::invalid_argument);
    EXPECT_THROW(mc::greeks(opt, mkt, cfg, mc::Bumps{1.0}), std::invalid_argument);
    EXPECT_THROW(mc::greeks(opt, mkt, cfg, mc::Bumps{0.01, -1e-3}), std::invalid_argument);
    EXPECT_THROW(mc::greeks(opt, mkt, mc::Config{1'001, 42, true}), std::invalid_argument);
}

}  // namespace
}  // namespace ope
