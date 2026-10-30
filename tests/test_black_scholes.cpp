#include "ope/BlackScholes.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ope {
namespace {

EuropeanOption call(double K, double T) { return {OptionType::Call, K, T}; }
EuropeanOption put(double K, double T) { return {OptionType::Put, K, T}; }

// A grid spanning ITM/ATM/OTM, short/long dated, low/high vol, with and
// without rates and dividends. Used by the identity-style tests below.
struct Case {
    double S, K, T, r, vol, q;
};

std::vector<Case> grid() {
    std::vector<Case> cases;
    for (double S : {70.0, 95.0, 100.0, 105.0, 140.0})
        for (double T : {0.05, 0.5, 2.0})
            for (double vol : {0.05, 0.25, 0.8})
                for (double r : {0.0, 0.05})
                    for (double q : {0.0, 0.03}) cases.push_back({S, 100.0, T, r, vol, q});
    return cases;
}

// ---- Reference values ------------------------------------------------------
// Published textbook numbers. Tolerances match the precision each source
// prints (Hull: 2-3 significant decimals; Haug: 4 decimals).

TEST(BlackScholesReference, HullStockOptionExample) {
    // Hull, Options, Futures, and Other Derivatives — Black-Scholes-Merton
    // chapter worked example: S=42, K=40, r=10%, sigma=20%, T=0.5 -> c=4.76, p=0.81.
    const Market mkt{42.0, 0.10, 0.20};
    EXPECT_NEAR(bs::price(call(40.0, 0.5), mkt), 4.76, 0.005);
    EXPECT_NEAR(bs::price(put(40.0, 0.5), mkt), 0.81, 0.005);
}

TEST(BlackScholesReference, HullGreekLettersExample) {
    // Hull — "Greek letters" chapter running example: a 20-week call,
    // S=49, K=50, r=5%, sigma=20%, T=20/52.
    // Hull reports: price ~2.40, delta 0.522, gamma 0.066, vega 12.1,
    // theta -4.31 (per year), rho 8.91.
    const Market mkt{49.0, 0.05, 0.20};
    const EuropeanOption opt = call(50.0, 20.0 / 52.0);
    const Greeks g = bs::greeks(opt, mkt);
    EXPECT_NEAR(bs::price(opt, mkt), 2.40, 0.005);
    EXPECT_NEAR(g.delta, 0.522, 0.0005);
    EXPECT_NEAR(g.gamma, 0.066, 0.0005);
    EXPECT_NEAR(g.vega, 12.1, 0.05);
    EXPECT_NEAR(g.theta, -4.31, 0.005);
    EXPECT_NEAR(g.rho, 8.91, 0.005);
}

TEST(BlackScholesReference, HaugGeneralizedBlackScholesCall) {
    // Haug, The Complete Guide to Option Pricing Formulas — generalized BSM
    // example: S=60, K=65, T=0.25, r=8%, b=r (no dividend), sigma=30% -> c=2.1334.
    EXPECT_NEAR(bs::price(call(65.0, 0.25), Market{60.0, 0.08, 0.30}), 2.1334, 5e-5);
}

TEST(BlackScholesReference, HaugMertonDividendPut) {
    // Haug — Merton (1973) continuous-dividend example:
    // S=100, K=95, T=0.5, r=10%, q=5%, sigma=20% -> p=2.4648.
    EXPECT_NEAR(bs::price(put(95.0, 0.5), Market{100.0, 0.10, 0.20, 0.05}), 2.4648, 5e-5);
}

// ---- No-arbitrage identities -----------------------------------------------

TEST(BlackScholesIdentities, PutCallParityHoldsAcrossGrid) {
    // C - P = S e^{-qT} - K e^{-rT}, model-free.
    for (const Case& c : grid()) {
        SCOPED_TRACE(::testing::Message() << "S=" << c.S << " T=" << c.T << " vol=" << c.vol
                                          << " r=" << c.r << " q=" << c.q);
        const Market mkt{c.S, c.r, c.vol, c.q};
        const double lhs = bs::price(call(c.K, c.T), mkt) - bs::price(put(c.K, c.T), mkt);
        const double rhs = c.S * std::exp(-c.q * c.T) - c.K * std::exp(-c.r * c.T);
        EXPECT_NEAR(lhs, rhs, 1e-10 * c.S);
    }
}

TEST(BlackScholesIdentities, PricesRespectStaticArbitrageBounds) {
    // max(S e^{-qT} - K e^{-rT}, 0) <= C <= S e^{-qT}
    // max(K e^{-rT} - S e^{-qT}, 0) <= P <= K e^{-rT}
    for (const Case& c : grid()) {
        const Market mkt{c.S, c.r, c.vol, c.q};
        const double fwd_s = c.S * std::exp(-c.q * c.T);
        const double pv_k = c.K * std::exp(-c.r * c.T);
        const double C = bs::price(call(c.K, c.T), mkt);
        const double P = bs::price(put(c.K, c.T), mkt);
        const double eps = 1e-12 * c.S;
        EXPECT_GE(C, std::max(fwd_s - pv_k, 0.0) - eps);
        EXPECT_LE(C, fwd_s + eps);
        EXPECT_GE(P, std::max(pv_k - fwd_s, 0.0) - eps);
        EXPECT_LE(P, pv_k + eps);
    }
}

TEST(BlackScholesIdentities, GreekParityRelations) {
    // Differentiating put-call parity: call and put share gamma and vega, and
    // delta_C - delta_P = e^{-qT}, rho_C - rho_P = T K e^{-rT}.
    for (const Case& c : grid()) {
        const Market mkt{c.S, c.r, c.vol, c.q};
        const Greeks gc = bs::greeks(call(c.K, c.T), mkt);
        const Greeks gp = bs::greeks(put(c.K, c.T), mkt);
        EXPECT_NEAR(gc.gamma, gp.gamma, 1e-12);
        EXPECT_NEAR(gc.vega, gp.vega, 1e-10);
        EXPECT_NEAR(gc.delta - gp.delta, std::exp(-c.q * c.T), 1e-12);
        EXPECT_NEAR(gc.rho - gp.rho, c.T * c.K * std::exp(-c.r * c.T), 1e-10);
        EXPECT_GT(gc.gamma, 0.0);
        EXPECT_GT(gc.vega, 0.0);
    }
}

// ---- Analytic Greeks vs numerical derivatives of the analytic price --------
// Independent of the Greek formulas: if a sign or a term is wrong in greeks(),
// these disagree. Central differences have O(h^2) truncation error, so the
// bump sizes and tolerances below are chosen with that in mind.

TEST(BlackScholesGreeks, MatchCentralFiniteDifferencesOfPrice) {
    for (const Case& c : grid()) {
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            SCOPED_TRACE(::testing::Message()
                         << (type == OptionType::Call ? "call" : "put") << " S=" << c.S
                         << " T=" << c.T << " vol=" << c.vol << " r=" << c.r << " q=" << c.q);
            const EuropeanOption opt{type, c.K, c.T};
            const Market mkt{c.S, c.r, c.vol, c.q};
            const Greeks g = bs::greeks(opt, mkt);

            auto px = [&](double S, double T, double r, double vol) {
                return bs::price(EuropeanOption{type, c.K, T}, Market{S, r, vol, c.q});
            };
            auto delta_at = [&](double S) {
                return bs::greeks(opt, Market{S, c.r, c.vol, c.q}).delta;
            };

            // Spot bump scaled to the terminal distribution's width, S sigma
            // sqrt(T): a fixed 1e-4 S bump is ~1% of the width for a 5%-vol,
            // 0.05y option, where gamma is so peaked that O(h^2) truncation
            // alone (~1.3e-5 rel) exceeds the tolerance.
            const double hS = 1e-3 * c.S * c.vol * std::sqrt(c.T);
            const double hT = 1e-5;
            const double hR = 1e-5;
            const double hV = 1e-5;

            const double fd_delta = (px(c.S + hS, c.T, c.r, c.vol) -
                                     px(c.S - hS, c.T, c.r, c.vol)) / (2 * hS);
            const double fd_gamma = (delta_at(c.S + hS) - delta_at(c.S - hS)) / (2 * hS);
            const double fd_vega = (px(c.S, c.T, c.r, c.vol + hV) -
                                    px(c.S, c.T, c.r, c.vol - hV)) / (2 * hV);
            // theta = dV/dt = -dV/dT
            const double fd_theta = -(px(c.S, c.T + hT, c.r, c.vol) -
                                      px(c.S, c.T - hT, c.r, c.vol)) / (2 * hT);
            const double fd_rho = (px(c.S, c.T, c.r + hR, c.vol) -
                                   px(c.S, c.T, c.r - hR, c.vol)) / (2 * hR);

            // Relative-plus-absolute tolerance: Greeks span many orders of
            // magnitude across the grid (far-OTM short-dated gamma ~ 1e-30).
            auto tol = [](double x) { return 1e-6 + 1e-5 * std::abs(x); };
            EXPECT_NEAR(g.delta, fd_delta, tol(g.delta));
            EXPECT_NEAR(g.gamma, fd_gamma, tol(g.gamma));
            EXPECT_NEAR(g.vega, fd_vega, tol(g.vega));
            EXPECT_NEAR(g.theta, fd_theta, tol(g.theta));
            EXPECT_NEAR(g.rho, fd_rho, tol(g.rho));
        }
    }
}

// ---- Degenerate inputs -----------------------------------------------------

TEST(BlackScholesDegenerate, ZeroExpiryIsIntrinsicValue) {
    const Market mkt{110.0, 0.05, 0.30, 0.02};
    EXPECT_DOUBLE_EQ(bs::price(call(100.0, 0.0), mkt), 10.0);
    EXPECT_DOUBLE_EQ(bs::price(put(100.0, 0.0), mkt), 0.0);
    EXPECT_DOUBLE_EQ(bs::price(call(120.0, 0.0), mkt), 0.0);
    EXPECT_DOUBLE_EQ(bs::price(put(120.0, 0.0), mkt), 10.0);

    const Greeks itm = bs::greeks(call(100.0, 0.0), mkt);
    EXPECT_DOUBLE_EQ(itm.delta, 1.0);
    EXPECT_DOUBLE_EQ(itm.gamma, 0.0);
    EXPECT_DOUBLE_EQ(itm.vega, 0.0);
    EXPECT_DOUBLE_EQ(itm.rho, 0.0);
    const Greeks otm = bs::greeks(call(120.0, 0.0), mkt);
    EXPECT_DOUBLE_EQ(otm.delta, 0.0);
}

TEST(BlackScholesDegenerate, ZeroVolIsDiscountedForwardIntrinsic) {
    // With sigma = 0 the stock grows deterministically at r - q, so the option
    // is worth its discounted payoff on the forward.
    const double S = 100.0, K = 100.0, T = 1.0, r = 0.05, q = 0.01;
    const Market mkt{S, r, 0.0, q};
    const double fwd_intrinsic = S * std::exp(-q * T) - K * std::exp(-r * T);
    ASSERT_GT(fwd_intrinsic, 0.0);
    EXPECT_NEAR(bs::price(call(K, T), mkt), fwd_intrinsic, 1e-12);
    EXPECT_DOUBLE_EQ(bs::price(put(K, T), mkt), 0.0);

    const Greeks g = bs::greeks(call(K, T), mkt);
    EXPECT_NEAR(g.delta, std::exp(-q * T), 1e-15);
    EXPECT_DOUBLE_EQ(g.gamma, 0.0);
    EXPECT_DOUBLE_EQ(g.vega, 0.0);
    EXPECT_NEAR(g.rho, T * K * std::exp(-r * T), 1e-12);
}

TEST(BlackScholesDegenerate, ClosedFormConvergesToDegenerateLimit) {
    // The special-cased branch should be the limit of the general formula,
    // not a discontinuity: shrinking sigma or T approaches it smoothly.
    const Market base{105.0, 0.04, 0.0, 0.01};
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const EuropeanOption opt{type, 100.0, 0.75};
        const double limit = bs::price(opt, base);
        const Greeks glimit = bs::greeks(opt, base);
        Market tiny = base;
        tiny.vol = 1e-6;
        EXPECT_NEAR(bs::price(opt, tiny), limit, 1e-9);
        const Greeks gt = bs::greeks(opt, tiny);
        EXPECT_NEAR(gt.delta, glimit.delta, 1e-9);
        EXPECT_NEAR(gt.theta, glimit.theta, 1e-6);
        EXPECT_NEAR(gt.rho, glimit.rho, 1e-6);
    }

    const Market vol_mkt{105.0, 0.04, 0.3, 0.01};
    EXPECT_NEAR(bs::price(call(100.0, 1e-10), vol_mkt), 5.0, 1e-6);
    EXPECT_NEAR(bs::price(put(100.0, 1e-10), vol_mkt), 0.0, 1e-6);
}

TEST(BlackScholesDegenerate, AtTheForwardKinkDeltaIsOneHalf) {
    // S e^{-qT} == K e^{-rT} with sigma = 0: delta takes the midpoint, which
    // is also what N(d1) tends to as sigma -> 0 at the forward.
    const Market mkt{100.0, 0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(bs::price(call(100.0, 1.0), mkt), 0.0);
    EXPECT_DOUBLE_EQ(bs::greeks(call(100.0, 1.0), mkt).delta, 0.5);
    EXPECT_DOUBLE_EQ(bs::greeks(put(100.0, 1.0), mkt).delta, -0.5);
}

TEST(BlackScholesDegenerate, FarOtmTailIsPositiveNotZero) {
    // The erfc-based N(x) keeps precision in the tail; a naive
    // 0.5 * (1 + erf(x)) version rounds both N(d1) and N(d2) to exactly 0 here.
    const double px = bs::price(call(300.0, 0.25), Market{100.0, 0.0, 0.2});
    EXPECT_GT(px, 0.0);
    EXPECT_LT(px, 1e-20);
}

TEST(BlackScholesValidation, RejectsInvalidInputs) {
    const Market ok{100.0, 0.05, 0.2};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(bs::price(call(100.0, 1.0), Market{0.0, 0.05, 0.2}), std::invalid_argument);
    EXPECT_THROW(bs::price(call(100.0, 1.0), Market{-1.0, 0.05, 0.2}), std::invalid_argument);
    EXPECT_THROW(bs::price(call(0.0, 1.0), ok), std::invalid_argument);
    EXPECT_THROW(bs::price(call(100.0, -0.1), ok), std::invalid_argument);
    EXPECT_THROW(bs::price(call(100.0, 1.0), Market{100.0, 0.05, -0.2}), std::invalid_argument);
    EXPECT_THROW(bs::price(call(100.0, 1.0), Market{nan, 0.05, 0.2}), std::invalid_argument);
    EXPECT_THROW(bs::greeks(call(nan, 1.0), ok), std::invalid_argument);
}

}  // namespace
}  // namespace ope
