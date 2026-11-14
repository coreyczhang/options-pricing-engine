#include "ope/ImpliedVol.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "ope/BlackScholes.hpp"

namespace ope {
namespace {

EuropeanOption call(double K, double T) { return {OptionType::Call, K, T}; }
EuropeanOption put(double K, double T) { return {OptionType::Put, K, T}; }

// Spot 100, so strikes are moneyness x 100. Spans 1% to 300% vol, 1 day to 5
// years, deep ITM to deep OTM: 2,016 contracts with call/put and q.
struct Case {
    EuropeanOption opt;
    Market mkt;  // mkt.vol is the true vol the quote is generated from
};

std::vector<Case> grid() {
    std::vector<Case> cases;
    for (OptionType type : {OptionType::Call, OptionType::Put})
        for (double vol : {0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3.0})
            for (double m : {0.5, 0.7, 0.85, 0.95, 1.0, 1.05, 1.2, 1.5, 2.0})
                for (double T : {1.0 / 365, 7.0 / 365, 0.1, 0.5, 1.0, 2.0, 5.0})
                    for (double q : {0.0, 0.03})
                        cases.push_back({{type, 100.0 * m, T}, {100.0, 0.05, vol, q}});
    return cases;
}

// Solve from the automatic starting guess (mkt.vol = 0).
iv::Result solve_auto(const EuropeanOption& opt, Market mkt, double price) {
    mkt.vol = 0.0;
    return iv::solve(opt, mkt, price);
}

// The best any solver can do: bs::price carries rounding error of order
// eps * scale, scale = max(S e^{-qT}, K e^{-rT}), and dividing a price error by
// vega gives the vol error it induces. The solver stops once
// |price residual| <= 8 eps (price + 2 * intrinsic) <= 24 eps * scale, plus a
// few eps * scale of rounding in the two prices compared: <= ~28 eps * scale /
// vega. 64x gives 2x headroom over that derived bound; measured peak is 10x.
double vol_error_bound(const Case& c) {
    const double scale = std::max(c.mkt.spot * std::exp(-c.mkt.div_yield * c.opt.expiry),
                                  c.opt.strike * std::exp(-c.mkt.rate * c.opt.expiry));
    return 1e-10 + 64.0 * DBL_EPSILON * scale / bs::greeks(c.opt, c.mkt).vega;
}

// ---- Reference value -------------------------------------------------------

TEST(ImpliedVolReference, HullExample) {
    // Hull — Black-Scholes-Merton chapter, implied volatility example:
    // S=21, K=20, r=10%, T=0.25, call price 1.875 -> implied vol 0.235.
    const iv::Result r = iv::solve(call(20.0, 0.25), {21.0, 0.10, 0.0}, 1.875);
    ASSERT_EQ(r.status, iv::Status::Converged);
    EXPECT_NEAR(r.vol, 0.235, 0.0005);
}

// ---- Round trip: price at sigma, solve, recover sigma ----------------------

TEST(ImpliedVolRoundTrip, RecoversVolToTheConditioningLimitOnGrid) {
    int converged = 0, at_bound = 0, max_iter = 0;
    for (const Case& c : grid()) {
        const double price = bs::price(c.opt, c.mkt);
        const iv::Result r = solve_auto(c.opt, c.mkt, price);
        SCOPED_TRACE(::testing::Message() << "K=" << c.opt.strike << " T=" << c.opt.expiry
                                          << " vol=" << c.mkt.vol << " q=" << c.mkt.div_yield
                                          << (c.opt.type == OptionType::Call ? " call" : " put"));
        max_iter = std::max(max_iter, r.iterations);
        if (r.status == iv::Status::AtLowerBound) {
            // Only legitimate when the time value is lost to rounding (or the
            // OTM price underflowed to 0): nothing left to infer a vol from.
            Market zero = c.mkt;
            zero.vol = 0.0;
            EXPECT_LE(price - bs::price(c.opt, zero), 1e-12 * c.mkt.spot);
            ++at_bound;
            continue;
        }
        ASSERT_EQ(r.status, iv::Status::Converged);
        EXPECT_NEAR(r.vol, c.mkt.vol, vol_error_bound(c));
        ++converged;
    }
    // Pin the split so a change in solver behaviour shows up here.
    EXPECT_EQ(converged + at_bound, 2016);
    EXPECT_GT(converged, 1600);
    EXPECT_LE(max_iter, 40);  // measured max: 30 (macOS and Linux)
}

TEST(ImpliedVolRoundTrip, WellConditionedQuotesRecoverToNearMachinePrecision) {
    // Near the money with non-trivial time value: vega is large, so the only
    // error left is the solver's own tolerance.
    for (double m : {0.9, 1.0, 1.1})
        for (double T : {0.25, 1.0})
            for (double vol : {0.1, 0.3, 0.6}) {
                const EuropeanOption opt = call(100.0 * m, T);
                const Market mkt{100.0, 0.03, vol, 0.01};
                const iv::Result r = solve_auto(opt, mkt, bs::price(opt, mkt));
                ASSERT_EQ(r.status, iv::Status::Converged);
                EXPECT_NEAR(r.vol, vol, 1e-12);
                // Measured: 4 near the money; up to 10 for low vol away from
                // it, where the vega-peak start (~0.9 at 10% OTM, 3 months)
                // sits far above the answer. A sharper initial guess would cut
                // that; convergence itself is guaranteed by the bracket.
                EXPECT_LE(r.iterations, 12);
            }
}

TEST(ImpliedVolRoundTrip, CallAndPutAtSameStrikeGiveSameVol) {
    // Put-call parity: in BS both are priced off one sigma, so a consistent
    // pair of quotes must imply the same vol. Holds where neither is lost to
    // rounding against intrinsic.
    for (double m : {0.8, 1.0, 1.25}) {
        const Market mkt{100.0, 0.05, 0.3, 0.02};
        const double vc = solve_auto(call(100.0 * m, 1.0), mkt,
                                     bs::price(call(100.0 * m, 1.0), mkt)).vol;
        const double vp = solve_auto(put(100.0 * m, 1.0), mkt,
                                     bs::price(put(100.0 * m, 1.0), mkt)).vol;
        EXPECT_NEAR(vc, vp, 1e-11);
    }
}

TEST(ImpliedVolRoundTrip, TinyOtmPriceStillHasAnImpliedVol) {
    // 1% vol, 5% OTM, one week: price ~1e-264. erfc keeps it accurate, and the
    // solver must not mistake it for "zero".
    const EuropeanOption opt = call(105.0, 7.0 / 365);
    const Market mkt{100.0, 0.05, 0.01};
    const double price = bs::price(opt, mkt);
    ASSERT_GT(price, 0.0);
    ASSERT_LT(price, 1e-200);
    const iv::Result r = solve_auto(opt, mkt, price);
    ASSERT_EQ(r.status, iv::Status::Converged);
    EXPECT_NEAR(r.vol, 0.01, 1e-12);
}

TEST(ImpliedVolRoundTrip, StartingGuessChangesIterationsNotTheAnswer) {
    const EuropeanOption opt = put(90.0, 0.5);
    const Market truth{100.0, 0.04, 0.35, 0.01};
    const double price = bs::price(opt, truth);
    const double automatic = solve_auto(opt, truth, price).vol;
    for (double guess : {1e-4, 0.05, 0.35, 2.0, 9.0}) {
        Market mkt = truth;
        mkt.vol = guess;
        const iv::Result r = iv::solve(opt, mkt, price);
        ASSERT_EQ(r.status, iv::Status::Converged) << "guess " << guess;
        EXPECT_NEAR(r.vol, automatic, 1e-12) << "guess " << guess;
    }
    // A guess at the answer converges immediately.
    Market exact = truth;
    EXPECT_LE(iv::solve(opt, exact, price).iterations, 2);
}

// ---- Quotes no vol can explain ----------------------------------------------

TEST(ImpliedVolBounds, PriceBelowIntrinsicIsArbitrage) {
    // Forward intrinsic of the call: 100 - 90 e^{-0.05}.
    const EuropeanOption opt = call(90.0, 1.0);
    const Market mkt{100.0, 0.05, 0.0};
    const double intrinsic = 100.0 - 90.0 * std::exp(-0.05);
    const iv::Result r = iv::solve(opt, mkt, intrinsic - 0.01);
    EXPECT_EQ(r.status, iv::Status::BelowLowerBound);
    EXPECT_TRUE(std::isnan(r.vol));
}

TEST(ImpliedVolBounds, PriceAtIntrinsicImpliesZeroVol) {
    const EuropeanOption opt = call(90.0, 1.0);
    const Market mkt{100.0, 0.05, 0.0};
    const double intrinsic = 100.0 - 90.0 * std::exp(-0.05);
    for (double p : {intrinsic, std::nextafter(intrinsic, 0.0)}) {  // 1 ulp under: rounding
        const iv::Result r = iv::solve(opt, mkt, p);
        EXPECT_EQ(r.status, iv::Status::AtLowerBound);
        EXPECT_EQ(r.vol, 0.0);
    }
    // OTM: a zero price is on the bound too.
    EXPECT_EQ(iv::solve(call(150.0, 0.1), mkt, 0.0).status, iv::Status::AtLowerBound);
}

TEST(ImpliedVolBounds, PriceAtOrAboveUpperBoundIsArbitrage) {
    // A call can never be worth more than the stock; a put never more than the
    // discounted strike. Both are only approached as sigma -> infinity.
    const Market mkt{100.0, 0.05, 0.0, 0.02};
    const double fwd_spot = 100.0 * std::exp(-0.02);
    EXPECT_EQ(iv::solve(call(100.0, 1.0), mkt, fwd_spot).status, iv::Status::AboveUpperBound);
    EXPECT_EQ(iv::solve(put(100.0, 1.0), mkt, 100.0 * std::exp(-0.05)).status,
              iv::Status::AboveUpperBound);
    EXPECT_TRUE(std::isnan(iv::solve(call(100.0, 1.0), mkt, 200.0).vol));
}

TEST(ImpliedVolBounds, VolBeyondTheCapIsReportedNotSilentlyClamped) {
    const EuropeanOption opt = call(100.0, 1.0);
    const Market truth{100.0, 0.0, 12.0};  // 1200% vol
    const double price = bs::price(opt, truth);
    EXPECT_EQ(solve_auto(opt, truth, price).status, iv::Status::VolAboveMax);

    iv::Settings wide;
    wide.max_vol = 20.0;
    const iv::Result r = iv::solve(opt, {100.0, 0.0, 0.0}, price, wide);
    ASSERT_EQ(r.status, iv::Status::Converged);
    EXPECT_NEAR(r.vol, 12.0, 1e-8);
}

// ---- Invalid inputs ---------------------------------------------------------

TEST(ImpliedVolValidation, RejectsUndefinedOrInvalidInputs) {
    const Market mkt{100.0, 0.05, 0.0};
    EXPECT_THROW(iv::solve(call(100.0, 0.0), mkt, 1.0), std::invalid_argument);  // T == 0
    EXPECT_THROW(iv::solve(call(100.0, 1.0), mkt, -1.0), std::invalid_argument);
    EXPECT_THROW(iv::solve(call(100.0, 1.0), mkt, std::nan("")), std::invalid_argument);
    EXPECT_THROW(iv::solve(call(100.0, 1.0), {-100.0, 0.05, 0.0}, 1.0),
                 std::invalid_argument);
    iv::Settings bad;
    bad.max_iter = 0;
    EXPECT_THROW(iv::solve(call(100.0, 1.0), mkt, 10.0, bad), std::invalid_argument);
}

}  // namespace
}  // namespace ope
