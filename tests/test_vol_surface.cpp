#include "ope/VolSurface.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "ope/BlackScholes.hpp"

namespace ope {
namespace {

using surface::Quote;
using surface::SkipReason;
using surface::Violation;
using surface::VolSurface;

constexpr double kSpot = 100.0;
constexpr double kRate = 0.03;
constexpr double kDiv = 0.01;

// ---- A known, arbitrage-free smile to generate quotes from ------------------
// SSVI (Gatheral & Jacquier 2014): total variance
//   w(k, T) = theta/2 * (1 + rho phi k + sqrt((phi k + rho)^2 + 1 - rho^2))
// with ATM variance theta(T) and phi(theta) = eta / (theta^g (1 + theta)^(1-g)).
// For g = 1/2 it is free of static arbitrage when eta (1 + |rho|) <= 2 and
// theta is increasing; here eta = 1, rho = -0.6 (equity-style skew).
double theta(double T) { return 0.04 * T + 0.02 * (1.0 - std::exp(-2.0 * T)); }

double ssvi_w(double k, double T) {
    const double th = theta(T);
    const double rho = -0.6;
    const double phi = 1.0 / (std::sqrt(th) * std::sqrt(1.0 + th));
    return th / 2.0 * (1.0 + rho * phi * k + std::sqrt((phi * k + rho) * (phi * k + rho) +
                                                       1.0 - rho * rho));
}

double forward(double T) { return kSpot * std::exp((kRate - kDiv) * T); }

double ssvi_vol(double K, double T) { return std::sqrt(ssvi_w(std::log(K / forward(T)), T) / T); }

double price(OptionType type, double K, double T, double vol) {
    return bs::price({type, K, T}, {kSpot, kRate, vol, kDiv});
}

// Calls and puts at every strike in [60, 160] spaced dK, at each expiry.
std::vector<Quote> ssvi_quotes(double dK, const std::vector<double>& expiries) {
    std::vector<Quote> out;
    for (double T : expiries)
        for (double K = 60.0; K <= 160.0 + 1e-9; K += dK)
            for (OptionType type : {OptionType::Call, OptionType::Put})
                out.push_back({type, K, T, price(type, K, T, ssvi_vol(K, T))});
    return out;
}

const std::vector<double> kExpiries = {1.0 / 12, 0.25, 0.5, 1.0, 2.0};

// Worst |surface vol - SSVI vol| halfway between quoted strikes.
double mid_strike_error(double dK) {
    const VolSurface s(ssvi_quotes(dK, kExpiries), kSpot, kRate, kDiv);
    double worst = 0.0;
    for (double T : kExpiries)
        for (double K = 60.0 + dK / 2; K < 160.0; K += dK)
            worst = std::max(worst, std::abs(s.vol(K, T) - ssvi_vol(K, T)));
    return worst;
}

// ---- Nodes: the quotes themselves -------------------------------------------

TEST(VolSurfaceNodes, ReproduceTheGeneratingVolsAndUseOtmQuotes) {
    const VolSurface s(ssvi_quotes(5.0, kExpiries), kSpot, kRate, kDiv);
    ASSERT_EQ(s.slices().size(), kExpiries.size());
    for (const surface::Slice& slice : s.slices()) {
        EXPECT_EQ(slice.nodes.size(), 21u);  // 60, 65, ..., 160
        for (const surface::Node& n : slice.nodes) {
            EXPECT_NEAR(n.vol, ssvi_vol(n.strike, slice.expiry), 1e-12);
            EXPECT_EQ(n.type, n.strike >= slice.forward ? OptionType::Call : OptionType::Put);
        }
    }
    // Every strike/expiry had a call and a put; exactly one of each pair is used.
    EXPECT_EQ(s.skipped().size(), 21u * kExpiries.size());
    for (const surface::Skipped& sk : s.skipped())
        EXPECT_EQ(sk.reason, SkipReason::ReplacedByOtm);
}

TEST(VolSurfaceNodes, CallsOnlyOrPutsOnlyGiveTheSameSurface) {
    // Put-call parity: consistent quotes imply one vol per strike, so which
    // side was quoted must not matter (beyond ITM conditioning).
    std::vector<Quote> calls, puts;
    for (const Quote& q : ssvi_quotes(10.0, {0.5, 1.0}))
        (q.type == OptionType::Call ? calls : puts).push_back(q);
    const VolSurface c(calls, kSpot, kRate, kDiv), p(puts, kSpot, kRate, kDiv);
    for (double K : {65.0, 90.0, 100.0, 117.0, 150.0})
        for (double T : {0.5, 0.75, 1.0}) EXPECT_NEAR(c.vol(K, T), p.vol(K, T), 1e-9);
}

// ---- Interpolation ------------------------------------------------------------

TEST(VolSurfaceInterpolation, StrikeErrorShrinksQuadraticallyWithSpacing) {
    // Linear interpolation in total variance has O(h^2) error; halving the
    // strike spacing should cut the worst mid-strike error ~4x.
    const double e5 = mid_strike_error(5.0);
    const double e25 = mid_strike_error(2.5);
    const double e125 = mid_strike_error(1.25);
    EXPECT_GT(e5 / e25, 3.0);
    EXPECT_LT(e5 / e25, 5.0);
    EXPECT_GT(e25 / e125, 3.0);
    EXPECT_LT(e25 / e125, 5.0);
    EXPECT_LT(e25, 1.5e-3);  // ~0.1 vol point at $2.50 strikes (measured 1.07e-3)
}

TEST(VolSurfaceInterpolation, ExpiryErrorShrinksAsExpiriesAreAdded) {
    auto between_expiry_error = [](const std::vector<double>& expiries) {
        const VolSurface s(ssvi_quotes(2.5, expiries), kSpot, kRate, kDiv);
        double worst = 0.0;
        for (double T : {0.5625, 0.625, 0.6875, 0.8125, 0.875, 0.9375})
            for (double K = 70.0; K <= 140.0; K += 2.5)
                worst = std::max(worst, std::abs(s.vol(K, T) - ssvi_vol(K, T)));
        return worst;
    };
    const double e2 = between_expiry_error({0.5, 1.0});
    const double e3 = between_expiry_error({0.5, 0.75, 1.0});
    const double e5 = between_expiry_error({0.5, 0.625, 0.75, 0.875, 1.0});
    // Measured 2.8e-3, 1.0e-3, 3.2e-4: approaching the O(h^2) ratio of 4.
    EXPECT_GT(e2 / e3, 2.0);
    EXPECT_GT(e3 / e5, 2.5);
    EXPECT_LT(e5, 5e-4);
}

TEST(VolSurfaceInterpolation, ExtrapolatesWithConstantVol) {
    const VolSurface s(ssvi_quotes(10.0, {0.5, 1.0}), kSpot, kRate, kDiv);
    const surface::Slice& first = s.slices().front();
    const surface::Slice& last = s.slices().back();
    // Beyond the outermost strikes: flat at the edge node's vol.
    EXPECT_NEAR(s.vol(400.0, 1.0), last.nodes.back().vol, 1e-14);
    EXPECT_NEAR(s.vol(20.0, 1.0), last.nodes.front().vol, 1e-14);
    // Before the first expiry and after the last: the slice's vol at the same
    // log-moneyness.
    for (double k : {-0.2, 0.0, 0.15}) {
        EXPECT_NEAR(s.vol(forward(0.1) * std::exp(k), 0.1),
                    s.vol(first.forward * std::exp(k), first.expiry), 1e-14);
        EXPECT_NEAR(s.vol(forward(3.0) * std::exp(k), 3.0),
                    s.vol(last.forward * std::exp(k), last.expiry), 1e-14);
    }
}

// ---- Bad quotes ------------------------------------------------------------------

TEST(VolSurfaceQuotes, BadQuotesAreSkippedWithReasonsAndDontMoveTheSurface) {
    const std::vector<Quote> clean = ssvi_quotes(10.0, {0.5, 1.0});
    std::vector<Quote> dirty = clean;
    const double T = 1.0;
    dirty.push_back({OptionType::Call, 100.0, T, 999.0});   // duplicate of a clean call
    dirty.push_back({OptionType::Put, 105.0, T, 200.0});    // > K e^{-rT}
    dirty.push_back({OptionType::Call, 75.0, T, 1.0});      // < forward intrinsic
    dirty.push_back({OptionType::Call, 300.0, T, 0.0});     // zero bid, far OTM
    dirty.push_back({OptionType::Call, -5.0, T, 1.0});      // invalid
    dirty.push_back({OptionType::Put, 100.0, 0.0, 1.0});    // expired

    // Bad quotes leave no node, so the surface is bit-for-bit the clean one.
    const VolSurface a(clean, kSpot, kRate, kDiv), b(dirty, kSpot, kRate, kDiv);
    for (double K : {62.0, 80.0, 100.0, 105.0, 133.0})
        for (double t : {0.5, 0.8, 1.0}) EXPECT_EQ(a.vol(K, t), b.vol(K, t));

    auto count = [&](SkipReason r) {
        return std::count_if(b.skipped().begin(), b.skipped().end(),
                             [&](const surface::Skipped& s) { return s.reason == r; });
    };
    EXPECT_EQ(count(SkipReason::Duplicate), 1);
    EXPECT_EQ(count(SkipReason::AboveUpperBound), 1);
    EXPECT_EQ(count(SkipReason::BelowLowerBound), 1);
    EXPECT_EQ(count(SkipReason::AtLowerBound), 1);
    EXPECT_EQ(count(SkipReason::InvalidQuote), 2);
    EXPECT_EQ(count(SkipReason::ReplacedByOtm), static_cast<long>(clean.size() / 2));
}

TEST(VolSurfaceQuotes, FallsBackToTheItmQuoteWhenTheOtmOneIsBroken) {
    const double T = 1.0, K = 90.0;  // below the forward: the put is OTM
    const double vol = ssvi_vol(K, T);
    const std::vector<Quote> quotes = {
        {OptionType::Put, K, T, 1e6},  // junk OTM quote
        {OptionType::Call, K, T, price(OptionType::Call, K, T, vol)},
    };
    const VolSurface s(quotes, kSpot, kRate, kDiv);
    ASSERT_EQ(s.slices().size(), 1u);
    ASSERT_EQ(s.slices()[0].nodes.size(), 1u);
    EXPECT_EQ(s.slices()[0].nodes[0].type, OptionType::Call);
    EXPECT_NEAR(s.slices()[0].nodes[0].vol, vol, 1e-10);
    ASSERT_EQ(s.skipped().size(), 1u);
    EXPECT_EQ(s.skipped()[0].reason, SkipReason::AboveUpperBound);
}

// ---- Static-arbitrage diagnostics ---------------------------------------------

TEST(VolSurfaceArbitrage, SsviQuotesHaveNone) {
    for (double dK : {20.0, 5.0, 1.25}) {
        const VolSurface s(ssvi_quotes(dK, kExpiries), kSpot, kRate, kDiv);
        EXPECT_TRUE(s.arbitrage_violations().empty()) << "dK=" << dK;
    }
}

TEST(VolSurfaceArbitrage, DetectsCalendarArbitrage) {
    // A 2y slice quoted at a flat 10% vol has less total variance than the 1y
    // SSVI slice (ATM ~21%): a calendar spread would be free money.
    std::vector<Quote> quotes = ssvi_quotes(10.0, {1.0});
    for (double K = 60.0; K <= 160.0; K += 10.0)
        quotes.push_back({OptionType::Call, K, 2.0, price(OptionType::Call, K, 2.0, 0.10)});
    const VolSurface s(quotes, kSpot, kRate, kDiv);
    const std::vector<Violation> v = s.arbitrage_violations();
    ASSERT_FALSE(v.empty());
    for (const Violation& x : v) {
        EXPECT_EQ(x.kind, Violation::Kind::Calendar);
        EXPECT_EQ(x.expiry, 2.0);
        EXPECT_GT(x.size, 0.0);
    }
}

TEST(VolSurfaceArbitrage, DetectsButterflyArbitrage) {
    // Spike one strike's vol: its call price jumps above the line between its
    // neighbours, so the call price is no longer convex in strike there.
    std::vector<Quote> quotes;
    const double T = 0.5;
    for (double K = 80.0; K <= 120.0; K += 5.0) {
        const double vol = K == 100.0 ? 0.35 : ssvi_vol(K, T);
        quotes.push_back({OptionType::Call, K, T, price(OptionType::Call, K, T, vol)});
    }
    const VolSurface s(quotes, kSpot, kRate, kDiv);
    const std::vector<Violation> v = s.arbitrage_violations();
    ASSERT_FALSE(v.empty());
    for (const Violation& x : v) EXPECT_EQ(x.kind, Violation::Kind::Butterfly);
    EXPECT_TRUE(std::any_of(v.begin(), v.end(),
                            [](const Violation& x) { return x.strike == 100.0; }));
}

// ---- Invalid inputs ---------------------------------------------------------------

TEST(VolSurfaceValidation, RejectsUnusableInputs) {
    const std::vector<Quote> good = ssvi_quotes(10.0, {1.0});
    EXPECT_THROW(VolSurface({}, kSpot, kRate), std::invalid_argument);
    EXPECT_THROW(VolSurface({{OptionType::Call, 100.0, 1.0, 1e6}}, kSpot, kRate),
                 std::invalid_argument);  // nothing invertible
    EXPECT_THROW(VolSurface(good, -1.0, kRate), std::invalid_argument);
    EXPECT_THROW(VolSurface(good, kSpot, std::nan("")), std::invalid_argument);
    const VolSurface s(good, kSpot, kRate, kDiv);
    EXPECT_THROW(s.vol(0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(s.vol(100.0, 0.0), std::invalid_argument);
}

}  // namespace
}  // namespace ope
