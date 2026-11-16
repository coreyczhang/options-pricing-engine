#include "ope/MonteCarlo.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "ope/BlackScholes.hpp"
#include "ope/NormalRng.hpp"

namespace ope {
namespace {

// These are statistical tests made deterministic by fixed seeds. Thresholds
// are set at ~4 standard errors (or the equivalent for the statistic under
// test), so they are not tuned to a lucky seed: for a random seed each check
// would fail with probability ~1e-4. A failure therefore points at a real
// bias, not noise.
constexpr double kSigmas = 4.0;

struct Case {
    OptionType type;
    double S, K, T, r, vol, q;
};

std::vector<Case> cases() {
    return {
        {OptionType::Call, 100.0, 100.0, 1.0, 0.05, 0.20, 0.00},   // ATM
        {OptionType::Put, 100.0, 100.0, 1.0, 0.05, 0.20, 0.00},
        {OptionType::Call, 100.0, 130.0, 0.5, 0.03, 0.25, 0.02},   // OTM call
        {OptionType::Put, 100.0, 130.0, 0.5, 0.03, 0.25, 0.02},    // ITM put
        {OptionType::Call, 100.0, 70.0, 2.0, 0.01, 0.40, 0.03},    // deep ITM, high vol
        {OptionType::Put, 100.0, 90.0, 0.1, 0.05, 0.15, 0.00},     // short-dated OTM put
        {OptionType::Call, 42.0, 40.0, 0.5, 0.10, 0.20, 0.00},     // Hull example
    };
}

double bs_price(const Case& c) {
    return bs::price(EuropeanOption{c.type, c.K, c.T}, Market{c.S, c.r, c.vol, c.q});
}

mc::Result mc_price(const Case& c, std::size_t paths, std::uint64_t seed,
                    bool antithetic = false) {
    return mc::price(EuropeanOption{c.type, c.K, c.T}, Market{c.S, c.r, c.vol, c.q},
                     mc::Config{paths, seed, antithetic});
}

// Variance reduction factor at equal payoff evaluations.
double vrf(const Case& c, std::size_t paths, std::uint64_t seed) {
    const double ratio =
        mc_price(c, paths, seed, false).std_error / mc_price(c, paths, seed, true).std_error;
    return ratio * ratio;
}

// ---- Agreement with the closed form ---------------------------------------

TEST(MonteCarlo, AgreesWithBlackScholesWithinStdErrors) {
    for (const Case& c : cases()) {
        SCOPED_TRACE(::testing::Message() << (c.type == OptionType::Call ? "call" : "put")
                                          << " S=" << c.S << " K=" << c.K << " T=" << c.T);
        const mc::Result r = mc_price(c, 200'000, 7);
        EXPECT_GT(r.std_error, 0.0);
        EXPECT_NEAR(r.price, bs_price(c), kSigmas * r.std_error);
    }
}

TEST(MonteCarlo, ErrorToleranceShrinksAsOneOverSqrtPaths) {
    // Each 10x more paths should shrink the standard error by sqrt(10), and
    // the estimate must stay inside the (shrinking) 4-sigma band at every N.
    const Case c = cases()[0];
    const double exact = bs_price(c);
    double prev_se = 0.0;
    for (std::size_t n : {1'000u, 10'000u, 100'000u, 1'000'000u}) {
        SCOPED_TRACE(::testing::Message() << "paths=" << n);
        const mc::Result r = mc_price(c, n, 11);
        EXPECT_NEAR(r.price, exact, kSigmas * r.std_error);
        if (prev_se > 0.0) {
            // The se estimate itself is noisy (~1/sqrt(2n) relative); 10% is
            // generous at n >= 1000.
            EXPECT_NEAR(prev_se / r.std_error, std::sqrt(10.0), 0.1 * std::sqrt(10.0));
        }
        prev_se = r.std_error;
    }
    // Concrete scale: ATM 1y 20%-vol call (~10.45) to within ~1 cent at 1M paths.
    EXPECT_LT(prev_se, 0.02);
}

// ---- The error bar is honest ----------------------------------------------
// A pricer that reports a wrong std_error would pass the tests above for the
// wrong reason, so check the reported number against the observed spread of
// many independent runs. Run for both modes: the antithetic se is the easy
// one to get wrong (it must be computed over pairs, not the 2N halves).

class MonteCarloStdError : public ::testing::TestWithParam<bool> {};

TEST_P(MonteCarloStdError, ReportedStdErrorMatchesSpreadAcrossSeeds) {
    const bool antithetic = GetParam();
    const Case c = cases()[0];
    const double exact = bs_price(c);
    constexpr int kRuns = 400;
    constexpr std::size_t kPaths = 2'000;

    double sum = 0.0, sum_sq_err = 0.0, sum_se = 0.0;
    int covered = 0;
    for (int seed = 1; seed <= kRuns; ++seed) {
        const mc::Result r =
            mc_price(c, kPaths, static_cast<std::uint64_t>(seed), antithetic);
        const double err = r.price - exact;
        sum += r.price;
        sum_sq_err += err * err;
        sum_se += r.std_error;
        if (std::abs(err) <= 1.96 * r.std_error) ++covered;
    }
    const double mean_price = sum / kRuns;
    const double rmse = std::sqrt(sum_sq_err / kRuns);
    const double mean_se = sum_se / kRuns;

    // Unbiased: the average of 400 estimates sits within 4 se/sqrt(400) of BS.
    EXPECT_NEAR(mean_price, exact, kSigmas * mean_se / std::sqrt(double(kRuns)));
    // Calibrated: actual RMS error ~= reported se. With 400 runs the RMSE is
    // known to ~1/sqrt(800) = 3.5%; allow 15% (~4 sigma).
    EXPECT_NEAR(rmse / mean_se, 1.0, 0.15);
    // 95% CIs cover the truth ~95% of the time: expected 380 of 400,
    // binomial sd ~4.4, so +/- 18 is ~4 sigma.
    EXPECT_NEAR(covered, 380, 18);
}

INSTANTIATE_TEST_SUITE_P(PlainAndAntithetic, MonteCarloStdError, ::testing::Bool(),
                         [](const auto& info) { return info.param ? "Antithetic" : "Plain"; });

// ---- Antithetic variates ---------------------------------------------------

TEST(MonteCarloAntithetic, AgreesWithBlackScholesWithinStdErrors) {
    for (const Case& c : cases()) {
        SCOPED_TRACE(::testing::Message() << (c.type == OptionType::Call ? "call" : "put")
                                          << " S=" << c.S << " K=" << c.K << " T=" << c.T);
        const mc::Result r = mc_price(c, 200'000, 7, true);
        EXPECT_GT(r.std_error, 0.0);
        EXPECT_NEAR(r.price, bs_price(c), kSigmas * r.std_error);
    }
}

TEST(MonteCarloAntithetic, ReducesVarianceMostWhenPayoffIsNearLinear) {
    // VRF = 1 / (1 + corr(f(Z), f(-Z))). Deep ITM a call is nearly linear in
    // S_T (rho -> -1, big win); ATM the kink makes it ~half linear (rho ~ -0.5,
    // VRF ~ 2); far OTM almost no path pair has both halves paying, so
    // rho ~ 0 and VRF ~ 1. With 1M paths each VRF estimate is good to a few %.
    const auto call = [](double K) {
        return Case{OptionType::Call, 100.0, K, 1.0, 0.05, 0.20, 0.0};
    };
    const double deep_itm = vrf(call(70.0), 1'000'000, 3);
    const double atm = vrf(call(100.0), 1'000'000, 3);
    const double far_otm = vrf(call(160.0), 1'000'000, 3);
    EXPECT_GT(deep_itm, 10.0);
    EXPECT_NEAR(atm, 2.0, 0.2);
    EXPECT_NEAR(far_otm, 1.0, 0.1);
    EXPECT_GT(deep_itm, atm);
    EXPECT_GT(atm, far_otm);
}

TEST(MonteCarloAntithetic, NeverMeaningfullyHurtsMonotonePayoffs) {
    // Calls and puts are monotone in Z, so f(Z) and f(-Z) are negatively
    // correlated (rho <= 0) and VRF >= 1 in theory. Measured VRF can dip
    // below 1 from estimation noise when the true value is ~1, hence the 0.9
    // floor rather than 1.0. That noise is worst for rare, heavy-tailed
    // payoffs: for the short-dated OTM put (true VRF ~1) the estimate's sd is
    // ~0.05 at 200k paths — 1.5% of seeds dipped below 0.9 — and ~0.016 at
    // 2M, which puts the floor ~6 sd away.
    for (const Case& c : cases()) {
        SCOPED_TRACE(::testing::Message() << (c.type == OptionType::Call ? "call" : "put")
                                          << " K=" << c.K << " T=" << c.T);
        EXPECT_GT(vrf(c, 2'000'000, 5), 0.9);
    }
}

TEST(MonteCarloAntithetic, ZeroVolIsExactWithZeroError) {
    const EuropeanOption opt{OptionType::Put, 110.0, 1.0};
    const Market mkt{100.0, 0.02, 0.0};
    const mc::Result r = mc::price(opt, mkt, mc::Config{1'000, 42, true});
    EXPECT_NEAR(r.price, bs::price(opt, mkt), 1e-12);
    EXPECT_EQ(r.std_error, 0.0);
}

TEST(MonteCarloAntithetic, RejectsOddOrSinglePairPathCounts) {
    const EuropeanOption opt{OptionType::Call, 100.0, 1.0};
    const Market mkt{100.0, 0.05, 0.2};
    EXPECT_THROW(mc::price(opt, mkt, mc::Config{1'001, 42, true}), std::invalid_argument);
    // One pair has no sample variance: it would report a misleading 0 +/- 0.
    EXPECT_THROW(mc::price(opt, mkt, mc::Config{2, 42, true}), std::invalid_argument);
    EXPECT_NO_THROW(mc::price(opt, mkt, mc::Config{4, 42, true}));
    // Plain MC has no parity requirement.
    EXPECT_NO_THROW(mc::price(opt, mkt, mc::Config{1'001, 42, false}));
}

// ---- Determinism -----------------------------------------------------------

TEST(MonteCarlo, PriceAndStdErrorAreTheSampleStatisticsOfTheDiscountedPayoffs) {
    // Recompute tiny MC runs by hand from the same normal draws: the price is
    // the mean of e^{-rT} payoff(S_T), the std error is the sample standard
    // deviation (n - 1 denominator) over sqrt(n), and for antithetic the
    // sampling unit is the pair average. Deterministic, so it pins down
    // details below any statistical test's resolution: dividing by n instead
    // of n - 1 changes the std error by ~1/(2n), 2.5e-6 relative at 200k
    // paths, and the mutation check (scripts/mutation_check.py) confirmed no
    // statistical test here could see it. At n = 5 it is an 11% change.
    const Case c = cases()[4];  // deep ITM call with q != 0: every path pays
    const double T = c.T;
    const double drift = (c.r - c.q - 0.5 * c.vol * c.vol) * T;
    const double disc = std::exp(-c.r * T);
    auto value = [&](double z) {
        return disc * std::max(c.S * std::exp(drift + c.vol * std::sqrt(T) * z) - c.K, 0.0);
    };
    for (bool antithetic : {false, true}) {
        const std::size_t paths = antithetic ? 10 : 5;
        NormalRng rng(123);
        std::vector<double> units;
        for (std::size_t i = 0; i < (antithetic ? paths / 2 : paths); ++i) {
            const double z = rng.next();
            units.push_back(antithetic ? 0.5 * (value(z) + value(-z)) : value(z));
        }
        double mean = 0.0;
        for (double u : units) mean += u / units.size();
        double ss = 0.0;
        for (double u : units) ss += (u - mean) * (u - mean);
        const double se = std::sqrt(ss / (units.size() - 1.0) / units.size());

        const mc::Result r = mc::price(EuropeanOption{c.type, c.K, c.T},
                                       Market{c.S, c.r, c.vol, c.q},
                                       mc::Config{paths, 123, antithetic});
        ASSERT_GT(se, 0.0);  // some paths finish in the money
        EXPECT_NEAR(r.price, mean, 1e-12 * std::max(1.0, mean)) << "antithetic=" << antithetic;
        EXPECT_NEAR(r.std_error, se, 1e-12 * se) << "antithetic=" << antithetic;
    }
}

TEST(MonteCarlo, SameSeedReproducesExactly) {
    const Case c = cases()[0];
    const mc::Result a = mc_price(c, 10'000, 123);
    const mc::Result b = mc_price(c, 10'000, 123);
    EXPECT_EQ(a.price, b.price);
    EXPECT_EQ(a.std_error, b.std_error);
    EXPECT_NE(mc_price(c, 10'000, 124).price, a.price);
}

// ---- Known limitation: rare events ----------------------------------------

TEST(MonteCarloLimitation, FarOtmWithNoHitsReportsZeroPlusMinusZero) {
    // Pins down a property of plain MC rather than a feature: a put ~4.7 sd
    // out of the money pays off on ~1 path in a million. With 200k paths none
    // hit, so the estimate is 0 and the sample-variance error bar is also 0 —
    // confidently wrong. The CLT-based std_error is only meaningful once
    // enough paths land in the payoff region. Importance sampling (shifting
    // the drift toward the strike) is the standard fix; out of scope for v1.
    const EuropeanOption opt{OptionType::Put, 80.0, 0.1};
    const Market mkt{100.0, 0.05, 0.15};
    const double exact = bs::price(opt, mkt);
    ASSERT_GT(exact, 1e-7);
    ASSERT_LT(exact, 1e-6);
    const mc::Result r = mc::price(opt, mkt, mc::Config{200'000, 7});
    EXPECT_EQ(r.price, 0.0);
    EXPECT_EQ(r.std_error, 0.0);
}

// ---- Degenerate inputs -----------------------------------------------------

TEST(MonteCarlo, ZeroVolIsExactWithZeroError) {
    // No randomness left: every path is the forward, so MC == BS limit exactly.
    const EuropeanOption opt{OptionType::Call, 100.0, 1.0};
    const Market mkt{100.0, 0.05, 0.0, 0.01};
    const mc::Result r = mc::price(opt, mkt, mc::Config{1'000});
    EXPECT_NEAR(r.price, bs::price(opt, mkt), 1e-12);
    EXPECT_EQ(r.std_error, 0.0);
}

TEST(MonteCarlo, ZeroExpiryIsIntrinsicWithZeroError) {
    const Market mkt{110.0, 0.05, 0.3};
    const mc::Result c = mc::price(EuropeanOption{OptionType::Call, 100.0, 0.0}, mkt, {1'000});
    const mc::Result p = mc::price(EuropeanOption{OptionType::Put, 100.0, 0.0}, mkt, {1'000});
    EXPECT_DOUBLE_EQ(c.price, 10.0);
    EXPECT_DOUBLE_EQ(p.price, 0.0);
    EXPECT_EQ(c.std_error, 0.0);
}

TEST(MonteCarlo, RejectsInvalidInputs) {
    const EuropeanOption opt{OptionType::Call, 100.0, 1.0};
    const Market ok{100.0, 0.05, 0.2};
    EXPECT_THROW(mc::price(opt, ok, mc::Config{1}), std::invalid_argument);
    EXPECT_THROW(mc::price(opt, ok, mc::Config{0}), std::invalid_argument);
    EXPECT_THROW(mc::price(opt, Market{-1.0, 0.05, 0.2}, mc::Config{100}),
                 std::invalid_argument);
}

}  // namespace
}  // namespace ope
