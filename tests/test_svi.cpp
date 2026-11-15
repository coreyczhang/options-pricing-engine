#include "ope/Svi.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "ope/BlackScholes.hpp"
#include "ope/NormalRng.hpp"
#include "ope/VolSurface.hpp"

namespace ope {
namespace {

// The SSVI surface of test_vol_surface.cpp. Every SSVI slice is exactly a raw
// SVI slice (Gatheral & Jacquier 2014, Lemma 3.2) with
//   a = theta (1 - rho^2) / 2,  b = theta phi / 2,  m = -rho / phi,
//   sigma = sqrt(1 - rho^2) / phi,
// so a correct fitter must recover these parameters from its quotes.
constexpr double kSpot = 100.0, kRate = 0.03, kDiv = 0.01, kRho = -0.6;

double theta(double T) { return 0.04 * T + 0.02 * (1.0 - std::exp(-2.0 * T)); }

svi::Params ssvi_slice(double T) {
    const double th = theta(T);
    const double phi = 1.0 / (std::sqrt(th) * std::sqrt(1.0 + th));
    return {th * (1.0 - kRho * kRho) / 2.0, th * phi / 2.0, kRho, -kRho / phi,
            std::sqrt(1.0 - kRho * kRho) / phi};
}

double forward(double T) { return kSpot * std::exp((kRate - kDiv) * T); }

// Log-moneyness of $5 strikes 60..160: 21 points.
std::vector<double> strikes_k(double T) {
    std::vector<double> k;
    for (double K = 60.0; K <= 160.0 + 1e-9; K += 5.0) k.push_back(std::log(K / forward(T)));
    return k;
}

const std::vector<double> kExpiries = {1.0 / 12, 0.25, 0.5, 1.0, 2.0};

// ---- Fitting a single slice ----------------------------------------------------

TEST(SviFit, RecoversTheParametersOfAnSsviSlice) {
    for (double T : kExpiries) {
        const svi::Params truth = ssvi_slice(T);
        const std::vector<double> k = strikes_k(T);
        std::vector<double> w;
        for (double x : k) w.push_back(truth.w(x));
        const svi::Fit f = svi::fit(k, w, T);
        SCOPED_TRACE(::testing::Message() << "T=" << T);
        EXPECT_NEAR(f.params.a, truth.a, 1e-7);
        EXPECT_NEAR(f.params.b, truth.b, 1e-7);
        EXPECT_NEAR(f.params.rho, truth.rho, 1e-6);
        EXPECT_NEAR(f.params.m, truth.m, 1e-6);
        EXPECT_NEAR(f.params.sigma, truth.sigma, 1e-6);
        EXPECT_LT(f.max_abs_vol, 1e-8);  // measured ~3e-9
        EXPECT_EQ(f.points, 21);
        EXPECT_GT(svi::min_durrleman_g(f.params, -1.5, 1.5).g, 0.0);
    }
}

TEST(SviFit, OnNoisyQuotesResidualAndErrorMatchLeastSquaresTheory) {
    // 20 bp of Gaussian noise on each vol. For a least-squares fit with p = 5
    // parameters to n = 21 points (near-linear near the optimum), theory says:
    //   residual vs quotes:  RMS ~ s sqrt((n - p) / n) = 0.873 s
    //   error vs the truth:  RMS ~ s sqrt(p / n)       = 0.488 s
    // The first shows the optimizer reaches the least-squares optimum (a stuck
    // fit leaves larger residuals); the second is the gain over linear
    // interpolation, whose error at the quotes is the noise itself (RMS s).
    // Averaged over 60 noise draws; the estimate's sd is ~2-3%.
    constexpr double s = 0.002;
    constexpr int kRuns = 60;
    const double T = 1.0;
    const svi::Params truth = ssvi_slice(T);
    const std::vector<double> k = strikes_k(T);
    double resid_sq = 0.0, truth_sq = 0.0;
    for (int run = 0; run < kRuns; ++run) {
        NormalRng rng(5000 + run);
        std::vector<double> w;
        for (double x : k) {
            const double vol = std::sqrt(truth.w(x) / T) + s * rng.next();
            w.push_back(vol * vol * T);
        }
        const svi::Fit f = svi::fit(k, w, T);
        resid_sq += f.rmse_vol * f.rmse_vol;
        for (double x : k) {
            const double e = std::sqrt(f.params.w(x) / T) - std::sqrt(truth.w(x) / T);
            truth_sq += e * e / k.size();
        }
        EXPECT_GT(svi::min_durrleman_g(f.params, k.front(), k.back()).g, 0.0) << "run " << run;
    }
    const double resid = std::sqrt(resid_sq / kRuns) / s;
    const double err = std::sqrt(truth_sq / kRuns) / s;
    EXPECT_NEAR(resid, std::sqrt(16.0 / 21.0), 0.08);
    EXPECT_GT(err, 0.38);
    EXPECT_LT(err, 0.62);
}

TEST(SviFit, LeeWingBoundBindsWhenTheDataAskForSteeperWings) {
    // w = 0.04 + 3 |k| has wing slope 3, above Lee's bound of 2: no
    // arbitrage-free smile can have it, so the fit must stop at the bound.
    std::vector<double> k, w;
    for (double x = -1.0; x <= 1.0 + 1e-9; x += 0.1) {
        k.push_back(x);
        w.push_back(0.04 + 3.0 * std::abs(x));
    }
    const svi::Params p = svi::fit(k, w, 1.0).params;
    EXPECT_LE(p.b * (1.0 + p.rho), 2.0 + 1e-9);
    EXPECT_LE(p.b * (1.0 - p.rho), 2.0 + 1e-9);
    EXPECT_NEAR(p.b * (1.0 + std::abs(p.rho)), 2.0, 1e-6);  // the bound is active
    EXPECT_GE(p.a, 0.0);
    EXPECT_LE(std::abs(p.rho), 1.0);
}

TEST(SviFit, DurrlemanConditionFlagsAKnownArbitrageableSlice) {
    // Axel Vogt's raw SVI parameters (Gatheral & Jacquier 2014, Example 3.1):
    // a well-behaved-looking smile whose implied density goes negative.
    const svi::Params vogt{-0.0410, 0.1331, 0.3060, 0.3586, 0.4153};
    const svi::DurrlemanMin g = svi::min_durrleman_g(vogt, -1.5, 1.5);
    EXPECT_LT(g.g, 0.0);
    EXPECT_GT(g.k, 0.0);  // the violation is on the right wing...
    EXPECT_LT(g.k, 1.0);  // ...not at the grid's edge
}

TEST(SviFit, DerivativesMatchFiniteDifferences) {
    const svi::Params p = ssvi_slice(0.5);
    const double h = 1e-5;
    for (double k : {-0.8, -0.1, 0.0, 0.11, 0.6}) {
        EXPECT_NEAR(p.dw(k), (p.w(k + h) - p.w(k - h)) / (2 * h), 1e-8);
        EXPECT_NEAR(p.d2w(k), (p.dw(k + h) - p.dw(k - h)) / (2 * h), 1e-6);
    }
}

TEST(SviFit, RejectsUnusableInput) {
    const std::vector<double> k = {-0.2, -0.1, 0.0, 0.1, 0.2};
    const std::vector<double> w = {0.05, 0.045, 0.04, 0.042, 0.046};
    EXPECT_NO_THROW(svi::fit(k, w, 1.0));
    EXPECT_THROW(svi::fit({-0.1, 0.0, 0.1, 0.2}, {0.05, 0.04, 0.04, 0.05}, 1.0),
                 std::invalid_argument);  // 4 points, 5 parameters
    EXPECT_THROW(svi::fit(k, {0.05, 0.04}, 1.0), std::invalid_argument);
    EXPECT_THROW(svi::fit(k, {0.05, 0.045, 0.0, 0.042, 0.046}, 1.0), std::invalid_argument);
    EXPECT_THROW(svi::fit(k, {0.05, 0.045, std::nan(""), 0.042, 0.046}, 1.0),
                 std::invalid_argument);
    EXPECT_THROW(svi::fit(k, w, 0.0), std::invalid_argument);
}

// ---- SVI as the surface's smile model -------------------------------------------

std::vector<surface::Quote> ssvi_quotes(double dK, double noise, std::uint64_t seed) {
    NormalRng rng(seed);
    std::vector<surface::Quote> out;
    for (double T : kExpiries) {
        const svi::Params p = ssvi_slice(T);
        for (double K = 60.0; K <= 160.0 + 1e-9; K += dK) {
            const double vol = std::sqrt(p.w(std::log(K / forward(T))) / T) + noise * rng.next();
            for (OptionType type : {OptionType::Call, OptionType::Put})
                out.push_back({type, K, T, bs::price({type, K, T}, {kSpot, kRate, vol, kDiv})});
        }
    }
    return out;
}

double true_vol(double K, double T) {
    return std::sqrt(ssvi_slice(T).w(std::log(K / forward(T))) / T);
}

TEST(VolSurfaceSvi, MatchesTheSmileBetweenAndBeyondTheQuotedStrikes) {
    // $10 strikes, no noise. Linear interpolation is off by up to ~1.5 vol
    // points between these strikes (results §9); SVI is exact because the
    // smile that generated the quotes is itself SVI. Beyond the quoted
    // strikes this is the best case for SVI, not evidence about real wings.
    const auto quotes = ssvi_quotes(10.0, 0.0, 1);
    const surface::VolSurface lin(quotes, kSpot, kRate, kDiv);
    const surface::VolSurface fit(quotes, kSpot, kRate, kDiv, surface::SmileModel::Svi);
    double worst_lin = 0.0, worst_svi = 0.0;
    for (double T : kExpiries)
        for (double K = 61.0; K < 160.0; K += 2.0) {
            worst_lin = std::max(worst_lin, std::abs(lin.vol(K, T) - true_vol(K, T)));
            worst_svi = std::max(worst_svi, std::abs(fit.vol(K, T) - true_vol(K, T)));
        }
    EXPECT_GT(worst_lin, 1e-3);
    EXPECT_LT(worst_svi, 1e-8);
    for (double K : {30.0, 250.0}) EXPECT_NEAR(fit.vol(K, 1.0), true_vol(K, 1.0), 1e-7);
    for (const surface::Slice& s : fit.slices()) EXPECT_TRUE(s.svi.has_value());
    EXPECT_TRUE(fit.arbitrage_violations().empty());
}

TEST(VolSurfaceSvi, NoisyQuotesBreakConvexityOfTheLinearSmileButNotTheFit) {
    // 20 bp of noise on each vol at $2.50 strikes: neighbouring quotes then
    // disagree by more than the smile's curvature, so call prices interpolated
    // through them are not convex — the linear surface reports butterfly
    // arbitrage that is pure noise. The SVI fit averages it out.
    int linear_butterflies = 0, svi_violations = 0;
    for (std::uint64_t seed = 1; seed <= 10; ++seed) {
        const auto quotes = ssvi_quotes(2.5, 0.002, seed);
        const surface::VolSurface lin(quotes, kSpot, kRate, kDiv);
        const surface::VolSurface fit(quotes, kSpot, kRate, kDiv, surface::SmileModel::Svi);
        for (const auto& v : lin.arbitrage_violations())
            linear_butterflies += v.kind == surface::Violation::Kind::Butterfly;
        svi_violations += static_cast<int>(fit.arbitrage_violations().size());
    }
    EXPECT_GT(linear_butterflies, 10);
    EXPECT_EQ(svi_violations, 0);
}

TEST(VolSurfaceSvi, SlicesWithTooFewQuotesStayLinear) {
    std::vector<surface::Quote> quotes = ssvi_quotes(10.0, 0.0, 1);
    for (double K : {90.0, 100.0, 110.0})  // a 3-strike expiry
        quotes.push_back({OptionType::Call, K, 3.0,
                          bs::price({OptionType::Call, K, 3.0}, {kSpot, kRate, 0.2, kDiv})});
    const surface::VolSurface s(quotes, kSpot, kRate, kDiv, surface::SmileModel::Svi);
    ASSERT_EQ(s.slices().size(), kExpiries.size() + 1);
    EXPECT_FALSE(s.slices().back().svi.has_value());
    EXPECT_NEAR(s.vol(100.0, 3.0), 0.2, 1e-10);
}

}  // namespace
}  // namespace ope
