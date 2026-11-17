#include "ope/Ssvi.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "ope/BlackScholes.hpp"
#include "ope/NormalRng.hpp"
#include "ope/VolSurface.hpp"

namespace ope {
namespace {

// The SSVI surface of test_vol_surface.cpp and test_svi.cpp:
// theta(T) = 0.04 T + 0.02 (1 - e^{-2T}), rho = -0.6, eta = 1, gamma = 1/2.
constexpr double kSpot = 100.0, kRate = 0.03, kDiv = 0.01;
const ssvi::Params kTruth{-0.6, 1.0, 0.5};

double theta(double T) { return 0.04 * T + 0.02 * (1.0 - std::exp(-2.0 * T)); }
double forward(double T) { return kSpot * std::exp((kRate - kDiv) * T); }
double true_vol(double K, double T) {
    return std::sqrt(kTruth.w(std::log(K / forward(T)), theta(T)) / T);
}

const std::vector<double> kExpiries = {1.0 / 12, 0.25, 0.5, 1.0, 2.0};

// $5 strikes 60..160 (21 points) with total variance from `vol`.
template <class Vol>
ssvi::SliceData slice(double T, const Vol& vol) {
    ssvi::SliceData d{T, {}, {}};
    for (double K = 60.0; K <= 160.0 + 1e-9; K += 5.0) {
        d.k.push_back(std::log(K / forward(T)));
        const double v = vol(K, T);
        d.w.push_back(v * v * T);
    }
    return d;
}

std::vector<surface::Quote> quotes_from(const std::vector<double>& Ts, double dK,
                                        double noise, std::uint64_t seed) {
    NormalRng rng(seed);
    std::vector<surface::Quote> q;
    for (double T : Ts)
        for (double K = 60.0; K <= 160.0 + 1e-9; K += dK) {
            const double v = true_vol(K, T) + noise * rng.next();
            for (OptionType t : {OptionType::Call, OptionType::Put})
                q.push_back({t, K, T, bs::price({t, K, T}, {kSpot, kRate, v, kDiv})});
        }
    return q;
}

int count(const std::vector<surface::Violation>& v, surface::Violation::Kind k) {
    return static_cast<int>(std::count_if(v.begin(), v.end(), [&](const auto& x) { return x.kind == k; }));
}


template <class Vol>
std::vector<ssvi::SliceData> slices_of(const Vol& vol) {
    std::vector<ssvi::SliceData> out;
    for (double T : kExpiries) out.push_back(slice(T, vol));
    return out;
}

auto vol_of(const ssvi::Params& p) {
    return [p](double K, double T) { return std::sqrt(p.w(std::log(K / forward(T)), theta(T)) / T); };
}

// Smallest Durrleman g over the fitted slices and the expiries between them.
double min_g(const ssvi::Fit& f) {
    double g = 1e300;
    for (double T = 0.02; T <= 2.5; T += 0.01)
        g = std::min(g, svi::min_durrleman_g(f.params.raw(f.theta(T)), -1.5, 1.5).g);
    return g;
}

// ---- The model ---------------------------------------------------------------

TEST(SsviParams, RawSviFormIsTheSameSmile) {
    for (double th : {0.005, 0.05, 0.5})
        for (const ssvi::Params& p :
             {kTruth, ssvi::Params{0.3, 1.2, 0.2}, ssvi::Params{-0.95, 0.9, 0.0}})
            for (double k = -2.0; k <= 2.0; k += 0.05)
                EXPECT_NEAR(p.raw(th).w(k), p.w(k, th), 1e-14 * p.w(k, th));
    const ssvi::Params flat{-0.6, 0.0, 0.5};  // eta = 0: w = theta everywhere
    EXPECT_DOUBLE_EQ(flat.raw(0.04).w(0.7), 0.04);
}

// ---- Fitting -----------------------------------------------------------------

TEST(SsviFit, RecoversTheGeneratingSurface) {
    const ssvi::Fit f = ssvi::fit(slices_of(true_vol));
    EXPECT_NEAR(f.params.rho, kTruth.rho, 1e-9);
    EXPECT_NEAR(f.params.eta, kTruth.eta, 1e-9);
    EXPECT_NEAR(f.params.gamma, kTruth.gamma, 1e-9);  // on its bound, 1/2
    ASSERT_EQ(f.thetas.size(), kExpiries.size());
    for (std::size_t i = 0; i < kExpiries.size(); ++i) {
        EXPECT_DOUBLE_EQ(f.expiries[i], kExpiries[i]);
        EXPECT_NEAR(f.thetas[i], theta(kExpiries[i]), 1e-9 * theta(kExpiries[i]));
    }
    EXPECT_LT(f.max_abs_vol, 1e-12);  // measured ~2e-14
    EXPECT_EQ(f.points, 105);
    EXPECT_GT(min_g(f), 0.0);
}

TEST(SsviFit, OnNoisyQuotesResidualAndErrorMatchLeastSquaresTheory) {
    // As the SVI version (test_svi.cpp), but jointly: p = 3 globals + 5 thetas
    // = 8 parameters for n = 105 points, so
    //   residual vs quotes:  RMS ~ s sqrt((n - p) / n) = 0.961 s
    //   error vs the truth:  RMS ~ s sqrt(p / n)       = 0.276 s
    // against 0.488 s for per-slice SVI. The error check also pins down the
    // weighting: weights 1/w instead of 1/(w T) count the 2y slice 24x the 1m
    // one, and measure 0.325 s here (correct weights: residual 0.954 s, error
    // 0.283 s). 60 noise draws; the error estimate's sd is ~3%.
    constexpr double s = 0.002;
    constexpr int kRuns = 60;
    double resid_sq = 0.0, err_sq = 0.0;
    int n = 0;
    for (int run = 0; run < kRuns; ++run) {
        NormalRng rng(7000 + run);
        const ssvi::Fit f = ssvi::fit(
            slices_of([&](double K, double T) { return true_vol(K, T) + s * rng.next(); }));
        resid_sq += f.rmse_vol * f.rmse_vol;
        for (double T : kExpiries)
            for (double k : slice(T, true_vol).k) {
                const double e = std::sqrt(f.w(k, T) / T) - std::sqrt(kTruth.w(k, theta(T)) / T);
                err_sq += e * e;
                ++n;
            }
    }
    const double resid = std::sqrt(resid_sq / kRuns) / s;
    const double err = std::sqrt(err_sq / n) / s;
    EXPECT_NEAR(resid, std::sqrt(97.0 / 105.0), 0.03);
    EXPECT_GT(err, 0.24);
    EXPECT_LT(err, 0.305);
}

TEST(SsviFit, StaysArbitrageFreeWhenTheDataAreNot) {
    // eta = 4, rho = -0.6: eta (1 + |rho|) = 6.4, and these slices really do
    // have butterfly arbitrage (negative Durrleman g).
    const ssvi::Params bad{-0.6, 4.0, 0.5};
    double truth_g = 1e300;
    for (double T : kExpiries)
        truth_g = std::min(truth_g, svi::min_durrleman_g(bad.raw(theta(T)), -1.5, 1.5).g);
    ASSERT_LT(truth_g, -0.1);
    const ssvi::Fit f = ssvi::fit(slices_of(vol_of(bad)));
    EXPECT_LE(f.params.eta * (1.0 + std::abs(f.params.rho)), 2.0 + 1e-12);
    EXPECT_GE(min_g(f), 0.0);

    // gamma = 0.9 is outside (0, 1/2]: the fit stops at the bound.
    const ssvi::Fit g = ssvi::fit(slices_of(vol_of({-0.6, 0.5, 0.9})));
    EXPECT_LE(g.params.gamma, 0.5);
    EXPECT_GT(g.params.gamma, 0.5 - 1e-9);
    EXPECT_GE(min_g(g), 0.0);
}

TEST(SsviFit, TheEtaBoundIsSufficientNotNecessary) {
    // eta = 1.8 breaks eta (1 + |rho|) <= 2, yet this surface has no
    // butterfly arbitrage (Durrleman g > 0 everywhere; arbitrage starts
    // between eta = 2 and 3). The fit still has to respect the bound, and pays
    // for it in fit error: the price of a sufficient condition.
    const ssvi::Params steep{-0.6, 1.8, 0.5};
    for (double T : kExpiries)
        EXPECT_GT(svi::min_durrleman_g(steep.raw(theta(T)), -1.5, 1.5).g, 0.2);
    const ssvi::Fit f = ssvi::fit(slices_of(vol_of(steep)));
    EXPECT_NEAR(f.params.eta * (1.0 + std::abs(f.params.rho)), 2.0, 1e-9);
    EXPECT_GT(f.rmse_vol, 0.01);  // measured 2.3 vol points
}

TEST(SsviFit, PoolsSlicesWhoseQuotesHaveCalendarArbitrage) {
    // The 0.6y quotes have *less* total variance than the 0.5y ones (3% lower
    // vol than the 0.5y smile rescaled), which no arbitrage-free surface can
    // match. The fit's theta(0.5) <= theta(0.6) constraint binds and PAV
    // pools the two slices into one theta.
    const std::vector<ssvi::SliceData> d = {
        slice(0.5, true_vol),
        slice(0.6, [](double K, double) { return 0.97 * true_vol(K, 0.5) * std::sqrt(0.5 / 0.6); })};
    const ssvi::Fit f = ssvi::fit(d);
    EXPECT_EQ(f.thetas[0], f.thetas[1]);
    EXPECT_GT(f.rmse_vol, 1e-3);  // the arbitrage cannot be fitted away
}

TEST(SsviFit, CannotFitASkewThatChangesWithExpiry) {
    // One rho for every expiry is SSVI's main restriction. Data with
    // rho(T) = -0.8 + 0.25 T (each slice still exactly raw SVI): per-slice
    // SVI fits it exactly, SSVI misses by ~2 vol points RMS.
    auto vol = [](double K, double T) {
        const ssvi::Params p{-0.8 + 0.25 * T, 1.0, 0.5};
        return std::sqrt(p.w(std::log(K / forward(T)), theta(T)) / T);
    };
    const std::vector<ssvi::SliceData> d = slices_of(vol);
    for (const ssvi::SliceData& s : d) EXPECT_LT(svi::fit(s.k, s.w, s.expiry).max_abs_vol, 1e-7);
    const ssvi::Fit f = ssvi::fit(d);
    EXPECT_GT(f.rmse_vol, 0.01);  // measured 202 bp
}

TEST(SsviFit, ThetaInterpolatesLinearlyAndExtrapolatesAtConstantAtmVol) {
    const ssvi::Fit f = ssvi::fit(slices_of(true_vol));
    const double t0 = f.thetas[0], t1 = f.thetas[1];
    EXPECT_NEAR(f.theta(kExpiries[0] / 2), t0 / 2, 1e-15);
    EXPECT_NEAR(f.theta((kExpiries[0] + kExpiries[1]) / 2), (t0 + t1) / 2, 1e-15);
    EXPECT_NEAR(f.theta(4.0), f.thetas.back() * 2.0, 1e-15);
    for (double T = 0.01; T < 5.0; T += 0.01) EXPECT_LE(f.theta(T), f.theta(T + 0.01));
}

TEST(SsviFit, RejectsBadInput) {
    using V = std::vector<double>;
    const ssvi::SliceData ok{0.5, {-0.1, 0.0, 0.1}, {0.02, 0.018, 0.017}};
    EXPECT_NO_THROW(ssvi::fit({ok}));
    EXPECT_THROW(ssvi::fit({}), std::invalid_argument);
    EXPECT_THROW(ssvi::fit({{0.5, {0.0, 0.1}, {0.02, 0.02}}}), std::invalid_argument);
    EXPECT_THROW(ssvi::fit({{0.5, V{-0.1, 0.0, 0.1}, V{0.02, 0.02}}}), std::invalid_argument);
    EXPECT_THROW(ssvi::fit({{0.5, V{-0.1, 0.0, 0.1}, V{0.02, 0.0, 0.02}}}), std::invalid_argument);
    EXPECT_THROW(ssvi::fit({{0.5, V{-0.1, NAN, 0.1}, V{0.02, 0.02, 0.02}}}), std::invalid_argument);
    EXPECT_THROW(ssvi::fit({{0.0, ok.k, ok.w}}), std::invalid_argument);
    EXPECT_THROW(ssvi::fit({ok, ok}), std::invalid_argument);  // repeated expiry
}

// ---- On the surface -------------------------------------------------------------

TEST(VolSurfaceSsvi, IsArbitrageFreeWherePerSliceSviFitsCross) {
    // Expiries 0.01y (~3.7 days) apart with 20 bp noise: the per-slice SVI
    // smiles cross each other (calendar arbitrage) in dozens of places; the
    // joint fit cannot produce any.
    for (std::uint64_t seed : {1, 2, 3}) {
        const auto q = quotes_from({0.5, 0.51, 0.52}, 5.0, 0.002, seed);
        const surface::VolSurface svi(q, kSpot, kRate, kDiv, surface::SmileModel::Svi);
        const surface::VolSurface ssvi(q, kSpot, kRate, kDiv, surface::SmileModel::Ssvi);
        EXPECT_GT(count(svi.arbitrage_violations(), surface::Violation::Kind::Calendar), 10)
            << "seed " << seed;  // measured 30, 13, 88
        EXPECT_TRUE(ssvi.arbitrage_violations().empty()) << "seed " << seed;
    }
}

TEST(VolSurfaceSsvi, MatchesTheTruthAtBetweenAndBeforeTheQuotedExpiries) {
    const auto q = quotes_from(kExpiries, 5.0, 0.0, 1);
    const surface::VolSurface lin(q, kSpot, kRate, kDiv);
    const surface::VolSurface ss(q, kSpot, kRate, kDiv, surface::SmileModel::Ssvi);
    ASSERT_TRUE(ss.ssvi().has_value());
    for (double T : kExpiries)
        for (double K = 40.0; K <= 220.0; K += 2.5)  // including beyond the quoted strikes
            EXPECT_NEAR(ss.vol(K, T), true_vol(K, T), 1e-8);
    // Between and before expiries the only error is theta linear in T (the
    // truth's theta is concave). Linear interpolation also errs in the smile.
    double worst_between_lin = 0.0, worst_between_ss = 0.0, worst_before_lin = 0.0,
           worst_before_ss = 0.0;
    for (double T = 0.02; T <= 2.0 + 1e-9; T += 0.01)
        for (double K = 60.0; K <= 160.0; K += 1.0) {
            const double el = std::abs(lin.vol(K, T) - true_vol(K, T));
            const double es = std::abs(ss.vol(K, T) - true_vol(K, T));
            double& wl = T < kExpiries[0] ? worst_before_lin : worst_between_lin;
            double& ws = T < kExpiries[0] ? worst_before_ss : worst_between_ss;
            wl = std::max(wl, el);
            ws = std::max(ws, es);
        }
    EXPECT_LT(worst_between_ss, 0.0035);  // measured 31 bp (linear: 113 bp)
    EXPECT_GT(worst_between_lin, 0.01);
    EXPECT_LT(worst_before_ss, 0.008);  // measured 76 bp (linear: 2625 bp)
    EXPECT_GT(worst_before_lin, 0.2);
    EXPECT_TRUE(ss.arbitrage_violations().empty());
}

TEST(VolSurfaceSsvi, SlicesWithTooFewQuotesTakeTheModelValue) {
    std::vector<surface::Quote> q = quotes_from(kExpiries, 5.0, 0.0, 1);
    for (double K : {95.0, 105.0})  // a 0.75y expiry with only two strikes
        q.push_back({OptionType::Call, K, 0.75,
                     bs::price({OptionType::Call, K, 0.75}, {kSpot, kRate, 0.5, kDiv})});
    const surface::VolSurface ss(q, kSpot, kRate, kDiv, surface::SmileModel::Ssvi);
    EXPECT_EQ(ss.slices().size(), 6u);
    EXPECT_EQ(ss.ssvi()->expiries.size(), 5u);
    EXPECT_NEAR(ss.vol(100.0, 0.75), std::sqrt(ss.ssvi()->w(std::log(100.0 / forward(0.75)), 0.75) / 0.75), 1e-15);
    EXPECT_LT(ss.vol(100.0, 0.75), 0.3);  // not the 50% vol of its quotes
}

}  // namespace
}  // namespace ope
