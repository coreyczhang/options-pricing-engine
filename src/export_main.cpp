// ope_export: write the CSVs behind the plots in docs/img/.
//
//   ope_export <output_dir>
//   python3 scripts/plot.py <output_dir> docs/img
//
// Plotting lives in Python so the C++ library stays dependency-free; this
// program only computes. Everything is deterministic (fixed MC seed, synthetic
// quotes), so the plots regenerate bit-for-bit.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "ope/BlackScholes.hpp"
#include "ope/MonteCarlo.hpp"
#include "ope/NormalRng.hpp"
#include "ope/VolSurface.hpp"

namespace {

using ope::EuropeanOption;
using ope::Market;
using ope::OptionType;

FILE* open_csv(const std::string& dir, const char* name, const char* header) {
    const std::string path = dir + "/" + name;
    FILE* f = std::fopen(path.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "error: cannot write %s\n", path.c_str());
        std::exit(1);
    }
    std::fprintf(f, "%s\n", header);
    return f;
}

// Analytic Greeks of a K=100 call across spot, at several expiries.
void greeks_vs_spot(const std::string& dir) {
    FILE* f = open_csv(dir, "greeks_vs_spot.csv", "expiry,spot,price,delta,gamma,vega,theta,rho");
    for (double T : {1.0 / 52, 1.0 / 12, 0.25, 1.0})
        for (double S = 60.0; S <= 140.0 + 1e-9; S += 0.5) {
            const EuropeanOption opt{OptionType::Call, 100.0, T};
            const Market mkt{S, 0.05, 0.20};
            const ope::Greeks g = ope::bs::greeks(opt, mkt);
            std::fprintf(f, "%.6f,%.2f,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", T, S,
                         ope::bs::price(opt, mkt), g.delta, g.gamma, g.vega, g.theta, g.rho);
        }
    std::fclose(f);
}

// Bump-and-reprice MC Greeks (with std errors) against the analytic curve.
// One independent seed per spot, so each point is its own draw. (With one
// shared seed every point carries the same sampling error: seed 7, used by
// several tests, happens to have a first-200k-draw mean of +3.95 standard
// errors, which put every point of an earlier version of this plot ~3 se high.)
void mc_vs_analytic(const std::string& dir) {
    FILE* f = open_csv(dir, "mc_greeks.csv",
                       "spot,greek,mc,mc_se,analytic,paths");
    constexpr std::size_t kPaths = 200'000;
    const EuropeanOption opt{OptionType::Call, 100.0, 0.25};
    std::uint64_t seed = 2'000;
    for (double S = 75.0; S <= 125.0 + 1e-9; S += 5.0) {
        const Market mkt{S, 0.05, 0.20};
        ++seed;
        const ope::Greeks a = ope::bs::greeks(opt, mkt);
        const ope::mc::Result p = ope::mc::price(opt, mkt, {kPaths, seed});
        const ope::mc::GreeksResult g = ope::mc::greeks(opt, mkt, {kPaths, seed});
        auto row = [&](const char* name, double mc, double se, double exact) {
            std::fprintf(f, "%.1f,%s,%.10g,%.10g,%.10g,%zu\n", S, name, mc, se, exact, kPaths);
        };
        row("price", p.price, p.std_error, ope::bs::price(opt, mkt));
        row("delta", g.value.delta, g.std_error.delta, a.delta);
        row("gamma", g.value.gamma, g.std_error.gamma, a.gamma);
        row("vega", g.value.vega, g.std_error.vega, a.vega);
        row("theta", g.value.theta, g.std_error.theta, a.theta);
        row("rho", g.value.rho, g.std_error.rho, a.rho);
    }
    std::fclose(f);
}

// ---- Implied-vol surface from synthetic SSVI quotes ---------------------------
// Same SSVI parameters as tests/test_vol_surface.cpp: theta(T) = 0.04 T +
// 0.02 (1 - e^{-2T}), rho = -0.6, eta = 1, gamma = 1/2 (arbitrage-free).
constexpr double kSpot = 100.0, kRate = 0.03, kDiv = 0.01;

double ssvi_vol(double K, double T) {
    const double th = 0.04 * T + 0.02 * (1.0 - std::exp(-2.0 * T));
    const double rho = -0.6;
    const double phi = 1.0 / (std::sqrt(th) * std::sqrt(1.0 + th));
    const double k = std::log(K / (kSpot * std::exp((kRate - kDiv) * T)));
    const double w = th / 2.0 * (1.0 + rho * phi * k +
                                 std::sqrt((phi * k + rho) * (phi * k + rho) + 1.0 - rho * rho));
    return std::sqrt(w / T);
}

void vol_surface(const std::string& dir) {
    // A chain-like grid: $5 strikes from 60 to 160, five listed expiries,
    // calls and puts at every strike.
    const std::vector<double> expiries = {1.0 / 12, 0.25, 0.5, 1.0, 2.0};
    std::vector<ope::surface::Quote> quotes;
    for (double T : expiries)
        for (double K = 60.0; K <= 160.0 + 1e-9; K += 5.0)
            for (OptionType type : {OptionType::Call, OptionType::Put})
                quotes.push_back({type, K, T,
                                  ope::bs::price({type, K, T}, {kSpot, kRate, ssvi_vol(K, T), kDiv})});
    const ope::surface::VolSurface s(quotes, kSpot, kRate, kDiv);

    FILE* nodes = open_csv(dir, "iv_nodes.csv", "expiry,strike,vol,type,true_vol");
    for (const ope::surface::Slice& slice : s.slices())
        for (const ope::surface::Node& n : slice.nodes)
            std::fprintf(nodes, "%.6f,%.2f,%.12g,%s,%.12g\n", slice.expiry, n.strike, n.vol,
                         n.type == OptionType::Call ? "call" : "put",
                         ssvi_vol(n.strike, slice.expiry));
    std::fclose(nodes);

    FILE* grid = open_csv(dir, "iv_surface.csv", "expiry,strike,vol,true_vol");
    for (double T = 0.05; T <= 2.0 + 1e-9; T += 0.05)
        for (double K = 60.0; K <= 160.0 + 1e-9; K += 1.0)
            std::fprintf(grid, "%.4f,%.2f,%.12g,%.12g\n", T, K, s.vol(K, T), ssvi_vol(K, T));
    // Each quoted expiry on a fine strike grid, for the smile/error panels.
    for (double T : expiries)
        for (double K = 60.0; K <= 160.0 + 1e-9; K += 0.5)
            std::fprintf(grid, "%.6f,%.2f,%.12g,%.12g\n", T, K, s.vol(K, T), ssvi_vol(K, T));
    std::fclose(grid);

    std::printf("surface: %zu quotes, %zu slices, %zu skipped, %zu arbitrage violations\n",
                quotes.size(), s.slices().size(), s.skipped().size(),
                s.arbitrage_violations().size());
}

// ---- SVI fit vs linear interpolation on noisy quotes -------------------------
// The same SSVI smile with 20 bp of Gaussian noise on every vol, $2.50 strikes,
// at 1 month and 1 year (NormalRng: identical draws on every platform).
void svi_fit(const std::string& dir) {
    const std::vector<double> expiries = {1.0 / 12, 1.0};
    ope::NormalRng rng(1);
    std::vector<ope::surface::Quote> quotes;
    for (double T : expiries)
        for (double K = 60.0; K <= 160.0 + 1e-9; K += 2.5) {
            const double vol = ssvi_vol(K, T) + 0.002 * rng.next();
            for (OptionType type : {OptionType::Call, OptionType::Put})
                quotes.push_back(
                    {type, K, T, ope::bs::price({type, K, T}, {kSpot, kRate, vol, kDiv})});
        }
    const ope::surface::VolSurface lin(quotes, kSpot, kRate, kDiv);
    const ope::surface::VolSurface fit(quotes, kSpot, kRate, kDiv,
                                       ope::surface::SmileModel::Svi);

    FILE* q = open_csv(dir, "svi_quotes.csv", "expiry,strike,quote_vol");
    for (const ope::surface::Slice& s : fit.slices())
        for (const ope::surface::Node& n : s.nodes)
            std::fprintf(q, "%.6f,%.2f,%.12g\n", s.expiry, n.strike, n.vol);
    std::fclose(q);

    FILE* c = open_csv(dir, "svi_curves.csv", "expiry,strike,svi_vol,linear_vol,true_vol");
    for (double T : expiries)
        for (double K = 40.0; K <= 220.0 + 1e-9; K += 0.5)
            std::fprintf(c, "%.6f,%.2f,%.12g,%.12g,%.12g\n", T, K, fit.vol(K, T), lin.vol(K, T),
                         ssvi_vol(K, T));
    std::fclose(c);

    int lin_bf = 0;
    for (const auto& v : lin.arbitrage_violations())
        lin_bf += v.kind == ope::surface::Violation::Kind::Butterfly;
    for (const ope::surface::Slice& s : fit.slices())
        std::printf("svi: T=%.4f rmse vs quotes %.1f bp (max %.1f bp)\n", s.expiry,
                    s.svi->rmse_vol * 1e4, s.svi->max_abs_vol * 1e4);
    std::printf("svi: butterfly violations: linear %d, svi %zu\n", lin_bf,
                fit.arbitrage_violations().size());
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: ope_export <output_dir>\n");
        return 2;
    }
    const std::string dir = argv[1];
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);  // failures surface in open_csv
    greeks_vs_spot(dir);
    mc_vs_analytic(dir);
    vol_surface(dir);
    svi_fit(dir);
    std::printf("wrote greeks_vs_spot.csv, mc_greeks.csv, iv_nodes.csv, iv_surface.csv, "
                "svi_quotes.csv, svi_curves.csv to %s\n",
                dir.c_str());
    return EXIT_SUCCESS;
}
