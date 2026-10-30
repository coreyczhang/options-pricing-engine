// ope_bench: validation & benchmark tables for docs/results.md.
//
// Prints GitHub-markdown tables so the numbers in the docs are regenerated,
// not hand-copied. Build Release for meaningful timings.
//
//   ./build/ope_bench

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <initializer_list>

#include "ope/BlackScholes.hpp"
#include "ope/MonteCarlo.hpp"

namespace {

using ope::EuropeanOption;
using ope::Market;
using ope::OptionType;

// Reference contract for every table: ATM 1y call, 20% vol, 5% rate.
const EuropeanOption kOpt{OptionType::Call, 100.0, 1.0};
const Market kMkt{100.0, 0.05, 0.20};
constexpr std::uint64_t kSeed = 42;

// RMS error vs Black-Scholes over `runs` independent seeds. One run's error
// is a single noisy draw — and the single-seed column reuses one seed for
// every N, so smaller runs are prefixes of larger ones and their errors are
// correlated row to row. The RMSE over independent seeds shows the 1/sqrt(N)
// law cleanly.
double rmse_over_seeds(std::size_t paths, int runs, double exact, bool antithetic = false) {
    double sum_sq = 0.0;
    for (int k = 0; k < runs; ++k) {
        const auto r =
            ope::mc::price(kOpt, kMkt, {paths, 1'000'000u + std::uint64_t(k), antithetic});
        sum_sq += (r.price - exact) * (r.price - exact);
    }
    return std::sqrt(sum_sq / runs);
}

void convergence_table() {
    const double exact = ope::bs::price(kOpt, kMkt);
    std::printf("## Plain MC convergence to Black-Scholes\n\n");
    std::printf("ATM call: S=100, K=100, T=1, r=5%%, sigma=20%%, q=0. "
                "Black-Scholes price = %.6f.\n\n", exact);
    std::printf("| paths | MC price (seed %llu) | std error | MC - BS | (MC - BS) / se "
                "| RMSE over runs | runs | RMSE x sqrt(N) |\n",
                static_cast<unsigned long long>(kSeed));
    std::printf("|------:|------:|------:|------:|------:|------:|------:|------:|\n");

    for (std::size_t n : {100u, 1'000u, 10'000u, 100'000u, 1'000'000u, 10'000'000u}) {
        const auto r = ope::mc::price(kOpt, kMkt, {n, kSeed});
        const double err = r.price - exact;
        // Fewer repeat runs at large N to keep the bench under a few seconds.
        const int runs = n <= 100'000 ? 200 : (n <= 1'000'000 ? 50 : 0);
        if (runs > 0) {
            const double rmse = rmse_over_seeds(n, runs, exact);
            std::printf("| %zu | %.6f | %.6f | %+.6f | %+.2f | %.6f | %d | %.3f |\n", n,
                        r.price, r.std_error, err, err / r.std_error, rmse, runs,
                        rmse * std::sqrt(double(n)));
        } else {
            std::printf("| %zu | %.6f | %.6f | %+.6f | %+.2f | - | - | - |\n", n, r.price,
                        r.std_error, err, err / r.std_error);
        }
    }
    std::printf("\nRMSE x sqrt(N) is constant at the payoff's discounted standard "
                "deviation iff error ~ 1/sqrt(N).\n\n");
}

// Best-of-`reps` wall time of one pricing call, in milliseconds. Min rather
// than mean: noise (scheduling, frequency scaling) only ever adds time.
template <typename F>
double best_ms(F&& f, int reps = 5) {
    double best = 1e300;
    for (int i = 0; i < reps; ++i) {
        const auto t0 = std::chrono::steady_clock::now();
        f();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

void variance_reduction_table() {
    constexpr std::size_t kPaths = 1'000'000;
    std::printf("## Antithetic variates: measured variance reduction\n\n");
    std::printf("S=100, T=1, r=5%%, sigma=20%%, q=0; %zu payoff evaluations per estimate "
                "(antithetic = %zu pairs), seed %llu.\n",
                kPaths, kPaths / 2, static_cast<unsigned long long>(kSeed));
    std::printf("VRF = (se_plain / se_anti)^2 at equal payoff evaluations; "
                "implied rho = corr(f(Z), f(-Z)) = 1/VRF - 1. "
                "Efficiency = VRF x t_plain / t_anti (variance reduction per unit "
                "of CPU time).\n\n");
    std::printf("| contract | BS | plain (± se) | antithetic (± se) | VRF | implied rho "
                "| t plain (ms) | t anti (ms) | efficiency |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    struct Row {
        const char* label;
        EuropeanOption opt;
    };
    const Row rows[] = {
        {"call K=70 (deep ITM)", {OptionType::Call, 70.0, 1.0}},
        {"call K=85", {OptionType::Call, 85.0, 1.0}},
        {"call K=100 (ATM)", {OptionType::Call, 100.0, 1.0}},
        {"call K=115", {OptionType::Call, 115.0, 1.0}},
        {"call K=130", {OptionType::Call, 130.0, 1.0}},
        {"call K=160 (far OTM)", {OptionType::Call, 160.0, 1.0}},
        {"put K=100 (ATM)", {OptionType::Put, 100.0, 1.0}},
        {"put K=130 (ITM)", {OptionType::Put, 130.0, 1.0}},
    };
    for (const Row& row : rows) {
        const double exact = ope::bs::price(row.opt, kMkt);
        ope::mc::Result plain{}, anti{};
        const double t_plain =
            best_ms([&] { plain = ope::mc::price(row.opt, kMkt, {kPaths, kSeed, false}); });
        const double t_anti =
            best_ms([&] { anti = ope::mc::price(row.opt, kMkt, {kPaths, kSeed, true}); });
        const double vrf = (plain.std_error / anti.std_error) *
                           (plain.std_error / anti.std_error);
        std::printf("| %s | %.4f | %.4f ± %.4f | %.4f ± %.4f | %.2f | %+.3f | %.1f | %.1f "
                    "| %.2f |\n",
                    row.label, exact, plain.price, plain.std_error, anti.price,
                    anti.std_error, vrf, 1.0 / vrf - 1.0, t_plain, t_anti,
                    vrf * t_plain / t_anti);
    }

    // Cross-check the se-based VRF against the actual spread of independent
    // estimates: if the antithetic std_error were computed wrongly (e.g. over
    // the 2N correlated halves), these two would disagree. Each MSE is known
    // to ~sqrt(2/runs), so runs must be large for the ratio to be sharp: at
    // 400 runs the ratio is only good to ~10%.
    constexpr std::size_t kSmall = 10'000;
    constexpr int kRuns = 2'000;
    const double exact = ope::bs::price(kOpt, kMkt);
    const double rmse_plain = rmse_over_seeds(kSmall, kRuns, exact, false);
    const double rmse_anti = rmse_over_seeds(kSmall, kRuns, exact, true);
    const auto one_plain = ope::mc::price(kOpt, kMkt, {kSmall, kSeed, false});
    const auto one_anti = ope::mc::price(kOpt, kMkt, {kSmall, kSeed, true});
    std::printf("\nCross-check, ATM call, %zu paths: VRF from reported se = %.2f; "
                "VRF from actual RMSE over %d independent seeds = %.2f.\n\n",
                kSmall,
                (one_plain.std_error / one_anti.std_error) *
                    (one_plain.std_error / one_anti.std_error),
                kRuns, (rmse_plain / rmse_anti) * (rmse_plain / rmse_anti));
}

// The finite differences mc::greeks computes, applied to the exact BS price
// (central differences; the reference contract is far from the domain edges).
ope::Greeks bs_fd(const EuropeanOption& opt, const Market& mkt, const ope::mc::Bumps& b) {
    auto px = [&](double S, double T, double r, double vol) {
        return ope::bs::price(EuropeanOption{opt.type, opt.strike, T},
                              Market{S, r, vol, mkt.div_yield});
    };
    const double S = mkt.spot, T = opt.expiry, r = mkt.rate, v = mkt.vol;
    const double h = b.spot_rel * S;
    const double up = px(S + h, T, r, v), mid = px(S, T, r, v), dn = px(S - h, T, r, v);
    return {(up - dn) / (2 * h), (up - 2 * mid + dn) / (h * h),
            (px(S, T, r, v + b.vol) - px(S, T, r, v - b.vol)) / (2 * b.vol),
            -(px(S, T + b.time, r, v) - px(S, T - b.time, r, v)) / (2 * b.time),
            (px(S, T, r + b.rate, v) - px(S, T, r - b.rate, v)) / (2 * b.rate)};
}

void greeks_table() {
    constexpr std::size_t kPaths = 1'000'000;
    const ope::mc::Bumps bumps{};
    std::printf("## Monte Carlo Greeks (bump-and-reprice, common random numbers)\n\n");
    std::printf("%zu paths, seed %llu, default bumps: spot 1%% of S, vol 0.001, expiry 1 "
                "day, rate 1 bp. Units as in ope::bs::greeks (vega per 1.00 vol, theta per "
                "year, rho per 1.00 rate).\n",
                kPaths, static_cast<unsigned long long>(kSeed));
    std::printf("\"bump bias\" = BS finite difference with the same bumps minus the analytic "
                "Greek: deterministic, no MC involved. z = (MC - BS finite difference) / se: "
                "pure sampling noise, should look like N(0,1).\n\n");

    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const EuropeanOption opt{type, kOpt.strike, kOpt.expiry};
        const ope::Greeks exact = ope::bs::greeks(opt, kMkt);
        const ope::Greeks fd = bs_fd(opt, kMkt, bumps);
        ope::mc::GreeksResult mc{};
        const double ms = best_ms([&] {
            mc = ope::mc::greeks(opt, kMkt, {kPaths, kSeed}, bumps);
        });

        std::printf("**ATM %s** (all five Greeks in one pass: %.0f ms)\n\n",
                    type == OptionType::Call ? "call" : "put", ms);
        std::printf("| Greek | analytic | bump bias | MC (± se) | rel. se | z |\n");
        std::printf("|---|---:|---:|---:|---:|---:|\n");
        auto row = [&](const char* name, double a, double f, double m, double se) {
            std::printf("| %s | %.5f | %+.1e | %.5f ± %.5f | %.2f%% | %+.2f |\n", name, a,
                        f - a, m, se, 100.0 * se / std::abs(a), (m - f) / se);
        };
        row("delta", exact.delta, fd.delta, mc.value.delta, mc.std_error.delta);
        row("gamma", exact.gamma, fd.gamma, mc.value.gamma, mc.std_error.gamma);
        row("vega", exact.vega, fd.vega, mc.value.vega, mc.std_error.vega);
        row("theta", exact.theta, fd.theta, mc.value.theta, mc.std_error.theta);
        row("rho", exact.rho, fd.rho, mc.value.rho, mc.std_error.rho);
        std::printf("\n");
    }
}

void bump_size_table() {
    constexpr std::size_t kPaths = 200'000;
    const ope::Greeks exact = ope::bs::greeks(kOpt, kMkt);
    std::printf("## Bump size: bias vs noise (ATM call, %zu paths)\n\n", kPaths);
    std::printf("Independent = up/down repriced with different seeds, se = "
                "sqrt(se_up^2 + se_dn^2) / 2h. Gamma bias = BS second difference minus "
                "analytic gamma (%.6f).\n\n", exact.gamma);
    std::printf("| spot bump h | delta se (CRN) | delta se (independent) | ratio "
                "| gamma (CRN, ± se) | gamma rel. se | gamma bump bias |\n");
    std::printf("|---:|---:|---:|---:|---:|---:|---:|\n");
    for (double rel : {1e-1, 3e-2, 1e-2, 1e-3, 1e-4}) {
        const ope::mc::Bumps bumps{rel};
        const auto crn = ope::mc::greeks(kOpt, kMkt, {kPaths, kSeed}, bumps);
        const double h = rel * kMkt.spot;
        Market up = kMkt, dn = kMkt;
        up.spot += h;
        dn.spot -= h;
        const double se_up = ope::mc::price(kOpt, up, {kPaths, 1001}).std_error;
        const double se_dn = ope::mc::price(kOpt, dn, {kPaths, 2002}).std_error;
        const double ind_se = std::sqrt(se_up * se_up + se_dn * se_dn) / (2 * h);
        const ope::Greeks fd = bs_fd(kOpt, kMkt, bumps);
        std::printf("| %g%% | %.5f | %.5f | %.0fx | %.5f ± %.5f | %.1f%% | %+.1e |\n",
                    100 * rel, crn.std_error.delta, ind_se, ind_se / crn.std_error.delta,
                    crn.value.gamma, crn.std_error.gamma,
                    100.0 * crn.std_error.gamma / exact.gamma, fd.gamma - exact.gamma);
    }
    std::printf("\n");
}

// Keeps benchmarked results observable so the optimizer cannot drop the work.
volatile double g_sink = 0.0;

void timing_table() {
    std::printf("## Timing per option\n\n");

    // Closed form: average over many calls with varying spot, so neither
    // constant folding nor a single cached input flatters the number.
    constexpr int kCalls = 2'000'000;
    double acc = 0.0;
    const double t_px = best_ms([&] {
        for (int i = 0; i < kCalls; ++i) {
            const Market m{80.0 + 40.0 * (i % 1000) / 1000.0, kMkt.rate, kMkt.vol};
            acc += ope::bs::price(kOpt, m);
        }
    });
    const double t_gk = best_ms([&] {
        for (int i = 0; i < kCalls; ++i) {
            const Market m{80.0 + 40.0 * (i % 1000) / 1000.0, kMkt.rate, kMkt.vol};
            acc += ope::bs::greeks(kOpt, m).delta;
        }
    });
    g_sink = acc;
    std::printf("Black-Scholes closed form: price %.0f ns/option, all five Greeks %.0f "
                "ns/option (average over %d calls).\n\n",
                t_px * 1e6 / kCalls, t_gk * 1e6 / kCalls, kCalls);

    // A fresh seed for every timed run. Re-running one seed lets the branch
    // predictor memorize the payoff's in/out-of-the-money pattern at small N:
    // measured 22.5 ns/path at 1k paths with a fixed seed vs 34.1 with varied
    // seeds (and 35.4 at 1M either way) — a ~1.5x flattering artifact.
    std::uint64_t timing_seed = 10'000;
    std::printf("Monte Carlo, ATM call, best of 5 runs (3 at 10M paths), fresh seed per run. "
                "\"Greeks\" = all five bump-and-reprice Greeks in one CRN pass.\n\n");
    std::printf("| paths | plain price | ns/path | antithetic price | ns/path "
                "| Greeks (plain) | ns/path | price se (plain) |\n");
    std::printf("|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for (std::size_t n : {1'000u, 10'000u, 100'000u, 1'000'000u, 10'000'000u}) {
        const int reps = n >= 10'000'000 ? 3 : 5;
        ope::mc::Result plain{};
        const double t_plain = best_ms(
            [&] { plain = ope::mc::price(kOpt, kMkt, {n, ++timing_seed, false}); }, reps);
        const double t_anti = best_ms(
            [&] { g_sink = ope::mc::price(kOpt, kMkt, {n, ++timing_seed, true}).price; },
            reps);
        const double t_greeks = best_ms(
            [&] { g_sink = ope::mc::greeks(kOpt, kMkt, {n, ++timing_seed}).value.delta; },
            reps);
        auto fmt_ms = [](double ms) {
            static char buf[32];
            if (ms < 1.0)
                std::snprintf(buf, sizeof buf, "%.0f µs", ms * 1e3);
            else
                std::snprintf(buf, sizeof buf, "%.1f ms", ms);
            return buf;
        };
        const double nd = static_cast<double>(n);
        std::printf("| %zu | %s | %.1f ", n, fmt_ms(t_plain), t_plain * 1e6 / nd);
        std::printf("| %s | %.1f ", fmt_ms(t_anti), t_anti * 1e6 / nd);
        std::printf("| %s | %.1f | %.4f |\n", fmt_ms(t_greeks), t_greeks * 1e6 / nd,
                    plain.std_error);
    }
    std::printf("\n");
}

}  // namespace

int main() {
    convergence_table();
    variance_reduction_table();
    greeks_table();
    bump_size_table();
    timing_table();
    return 0;
}
