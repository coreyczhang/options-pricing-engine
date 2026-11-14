# Options Pricing Engine

[![CI](https://github.com/coreyczhang/options-pricing-engine/actions/workflows/ci.yml/badge.svg)](https://github.com/coreyczhang/options-pricing-engine/actions/workflows/ci.yml)

A European option pricing and Greeks engine in C++17: a **Black-Scholes closed form** with
analytic Greeks and an **implied-volatility solver**, and a **Monte Carlo** pricer — plain, antithetic, and bump-and-reprice
Greeks with common random numbers — that is validated *against* the closed form. The
emphasis is on correctness you can check: every Monte Carlo number carries a standard
error, the standard error itself is tested, and the variance reduction is measured rather
than asserted.

See **[docs/results.md](docs/results.md)** for the validation write-up and
**[docs/benchmarks.md](docs/benchmarks.md)** for timings.

---

## Scope

**What this is.** A pricing and Greeks engine for European calls and puts under
Black-Scholes-Merton (constant rate, constant volatility, continuous dividend yield), with
two independent implementations that are checked against each other and against published
textbook values.

**What this isn't.**

- **No American or path-dependent options.** European payoffs only; no early exercise, no
  Asians or barriers, so no time-stepping or PDE/lattice methods.
- **No volatility surface or calibration (yet).** A single quote can be inverted to its
  implied vol (M6), but there is no smile, term structure, or surface fitting; that is M7.
- **No market data.** Inputs are numbers on the command line or in code; there is no feed,
  no day-count or calendar handling, and no discrete dividends.
- **Not optimized for throughput.** Single-threaded and scalar; see
  [benchmarks.md](docs/benchmarks.md#not-measured).

## Roadmap

| Milestone | Scope | State |
|-----------|-------|-------|
| M1 | Black-Scholes price + analytic Greeks, degenerate cases, reference-value tests | ✅ done |
| M2 | Plain Monte Carlo with standard errors; convergence to BS at 1/√N | ✅ done |
| M3 | Antithetic variates; measured variance reduction across moneyness | ✅ done |
| M4 | Bump-and-reprice MC Greeks with common random numbers vs analytic Greeks | ✅ done |
| M5 | Benchmarks, results write-up, CI | ✅ done |
| M6 | Implied-vol solver: safeguarded Newton, arbitrage-bound handling, conditioning-aware tests | ✅ done |
| M7 | Implied-vol surface from a quote grid; Greeks and surface plots | planned |

**Not done (see [Next steps](docs/results.md#next-steps)):** importance sampling for
deep-OTM options, a sharper implied-vol starting guess, pathwise / likelihood-ratio Greeks, control variates, quasi-Monte Carlo,
path-dependent payoffs.

## Build & test

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires a C++17 compiler and CMake ≥ 3.14. GoogleTest is fetched automatically via CMake
`FetchContent`. Builds Release by default, with `-Wall -Wextra -Wpedantic` clean.

## Usage

```bash
# pricer <call|put> <spot> <strike> <expiry_years> <rate> <vol> [div_yield]
./build/pricer call 49 50 0.3846 0.05 0.20

# pricer iv <call|put> <spot> <strike> <expiry_years> <rate> <price> [div_yield]
./build/pricer iv call 21 20 0.25 0.10 1.875      # Hull's example: 23.45%

# Regenerate every table in docs/ (Release build for meaningful timings)
./build/ope_bench
```

As a library:

```cpp
#include "ope/BlackScholes.hpp"
#include "ope/ImpliedVol.hpp"
#include "ope/MonteCarlo.hpp"

const ope::EuropeanOption opt{ope::OptionType::Call, /*strike=*/100.0, /*expiry=*/1.0};
const ope::Market mkt{/*spot=*/100.0, /*rate=*/0.05, /*vol=*/0.20};

double px = ope::bs::price(opt, mkt);                     // 10.450584
ope::Greeks g = ope::bs::greeks(opt, mkt);

auto mc = ope::mc::price(opt, mkt, {1'000'000, /*seed=*/42, /*antithetic=*/true});
// mc.price ± mc.std_error
auto mg = ope::mc::greeks(opt, mkt, {1'000'000});         // mg.value, mg.std_error

ope::iv::Result r = ope::iv::solve(opt, mkt, /*price=*/10.450584);  // r.vol ~ 0.20
// mkt.vol is only the starting guess; r.status flags quotes outside no-arbitrage bounds
```

## Results at a glance

- **Black-Scholes** matches Hull's and Haug's published prices and Greeks to every printed
  digit; put-call parity holds to 1e-10 and every analytic Greek matches a finite
  difference of the price on a 180-point grid.
- **Plain MC** converges at 1/√N across five orders of magnitude (RMSE × √N ≈ 14.5,
  constant), and its reported std error matches the actual spread of independent runs.
- **Antithetic variates** cut variance **2.0× ATM, 17× deep ITM, 50× for an ITM put, and
  ~1× far OTM** — exactly the pattern VRF = 1/(1 + ρ) predicts — and run ~1.7× faster, so
  the time-adjusted gain is 1.7–86×.
- **MC Greeks** (1M paths) match the analytic Greeks with all |z| ≤ 0.74; CRN keeps
  delta's noise flat as the bump shrinks, while independent seeds are 18× worse at a 1%
  bump and 1,809× worse at 0.01%.
- **Implied vol** recovers σ to ~1e-12 wherever the quote determines it, and otherwise to
  within the limit set by price rounding / vega (tested on 2,016 contracts, 1%–300% vol);
  median 5 iterations, 661 ns per solve.
- **Speed (Apple M3):** closed form 57 ns; MC 35 ns/path plain, 20 ns/path antithetic,
  72 ns/path for all five Greeks.

## Why these design choices

**Two pricers, one model.** The Monte Carlo engine prices the same GBM model as the closed
form, so any disagreement beyond sampling error is a bug. That turns Black-Scholes from a
feature into a test oracle — the MC machinery (error bars, variance reduction, Greeks) is
what would carry over to payoffs with no closed form.

**Exact terminal sampling, no time-stepping.** A European payoff depends only on S_T, which
is exactly lognormal under GBM, so each path is one normal draw. That removes
discretization bias and leaves sampling error as the *only* MC error — which is what makes
a "within 4 standard errors" test meaningful.

**Every MC result returns its standard error.** A Monte Carlo price without an error bar
can't be checked. The std error is computed over the independent sampling unit — a single
path, or for antithetic the pair average (f(Z) + f(−Z))/2, since the two halves are
correlated. Treating the halves as independent is a classic way to get the error bar wrong.

**Test the error bar, not just the price.** A pricer that under-reports its std error
passes every "within k·se" test for the wrong reason, so the suite also checks calibration:
over 400 independent runs the actual RMS error must match the reported se, and 95%
intervals must cover the truth ~95% of the time.

**Statistical tests that are deterministic but not seed-tuned.** Fixed seeds make every
test pass or fail reproducibly; ~4σ thresholds make a pass meaningful. The thresholds were
then re-checked under thousands of other seeds — which caught one that wasn't safe (see
[results §7](docs/results.md#7-how-the-statistical-tests-stay-deterministic)).

**Own normal transform instead of `std::normal_distribution`.** The standard fixes
`mt19937_64`'s output bit-for-bit but leaves the normal transform to the library, so
libstdc++ and libc++ give different draws from one seed. Box-Muller on top of
`mt19937_64` makes a seed mean the same paths on Linux CI and macOS, so the numbers in the
docs reproduce.

**Separate bump bias from sampling noise.** A finite-difference Greek has O(h²) bias *and*
MC noise. Testing MC against the analytic Greek with one loose tolerance blurs the two, so
MC is tested tightly against the *same finite difference of the BS price* (pure noise), and
that finite difference is tested against the analytic Greek (pure bias).

**Common random numbers inside the path loop.** Each bumped payoff is evaluated on the same
draw as the base payoff, and the Greek is the mean of per-path differences. This keeps
delta's noise independent of the bump size (independent seeds make it grow as 1/h) and
gives each Greek its own std error from one pass.

**Degenerate inputs return the exact limit.** At T = 0 or σ = 0 the BS formula divides by
zero; instead the engine returns the discounted forward intrinsic value, which is the
formula's limit, and tests that the general formula converges to it smoothly.

**Implied vol: bracketed Newton, statuses not exceptions, accuracy judged by conditioning.**
Newton with the analytic vega is quadratic near the root but can diverge or crawl on the
flat tails, so it runs inside a bracket and bisects whenever a step leaves it or fails to
halve the last one — convergence is guaranteed, not hoped for. Quotes outside the
no-arbitrage bounds are returned as a status because real chains contain them and a surface
builder must skip, not abort. And the round-trip test doesn't demand 1e-12 everywhere —
that would be unattainable for a deep-ITM option whose vega is 1e-10 — it demands accuracy
to within what price rounding divided by vega allows.

**Numerics where they matter.** N(x) uses `erfc`, not `1 + erf`, so far-OTM prices keep
precision instead of rounding to 0; variance uses Welford's one-pass algorithm instead of
Σx² − n·x̄², which cancels catastrophically when the spread is small relative to the level.

## Known limitations

- **Rare events break plain MC's error bar.** A put ~4.7σ out of the money pays on about
  one path in a million; at 200k paths none hit and MC reports **0 ± 0**. This is pinned
  down by a test rather than hidden. Fix: importance sampling (not implemented).
- **Bump-and-reprice gamma is noisy.** Its variance grows like 1/h, so gamma needs a
  coarser bump than delta (bias vs noise in [results §5](docs/results.md#5-bump-size-bias-vs-noise-and-why-common-random-numbers)).
  Default bumps are shared across Greeks; at 1M paths delta's 1%-bump bias is ~0.15 se,
  so for much higher path counts delta should use a smaller bump.
- **Model risk is out of scope.** Everything is validated *within* Black-Scholes; nothing
  here says the model is right for any real market.
- **Timings are from one machine** and single-threaded.

## Layout

```
include/ope/   Option.hpp (types, payoff, validation), BlackScholes.hpp,
               ImpliedVol.hpp, MonteCarlo.hpp, NormalRng.hpp
src/           Option.cpp, BlackScholes.cpp, ImpliedVol.cpp, MonteCarlo.cpp,
               main.cpp (pricer CLI), bench_main.cpp (ope_bench)
tests/         test_black_scholes.cpp, test_implied_vol.cpp, test_monte_carlo.cpp,
               test_mc_greeks.cpp
docs/          results.md, benchmarks.md
```

## License

MIT — see [LICENSE](LICENSE).
