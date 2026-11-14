# Benchmarks

Run with a Release build (`ope_bench` target):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
./build/ope_bench
```

Numbers below are from an **Apple M3** laptop (AppleClang 21, `-O3`, single thread); they
are **machine-dependent** — treat them as relative, not absolute. Each MC figure is the
best of 5 runs (3 at 10M paths) with a fresh seed per run.

## Results

**Black-Scholes closed form:** price **57 ns/option**; price + all five Greeks
**66 ns/option** (average over 2,000,000 calls with varying spot).

**Implied vol:** **661 ns/solve**, 4.6 Newton/bisection iterations on average (1y calls,
K 80–120, true vol 10%–46%, automatic start). Each iteration is one price and one vega, so
this is ~5× the cost of a price; iteration counts across a wider grid are in
[results §8](results.md#8-implied-volatility-inverting-black-scholes).

**Monte Carlo, ATM call:**

| paths | plain price | ns/path | antithetic price | ns/path | Greeks (plain) | ns/path | price se (plain) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 35 µs | 35.3 | 21 µs | 21.2 | 72 µs | 72.0 | 0.4767 |
| 10000 | 357 µs | 35.7 | 205 µs | 20.5 | 706 µs | 70.6 | 0.1472 |
| 100000 | 3.5 ms | 35.1 | 2.0 ms | 20.4 | 7.1 ms | 71.3 | 0.0465 |
| 1000000 | 35.5 ms | 35.5 | 20.5 ms | 20.5 | 71.5 ms | 71.5 | 0.0147 |
| 10000000 | 355.1 ms | 35.5 | 205.2 ms | 20.5 | 717.7 ms | 71.8 | 0.0047 |

## Reading these honestly

- **Cost per path is flat in N.** Each path is one normal draw, an `exp`, and a payoff;
  there is no per-path memory, so nothing changes as N grows. Total time is linear in N
  while the error falls as 1/√N: a penny of accuracy on this contract (se ≈ 0.01) costs
  ~2.2M paths ≈ 75 ms plain, or ~1.1M paths ≈ 22 ms antithetic.
- **Per-unit work, not the payoff, is the dominant cost.** Antithetic does the same number
  of `exp` + payoff evaluations as plain, but half as many normal draws and half as many
  accumulator updates, and runs ~1.7× faster. Solving the two rows: ≈ 30 ns per sampling
  unit (one normal from Box-Muller's `log`/`sqrt`/`sin`/`cos`, plus a Welford update whose
  divide sits on a serial dependency chain) vs ≈ 5.5 ns for `exp` + payoff. These numbers
  can't split the 30 ns between the RNG and the accumulator; a profile would. Likely next
  optimizations are a faster normal generator (ziggurat — not done because it would change
  the draws behind every published number) and a cheaper accumulator (e.g. sum and sum of
  squares around a shifted origin).
- **All five Greeks cost ~2× one price**, not the ~9× of repricing separately. One CRN
  pass evaluates 1 base `exp`, 6 bumped `exp`s (vol, expiry, rate up/down; spot bumps are
  free because S_T is linear in S) and 9 payoffs, but shares one normal draw per path.
  Hoisting the bumped scenarios' drift/√T/discount constants out of the path loop took this
  from 95 ms to 71 ms per 1M paths with bit-identical output.
- **Closed form vs Monte Carlo:** 57 ns vs ~75 ms for a penny-accurate plain MC price —
  about six orders of magnitude. That gap is the point of the comparison: for a European option
  MC is a *validation tool* against the closed form, and the machinery (error bars,
  variance reduction, CRN Greeks) is what carries over to payoffs that have no closed form.

## A benchmarking artifact that was caught

The first version of this table re-ran the same seed for each of the 5 timing repetitions.
At small N that let the CPU's branch predictor memorize the payoff's in/out-of-the-money
pattern: **22.5 ns/path at 1k paths with a fixed seed vs 34.1 with varied seeds** (35.4 at
1M paths either way). That would have reported small-N pricing as ~1.5× faster than it is
for any real, never-repeated input. The timing loop now uses a fresh seed per run, which
is why ns/path is flat down to 1k paths.

## Not measured

- **Multi-threading.** Everything is single-threaded. Paths are independent, so MC is
  embarrassingly parallel, but a parallel version needs per-thread RNG streams with a
  deterministic split to keep the "seed ⇒ exact result" guarantee the tests rely on.
- **SIMD / batched transcendental functions.** The path loop is scalar.
