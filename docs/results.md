# Validation Results

Every table below is printed by `ope_bench` (except §1, which is from `pricer` and the
test suite) and pasted unedited; rerun it to regenerate them:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
./build/ope_bench
```

Unless stated otherwise the reference contract is an **ATM European call: S = 100, K = 100,
T = 1y, r = 5%, σ = 20%, q = 0**, Black-Scholes price **10.450584**. Monte Carlo runs use
the exact lognormal terminal draw (no time-stepping), so the only MC error is sampling
error. Seeds are fixed; the RNG is `mt19937_64` + our own Box-Muller, so a seed means the
same paths on Linux and macOS.

## 1. Black-Scholes against published values

| Source | Contract | Published | This engine |
|---|---|---:|---:|
| Hull, BSM chapter example | call S=42 K=40 r=10% σ=20% T=0.5 | 4.76 | 4.759422 |
| | put, same | 0.81 | 0.808599 |
| Haug, generalized BSM example | call S=60 K=65 r=b=8% σ=30% T=0.25 | 2.1334 | 2.133368 |
| Haug, Merton (1973) dividend example | put S=100 K=95 r=10% q=5% σ=20% T=0.5 | 2.4648 | 2.464788 |
| Hull, "Greek letters" example (20-week call, S=49 K=50 r=5% σ=20%) | price | ≈2.40 | 2.400527 |
| | delta | 0.522 | 0.521605 |
| | gamma | 0.066 | 0.065544 |
| | vega (per 1.00 vol) | 12.1 | 12.105480 |
| | theta (per year) | −4.31 | −4.305330 |
| | rho (per 1.00 rate) | 8.91 | 8.906962 |

Every value agrees to the precision the source prints. Beyond these spot checks the test
suite verifies, on a 180-point grid (ITM/ATM/OTM, T from 0.05 to 2y, σ from 5% to 80%,
with and without rates and dividends): put-call parity to 1e-10, static no-arbitrage
bounds, Greek parity relations, and **every analytic Greek against a central finite
difference of the analytic price** — an independent check that catches a wrong sign or a
missing term in a Greek formula.

## 2. Plain Monte Carlo converges to Black-Scholes

| paths | MC price (seed 42) | std error | MC - BS | (MC - BS) / se | RMSE over runs | runs | RMSE x sqrt(N) |
|------:|------:|------:|------:|------:|------:|------:|------:|
| 100 | 9.100140 | 1.181609 | -1.350443 | -1.14 | 1.471632 | 200 | 14.716 |
| 1000 | 10.117523 | 0.433487 | -0.333060 | -0.77 | 0.447845 | 200 | 14.162 |
| 10000 | 10.323103 | 0.145967 | -0.127481 | -0.87 | 0.141440 | 200 | 14.144 |
| 100000 | 10.397135 | 0.046227 | -0.053449 | -1.16 | 0.047037 | 200 | 14.874 |
| 1000000 | 10.441925 | 0.014695 | -0.008659 | -0.59 | 0.014184 | 50 | 14.184 |
| 10000000 | 10.455111 | 0.004657 | +0.004528 | +0.97 | - | - | - |

**Reading it.** The last column is the convergence check: if error ∝ 1/√N, RMSE × √N is
constant (at the discounted payoff's standard deviation, ≈ 14.5 here). It holds across
five orders of magnitude. Each extra decimal digit costs 100× the paths.

The single-seed errors are all negative up to 1M paths. That is **not bias**: every row
reuses seed 42, so the 1,000-path run starts with the 100-path run's paths and their errors
are correlated row to row. The RMSE column uses independent seeds for exactly this reason.

**The error bar is checked, not assumed.** A pricer that under-reports its std error would
pass every "within 4 se" test for the wrong reason, so the tests also run 400 independent
pricings and check that the average is unbiased, that the actual RMS error matches the
reported se to within 15%, and that 95% intervals cover the true price ~380/400 times.

## 3. Antithetic variates: measured variance reduction

1,000,000 payoff evaluations per estimate (antithetic = 500,000 pairs Z, −Z), seed 42.
VRF = (se_plain / se_anti)² at **equal payoff evaluations**; implied ρ = corr(f(Z), f(−Z))
= 1/VRF − 1. Efficiency = VRF × t_plain / t_anti, i.e. variance reduction per unit of CPU
time.

| contract | BS | plain (± se) | antithetic (± se) | VRF | implied rho | t plain (ms) | t anti (ms) | efficiency |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| call K=70 (deep ITM) | 33.5401 | 33.5248 ± 0.0199 | 33.5346 ± 0.0048 | 17.36 | -0.942 | 35.5 | 20.5 | 30.10 |
| call K=85 | 20.4693 | 20.4558 ± 0.0184 | 20.4616 ± 0.0076 | 5.87 | -0.830 | 35.5 | 20.5 | 10.16 |
| call K=100 (ATM) | 10.4506 | 10.4419 ± 0.0147 | 10.4439 ± 0.0104 | 2.01 | -0.501 | 35.3 | 20.4 | 3.46 |
| call K=115 | 4.4666 | 4.4593 ± 0.0101 | 4.4602 ± 0.0090 | 1.24 | -0.194 | 35.4 | 20.4 | 2.15 |
| call K=130 | 1.6396 | 1.6318 ± 0.0061 | 1.6312 ± 0.0059 | 1.07 | -0.067 | 35.4 | 20.4 | 1.85 |
| call K=160 (far OTM) | 0.1590 | 0.1556 ± 0.0018 | 0.1576 ± 0.0019 | 0.98 | +0.022 | 35.5 | 20.4 | 1.70 |
| put K=100 (ATM) | 5.5735 | 5.5790 ± 0.0087 | 5.5707 ± 0.0066 | 1.71 | -0.416 | 35.1 | 20.5 | 2.92 |
| put K=130 (ITM) | 25.2994 | 25.3058 ± 0.0169 | 25.2948 ± 0.0024 | 50.25 | -0.980 | 35.1 | 20.4 | 86.28 |

Cross-check, ATM call, 10,000 paths: VRF from reported se = 2.02; VRF from the actual RMSE
over 2,000 independent seeds = 1.91 (each MSE is known to ~±3% at 2,000 runs, so the
ratio to ~±4.5%).

**One formula explains the table.** The pair average has variance Var(f)(1 + ρ)/2, so at
equal work VRF = 1/(1 + ρ):

- **Deep ITM**, a call is nearly linear in S_T, so f(Z) and f(−Z) move almost exactly
  opposite: ρ → −1, VRF 17–50×.
- **ATM**, the payoff is linear on one side of the kink and flat on the other: ρ ≈ −0.5,
  VRF ≈ 2.
- **Far OTM**, Z and −Z almost never both pay: ρ ≈ 0, VRF ≈ 1. Theory says ρ = −μ²/σ² ≈
  −0.008 here (VRF ≈ 1.01); the measured 0.98 is estimation noise — over 20 seeds the VRF
  averaged 1.001 (range 0.95–1.03).

For any monotone payoff (a vanilla call or put) ρ ≤ 0, so antithetic never hurts. A
symmetric payoff such as a straddle would have ρ > 0 and antithetic would *increase*
variance — stated from theory, not measured (the engine prices calls and puts only).

**Antithetic is also ~1.7× faster** at equal payoff evaluations: it does half as many
normal draws and accumulator updates, and that per-unit work costs more than the exp and
the payoff (≈ 30 ns vs ≈ 5.5 ns; see [benchmarks.md](benchmarks.md)). That is why
efficiency exceeds VRF on every row, even the far-OTM call.

**The std error is computed over pairs.** The independent unit is (f(Z) + f(−Z))/2, not
each half; the calibration test from §2 is run for both modes. Over 4,000 independent
seeds, actual MSE / mean reported se² = 1.001 (antithetic) and 0.990 (plain).

## 4. Monte Carlo Greeks against analytic Greeks

1,000,000 paths, seed 42. Bump-and-reprice with **common random numbers**: every bumped
payoff is evaluated on the same draw Z as the base payoff, and each Greek is the sample
mean of the per-path finite differences, which also gives it a std error. Default bumps:
spot 1% of S, vol 0.001, expiry 1 day, rate 1 bp.

Two error sources are reported separately:

- **bump bias** = (Black-Scholes finite difference with the same bumps) − (analytic
  Greek). Deterministic; no Monte Carlo involved.
- **z** = (MC − BS finite difference) / se. Pure sampling noise, should look like N(0, 1).

**ATM call** (all five Greeks in one pass: 71 ms)

| Greek | analytic | bump bias | MC (± se) | rel. se | z |
|---|---:|---:|---:|---:|---:|
| delta | 0.63683 | -8.6e-05 | 0.63649 ± 0.00057 | 0.09% | -0.44 |
| gamma | 0.01876 | -2.3e-06 | 0.01870 ± 0.00011 | 0.57% | -0.57 |
| vega | 37.52403 | -3.0e-05 | 37.47228 ± 0.07548 | 0.20% | -0.69 |
| theta | -6.41403 | -3.7e-06 | -6.40786 ± 0.00885 | 0.14% | +0.70 |
| rho | 53.23248 | -7.7e-07 | 53.21212 ± 0.04722 | 0.09% | -0.43 |

**ATM put** (all five Greeks in one pass: 71 ms)

| Greek | analytic | bump bias | MC (± se) | rel. se | z |
|---|---:|---:|---:|---:|---:|
| delta | -0.36317 | -8.6e-05 | -0.36337 ± 0.00041 | 0.11% | -0.28 |
| gamma | 0.01876 | -2.3e-06 | 0.01870 ± 0.00011 | 0.57% | -0.57 |
| vega | 37.52403 | -3.0e-05 | 37.55537 ± 0.04775 | 0.13% | +0.66 |
| theta | -1.65788 | -3.7e-06 | -1.66002 ± 0.00290 | 0.17% | -0.74 |
| rho | -41.89046 | -9.3e-07 | -41.91082 ± 0.04722 | 0.11% | -0.43 |

All ten |z| ≤ 0.74. The bump bias is below the 1M-path sampling error for every Greek, but
by very different margins: bias / se ≈ 0.15 for delta, 0.02 for gamma, 4e-4 for vega and
theta, 2e-5 for rho. Delta is the tight one — a 1% spot bump still picks up delta's
curvature — so its bias would become visible at roughly 50× more paths (se shrinks √50 ≈
7×). For high-precision delta, shrink `Bumps::spot_rel`: with CRN this costs nothing in
delta noise (§5), only in gamma noise, so in practice delta and gamma would use different
bumps. Call and put gamma agree to every digit because C − P = S·e^{−qT} − K·e^{−rT}
holds *path by path*: the payoffs differ by a linear function of S_T, whose second
difference is zero on every draw.

## 5. Bump size: bias vs. noise, and why common random numbers

ATM call, 200,000 paths. "Independent" = the up and down prices repriced with different
seeds, se = √(se_up² + se_dn²) / 2h.

| spot bump h | delta se (CRN) | delta se (independent) | ratio | gamma (CRN, ± se) | gamma rel. se | gamma bump bias |
|---:|---:|---:|---:|---:|---:|---:|
| 10% | 0.00117 | 0.00240 | 2x | 0.01853 ± 0.00006 | 0.3% | -2.3e-04 |
| 3% | 0.00125 | 0.00779 | 6x | 0.01871 ± 0.00013 | 0.7% | -2.1e-05 |
| 1% | 0.00128 | 0.02330 | 18x | 0.01904 ± 0.00024 | 1.3% | -2.3e-06 |
| 0.1% | 0.00129 | 0.23289 | 181x | 0.01974 ± 0.00079 | 4.2% | -2.3e-08 |
| 0.01% | 0.00129 | 2.32893 | 1809x | 0.02235 ± 0.00267 | 14.2% | -1.9e-11 |

- **Delta, independent seeds:** both prices' full noise survives into (V_up − V_dn)/2h,
  so the std error grows as 1/h — 1,809× worse than CRN at a 0.01% bump.
- **Delta, CRN:** the per-path difference is bounded (the payoff is Lipschitz in S), so
  the noise is flat in h and small bumps are free. At h = 1% the predicted ratio is
  √2·sd(payoff)/(2h) ÷ sd(per-path delta) ≈ 10.3 / 0.55 ≈ 18.5; measured 18.
- **Gamma, even with CRN:** the second difference is nonzero only on paths ending within
  ~h of the strike, where it is ~1/h, so variance ∝ 1/h and se ∝ 1/√h (0.3% → 14% as h
  shrinks 1000×), while the O(h²) bias falls 10⁷×. At 10% bias dominates (1.2% bias vs
  0.3% noise); below ~0.1% noise dominates. The 1% default sits between. Pathwise or
  likelihood-ratio estimators avoid this; they are out of scope for v1.

## 6. Known limitation: rare events

A put ~4.7 standard deviations out of the money (S=100, K=80, T=0.1, σ=15%) has a true
price of 6.2e-7 — roughly one path in a million pays. With 200,000 paths, none did, and
plain MC returned **0 ± 0**: wrong, with an error bar claiming certainty. The CLT-based std
error only means something once enough paths land where the payoff is nonzero. This is
pinned down by a test (`MonteCarloLimitation.FarOtmWithNoHitsReportsZeroPlusMinusZero`)
rather than hidden; the standard fix is importance sampling (shift the drift toward the
strike and reweight), which is not implemented.

## 7. How the statistical tests stay deterministic

Every MC test uses a fixed seed, so it passes or fails deterministically, with thresholds
at ~4 standard errors so they do not depend on a lucky seed. That claim was checked by
re-running each statistical check under thousands of other seeds (agreement within 4 se:
0/700 plain, 0/700 antithetic, 0/2,400 Greek checks; std-error calibration: 0/100 blocks
per mode; convergence-rate and CRN/gamma scaling checks: 0/200 and 0/100). The sweep also
caught a threshold that was *not* safe: "antithetic VRF > 0.9 for monotone payoffs" failed
1.5% of seeds at 200k paths on a short-dated OTM put, whose rare, heavy-tailed payoff makes
the VRF estimate noisy (sd ≈ 0.05). That test now uses 2M paths (sd ≈ 0.016; minimum over
100 seeds 0.970).

## Next steps

- **Importance sampling** for deep-OTM options (fixes §6).
- **Pathwise and likelihood-ratio Greeks**: pathwise delta/vega for lower variance than
  bump-and-reprice; likelihood-ratio gamma to avoid the 1/√h noise of §5.
- **Control variates** (e.g. the discounted terminal stock, whose expectation is known) —
  a variance reduction that, unlike antithetic, also helps away from the money.
- **Quasi-Monte Carlo** (Sobol) for error closer to 1/N than 1/√N on smooth problems.
- **Path-dependent payoffs** (Asian, barrier), which need time-stepping and would
  introduce discretization bias to measure alongside sampling error.
