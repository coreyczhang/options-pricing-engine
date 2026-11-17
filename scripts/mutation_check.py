"""Mutation check: does the test suite catch deliberately injected bugs?

    python3 scripts/mutation_check.py [--keep] [MUTANT_NAME ...]

Copies the repository (without build directories) to a temporary directory,
then for each mutant applies its edits to the copy, rebuilds ope_tests, runs
the full suite and records which tests fail. The working tree is never
modified. A mutant no test catches is a hole in the validation; see
docs/results.md for the write-up.

The point is the agreement between Monte Carlo and Black-Scholes: a bug in
only one pricer should break agreement, while a convention shared wrongly by
both cannot — those "shared" mutants must be caught by something independent
of Monte Carlo (the textbook reference values).
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import namedtuple
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
Edit = namedtuple("Edit", "file old new")
Mutant = namedtuple("Mutant", "name group description edits")

MC = "src/MonteCarlo.cpp"
BS = "src/BlackScholes.cpp"
RNG = "include/ope/NormalRng.hpp"
SSVI = "src/Ssvi.cpp"

MC_PRICE_DRIFT = "const double drift = (mkt.rate - mkt.div_yield - 0.5 * mkt.vol * mkt.vol) * T;"
MC_PRICE_DISCOUNT = "const double discount = std::exp(-mkt.rate * T);"
MC_PRICE_ADD = ("stats.add(cfg.antithetic ? 0.5 * (payoff_at(z) + payoff_at(-z)) : "
                "payoff_at(z));")
MC_RESULT = "return Result{discount * stats.mean(), discount * stats.std_error(), cfg.paths};"
MC_SCENARIO = "return Scenario{(rate - q - 0.5 * vol * vol) * t, vol * std::sqrt(t),"
MC_SCENARIO_DISC = "                        std::exp(-rate * t)};"
BS_FWD = "t.fwd_spot = mkt.spot * std::exp(-mkt.div_yield * T);"
BS_PVK = "t.pv_strike = opt.strike * std::exp(-mkt.rate * T);"
BS_DISC_Q = "const double disc_q = std::exp(-q * T);"

MUTANTS = [
    # ---- Monte Carlo pricer only: should break MC-vs-BS agreement ----------
    Mutant("mc-drift-no-ito", "MC only", "MC drift drops the -sigma^2/2 Ito term",
           [Edit(MC, MC_PRICE_DRIFT,
                 "const double drift = (mkt.rate - mkt.div_yield) * T;")]),
    Mutant("mc-drift-no-div", "MC only", "MC drift ignores the dividend yield",
           [Edit(MC, MC_PRICE_DRIFT,
                 "const double drift = (mkt.rate - 0.5 * mkt.vol * mkt.vol) * T;")]),
    Mutant("mc-discount-r-minus-q", "MC only", "MC discounts at r - q instead of r",
           [Edit(MC, MC_PRICE_DISCOUNT,
                 "const double discount = std::exp(-(mkt.rate - mkt.div_yield) * T);")]),
    Mutant("mc-vol-sqrt2", "MC only",
           "Box-Muller radius sqrt(-log u) instead of sqrt(-2 log u): normals have variance 1/2",
           [Edit(RNG, "std::sqrt(-2.0 * std::log(u1))", "std::sqrt(-std::log(u1))")]),
    # ---- Error bar: agreement tests alone pass more easily, not less ------
    Mutant("se-inflated-2x", "error bar", "reported std error is 2x too large",
           [Edit(MC, MC_RESULT,
                 "return Result{discount * stats.mean(), 2.0 * discount * stats.std_error(), "
                 "cfg.paths};")]),
    Mutant("se-understated-half", "error bar", "reported std error is half what it should be",
           [Edit(MC, MC_RESULT,
                 "return Result{discount * stats.mean(), 0.5 * discount * stats.std_error(), "
                 "cfg.paths};")]),
    Mutant("antithetic-se-over-halves", "error bar",
           "antithetic std error computed treating Z and -Z halves as independent",
           [Edit(MC, MC_PRICE_ADD,
                 "if (cfg.antithetic) { stats.add(payoff_at(z)); stats.add(payoff_at(-z)); } "
                 "else { stats.add(payoff_at(z)); }")]),
    Mutant("rng-duplicate-draws", "error bar",
           "Box-Muller spare uses cos again: every normal is drawn twice",
           [Edit(RNG, "spare_ = radius * std::sin(angle);", "spare_ = radius * std::cos(angle);")]),
    Mutant("se-n-not-n-minus-1", "error bar",
           "variance divides by n instead of n - 1: a 1/(2n) change, below any statistical test",
           [Edit(MC, "return std::sqrt(m2_ / (n - 1.0) / n);", "return std::sqrt(m2_ / n / n);")]),
    # ---- Antithetic --------------------------------------------------------
    Mutant("antithetic-no-negation", "antithetic",
           "antithetic pair evaluates f(Z) twice instead of f(Z) and f(-Z)",
           [Edit(MC, MC_PRICE_ADD,
                 "stats.add(cfg.antithetic ? 0.5 * (payoff_at(z) + payoff_at(z)) : "
                 "payoff_at(z));")]),
    # ---- Monte Carlo Greeks ------------------------------------------------
    Mutant("greeks-delta-over-h", "MC Greeks", "delta central difference divided by h, not 2h",
           [Edit(MC, "c.delta = (up - dn) / (2.0 * h_s);", "c.delta = (up - dn) / h_s;")]),
    Mutant("greeks-theta-sign", "MC Greeks", "theta returned as dV/dT instead of dV/dt",
           [Edit(MC, "c.theta = -(s_t_up.value", "c.theta = (s_t_up.value")]),
    Mutant("greeks-vega-no-crn", "MC Greeks",
           "vega's down-bump uses -Z: unbiased, but common random numbers are lost",
           [Edit(MC, "s_vol_dn.value(opt, S, z)", "s_vol_dn.value(opt, S, -z)")]),
    Mutant("greeks-theta-no-crn", "MC Greeks",
           "theta's down-bump uses -Z: unbiased, but common random numbers are lost",
           [Edit(MC, "s_t_dn.value(opt, S, z)", "s_t_dn.value(opt, S, -z)")]),
    Mutant("greeks-rho-no-crn", "MC Greeks",
           "rho's down-bump uses -Z: unbiased, but common random numbers are lost",
           [Edit(MC, "s_r_dn.value(opt, S, z)", "s_r_dn.value(opt, S, -z)")]),
    Mutant("greeks-delta-no-crn", "MC Greeks",
           "delta's down-bump uses an independent-looking draw (-Z): CRN lost for delta",
           [Edit(MC, "const double dn = base.discount * payoff(opt.type, K, (S - h_s) * g);",
                 "const double dn = base.discount * payoff(opt.type, K, (S - h_s) * "
                 "std::exp(base.drift - base.vol_sqrt_t * z));")]),
    # ---- Black-Scholes only -----------------------------------------------
    Mutant("bs-swap-d1-d2", "BS only", "call price uses N(d2) on the stock leg and N(d1) on the strike leg",
           [Edit(BS, "return t.fwd_spot * norm_cdf(t.d1) - t.pv_strike * norm_cdf(t.d2);",
                 "return t.fwd_spot * norm_cdf(t.d2) - t.pv_strike * norm_cdf(t.d1);")]),
    Mutant("bs-d1-minus-half-var", "BS only", "d1 uses -sigma sqrt(T)/2 instead of +",
           [Edit(BS, "/ sig_sqrt_t + 0.5 * sig_sqrt_t;", "/ sig_sqrt_t - 0.5 * sig_sqrt_t;")]),
    # ---- Shared: the same wrong convention in both pricers ----------------
    Mutant("shared-ignore-dividends", "shared",
           "both pricers ignore the dividend yield",
           [Edit(BS, BS_FWD, "t.fwd_spot = mkt.spot;"),
            Edit(BS, BS_DISC_Q, "const double disc_q = 1.0;"),
            Edit(MC, MC_PRICE_DRIFT,
                 "const double drift = (mkt.rate - 0.5 * mkt.vol * mkt.vol) * T;"),
            Edit(MC, MC_SCENARIO, "return Scenario{(rate - 0.5 * vol * vol) * t, vol * std::sqrt(t),")]),
    Mutant("shared-no-discounting", "shared",
           "both pricers return the undiscounted (forward) value",
           [Edit(BS, BS_FWD, "t.fwd_spot = mkt.spot * std::exp((mkt.rate - mkt.div_yield) * T);"),
            Edit(BS, BS_PVK, "t.pv_strike = opt.strike;"),
            Edit(MC, MC_PRICE_DISCOUNT, "const double discount = 1.0;"),
            Edit(MC, MC_SCENARIO_DISC, "                        1.0};")]),
    Mutant("shared-div-as-rate", "shared",
           "both pricers grow the stock at r + q (dividend sign flipped)",
           [Edit(BS, BS_FWD, "t.fwd_spot = mkt.spot * std::exp(+mkt.div_yield * T);"),
            Edit(BS, BS_DISC_Q, "const double disc_q = std::exp(+q * T);"),
            Edit(MC, MC_PRICE_DRIFT,
                 "const double drift = (mkt.rate + mkt.div_yield - 0.5 * mkt.vol * mkt.vol) * T;"),
            Edit(MC, MC_SCENARIO, "return Scenario{(rate + q - 0.5 * vol * vol) * t, vol * std::sqrt(t),")]),
    # ---- Joint SSVI fit ----------------------------------------------------
    Mutant("ssvi-weight-1-over-w", "SSVI",
           "weights 1/w instead of 1/(w T): long expiries count T times more",
           [Edit(SSVI, "s.wt.push_back(1.0 / (d.w[i] * d.expiry));", "s.wt.push_back(1.0 / d.w[i]);")]),
    Mutant("ssvi-no-pav", "SSVI", "thetas fitted per slice, never pooled: calendar arbitrage allowed",
           [Edit(SSVI, "while (stack.size() >= 2 && stack[stack.size() - 2].theta > stack.back().theta)",
                 "while (false && stack.size() >= 2)")]),
    Mutant("ssvi-eta-bound-2x", "SSVI", "eta allowed up to 4 / (1 + |rho|) instead of 2 / (1 + |rho|)",
           [Edit(SSVI, "return {rho, (1.0 + std::sin(x[2])) / (1.0 + std::abs(rho)),",
                 "return {rho, 2.0 * (1.0 + std::sin(x[2])) / (1.0 + std::abs(rho)),")]),
    Mutant("ssvi-flat-theta-before-first", "SSVI",
           "theta held at theta_1 before the first expiry instead of scaled by T / T_1",
           [Edit(SSVI, "if (T <= expiries.front()) return thetas.front() * T / expiries.front();",
                 "if (T <= expiries.front()) return thetas.front();")]),
]

# Tests that compare Monte Carlo against Black-Scholes directly (substrings, so
# parameterized instances like .../Plain and .../Antithetic both match).
AGREEMENT = ("AgreesWithBlackScholesWithinStdErrors",
             "MatchBlackScholesFiniteDifferenceWithinStdErrors")


def run(cmd, cwd, **kw):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, **kw)


def build_and_test(src, build):
    b = run(["cmake", "--build", str(build), "--target", "ope_tests", "--parallel"], src)
    if b.returncode != 0:
        return None, b.stdout[-2000:] + b.stderr[-2000:]
    out = build / "result.json"
    if out.exists():
        out.unlink()
    run([str(build / "tests" / "ope_tests"), f"--gtest_output=json:{out}"], src)
    data = json.loads(out.read_text())
    failed = []
    for suite in data["testsuites"]:
        for t in suite["testsuite"]:
            if t.get("failures"):
                failed.append(f"{suite['name']}.{t['name']}")
    return failed, None


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    keep = "--keep" in sys.argv
    selected = [m for m in MUTANTS if not args or m.name in args]

    work = Path(tempfile.mkdtemp(prefix="ope-mutation-"))
    src, build = work / "src", work / "build"
    shutil.copytree(ROOT, src, ignore=shutil.ignore_patterns("build*", ".git", "cmake-build-*"))
    gtest = ROOT / "build" / "_deps" / "googletest-src"
    cfg = ["cmake", "-S", str(src), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release"]
    if gtest.exists():  # reuse the already-fetched GoogleTest: no network needed
        cfg.append(f"-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST={gtest}")
    c = run(cfg, src)
    if c.returncode != 0:
        sys.exit(c.stderr)

    baseline, err = build_and_test(src, build)
    if err or baseline:
        sys.exit(f"baseline must pass cleanly: {err or baseline}")
    print(f"baseline: all tests pass ({work})\n", flush=True)

    results = []
    for m in selected:
        originals = {}
        for e in m.edits:
            path = src / e.file
            if e.file not in originals:
                originals[e.file] = path.read_text()
            text = path.read_text()
            if text.count(e.old) != 1:
                sys.exit(f"{m.name}: edit target not found exactly once in {e.file}: {e.old!r}")
            path.write_text(text.replace(e.old, e.new))
        t0 = time.time()
        failed, err = build_and_test(src, build)
        for f, text in originals.items():
            (src / f).write_text(text)
        if err:
            print(f"{m.name}: BUILD FAILED\n{err}")
            sys.exit(1)
        caught_by_agreement = [f for f in failed if any(a in f for a in AGREEMENT)]
        results.append({"name": m.name, "group": m.group, "description": m.description,
                        "failed": failed, "agreement": caught_by_agreement})
        status = "CAUGHT" if failed else "SURVIVED"
        print(f"[{status:8}] {m.name:28} {len(failed):3} failing "
              f"({len(caught_by_agreement)} MC-vs-BS)  {time.time() - t0:.0f}s", flush=True)
        for f in failed[:8]:
            print(f"             - {f}")
        if len(failed) > 8:
            print(f"             ... {len(failed) - 8} more")

    # Restore a clean build so --keep leaves a usable tree.
    build_and_test(src, build)
    (work / "mutation_results.json").write_text(json.dumps(results, indent=2))
    survived = [r["name"] for r in results if not r["failed"]]
    print(f"\n{len(results) - len(survived)}/{len(results)} mutants caught; "
          f"survived: {', '.join(survived) or 'none'}")
    print(f"results: {work / 'mutation_results.json'}")
    if not keep:
        shutil.rmtree(src)
        shutil.rmtree(build)
    sys.exit(1 if survived else 0)


if __name__ == "__main__":
    main()
