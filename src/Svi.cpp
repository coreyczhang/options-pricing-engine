#include "ope/Svi.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ope::svi {

double Params::w(double k) const {
    const double x = k - m;
    return a + b * (rho * x + std::sqrt(x * x + sigma * sigma));
}

double Params::dw(double k) const {
    const double x = k - m;
    return b * (rho + x / std::sqrt(x * x + sigma * sigma));
}

double Params::d2w(double k) const {
    const double x = k - m;
    const double r = std::sqrt(x * x + sigma * sigma);
    return b * sigma * sigma / (r * r * r);
}

DurrlemanMin min_durrleman_g(const Params& p, double k_lo, double k_hi, int n) {
    DurrlemanMin worst{std::numeric_limits<double>::infinity(), k_lo};
    for (int i = 0; i < n; ++i) {
        const double k = k_lo + (k_hi - k_lo) * i / (n - 1);
        const double w = p.w(k);
        if (!(w > 0.0)) return {-std::numeric_limits<double>::infinity(), k};
        const double w1 = p.dw(k);
        const double t = 1.0 - k * w1 / (2.0 * w);
        const double g = t * t - w1 * w1 / 4.0 * (1.0 / w + 0.25) + p.d2w(k) / 2.0;
        if (g < worst.g) worst = {g, k};
    }
    return worst;
}

namespace {

// ---- Inner problem: weighted least squares in (a, d, c) over a polytope ----

struct Inner {
    double sse;
    double a, d, c;
};

// Solve the n x n system in M (row-major, augmented with the rhs as column n)
// by Gaussian elimination with partial pivoting. False if (near-)singular.
template <std::size_t N>
bool solve(std::array<std::array<double, N + 1>, N>& M, int n, std::array<double, N>& x) {
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::abs(M[r][col]) > std::abs(M[piv][col])) piv = r;
        if (std::abs(M[piv][col]) < 1e-300) return false;
        std::swap(M[col], M[piv]);
        for (int r = col + 1; r < n; ++r) {
            const double f = M[r][col] / M[col][col];
            for (int c = col; c <= n; ++c) M[r][c] -= f * M[col][c];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = M[r][n];
        for (int c = r + 1; c < n; ++c) s -= M[r][c] * x[c];
        x[r] = s / M[r][r];
    }
    return true;
}

// Minimize sum_i wt_i (a + d y_i + c sqrt(y_i^2 + 1) - w_i)^2 subject to
// G (a, d, c) <= h. The problem is a convex QP in 3 variables with 6
// constraints; its optimum is the equality-constrained optimum of some active
// set of at most 3 constraints, so trying all 42 such sets and keeping the
// best feasible one is exact.
Inner inner_fit(const std::vector<double>& k, const std::vector<double>& w,
                const std::vector<double>& wt, double m, double sigma, double amax) {
    // Normal equations: objective = th' N th - 2 r' th + C.
    double N[3][3] = {}, r[3] = {}, C = 0.0;
    for (std::size_t i = 0; i < k.size(); ++i) {
        const double y = (k[i] - m) / sigma;
        const double x[3] = {1.0, y, std::sqrt(y * y + 1.0)};
        for (int p = 0; p < 3; ++p) {
            for (int q = 0; q < 3; ++q) N[p][q] += wt[i] * x[p] * x[q];
            r[p] += wt[i] * x[p] * w[i];
        }
        C += wt[i] * w[i] * w[i];
    }
    const double G[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, 1, -1}, {0, -1, -1}, {0, 1, 1}, {0, -1, 1}};
    const double h[6] = {0.0, amax, 0.0, 0.0, 2.0 * sigma, 2.0 * sigma};
    const double tol = 1e-12 * std::max({1.0, amax, 2.0 * sigma});

    Inner best{std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0};
    for (int set = 0; set < 64; ++set) {
        int active[3], na = 0;
        for (int j = 0; j < 6 && na <= 3; ++j)
            if (set & (1 << j)) {
                if (na == 3) { na = 4; break; }
                active[na++] = j;
            }
        if (na > 3) continue;

        // KKT system: [N G_S'; G_S 0] [th; lambda] = [r; h_S].
        constexpr std::size_t kMax = 6;
        const int n = 3 + na;
        std::array<std::array<double, kMax + 1>, kMax> M{};
        for (int p = 0; p < 3; ++p) {
            for (int q = 0; q < 3; ++q) M[p][q] = N[p][q];
            for (int s = 0; s < na; ++s) M[p][3 + s] = G[active[s]][p];
            M[p][n] = r[p];
        }
        for (int s = 0; s < na; ++s) {
            for (int q = 0; q < 3; ++q) M[3 + s][q] = G[active[s]][q];
            M[3 + s][n] = h[active[s]];
        }
        std::array<double, kMax> sol{};
        if (!solve<kMax>(M, n, sol)) continue;

        bool feasible = true;
        for (int j = 0; j < 6 && feasible; ++j)
            feasible = G[j][0] * sol[0] + G[j][1] * sol[1] + G[j][2] * sol[2] <= h[j] + tol;
        if (!feasible) continue;

        double sse = C;
        for (int p = 0; p < 3; ++p) {
            sse -= 2.0 * r[p] * sol[p];
            for (int q = 0; q < 3; ++q) sse += sol[p] * N[p][q] * sol[q];
        }
        if (sse < best.sse) best = {sse, sol[0], sol[1], sol[2]};
    }
    return best;
}

// ---- Outer problem: (m, ln sigma) by grid search, then Nelder-Mead ---------

struct Point {
    double m, s, f;
};

}  // namespace

Fit fit(const std::vector<double>& k, const std::vector<double>& w, double expiry) {
    if (k.size() != w.size()) throw std::invalid_argument("k and w sizes differ");
    if (k.size() < 5) throw std::invalid_argument("SVI has 5 parameters: need >= 5 points");
    if (!(std::isfinite(expiry) && expiry > 0.0)) throw std::invalid_argument("need expiry > 0");
    for (std::size_t i = 0; i < k.size(); ++i)
        if (!std::isfinite(k[i]) || !(std::isfinite(w[i]) && w[i] > 0.0))
            throw std::invalid_argument("need finite k and finite w > 0");

    std::vector<double> wt(w.size());
    for (std::size_t i = 0; i < w.size(); ++i) wt[i] = 1.0 / w[i];
    const double amax = *std::max_element(w.begin(), w.end());
    const auto [kmin_it, kmax_it] = std::minmax_element(k.begin(), k.end());
    const double kmin = *kmin_it, kmax = *kmax_it;
    const double span = std::max(kmax - kmin, 1e-3);

    auto objective = [&](double m, double s) {
        if (!std::isfinite(m) || !std::isfinite(s) || s < std::log(1e-4) || s > std::log(10.0))
            return std::numeric_limits<double>::infinity();
        return inner_fit(k, w, wt, m, std::exp(s), amax).sse;
    };

    // Coarse grid: m over the data range widened by half a span each side,
    // sigma log-spaced over [1e-3, 2].
    Point best{0.0, 0.0, std::numeric_limits<double>::infinity()};
    constexpr int kGridM = 41, kGridS = 25;
    for (int i = 0; i < kGridM; ++i) {
        const double m = kmin - 0.5 * span + 2.0 * span * i / (kGridM - 1);
        for (int j = 0; j < kGridS; ++j) {
            const double s = std::log(1e-3) + (std::log(2.0) - std::log(1e-3)) * j / (kGridS - 1);
            const double f = objective(m, s);
            if (f < best.f) best = {m, s, f};
        }
    }

    // Nelder-Mead from the best grid point (standard coefficients).
    std::array<Point, 3> simplex = {best,
                                    Point{best.m + 0.05 * span, best.s, 0.0},
                                    Point{best.m, best.s + 0.3, 0.0}};
    for (int v = 1; v < 3; ++v) simplex[v].f = objective(simplex[v].m, simplex[v].s);
    for (int iter = 0; iter < 2000; ++iter) {
        std::sort(simplex.begin(), simplex.end(),
                  [](const Point& p, const Point& q) { return p.f < q.f; });
        const double size = std::max(std::abs(simplex[2].m - simplex[0].m) / span,
                                     std::abs(simplex[2].s - simplex[0].s));
        if (size < 1e-11) break;
        const Point c{(simplex[0].m + simplex[1].m) / 2, (simplex[0].s + simplex[1].s) / 2, 0.0};
        auto at = [&](double t) {
            Point p{c.m + t * (simplex[2].m - c.m), c.s + t * (simplex[2].s - c.s), 0.0};
            p.f = objective(p.m, p.s);
            return p;
        };
        const Point refl = at(-1.0);
        if (refl.f < simplex[0].f) {
            const Point exp = at(-2.0);
            simplex[2] = exp.f < refl.f ? exp : refl;
        } else if (refl.f < simplex[1].f) {
            simplex[2] = refl;
        } else {
            const Point con = refl.f < simplex[2].f ? at(-0.5) : at(0.5);
            if (con.f < std::min(refl.f, simplex[2].f)) {
                simplex[2] = con;
            } else {  // shrink toward the best vertex
                for (int v = 1; v < 3; ++v) {
                    simplex[v].m = (simplex[v].m + simplex[0].m) / 2;
                    simplex[v].s = (simplex[v].s + simplex[0].s) / 2;
                    simplex[v].f = objective(simplex[v].m, simplex[v].s);
                }
            }
        }
    }
    const Point& opt = *std::min_element(simplex.begin(), simplex.end(),
                                         [](const Point& p, const Point& q) { return p.f < q.f; });

    const double sigma = std::exp(opt.s);
    const Inner in = inner_fit(k, w, wt, opt.m, sigma, amax);
    Fit out{};
    out.params = {in.a, in.c / sigma, in.c > 0.0 ? in.d / in.c : 0.0, opt.m, sigma};
    out.expiry = expiry;
    out.points = static_cast<int>(k.size());
    double sq = 0.0;
    for (std::size_t i = 0; i < k.size(); ++i) {
        const double e = std::sqrt(std::max(out.params.w(k[i]), 0.0) / expiry) -
                         std::sqrt(w[i] / expiry);
        sq += e * e;
        out.max_abs_vol = std::max(out.max_abs_vol, std::abs(e));
    }
    out.rmse_vol = std::sqrt(sq / k.size());
    return out;
}

}  // namespace ope::svi
