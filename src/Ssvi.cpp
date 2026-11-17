#include "ope/Ssvi.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ope::ssvi {

double Params::phi(double theta) const {
    return eta / (std::pow(theta, gamma) * std::pow(1.0 + theta, 1.0 - gamma));
}

double Params::w(double k, double theta) const {
    const double p = phi(theta);
    const double x = p * k + rho;
    return theta / 2.0 * (1.0 + rho * p * k + std::sqrt(x * x + 1.0 - rho * rho));
}

svi::Params Params::raw(double theta) const {
    const double p = phi(theta);
    if (!(p > 0.0)) return {theta, 0.0, rho, 0.0, 1.0};  // eta = 0: flat w = theta
    return {theta * (1.0 - rho * rho) / 2.0, theta * p / 2.0, rho, -rho / p,
            std::sqrt(1.0 - rho * rho) / p};
}

double Fit::theta(double T) const {
    if (T <= expiries.front()) return thetas.front() * T / expiries.front();
    if (T >= expiries.back()) return thetas.back() * T / expiries.back();
    const auto hi = std::upper_bound(expiries.begin(), expiries.end(), T);
    const std::size_t j = static_cast<std::size_t>(hi - expiries.begin());
    const double t = (T - expiries[j - 1]) / (expiries[j] - expiries[j - 1]);
    return thetas[j - 1] + t * (thetas[j] - thetas[j - 1]);
}

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

struct Slice {
    const SliceData* data;
    std::vector<double> wt;  // 1 / (w T)
    double log_lo, log_hi;   // search range for ln theta
};

double slice_sse(const Slice& s, const Params& p, double theta) {
    const double ph = p.phi(theta);
    const double c = 1.0 - p.rho * p.rho;
    double sse = 0.0;
    for (std::size_t i = 0; i < s.data->k.size(); ++i) {
        const double k = s.data->k[i];
        const double x = ph * k + p.rho;
        const double e = theta / 2.0 * (1.0 + p.rho * ph * k + std::sqrt(x * x + c)) - s.data->w[i];
        sse += s.wt[i] * e * e;
    }
    return sse;
}

// One theta shared by slices [b, e): scan ln theta, then golden-section search
// between the neighbours of the best scan point.
std::pair<double, double> profile_theta(const std::vector<Slice>& sl, std::size_t b,
                                        std::size_t e, const Params& p) {
    double lo = kInf, hi = -kInf;
    for (std::size_t i = b; i < e; ++i) {
        lo = std::min(lo, sl[i].log_lo);
        hi = std::max(hi, sl[i].log_hi);
    }
    auto f = [&](double x) {
        double s = 0.0;
        for (std::size_t i = b; i < e; ++i) s += slice_sse(sl[i], p, std::exp(x));
        return s;
    };
    constexpr int kScan = 60;
    const double step = (hi - lo) / (kScan - 1);
    int best = 0;
    double fbest = kInf;
    for (int j = 0; j < kScan; ++j) {
        const double v = f(lo + step * j);
        if (v < fbest) fbest = v, best = j;
    }
    double a = lo + step * std::max(best - 1, 0);
    double c = lo + step * std::min(best + 1, kScan - 1);
    constexpr double kGold = 0.6180339887498949;
    double x1 = c - kGold * (c - a), x2 = a + kGold * (c - a);
    double f1 = f(x1), f2 = f(x2);
    while (c - a > 1e-12) {
        if (f1 < f2) {
            c = x2, x2 = x1, f2 = f1;
            x1 = c - kGold * (c - a), f1 = f(x1);
        } else {
            a = x1, x1 = x2, f1 = f2;
            x2 = a + kGold * (c - a), f2 = f(x2);
        }
    }
    const double x = f1 < f2 ? x1 : x2;
    return {std::exp(x), std::min(f1, f2)};
}

// Best non-decreasing thetas for fixed globals, by pool-adjacent-violators.
double profile_thetas(const std::vector<Slice>& sl, const Params& p, std::vector<double>* out) {
    struct Block {
        std::size_t b, e;
        double theta, sse;
    };
    std::vector<Block> stack;
    for (std::size_t i = 0; i < sl.size(); ++i) {
        const auto [th, f] = profile_theta(sl, i, i + 1, p);
        stack.push_back({i, i + 1, th, f});
        while (stack.size() >= 2 && stack[stack.size() - 2].theta > stack.back().theta) {
            const std::size_t b = stack[stack.size() - 2].b, e = stack.back().e;
            stack.pop_back();
            stack.pop_back();
            const auto [pth, pf] = profile_theta(sl, b, e, p);
            stack.push_back({b, e, pth, pf});
        }
    }
    double total = 0.0;
    if (out) out->assign(sl.size(), 0.0);
    for (const Block& blk : stack) {
        total += blk.sse;
        if (out)
            for (std::size_t i = blk.b; i < blk.e; ++i) (*out)[i] = blk.theta;
    }
    return total;
}

// R^3 -> the arbitrage-free box: rho in [-1, 1], gamma in [0, 1/2],
// eta in [0, 2 / (1 + |rho|)].
using Vec = std::array<double, 3>;

Params to_params(const Vec& x) {
    const double rho = std::sin(x[0]);
    return {rho, (1.0 + std::sin(x[2])) / (1.0 + std::abs(rho)), 0.25 * (1.0 + std::sin(x[1]))};
}

Vec from_params(double rho, double gamma, double eta_frac) {
    return {std::asin(rho), std::asin(4.0 * gamma - 1.0), std::asin(2.0 * eta_frac - 1.0)};
}

struct Vertex {
    Vec x;
    double f;
};

// Nelder-Mead (standard coefficients) in 3 dimensions.
template <class F>
Vertex nelder_mead(const F& f, const Vec& x0, double step) {
    std::array<Vertex, 4> s;
    s[0] = {x0, f(x0)};
    for (int d = 0; d < 3; ++d) {
        Vec x = x0;
        x[d] += step;
        s[d + 1] = {x, f(x)};
    }
    for (int iter = 0; iter < 5000; ++iter) {
        std::sort(s.begin(), s.end(), [](const Vertex& a, const Vertex& b) { return a.f < b.f; });
        double size = 0.0;
        for (int v = 1; v < 4; ++v)
            for (int d = 0; d < 3; ++d) size = std::max(size, std::abs(s[v].x[d] - s[0].x[d]));
        if (size < 1e-10) break;
        Vec c{};
        for (int v = 0; v < 3; ++v)
            for (int d = 0; d < 3; ++d) c[d] += s[v].x[d] / 3.0;
        auto at = [&](double t) {
            Vec x;
            for (int d = 0; d < 3; ++d) x[d] = c[d] + t * (s[3].x[d] - c[d]);
            return Vertex{x, f(x)};
        };
        const Vertex refl = at(-1.0);
        if (refl.f < s[0].f) {
            const Vertex exp = at(-2.0);
            s[3] = exp.f < refl.f ? exp : refl;
        } else if (refl.f < s[2].f) {
            s[3] = refl;
        } else {
            const Vertex con = refl.f < s[3].f ? at(-0.5) : at(0.5);
            if (con.f < std::min(refl.f, s[3].f)) {
                s[3] = con;
            } else {  // shrink toward the best vertex
                for (int v = 1; v < 4; ++v) {
                    for (int d = 0; d < 3; ++d) s[v].x[d] = (s[v].x[d] + s[0].x[d]) / 2.0;
                    s[v].f = f(s[v].x);
                }
            }
        }
    }
    return *std::min_element(s.begin(), s.end(),
                             [](const Vertex& a, const Vertex& b) { return a.f < b.f; });
}

}  // namespace

Fit fit(std::vector<SliceData> data) {
    if (data.empty()) throw std::invalid_argument("need at least one slice");
    std::sort(data.begin(), data.end(),
              [](const SliceData& a, const SliceData& b) { return a.expiry < b.expiry; });
    std::vector<Slice> sl;
    int points = 0;
    for (std::size_t j = 0; j < data.size(); ++j) {
        const SliceData& d = data[j];
        if (!(std::isfinite(d.expiry) && d.expiry > 0.0))
            throw std::invalid_argument("need finite expiry > 0");
        if (j > 0 && d.expiry == data[j - 1].expiry)
            throw std::invalid_argument("repeated expiry");
        if (d.k.size() != d.w.size()) throw std::invalid_argument("k and w sizes differ");
        if (d.k.size() < 3) throw std::invalid_argument("need >= 3 points per slice");
        Slice s{&d, {}, 0.0, 0.0};
        double wmin = kInf, wmax = 0.0;
        for (std::size_t i = 0; i < d.k.size(); ++i) {
            if (!std::isfinite(d.k[i]) || !(std::isfinite(d.w[i]) && d.w[i] > 0.0))
                throw std::invalid_argument("need finite k and finite w > 0");
            s.wt.push_back(1.0 / (d.w[i] * d.expiry));
            wmin = std::min(wmin, d.w[i]);
            wmax = std::max(wmax, d.w[i]);
        }
        // ATM variance theta is bounded by the smile: min w = theta (1 - rho^2)
        // for an SSVI slice, so this range is generous unless |rho| ~ 1.
        s.log_lo = std::log(0.01 * wmin);
        s.log_hi = std::log(20.0 * wmax);
        points += static_cast<int>(d.k.size());
        sl.push_back(std::move(s));
    }

    auto objective = [&](const Vec& x) {
        for (double v : x)
            if (!std::isfinite(v)) return kInf;
        return profile_thetas(sl, to_params(x), nullptr);
    };

    // Grid of starts, then Nelder-Mead from the best few.
    std::vector<Vertex> starts;
    for (double rho : {-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9})
        for (double gamma : {0.1, 0.3, 0.5})
            for (double eta_frac : {0.2, 0.5, 0.9}) {
                const Vec x = from_params(rho, gamma, eta_frac);
                starts.push_back({x, objective(x)});
            }
    std::sort(starts.begin(), starts.end(),
              [](const Vertex& a, const Vertex& b) { return a.f < b.f; });
    Vertex best{{}, kInf};
    for (std::size_t i = 0; i < 3 && i < starts.size(); ++i) {
        const Vertex v = nelder_mead(objective, starts[i].x, 0.2);
        if (v.f < best.f) best = v;
    }

    Fit out{};
    out.params = to_params(best.x);
    profile_thetas(sl, out.params, &out.thetas);
    out.points = points;
    double sq = 0.0;
    for (std::size_t j = 0; j < sl.size(); ++j) {
        const SliceData& d = *sl[j].data;
        out.expiries.push_back(d.expiry);
        for (std::size_t i = 0; i < d.k.size(); ++i) {
            const double e = std::sqrt(out.params.w(d.k[i], out.thetas[j]) / d.expiry) -
                             std::sqrt(d.w[i] / d.expiry);
            sq += e * e;
            out.max_abs_vol = std::max(out.max_abs_vol, std::abs(e));
        }
    }
    out.rmse_vol = std::sqrt(sq / points);
    return out;
}

}  // namespace ope::ssvi
