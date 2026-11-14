#include "ope/VolSurface.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

#include "ope/BlackScholes.hpp"
#include "ope/ImpliedVol.hpp"

namespace ope::surface {
namespace {

SkipReason skip_reason(iv::Status s) {
    switch (s) {
        case iv::Status::BelowLowerBound: return SkipReason::BelowLowerBound;
        case iv::Status::AboveUpperBound: return SkipReason::AboveUpperBound;
        case iv::Status::AtLowerBound: return SkipReason::AtLowerBound;
        case iv::Status::VolAboveMax: return SkipReason::VolAboveMax;
        default: return SkipReason::NotConverged;
    }
}

bool valid_quote(const Quote& q) {
    return std::isfinite(q.strike) && q.strike > 0.0 && std::isfinite(q.expiry) &&
           q.expiry > 0.0 && std::isfinite(q.price) && q.price >= 0.0;
}

}  // namespace

VolSurface::VolSurface(const std::vector<Quote>& quotes, double spot, double rate,
                       double div_yield)
    : spot_(spot), rate_(rate), div_yield_(div_yield) {
    if (!(std::isfinite(spot) && spot > 0.0) || !std::isfinite(rate) ||
        !std::isfinite(div_yield))
        throw std::invalid_argument("invalid spot, rate or dividend yield");

    // Group by (expiry, strike); within a group, OTM quote first.
    std::map<std::pair<double, double>, std::vector<Quote>> groups;
    for (const Quote& q : quotes) {
        if (!valid_quote(q)) {
            skipped_.push_back({q, SkipReason::InvalidQuote});
            continue;
        }
        std::vector<Quote>& g = groups[{q.expiry, q.strike}];
        const bool dup = std::any_of(g.begin(), g.end(),
                                     [&](const Quote& o) { return o.type == q.type; });
        if (dup) skipped_.push_back({q, SkipReason::Duplicate});
        else g.push_back(q);
    }

    // std::map iterates in (expiry, strike) order, so each slice is built in
    // increasing strike and the previous strike's vol warm-starts the next solve.
    Slice* slice = nullptr;
    double guess = 0.0;
    for (auto& [key, group] : groups) {
        const double T = key.first;
        const double K = key.second;
        if (!slice || slice->expiry != T) {
            slices_.push_back({T, forward(T), {}});
            slice = &slices_.back();
            guess = 0.0;
        }
        const OptionType otm = K >= slice->forward ? OptionType::Call : OptionType::Put;
        std::stable_sort(group.begin(), group.end(), [&](const Quote& a, const Quote& b) {
            return (a.type == otm) > (b.type == otm);
        });

        bool used = false;
        for (const Quote& q : group) {
            if (used) {
                skipped_.push_back({q, SkipReason::ReplacedByOtm});
                continue;
            }
            const iv::Result r = iv::solve({q.type, K, T}, {spot_, rate_, guess, div_yield_},
                                           q.price);
            if (r.status != iv::Status::Converged) {
                skipped_.push_back({q, skip_reason(r.status)});
                continue;
            }
            slice->nodes.push_back(
                {K, std::log(K / slice->forward), r.vol, r.vol * r.vol * T, q.type});
            guess = r.vol;
            used = true;
        }
    }
    slices_.erase(std::remove_if(slices_.begin(), slices_.end(),
                                 [](const Slice& s) { return s.nodes.empty(); }),
                  slices_.end());
    if (slices_.empty()) throw std::invalid_argument("no quote could be inverted");
}

double VolSurface::forward(double expiry) const {
    return spot_ * std::exp((rate_ - div_yield_) * expiry);
}

double VolSurface::slice_variance(const Slice& s, double k) const {
    const std::vector<Node>& n = s.nodes;
    if (k <= n.front().log_moneyness) return n.front().total_variance;
    if (k >= n.back().log_moneyness) return n.back().total_variance;
    const auto hi = std::upper_bound(n.begin(), n.end(), k, [](double x, const Node& node) {
        return x < node.log_moneyness;
    });
    const auto lo = hi - 1;
    const double t = (k - lo->log_moneyness) / (hi->log_moneyness - lo->log_moneyness);
    return lo->total_variance + t * (hi->total_variance - lo->total_variance);
}

double VolSurface::total_variance(double k, double T) const {
    if (!(T > 0.0) || !std::isfinite(k)) throw std::invalid_argument("need T > 0, finite k");
    const Slice& first = slices_.front();
    const Slice& last = slices_.back();
    if (T <= first.expiry) return slice_variance(first, k) * T / first.expiry;
    if (T >= last.expiry) return slice_variance(last, k) * T / last.expiry;
    const auto hi = std::upper_bound(slices_.begin(), slices_.end(), T,
                                     [](double x, const Slice& s) { return x < s.expiry; });
    const auto lo = hi - 1;
    const double t = (T - lo->expiry) / (hi->expiry - lo->expiry);
    const double w_lo = slice_variance(*lo, k);
    return w_lo + t * (slice_variance(*hi, k) - w_lo);
}

double VolSurface::vol(double strike, double expiry) const {
    if (!(strike > 0.0) || !(expiry > 0.0))
        throw std::invalid_argument("need strike > 0 and expiry > 0");
    return std::sqrt(total_variance(std::log(strike / forward(expiry)), expiry) / expiry);
}

std::vector<Violation> VolSurface::arbitrage_violations() const {
    std::vector<Violation> out;

    // Calendar: compare consecutive slices at every node of either one.
    for (std::size_t i = 0; i + 1 < slices_.size(); ++i) {
        const Slice& a = slices_[i];
        const Slice& b = slices_[i + 1];
        auto check = [&](const Node& node) {
            const double wa = slice_variance(a, node.log_moneyness);
            const double wb = slice_variance(b, node.log_moneyness);
            if (wb < wa * (1.0 - 1e-10))
                out.push_back({Violation::Kind::Calendar, b.expiry,
                               b.forward * std::exp(node.log_moneyness), wa - wb});
        };
        for (const Node& node : a.nodes) check(node);
        for (const Node& node : b.nodes) check(node);
    }

    // Butterfly: call prices at the node strikes, slopes between neighbours
    // must lie in [-e^{-rT}, 0] and be non-decreasing.
    for (const Slice& s : slices_) {
        const std::size_t n = s.nodes.size();
        if (n < 2) continue;
        std::vector<double> c(n), slope(n - 1);
        for (std::size_t j = 0; j < n; ++j)
            c[j] = bs::price({OptionType::Call, s.nodes[j].strike, s.expiry},
                             {spot_, rate_, s.nodes[j].vol, div_yield_});
        const double tol = 1e-10 * spot_;
        const double min_slope = -std::exp(-rate_ * s.expiry);
        for (std::size_t j = 0; j + 1 < n; ++j) {
            slope[j] = (c[j + 1] - c[j]) / (s.nodes[j + 1].strike - s.nodes[j].strike);
            if (slope[j] > tol || slope[j] < min_slope - tol)
                out.push_back({Violation::Kind::Butterfly, s.expiry, s.nodes[j + 1].strike,
                               slope[j] > 0.0 ? slope[j] : min_slope - slope[j]});
        }
        for (std::size_t j = 0; j + 2 < n; ++j)
            if (slope[j] > slope[j + 1] + tol)
                out.push_back({Violation::Kind::Butterfly, s.expiry, s.nodes[j + 1].strike,
                               slope[j] - slope[j + 1]});
    }
    return out;
}

}  // namespace ope::surface
