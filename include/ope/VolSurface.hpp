#pragma once

#include <vector>

#include "ope/Option.hpp"

namespace ope::surface {

/// One option quote: a contract and its price (e.g. a mid).
struct Quote {
    OptionType type;
    double strike;
    double expiry;  ///< years
    double price;
};

/// Why a quote did not become a surface node.
enum class SkipReason {
    ReplacedByOtm,    ///< the out-of-the-money quote at the same strike/expiry was used instead
    Duplicate,        ///< same type, strike and expiry as an earlier quote
    InvalidQuote,     ///< non-positive strike or expiry, negative or non-finite price
    BelowLowerBound,  ///< price below forward intrinsic: arbitrage
    AboveUpperBound,  ///< price at or above S e^{-qT} (call) / K e^{-rT} (put): arbitrage
    AtLowerBound,     ///< no time value left (e.g. a zero bid); says nothing about vol
    VolAboveMax,
    NotConverged,
};

struct Skipped {
    Quote quote;
    SkipReason reason;
};

/// A quote that made it onto the surface.
struct Node {
    double strike;
    double log_moneyness;   ///< k = ln(K / F(T)), F(T) = S e^{(r-q)T}
    double vol;             ///< implied vol of the quote
    double total_variance;  ///< w = vol^2 T
    OptionType type;        ///< which quote it came from
};

/// All nodes at one expiry, sorted by strike.
struct Slice {
    double expiry;
    double forward;
    std::vector<Node> nodes;
};

/// A static-arbitrage violation found on the quoted nodes.
struct Violation {
    enum class Kind {
        Calendar,   ///< total variance falls with expiry at fixed log-moneyness
        Butterfly,  ///< call price not convex (or not monotone) in strike
    };
    Kind kind;
    double expiry;  ///< the later expiry (Calendar) or the slice's expiry (Butterfly)
    double strike;  ///< where it was detected
    double size;    ///< how far past the bound: variance for Calendar, slope for Butterfly
};

/// An implied-volatility surface built from a set of option quotes.
///
/// Construction inverts every quote with iv::solve and keeps, at each
/// strike/expiry, the out-of-the-money quote (put below the forward, call at or
/// above it) when both are given. Out of the money the price is all time value,
/// so vega/price is large and the implied vol is well determined; the ITM quote
/// is mostly intrinsic value and pins vol down less tightly (results §8). The
/// ITM quote is used only if the OTM one is missing or fails. Quotes that
/// cannot be inverted are recorded in skipped() rather than thrown, because
/// real chains contain stale, crossed and zero-bid quotes.
///
/// Interpolation is in total variance w = sigma^2 T against log-moneyness
/// k = ln(K / F(T)), the coordinates in which no-arbitrage conditions are
/// simplest:
///   - within an expiry: linear in w between nodes, flat w beyond the
///     outermost strikes (constant vol along the slice);
///   - across expiries: linear in w in T at fixed k, which keeps the surface
///     free of calendar arbitrage wherever the quoted slices are; before the
///     first expiry and after the last, w scales with T (constant vol).
///
/// It interpolates; it does not fit or smooth. Noisy quotes give a noisy
/// surface, and arbitrage in the quotes is reported (arbitrage_violations()),
/// not repaired. Fitting a parametric model (e.g. SVI) is the next step.
class VolSurface {
public:
    /// Throws std::invalid_argument on an invalid spot, rate or dividend
    /// yield, or if no quote can be inverted.
    VolSurface(const std::vector<Quote>& quotes, double spot, double rate,
               double div_yield = 0.0);

    /// Interpolated implied vol at (strike, expiry). Throws on K <= 0 or T <= 0.
    double vol(double strike, double expiry) const;

    /// Interpolated total variance w(k, T), k = ln(K / F(T)).
    double total_variance(double log_moneyness, double expiry) const;

    double forward(double expiry) const;

    const std::vector<Slice>& slices() const { return slices_; }
    const std::vector<Skipped>& skipped() const { return skipped_; }

    /// Static-arbitrage checks on the quoted nodes:
    ///   - calendar: at each node's log-moneyness, w must not fall from one
    ///     expiry to the next;
    ///   - butterfly: on each slice, call prices at the node strikes (from the
    ///     node vols) must be decreasing and convex in strike.
    /// Only the nodes are checked; linear interpolation between them can in
    /// principle create arbitrage the node checks don't see.
    std::vector<Violation> arbitrage_violations() const;

private:
    double slice_variance(const Slice& s, double log_moneyness) const;

    double spot_;
    double rate_;
    double div_yield_;
    std::vector<Slice> slices_;  // sorted by expiry
    std::vector<Skipped> skipped_;
};

}  // namespace ope::surface
