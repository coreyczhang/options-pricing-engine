#pragma once

#include <optional>
#include <vector>

#include "ope/Option.hpp"
#include "ope/Ssvi.hpp"
#include "ope/Svi.hpp"

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
    std::optional<svi::Fit> svi;  ///< set when the surface uses SVI and the slice has >= 5 nodes
};

/// How each expiry's smile is drawn through its nodes.
enum class SmileModel {
    Linear,  ///< linear in total variance between nodes: exact at quotes, follows their noise
    Svi,     ///< least-squares raw-SVI fit per expiry (svi::fit): smooth, wings by model
    Ssvi,    ///< one joint SSVI fit to all expiries (ssvi::fit): arbitrage-free by construction
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
/// not repaired, unless a fitted SmileModel is chosen.
class VolSurface {
public:
    /// Throws std::invalid_argument on an invalid spot, rate or dividend
    /// yield, or if no quote can be inverted.
    ///
    /// With SmileModel::Svi each expiry with >= 5 nodes is fitted with raw SVI
    /// and the fit replaces the linear smile at that expiry, including beyond
    /// the outermost strikes (SVI's linear-in-k wings instead of flat vol);
    /// slices with fewer nodes stay linear. Interpolation across expiries is
    /// the same as for Linear.
    ///
    /// With SmileModel::Ssvi all slices with >= 3 nodes are fitted jointly
    /// (throws if there are none) and the SSVI surface replaces the
    /// interpolation everywhere, across expiries too (theta linear in T, see
    /// ssvi::Fit::theta). Slices with fewer nodes do not enter the fit; the
    /// surface at their expiry is the model's.
    VolSurface(const std::vector<Quote>& quotes, double spot, double rate,
               double div_yield = 0.0, SmileModel model = SmileModel::Linear);

    /// Interpolated implied vol at (strike, expiry). Throws on K <= 0 or T <= 0.
    double vol(double strike, double expiry) const;

    /// Interpolated total variance w(k, T), k = ln(K / F(T)).
    double total_variance(double log_moneyness, double expiry) const;

    double forward(double expiry) const;

    const std::vector<Slice>& slices() const { return slices_; }
    const std::vector<Skipped>& skipped() const { return skipped_; }

    /// Static-arbitrage checks, within each slice's quoted strike range:
    ///   - calendar: w must not fall from one expiry to the next at fixed k.
    ///     Linear smiles are checked at every node of both slices, which is
    ///     exact (the difference of two piecewise-linear curves is extremal at
    ///     a breakpoint); SVI smiles on a 201-point k grid plus the nodes.
    ///   - butterfly, linear smiles: call prices at the node strikes must be
    ///     decreasing and convex in strike. Only the nodes are checked; linear
    ///     interpolation between them can create arbitrage this doesn't see.
    ///   - butterfly, SVI smiles: Durrleman's g(k) >= 0 (svi::min_durrleman_g),
    ///     reported at the worst k.
    ///   - SSVI: the same calendar grid check, and Durrleman on each slice.
    ///     Both should come back empty: SSVI's constraints rule them out.
    std::vector<Violation> arbitrage_violations() const;

    /// The joint fit when the surface was built with SmileModel::Ssvi.
    const std::optional<ssvi::Fit>& ssvi() const { return ssvi_; }

private:
    double slice_variance(const Slice& s, double log_moneyness) const;

    double spot_;
    double rate_;
    double div_yield_;
    std::vector<Slice> slices_;  // sorted by expiry
    std::vector<Skipped> skipped_;
    std::optional<ssvi::Fit> ssvi_;
};

}  // namespace ope::surface
