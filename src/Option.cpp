#include "ope/Option.hpp"

#include <cmath>
#include <stdexcept>

namespace ope {

void validate(const EuropeanOption& opt, const Market& mkt) {
    const bool finite = std::isfinite(opt.strike) && std::isfinite(opt.expiry) &&
                        std::isfinite(mkt.spot) && std::isfinite(mkt.rate) &&
                        std::isfinite(mkt.vol) && std::isfinite(mkt.div_yield);
    if (!finite) throw std::invalid_argument("non-finite input");
    if (mkt.spot <= 0.0) throw std::invalid_argument("spot must be > 0");
    if (opt.strike <= 0.0) throw std::invalid_argument("strike must be > 0");
    if (opt.expiry < 0.0) throw std::invalid_argument("expiry must be >= 0");
    if (mkt.vol < 0.0) throw std::invalid_argument("vol must be >= 0");
}

}  // namespace ope
