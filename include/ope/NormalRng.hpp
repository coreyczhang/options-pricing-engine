#pragma once

#include <cmath>
#include <cstdint>
#include <random>

namespace ope {

/// Standard-normal generator: std::mt19937_64 + Box-Muller.
///
/// Why not std::normal_distribution: the standard fixes mt19937_64's output
/// sequence bit-for-bit, but leaves the normal transform to the library —
/// libstdc++ and libc++ produce different draws from the same seed. Doing the
/// transform ourselves makes a seed mean the same paths on every platform (up
/// to last-ulp differences in libm's log/sin/cos), so test results and the
/// numbers in docs/ reproduce on Linux CI and on macOS.
class NormalRng {
public:
    explicit NormalRng(std::uint64_t seed) : engine_(seed) {}

    double next() {
        if (has_spare_) {
            has_spare_ = false;
            return spare_;
        }
        // u1 in (0, 1] so log(u1) is finite; u2 in [0, 1).
        const double u1 = 1.0 - uniform();
        const double u2 = uniform();
        const double radius = std::sqrt(-2.0 * std::log(u1));
        const double angle = kTwoPi * u2;
        spare_ = radius * std::sin(angle);
        has_spare_ = true;
        return radius * std::cos(angle);
    }

private:
    static constexpr double kTwoPi = 6.28318530717958647693;

    // Top 53 bits -> a double uniformly spaced on [0, 1).
    double uniform() { return static_cast<double>(engine_() >> 11) * 0x1.0p-53; }

    std::mt19937_64 engine_;
    double spare_ = 0.0;
    bool has_spare_ = false;
};

}  // namespace ope
