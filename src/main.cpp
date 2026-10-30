// pricer: price one European option and print its Greeks.
//
//   pricer <call|put> <spot> <strike> <expiry_years> <rate> <vol> [div_yield]
//
// Example (Hull, "Greek letters" chapter worked example):
//   pricer call 49 50 0.3846 0.05 0.20

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

#include "ope/BlackScholes.hpp"

namespace {

int usage() {
    std::fprintf(stderr,
                 "usage: pricer <call|put> <spot> <strike> <expiry_years> <rate> <vol> "
                 "[div_yield]\n");
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 7 && argc != 8) return usage();

    const std::string kind = argv[1];
    if (kind != "call" && kind != "put") return usage();

    try {
        const ope::EuropeanOption opt{
            kind == "call" ? ope::OptionType::Call : ope::OptionType::Put,
            std::stod(argv[3]), std::stod(argv[4])};
        const ope::Market mkt{std::stod(argv[2]), std::stod(argv[5]), std::stod(argv[6]),
                              argc == 8 ? std::stod(argv[7]) : 0.0};

        const double px = ope::bs::price(opt, mkt);
        const ope::Greeks g = ope::bs::greeks(opt, mkt);

        std::printf("Black-Scholes %s  S=%g K=%g T=%g r=%g sigma=%g q=%g\n", kind.c_str(),
                    mkt.spot, opt.strike, opt.expiry, mkt.rate, mkt.vol, mkt.div_yield);
        std::printf("  price  %12.6f\n", px);
        std::printf("  delta  %12.6f\n", g.delta);
        std::printf("  gamma  %12.6f\n", g.gamma);
        std::printf("  vega   %12.6f   (per 1.00 vol; %.6f per vol-point)\n", g.vega,
                    g.vega / 100.0);
        std::printf("  theta  %12.6f   (per year; %.6f per calendar day)\n", g.theta,
                    g.theta / 365.0);
        std::printf("  rho    %12.6f   (per 1.00 rate; %.6f per 1%%)\n", g.rho, g.rho / 100.0);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return EXIT_SUCCESS;
}
