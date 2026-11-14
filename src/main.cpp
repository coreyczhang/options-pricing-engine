// pricer: price one European option and print its Greeks, or back out the
// implied volatility of a quoted price.
//
//   pricer <call|put> <spot> <strike> <expiry_years> <rate> <vol> [div_yield]
//   pricer iv <call|put> <spot> <strike> <expiry_years> <rate> <price> [div_yield]
//
// Examples (Hull: "Greek letters" worked example; implied-vol example):
//   pricer call 49 50 0.3846 0.05 0.20
//   pricer iv call 21 20 0.25 0.10 1.875

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

#include "ope/BlackScholes.hpp"
#include "ope/ImpliedVol.hpp"

namespace {

int usage() {
    std::fprintf(stderr,
                 "usage: pricer <call|put> <spot> <strike> <expiry_years> <rate> <vol> "
                 "[div_yield]\n"
                 "       pricer iv <call|put> <spot> <strike> <expiry_years> <rate> <price> "
                 "[div_yield]\n");
    return 2;
}

const char* status_text(ope::iv::Status s) {
    switch (s) {
        case ope::iv::Status::Converged: return "converged";
        case ope::iv::Status::AtLowerBound: return "price at no-arbitrage lower bound";
        case ope::iv::Status::BelowLowerBound: return "price below intrinsic (arbitrage)";
        case ope::iv::Status::AboveUpperBound: return "price above upper bound (arbitrage)";
        case ope::iv::Status::VolAboveMax: return "implied vol above cap";
        case ope::iv::Status::NotConverged: return "did not converge";
    }
    return "unknown";
}

int implied_vol_main(int argc, char** argv) {
    if (argc != 8 && argc != 9) return usage();
    const std::string kind = argv[2];
    if (kind != "call" && kind != "put") return usage();

    const ope::EuropeanOption opt{kind == "call" ? ope::OptionType::Call : ope::OptionType::Put,
                                  std::stod(argv[4]), std::stod(argv[5])};
    const ope::Market mkt{std::stod(argv[3]), std::stod(argv[6]), 0.0,
                          argc == 9 ? std::stod(argv[8]) : 0.0};
    const double quote = std::stod(argv[7]);
    const ope::iv::Result r = ope::iv::solve(opt, mkt, quote);

    std::printf("Implied vol %s  S=%g K=%g T=%g r=%g q=%g price=%g\n", kind.c_str(), mkt.spot,
                opt.strike, opt.expiry, mkt.rate, mkt.div_yield, quote);
    std::printf("  status  %s (%d iterations)\n", status_text(r.status), r.iterations);
    if (r.status != ope::iv::Status::Converged && r.status != ope::iv::Status::AtLowerBound)
        return 1;

    ope::Market solved = mkt;
    solved.vol = r.vol;
    std::printf("  vol     %12.8f   (%.4f%%)\n", r.vol, 100.0 * r.vol);
    std::printf("  reprice %12.8f\n", ope::bs::price(opt, solved));
    std::printf("  vega    %12.6f   (per 1.00 vol; small vega => vol loosely determined)\n",
                ope::bs::greeks(opt, solved).vega);
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "iv") {
        try {
            return implied_vol_main(argc, argv);
        } catch (const std::exception& e) {
            std::fprintf(stderr, "error: %s\n", e.what());
            return 1;
        }
    }
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
