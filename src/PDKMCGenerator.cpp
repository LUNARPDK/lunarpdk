// Quick check of the proton Fermi-momentum sampler: loads the PDF from
// config/params.dat and prints a few sampled momenta.
//
// Usage: PDKMCGenerator [n_samples] [config] [seed]

#include <cstdint>
#include <iostream>
#include <optional>
#include <string>

#include "PDKConfig.h"

int main(int argc, char* argv[]) {
    int n_samples = 5;
    std::string config_path = "config/params.dat";
    std::optional<std::uint32_t> seed;

    try {
        if (argc > 1) n_samples = std::stoi(argv[1]);
        if (argc > 2) config_path = argv[2];
        if (argc > 3) seed = static_cast<std::uint32_t>(std::stoul(argv[3]));
    } catch (const std::exception&) {
        std::cerr << "Usage: " << argv[0] << " [n_samples] [config] [seed]\n";
        return 1;
    }

    pdk::Config cfg = pdk::load_config(config_path);
    std::mt19937 gen = pdk::make_rng(seed);

    std::cout << "Loaded model: f(p) = " << cfg.coeffs[0] << " + "
              << cfg.coeffs[1] << "p + " << cfg.coeffs[2] << "p^2 + "
              << cfg.coeffs[3] << "p^3,  p in [" << cfg.p_min << ", "
              << cfg.p_max << "]\n";

    for (int i = 0; i < n_samples; ++i) {
        double p = pdk::sample_momentum(cfg, gen);
        std::cout << "Sample " << i + 1
                  << " | proton momentum: " << p << " GeV/c\n";
    }

    return 0;
}
