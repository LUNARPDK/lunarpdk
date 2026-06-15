#ifndef PDK_CONFIG_H
#define PDK_CONFIG_H

// Shared configuration + Fermi-momentum sampling for the PDK Monte Carlo.
//
// Both generators (PDKMCGenerator, LunarPDKGenerator) include this header so
// that the proton Fermi-momentum PDF is defined in exactly one place and is
// driven by config/params.dat rather than being hard-coded.

#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace pdk {

struct Config {
    std::vector<double> coeffs;  // polynomial coefficients c0 + c1*p + c2*p^2 + ...
    double p_min;
    double p_max;
    double f_max;  // upper bound on f(p) over [p_min, p_max] for rejection sampling
};

inline Config default_config() {
    return Config{{0.0, 4.0, -16.0, 0.0}, 0.0, 0.25, 0.5};
}

// Evaluate the PDF polynomial via Horner's method. Returns 0 outside the range.
inline double evaluate_poly(double p, const Config& cfg) {
    if (p < cfg.p_min || p > cfg.p_max) return 0.0;
    double result = 0.0;
    for (std::size_t i = cfg.coeffs.size(); i-- > 0;) {
        result = result * p + cfg.coeffs[i];
    }
    return result;
}

// Load a Config from a params file. Lines beginning with '#' are comments.
// Expected tokens: c0 c1 c2 c3  p_min p_max f_max. Falls back to defaults on
// any error (missing file / malformed contents) and reports it on stderr.
inline Config load_config(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: could not open " << filename
                  << "; using built-in defaults.\n";
        return default_config();
    }

    std::ostringstream data;
    std::string line;
    while (std::getline(file, line)) {
        std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') continue;
        data << line << ' ';
    }

    std::istringstream tokens(data.str());
    Config cfg;
    double val;
    while (cfg.coeffs.size() < 4 && tokens >> val) {
        cfg.coeffs.push_back(val);
    }

    if (cfg.coeffs.size() < 4 ||
        !(tokens >> cfg.p_min >> cfg.p_max >> cfg.f_max)) {
        std::cerr << "Warning: malformed " << filename
                  << "; using built-in defaults.\n";
        return default_config();
    }
    return cfg;
}

// Construct a Mersenne-Twister RNG. A given seed yields a reproducible stream;
// std::nullopt draws a non-deterministic seed from std::random_device.
inline std::mt19937 make_rng(std::optional<std::uint32_t> seed) {
    if (seed) return std::mt19937(*seed);
    std::random_device rd;
    return std::mt19937(rd());
}

// Draw a proton Fermi momentum from f(p) via rejection sampling. Guards against
// a too-small f_max in the config, which would otherwise silently bias the
// sample by truncating the PDF.
inline double sample_momentum(const Config& cfg, std::mt19937& gen) {
    std::uniform_real_distribution<double> distP(cfg.p_min, cfg.p_max);
    std::uniform_real_distribution<double> distH(0.0, cfg.f_max);

    while (true) {
        double p = distP(gen);
        double fp = evaluate_poly(p, cfg);
        if (fp > cfg.f_max) {
            std::cerr << "Error: f(" << p << ") = " << fp
                      << " exceeds f_max = " << cfg.f_max
                      << "; rejection sampling is biased. Fix params.dat.\n";
            std::exit(1);
        }
        if (distH(gen) < fp) return p;
    }
}

}  // namespace pdk

#endif  // PDK_CONFIG_H
