#ifndef PDK_SPECTRAL_H
#define PDK_SPECTRAL_H

// Tabulated nuclear spectral function S(p, E) for the PDK generator.
//
// Loads a NuWro-style grid spectral function (data/sf/*.grid) and samples the
// joint (momentum, removal-energy) of a bound nucleon directly from it, so the
// `benhar` / `ankowski` momentum models are driven by real tabulated data rather
// than a parametrization. The table already encodes the full mean-field +
// correlated structure and the removal-energy distribution, so no separate SRC
// tail or binding model is applied on top.
//
// NuWro grid format (src/sf/gridfun2d.cc), units MeV:
//   eRes pRes        # of energy points, # of momentum points
//   eMin pMin        axis minima
//   eMax pMax        axis maxima
//   then pRes momentum blocks, each: p  followed by eRes "(e  S(p,e))" pairs.
// Stored momentum-outer / energy-inner. Axes are converted MeV -> GeV on load;
// only relative S values matter for sampling.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace pdk {

struct SpectralTable {
    int nP = 0, nE = 0;          // momentum / energy grid counts
    std::vector<double> p;       // momentum axis [GeV/c]  (size nP)
    std::vector<double> E;       // removal-energy axis [GeV] (size nE)
    std::vector<double> S;       // S(p,E) row-major [iP*nE + jE]
};

// Parse a NuWro grid spectral function. Returns false (with a message on stderr)
// on a missing or malformed file.
inline bool load_nuwro_grid(const std::string& path, SpectralTable& out) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Error: could not open spectral-function file " << path
                  << "\n";
        return false;
    }
    int eRes = 0, pRes = 0;
    double eMin, pMin, eMax, pMax;
    if (!(f >> eRes >> pRes >> eMin >> pMin >> eMax >> pMax) || eRes <= 0 ||
        pRes <= 0) {
        std::cerr << "Error: malformed header in " << path << "\n";
        return false;
    }
    out.nE = eRes;
    out.nP = pRes;
    out.p.assign(pRes, 0.0);
    out.E.assign(eRes, 0.0);
    out.S.assign(static_cast<std::size_t>(pRes) * eRes, 0.0);

    constexpr double kMeVtoGeV = 1.0e-3;
    for (int i = 0; i < pRes; ++i) {
        double p_mev;
        if (!(f >> p_mev)) {
            std::cerr << "Error: truncated momentum block " << i << " in "
                      << path << "\n";
            return false;
        }
        out.p[i] = p_mev * kMeVtoGeV;
        for (int j = 0; j < eRes; ++j) {
            double e_mev, val;
            if (!(f >> e_mev >> val)) {
                std::cerr << "Error: truncated S(p,E) data in " << path << "\n";
                return false;
            }
            if (i == 0) out.E[j] = e_mev * kMeVtoGeV;
            out.S[static_cast<std::size_t>(i) * eRes + j] = val;
        }
    }
    return true;
}

// Samples (p, E) from a tabulated S(p,E). The physical sampling weight is
// proportional to p^2 S(p,E) (since the spectral function is normalized as
// integral of 4*pi*p^2 S dp dE = 1). Momentum is drawn from the p^2-weighted
// marginal, then the removal energy from S(p,.) conditioned on that momentum;
// both are smeared uniformly within their grid cells.
class SpectralSampler {
public:
    SpectralSampler() = default;

    explicit SpectralSampler(SpectralTable table) : t_(std::move(table)) {
        const int nP = t_.nP, nE = t_.nE;
        dp_ = nP > 1 ? t_.p[1] - t_.p[0] : 0.0;
        dE_ = nE > 1 ? t_.E[1] - t_.E[0] : 0.0;

        // Per-row cumulative over E, and the p^2-weighted momentum marginal.
        cumE_.assign(static_cast<std::size_t>(nP) * nE, 0.0);
        rowSum_.assign(nP, 0.0);
        cumP_.assign(nP, 0.0);
        double running_p = 0.0;
        for (int i = 0; i < nP; ++i) {
            double running_e = 0.0;
            for (int j = 0; j < nE; ++j) {
                double s = t_.S[static_cast<std::size_t>(i) * nE + j];
                if (s < 0.0) s = 0.0;
                running_e += s;
                cumE_[static_cast<std::size_t>(i) * nE + j] = running_e;
            }
            rowSum_[i] = running_e;
            running_p += t_.p[i] * t_.p[i] * running_e;
            cumP_[i] = running_p;
        }
        total_ = running_p;
    }

    bool valid() const { return total_ > 0.0; }

    // Draw a momentum [GeV/c] and removal energy [GeV] into p_out, e_out.
    void sample(std::mt19937& gen, double& p_out, double& e_out) const {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        const int nE = t_.nE;

        // Momentum index from the cumulative p^2-weighted marginal.
        double up = u(gen) * total_;
        int i = static_cast<int>(
            std::upper_bound(cumP_.begin(), cumP_.end(), up) - cumP_.begin());
        if (i >= t_.nP) i = t_.nP - 1;

        // Energy index from the conditional cumulative of row i.
        const double* rowCum = &cumE_[static_cast<std::size_t>(i) * nE];
        double ue = u(gen) * rowSum_[i];
        int j = static_cast<int>(std::upper_bound(rowCum, rowCum + nE, ue) -
                                 rowCum);
        if (j >= nE) j = nE - 1;

        // Smear uniformly within the grid cell; keep both non-negative.
        p_out = t_.p[i] + (u(gen) - 0.5) * dp_;
        e_out = t_.E[j] + (u(gen) - 0.5) * dE_;
        if (p_out < 0.0) p_out = 0.0;
        if (e_out < 0.0) e_out = 0.0;
    }

private:
    SpectralTable t_;
    std::vector<double> cumE_;    // per-row cumulative over E [nP*nE]
    std::vector<double> rowSum_;  // sum_j S(p_i,E_j) [nP]
    std::vector<double> cumP_;    // cumulative p^2-weighted momentum marginal [nP]
    double total_ = 0.0;
    double dp_ = 0.0, dE_ = 0.0;
};

}  // namespace pdk

#endif  // PDK_SPECTRAL_H
