#ifndef PDK_MOMENTUM_H
#define PDK_MOMENTUM_H
// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//

// Initial-state proton momentum models for the PDK generator.
//
// Several nuclear models for the momentum distribution n(p) of a bound proton
// in argon-40 are provided. All values are toy/parametrized approximations
// tuned to reasonable argon numbers, suitable for spectrum studies but not for
// a precision oscillation/PDK analysis.
//
//   Polynomial       Toy PDF read from config/params.dat (rejection sampling).
//   GlobalFermiGas   n(p) ~ p^2 for p < k_F (single Fermi sphere).
//   LocalFermiGas    Superposition of Fermi spheres with k_F(r) set by the
//                    local nuclear density rho(r) (Woods-Saxon).
//   SRC              Mean-field (Fermi gas) bulk + a 1/p^4 short-range
//                    correlation tail above k_F.
//   SpectralFunction Local-Fermi-gas mean field + SRC correlation tail
//                    (a toy stand-in for a Benhar-style spectral function).
//   HarmonicOscillator Analytic shell-model n(p): a sum over the argon shells of
//                    the harmonic-oscillator momentum density |phi_{nl}(p)|^2,
//                    weighted by shell occupancy. Uses the (n,l) quantum numbers
//                    of the shell tables and shows genuine shell structure.
//   BodekRitchie     Global Fermi sphere bulk + a high-momentum correlation tail
//                    in the Bodek-Ritchie spirit (PRD 23 (1981) 1070): a larger
//                    tail fraction reaching further in p than the SRC model.
//   Gaussian         Smooth analytic baseline n(p) ~ p^2 exp(-p^2 / 2 sigma^2)
//                    (no sharp Fermi edge); sigma tuned to the argon k_F.
//   CorrelatedFermiGas Sharp Fermi sphere bulk + a 1/p^4 contact tail to ~2 k_F
//                    with an explicit correlated fraction (modern CFG).
//   Benhar           Tabulated spectral function S(p,E) read from a NuWro grid
//                    file (see PDKSpectral.h): draws the joint (momentum,
//                    removal-energy) directly from real argon data (JLab
//                    E12-14-012, proton and neutron). Supply the grid with
//                    --sf-file (defaults to config/sf/gsf_Ar40{P,N}.grid).
//   Ankowski         Effective Ankowski-Sobczyk spectral function built
//                    analytically (no grid): a shell mean field (harmonic-
//                    oscillator momentum profile paired with each shell's
//                    separation energy) plus a correlated 1/p^4 tail carrying
//                    the two-nucleon removal energy.
//
// For both benhar and ankowski the spectral function supersedes the analytic
// momentum sampling AND the binding model, so the BindingModel below does NOT
// apply to them (their removal energy is intrinsic to S(p,E)).
//
// The removal (separation) energy of a *mean-field* proton (analytic models
// only) is set by a separate, independently selectable BindingModel (see below):
//
//   Potential  Momentum-dependent optical potential V(k_F, p) of Juszczak,
//              Nowak & Sobczyk, Eur. Phys. J. C 39 (2005) 195 [nucl-th/0311051];
//              see nuclear_potential() below (default).
//   Constant   Fixed average separation energy e_sep_const.
//   Shell      Argon proton shell levels, Gaussian-smeared about their
//              separation energies (argon_proton_shells()).
//
// Correlated (SRC-tail) protons always take the two-nucleon removal energy
// E ~ E_offset + p^2/2M, independent of the binding model.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <random>
#include <string>

#include <cstdlib>

#include "PDKChannels.h"    // Nucleon, nucleon masses
#include "PDKConfig.h"
#include "PDKKinematics.h"
#include "PDKSpectral.h"    // tabulated S(p,E) loader + sampler

namespace pdk {

enum class MomentumModel {
    Polynomial,
    GlobalFermiGas,
    LocalFermiGas,
    SRC,
    SpectralFunction,
    HarmonicOscillator,
    BodekRitchie,
    Gaussian,
    CorrelatedFermiGas,
    Benhar,
    Ankowski,
};

// Treatment of the removal (separation) energy of a mean-field bound proton.
// Selected independently of the momentum model; correlated SRC-tail protons
// always use the two-nucleon removal energy regardless of this choice.
enum class BindingModel {
    Potential,  // momentum-dependent optical potential V(k_F, p) [nucl-th/0311051]
    Constant,   // fixed separation energy e_sep_const
    Shell,      // argon proton shell levels, Gaussian-smeared
};

// Nuclear parameters for argon-40 (GeV/c and fm units). The density profile and
// the binding potential follow Juszczak, Nowak & Sobczyk, Eur. Phys. J. C 39
// (2005) 195 [nucl-th/0311051]: a two-parameter Fermi (Woods-Saxon) density,
// Eq. (4), and the momentum-dependent optical potential V(k_F, p), Eq. (8).
struct NuclearParams {
    double k_fermi = 0.217;   // proton global Fermi momentum [GeV/c]: <k_F^Ar> = 217 MeV, Eq. (7)
    double k_fermi_n = 0.230; // neutron global Fermi momentum [GeV/c] (N>Z in Ar-40)
    double frac_p = 0.45;     // proton number fraction Z/A for argon-40 (18/40)
    double frac_n = 0.55;     // neutron number fraction N/A for argon-40 (22/40)
    double ws_c = 3.530;      // two-parameter Fermi half-density radius C [fm]
    double ws_a = 0.541;      // two-parameter Fermi diffuseness C1 [fm]
    double rho0 = 0.176;      // central nucleon density [fm^-3]
    double src_fraction = 0.20;  // fraction of nucleons in the SRC tail
    double k_max = 0.65;      // upper edge of the SRC tail [GeV/c]
    double r_max = 12.0;      // radial integration cutoff [fm]

    // Harmonic-oscillator shell model: oscillator length b [fm]. For argon
    // hbar*omega ~ 41 A^(-1/3) ~ 12 MeV, so b = sqrt(hbar / (M omega)) ~ 1.9 fm.
    // The momentum scale of |phi_{nl}(p)|^2 is set by 1/b; b is tunable.
    double ho_b = 1.9;

    // Gaussian model: width sigma [GeV/c]. With n(p) ~ p^2 exp(-p^2/2 sigma^2)
    // the mean square momentum is 3 sigma^2; choosing sigma = k_F / sqrt(5)
    // matches <p^2> = (3/5) k_F^2 of a Fermi sphere of radius k_F.
    double gauss_sigma = 0.217 / 2.2360679775;  // k_F / sqrt(5)

    // Bodek-Ritchie Fermi gas: high-momentum tail fraction and upper edge
    // [GeV/c]. Larger and longer than the SRC tail (Bodek & Ritchie, PRD 23
    // (1981) 1070).
    double br_fraction = 0.25;
    double br_kmax = 1.00;

    // Correlated Fermi gas: correlated (contact-tail) fraction. The 1/p^4 tail
    // runs from k_F up to 2 k_F (derived from the global Fermi momentum).
    double cfg_fraction = 0.20;

    // Ankowski-Sobczyk effective spectral function: fraction of nucleons in the
    // correlated (1/p^4 tail, two-nucleon removal energy) part; the rest form the
    // shell mean field.
    double as_corr_fraction = 0.20;

    // SRC pairs are not described by the mean-field potential; their removal
    // energy comes from the recoiling partner, E ~ offset + p^2 / 2M_p.
    double removal_src_offset = 0.020;  // SRC removal-energy offset [GeV]

    // Binding option.
    double e_sep_const = 0.030;  // fixed separation energy, 'constant' binding [GeV]
};

// Nucleon shell structure of argon-40, used by the 'shell' binding model and the
// Ankowski-Sobczyk spectral function. Separation energies are representative
// values for Ar-40 from (e,e'p) / electron-scattering systematics; the
// harmonic-oscillator quantum numbers (n, l) fix the momentum profile of each
// shell.
struct Shell {
    const char* name;
    int occ;         // nucleon occupancy
    double e_sep;    // mean separation (removal) energy [GeV]
    double e_width;  // Gaussian spread of the level [GeV]
    int n, l;        // radial / orbital quantum numbers (harmonic oscillator)
};

inline const std::vector<Shell>& argon_proton_shells() {
    static const std::vector<Shell> shells = {
        {"1s1/2", 2, 0.0520, 0.0090, 0, 0},
        {"1p3/2", 4, 0.0360, 0.0070, 0, 1},
        {"1p1/2", 2, 0.0320, 0.0060, 0, 1},
        {"1d5/2", 6, 0.0180, 0.0040, 0, 2},
        {"2s1/2", 2, 0.0130, 0.0030, 1, 0},
        {"1d3/2", 2, 0.0125, 0.0030, 0, 2},
    };  // occupancies sum to Z = 18
    return shells;
}

// Neutron shells of argon-40 (N = 22). Neutron separation energies are smaller
// than the proton ones (no Coulomb barrier), and the extra two neutrons start
// filling the 1f7/2 level. Representative values, not a fit.
inline const std::vector<Shell>& argon_neutron_shells() {
    static const std::vector<Shell> shells = {
        {"1s1/2", 2, 0.0480, 0.0090, 0, 0},
        {"1p3/2", 4, 0.0320, 0.0070, 0, 1},
        {"1p1/2", 2, 0.0280, 0.0060, 0, 1},
        {"1d5/2", 6, 0.0150, 0.0040, 0, 2},
        {"2s1/2", 2, 0.0110, 0.0030, 1, 0},
        {"1d3/2", 4, 0.0100, 0.0030, 0, 2},
        {"1f7/2", 2, 0.0099, 0.0030, 0, 3},
    };  // occupancies sum to N = 22
    return shells;
}

// Real part of the momentum-dependent nuclear optical potential V(k_F, p) felt
// by a bound nucleon, from Juszczak, Nowak & Sobczyk, Eur. Phys. J. C 39 (2005)
// 195 [nucl-th/0311051], Eq. (8):
//
//                       (a k_F)^2 (k_F + b)
//   V(k_F, p) = - ------------------------------------
//                  c^4 + d^3 k_F + e^3 p^2 / k_F + p^4
//
// V is negative (binding), deepest at low momentum, and rises smoothly to zero
// at high momentum. The paper fits a..e in MeV; the expression is homogeneous
// of degree 1, so passing k_F and p in GeV with the parameters in GeV returns V
// in GeV. Reproduces Fig. 2: V(p=0) = -59 MeV at k_F = 217 MeV.
inline double nuclear_potential(double k_F, double p) {
    if (k_F <= 0.0) return 0.0;
    constexpr double a = 0.206, b = 0.582, c = -0.322, d = 0.422, e = 0.289;
    const double c4 = c * c * c * c;
    const double d3 = d * d * d;
    const double e3 = e * e * e;
    const double num = (a * k_F) * (a * k_F) * (k_F + b);
    const double den = c4 + d3 * k_F + e3 * p * p / k_F + p * p * p * p;
    return -num / den;
}

// Woods-Saxon (two-parameter Fermi) nucleon density [fm^-3] of argon-40 at
// radius r [fm], rho(r) = rho0 / (1 + exp((r - C)/a)), from Juszczak, Nowak &
// Sobczyk [nucl-th/0311051], Eq. (4). Exposed as a free function so the
// initial-state momentum sampler (NucleonMomentumSampler) and the FSI cascade
// (PDKCascade.h) share one nuclear geometry instead of duplicating it.
inline double nucleon_density(double r, const NuclearParams& np) {
    return np.rho0 / (1.0 + std::exp((r - np.ws_c) / np.ws_a));
}

// Local Fermi momentum [GeV/c] of the nucleon species with number fraction
// `frac` (Z/A for protons, N/A for neutrons; 1/2 for the symmetric case) at
// radius r: k_F(r) = hbar*c (3 pi^2 frac rho(r))^(1/3).
inline double local_fermi_momentum(double r, const NuclearParams& np,
                                   double frac) {
    constexpr double kHbarC = 0.197327;  // GeV*fm
    return kHbarC *
           std::cbrt(3.0 * frac * M_PI * M_PI * nucleon_density(r, np));
}

inline const char* model_name(MomentumModel m) {
    switch (m) {
        case MomentumModel::Polynomial: return "polynomial";
        case MomentumModel::GlobalFermiGas: return "global Fermi gas";
        case MomentumModel::LocalFermiGas: return "local Fermi gas";
        case MomentumModel::SRC: return "short-range correlations";
        case MomentumModel::SpectralFunction: return "spectral function";
        case MomentumModel::HarmonicOscillator: return "harmonic-oscillator shell model";
        case MomentumModel::BodekRitchie: return "Bodek-Ritchie Fermi gas";
        case MomentumModel::Gaussian: return "Gaussian";
        case MomentumModel::CorrelatedFermiGas: return "correlated Fermi gas";
        case MomentumModel::Benhar: return "tabulated spectral function (benhar)";
        case MomentumModel::Ankowski: return "analytic spectral function (ankowski)";
    }
    return "unknown";
}

// Parse a model name (case-insensitive, with aliases). Returns false on no match.
inline bool parse_model(std::string s, MomentumModel& out) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (s == "polynomial" || s == "poly") out = MomentumModel::Polynomial;
    else if (s == "gfg" || s == "global") out = MomentumModel::GlobalFermiGas;
    else if (s == "lfg" || s == "local") out = MomentumModel::LocalFermiGas;
    else if (s == "src") out = MomentumModel::SRC;
    else if (s == "sf" || s == "spectral") out = MomentumModel::SpectralFunction;
    else if (s == "hosm" || s == "ho" || s == "shellmom") out = MomentumModel::HarmonicOscillator;
    else if (s == "br" || s == "bodek") out = MomentumModel::BodekRitchie;
    else if (s == "gauss" || s == "gaussian") out = MomentumModel::Gaussian;
    else if (s == "cfg") out = MomentumModel::CorrelatedFermiGas;
    else if (s == "benhar" || s == "bfr") out = MomentumModel::Benhar;
    else if (s == "ankowski" || s == "as") out = MomentumModel::Ankowski;
    else return false;
    return true;
}

inline const char* binding_name(BindingModel b) {
    switch (b) {
        case BindingModel::Potential: return "optical potential V(k_F, p)";
        case BindingModel::Constant: return "constant separation energy";
        case BindingModel::Shell: return "shell-model separation energies";
    }
    return "unknown";
}

// Parse a binding-model name (case-insensitive, with aliases). False on no match.
inline bool parse_binding(std::string s, BindingModel& out) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (s == "potential" || s == "pot" || s == "optical") out = BindingModel::Potential;
    else if (s == "constant" || s == "const" || s == "fixed") out = BindingModel::Constant;
    else if (s == "shell" || s == "shells") out = BindingModel::Shell;
    else return false;
    return true;
}

// Sampled initial state of the struck nucleon.
struct NucleonState {
    double p;      // momentum magnitude [GeV/c]
    double e_rem;  // removal (separation) energy [GeV]
};

// Samples the initial nucleon state under the chosen model for the given parent
// nucleon (proton or neutron). Construct once (it precomputes a sampling
// envelope) and call sample().
class NucleonMomentumSampler {
public:
    NucleonMomentumSampler(MomentumModel model, Config cfg, NuclearParams np,
                           BindingModel binding = BindingModel::Potential,
                           Nucleon nucleon = Nucleon::Proton,
                           const std::string& sf_path = "")
        : model_(model),
          binding_(binding),
          nucleon_(nucleon),
          cfg_(std::move(cfg)),
          np_(np) {
        m_nucleon_ = nucleon_mass(nucleon_);

        // The tabulated benhar model reads its joint (p, E) from a NuWro grid
        // file; the table supersedes the analytic momentum/binding treatment.
        // (ankowski builds its S(p,E) analytically and needs no grid.)
        if (model_ == MomentumModel::Benhar) {
            SpectralTable table;
            if (sf_path.empty() || !load_nuwro_grid(sf_path, table)) {
                std::cerr << "Error: model '" << model_name(model_)
                          << "' requires a spectral-function grid (--sf-file).\n";
                std::exit(1);
            }
            spectral_ = SpectralSampler(std::move(table));
            if (!spectral_.valid()) {
                std::cerr << "Error: spectral-function grid " << sf_path
                          << " has no positive weight.\n";
                std::exit(1);
            }
        }
        k_fermi_ = (nucleon_ == Nucleon::Proton) ? np_.k_fermi : np_.k_fermi_n;
        frac_ = (nucleon_ == Nucleon::Proton) ? np_.frac_p : np_.frac_n;

        // Precompute the maximum of the radial weight r^2 * rho(r) for the
        // Woods-Saxon rejection sampling used by the local-density models.
        const int steps = 2000;
        for (int i = 0; i <= steps; ++i) {
            double r = np_.r_max * i / steps;
            r2rho_max_ = std::max(r2rho_max_, r * r * density(r));
        }

        // Total occupancy for occupancy-weighted shell selection (the 'shell'
        // binding option), using this nucleon's shell table.
        for (const auto& s : shells()) total_occ_ += s.occ;

        // Precompute the envelope peak of each shell's harmonic-oscillator
        // momentum profile for the rejection sampling used by the 'hosm' model.
        const double ho_pmax = 6.0 * kHbarC / np_.ho_b;  // GeV/c
        for (const auto& s : shells()) {
            double wmax = 0.0;
            const int hsteps = 400;
            for (int i = 1; i <= hsteps; ++i) {
                wmax = std::max(wmax, ho_weight(s.n, s.l, ho_pmax * i / hsteps));
            }
            ho_wmax_.push_back(wmax);
        }
    }

    NucleonState sample(std::mt19937& gen) const {
        switch (model_) {
            case MomentumModel::Polynomial: {
                // Toy PDF; bind it at the global k_F.
                double p = sample_momentum(cfg_, gen);
                return {p, mean_field_removal(k_fermi_, p, gen)};
            }
            case MomentumModel::GlobalFermiGas: {
                double p = sample_fermi_sphere(k_fermi_, gen);
                return {p, mean_field_removal(k_fermi_, p, gen)};
            }
            case MomentumModel::LocalFermiGas: {
                double kf;
                double p = sample_lfg(gen, kf);
                return {p, mean_field_removal(kf, p, gen)};
            }
            case MomentumModel::SRC:
                if (draw_src(gen)) {
                    double p = sample_tail(k_fermi_, gen);
                    return {p, src_removal(p)};
                } else {
                    double p = sample_fermi_sphere(k_fermi_, gen);
                    return {p, mean_field_removal(k_fermi_, p, gen)};
                }
            case MomentumModel::SpectralFunction:
                if (draw_src(gen)) {
                    double p = sample_tail(k_fermi_, gen);
                    return {p, src_removal(p)};
                } else {
                    double kf;
                    double p = sample_lfg(gen, kf);
                    return {p, mean_field_removal(kf, p, gen)};
                }
            case MomentumModel::HarmonicOscillator: {
                // Pick a shell (occupancy-weighted), draw p from its HO momentum
                // profile, and bind it. Passing the shell index keeps the 'shell'
                // binding coherent (n,l <-> e_sep); 'potential'/'constant' ignore it.
                int shell = pick_shell(gen);
                double p = sample_ho_shell(shell, gen);
                return {p, mean_field_removal(k_fermi_, p, gen, shell)};
            }
            case MomentumModel::BodekRitchie:
                if (draw_frac(gen, np_.br_fraction)) {
                    double p = sample_tail_to(k_fermi_, np_.br_kmax, gen);
                    return {p, src_removal(p)};
                } else {
                    double p = sample_fermi_sphere(k_fermi_, gen);
                    return {p, mean_field_removal(k_fermi_, p, gen)};
                }
            case MomentumModel::Gaussian: {
                double p = sample_gaussian(gen);
                return {p, mean_field_removal(k_fermi_, p, gen)};
            }
            case MomentumModel::CorrelatedFermiGas:
                if (draw_frac(gen, np_.cfg_fraction)) {
                    double p = sample_tail_to(k_fermi_, 2.0 * k_fermi_, gen);
                    return {p, src_removal(p)};
                } else {
                    double p = sample_fermi_sphere(k_fermi_, gen);
                    return {p, mean_field_removal(k_fermi_, p, gen)};
                }
            case MomentumModel::Benhar: {
                // Tabulated spectral function: draw (p, E_rem) jointly from the
                // loaded NuWro grid (no separate SRC tail or binding model).
                double p, e_rem;
                spectral_.sample(gen, p, e_rem);
                return {p, e_rem};
            }
            case MomentumModel::Ankowski: {
                // Effective (Ankowski-Sobczyk) S(p,E) built analytically: a
                // correlated 1/p^4 tail with two-nucleon removal energy, plus a
                // shell mean field (HO momentum profile paired with the shell's
                // separation energy). The shell energies are intrinsic, so the
                // --binding choice does not apply (as for benhar).
                if (draw_frac(gen, np_.as_corr_fraction)) {
                    double p = sample_tail(k_fermi_, gen);
                    return {p, src_removal(p)};
                }
                int shell = pick_shell(gen);
                double p = sample_ho_shell(shell, gen);
                return {p, sample_shell_energy(gen, shell)};
            }
        }
        return {0.0, 0.0};
    }

private:
    // Shell table for this nucleon type (proton or neutron).
    const std::vector<Shell>& shells() const {
        return nucleon_ == Nucleon::Proton ? argon_proton_shells()
                                            : argon_neutron_shells();
    }

    // Effective removal (separation) energy of a mean-field nucleon implied by
    // the optical potential. The bound nucleon energy in the nuclear rest frame
    // is E = sqrt(p^2 + M_N^2) + V(k_F, p) (V < 0 binds it), so relative to a
    // free nucleon at rest the removed energy is
    //   E_rem = M_N - E = M_N - sqrt(p^2 + M_N^2) - V(k_F, p).
    // PDKKinematics then forms the off-shell invariant mass W^2 = E^2 - p^2.
    double potential_removal(double k_F, double p) const {
        return m_nucleon_ - std::sqrt(p * p + m_nucleon_ * m_nucleon_) -
               nuclear_potential(k_F, p);
    }

    // An SRC nucleon's partner recoils with momentum -p, so the removal energy
    // grows with momentum: E ~ E_offset + p^2 / (2 M_N).
    double src_removal(double p) const {
        return np_.removal_src_offset + p * p / (2.0 * m_nucleon_);
    }

    // Removal energy of a mean-field proton under the selected binding model.
    // For the shell model, pass the shell that was sampled for the momentum
    // (Ankowski); shell < 0 draws a shell weighted by occupancy.
    double mean_field_removal(double k_F, double p, std::mt19937& gen,
                              int shell = -1) const {
        switch (binding_) {
            case BindingModel::Potential: return potential_removal(k_F, p);
            case BindingModel::Constant: return np_.e_sep_const;
            case BindingModel::Shell: return sample_shell_energy(gen, shell);
        }
        return potential_removal(k_F, p);
    }

    // Pick a nucleon shell at random, weighted by occupancy.
    int pick_shell(std::mt19937& gen) const {
        std::uniform_int_distribution<int> u(0, total_occ_ - 1);
        int x = u(gen);
        const auto& sh = shells();
        for (std::size_t i = 0; i < sh.size(); ++i) {
            x -= sh[i].occ;
            if (x < 0) return static_cast<int>(i);
        }
        return 0;
    }

    // Gaussian-smeared separation energy of a shell (drawing one if shell < 0).
    // Resamples to keep the removal energy positive.
    double sample_shell_energy(std::mt19937& gen, int shell) const {
        if (shell < 0) shell = pick_shell(gen);
        const Shell& s = shells()[shell];
        std::normal_distribution<double> g(s.e_sep, s.e_width);
        double e = g(gen);
        return e > 0.0 ? e : s.e_sep;
    }

    static constexpr double kHbarC = 0.197327;  // GeV*fm

    // Woods-Saxon nucleon density [fm^-3] (shared free function).
    double density(double r) const { return nucleon_density(r, np_); }

    // Local Fermi momentum [GeV/c] from this nucleon's density rho_q = frac*rho
    // (shared free function; frac = Z/A for protons, N/A for neutrons).
    double local_kf(double r) const {
        return local_fermi_momentum(r, np_, frac_);
    }

    bool draw_src(std::mt19937& gen) const { return draw_frac(gen, np_.src_fraction); }

    // True with probability frac (Bernoulli draw), used to split a model into a
    // mean-field bulk and a high-momentum correlation tail.
    bool draw_frac(std::mt19937& gen, double frac) const {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        return u(gen) < frac;
    }

    // Sample p inside a Fermi sphere: n(p) ~ p^2 on [0, kF] -> p = kF * u^(1/3).
    double sample_fermi_sphere(double kf, std::mt19937& gen) const {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        return kf * std::cbrt(u(gen));
    }

    // Local Fermi gas: pick a radius weighted by r^2 rho(r), then sample the
    // local Fermi sphere with k_F(r). Reports the local k_F via kf_out so the
    // caller can evaluate the binding potential at the same density.
    double sample_lfg(std::mt19937& gen, double& kf_out) const {
        std::uniform_real_distribution<double> dr(0.0, np_.r_max);
        std::uniform_real_distribution<double> dy(0.0, r2rho_max_);
        double r;
        while (true) {
            r = dr(gen);
            if (dy(gen) < r * r * density(r)) break;
        }
        kf_out = local_kf(r);
        return sample_fermi_sphere(kf_out, gen);
    }

    // SRC tail: n(p) ~ 1/p^4 on [klow, k_max], so the radial weight p^2 n(p)
    // ~ 1/p^2. Inverse-CDF sampling of 1/p^2.
    double sample_tail(double klow, std::mt19937& gen) const {
        return sample_tail_to(klow, np_.k_max, gen);
    }

    // As sample_tail() but with an explicit upper edge kmax, so the Bodek-Ritchie
    // and correlated-Fermi-gas models can reuse the same 1/p^4 inverse-CDF with
    // their own tail reach.
    double sample_tail_to(double klow, double kmax, std::mt19937& gen) const {
        std::uniform_real_distribution<double> u(0.0, 1.0);
        double inv = 1.0 / klow - u(gen) * (1.0 / klow - 1.0 / kmax);
        return 1.0 / inv;
    }

    // Gaussian model: n(p) ~ p^2 exp(-p^2 / 2 sigma^2) is the speed distribution
    // of a 3D Gaussian momentum vector with per-component width sigma, so draw
    // three N(0, sigma) components and return the magnitude.
    double sample_gaussian(std::mt19937& gen) const {
        std::normal_distribution<double> g(0.0, np_.gauss_sigma);
        double px = g(gen), py = g(gen), pz = g(gen);
        return std::sqrt(px * px + py * py + pz * pz);
    }

    // Generalized Laguerre polynomial L_n^alpha(x) by the standard three-term
    // recurrence. Used by the harmonic-oscillator momentum profile.
    static double laguerre(int n, double alpha, double x) {
        double lkm1 = 0.0, lk = 1.0;  // L_0 = 1
        for (int k = 0; k < n; ++k) {
            double lkp1 = ((2 * k + 1 + alpha - x) * lk - (k + alpha) * lkm1) / (k + 1);
            lkm1 = lk;
            lk = lkp1;
        }
        return lk;
    }

    // Unnormalized radial momentum density p^2 |phi_{nl}(p)|^2 of a 3D isotropic
    // harmonic oscillator with length b [fm], evaluated at p [GeV/c]:
    //   w(p) ~ (p b/hbarc)^(2l+2) [L_n^{l+1/2}(x)]^2 exp(-x),  x = (p b/hbarc)^2.
    double ho_weight(int n, int l, double p) const {
        double x = (p * np_.ho_b / kHbarC);
        x = x * x;
        double lag = laguerre(n, l + 0.5, x);
        return std::pow(x, l + 1) * lag * lag * std::exp(-x);
    }

    // Sample p from the HO momentum profile of shell `shell` by rejection on a
    // bounded p-grid. The momentum scale is ~ hbarc/b, so an upper edge of a few
    // hbarc/b safely covers the distribution.
    double sample_ho_shell(int shell, std::mt19937& gen) const {
        const Shell& s = shells()[shell];
        const double pmax = 6.0 * kHbarC / np_.ho_b;  // GeV/c, well past the bulk
        std::uniform_real_distribution<double> dp(0.0, pmax);
        std::uniform_real_distribution<double> dy(0.0, ho_wmax_[shell]);
        while (true) {
            double p = dp(gen);
            if (dy(gen) < ho_weight(s.n, s.l, p)) return p;
        }
    }

    MomentumModel model_;
    BindingModel binding_;
    Nucleon nucleon_;
    Config cfg_;
    NuclearParams np_;
    double m_nucleon_ = kProtonMass;  // parent nucleon rest mass [GeV]
    double k_fermi_ = 0.217;          // global Fermi momentum for this nucleon
    double frac_ = 0.5;               // number fraction (Z/A or N/A)
    double r2rho_max_ = 0.0;
    int total_occ_ = 0;               // total occupancy (Z or N), for shell binding
    std::vector<double> ho_wmax_;     // per-shell HO envelope peak (hosm model)
    SpectralSampler spectral_;        // tabulated S(p,E) for benhar/ankowski
};

}  // namespace pdk

#endif  // PDK_MOMENTUM_H
