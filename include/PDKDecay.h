#ifndef PDK_DECAY_H
#define PDK_DECAY_H

// Optional secondary decays of the escaped unstable mesons, applied AFTER the
// FSI cascade and gated by the --decay-mesons flag. This is a spectrum-level
// treatment of the observable final state, not a full decay package:
//
//   pi0 -> gamma gamma                       (dominant mode)
//   eta -> gamma gamma                       (dominant mode used for all eta)
//   K0  -> K0_S (50%) / K0_L (50%)
//   K0_S -> pi+ pi- (69%) / pi0 pi0 (31%)    (then the pi0 decay further)
//   K0_L                                     (long-lived: treated as escaping)
//
// Stable hadrons (K+, K-, p, n, Lambda, ...), leptons and photons pass through
// unchanged. Decays recurse into their products (e.g. K0_S -> pi0 pi0 -> 4 gamma)
// and conserve four-momentum exactly via the same isotropic two-body kinematics
// used by the cascade (decay_into in PDKCascade.h).

#include <random>
#include <vector>

#include "PDKCascade.h"     // Particle, decay_into, pdg_mass, fsi PDG constants
#include "PDKKinematics.h"  // LorentzVector

namespace pdk {

constexpr int kPdgGamma = 22;
constexpr int kPdgK0S = 310;
constexpr int kPdgK0L = 130;

// True if the meson has a secondary decay modelled here.
inline bool is_unstable_meson(int pdg) {
    return pdg == fsi::kPdgPi0 || pdg == fsi::kPdgEta ||
           pdg == fsi::kPdgK0 || pdg == fsi::kPdgK0bar;
}

// Decay every unstable meson in `in` into its observable daughters, recursing
// into the products. K0_L and all stable particles pass through. Returns the
// fully-decayed final-state list.
inline std::vector<Particle> decay_final_mesons(const std::vector<Particle>& in,
                                                std::mt19937& gen) {
    std::vector<Particle> work(in.rbegin(), in.rend());  // used as a LIFO stack
    std::vector<Particle> out;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::size_t guard = 0;
    while (!work.empty() && guard++ < 100000) {
        Particle q = work.back();
        work.pop_back();

        if (q.pdg == fsi::kPdgPi0 || q.pdg == fsi::kPdgEta) {  // -> gamma gamma
            LorentzVector g1, g2;
            if (decay_into(q.p4, 0.0, 0.0, gen, g1, g2)) {
                work.push_back(Particle{kPdgGamma, g1});
                work.push_back(Particle{kPdgGamma, g2});
            } else {
                out.push_back(q);
            }
            continue;
        }
        if (q.pdg == fsi::kPdgK0 || q.pdg == fsi::kPdgK0bar) {
            if (u(gen) < 0.5) {  // K0_L: long-lived, escapes
                out.push_back(Particle{kPdgK0L, q.p4});
            } else {             // K0_S -> pi+ pi- (69.2%) / pi0 pi0 (30.8%)
                bool charged = u(gen) < 0.692;
                int d1 = charged ? fsi::kPdgPiPlus : fsi::kPdgPi0;
                int d2 = charged ? fsi::kPdgPiMinus : fsi::kPdgPi0;
                LorentzVector a, b;
                if (decay_into(q.p4, pdg_mass(d1), pdg_mass(d2), gen, a, b)) {
                    work.push_back(Particle{d1, a});
                    work.push_back(Particle{d2, b});
                } else {
                    out.push_back(Particle{kPdgK0S, q.p4});
                }
            }
            continue;
        }
        out.push_back(q);  // stable hadron / lepton / photon
    }
    while (!work.empty()) { out.push_back(work.back()); work.pop_back(); }
    return out;
}

}  // namespace pdk

#endif  // PDK_DECAY_H
