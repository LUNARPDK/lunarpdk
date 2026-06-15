#ifndef PDK_CHANNELS_H
#define PDK_CHANNELS_H

// Nucleon-decay channels for the PDK generator.
//
// Each supported mode is a two-body decay of a bound nucleon (proton or neutron)
// into a lepton-side daughter (e+, mu+, or an invisible nu-bar) and a hadron-side
// daughter (a pseudoscalar/charged meson). The generator samples the initial
// nucleon momentum and removal energy (PDKMomentum.h), then performs the off-shell
// two-body decay (PDKKinematics.h) for the selected channel.
//
// Secondary meson decays (pi0 -> gamma gamma, eta, K0 -> ...) are NOT modelled;
// the primary daughter four-momenta are reported as-is.

#include <algorithm>
#include <cctype>
#include <string>

namespace pdk {

// Parent nucleon of a decay channel.
enum class Nucleon { Proton, Neutron };

// Free-nucleon rest masses [GeV].
constexpr double kProtonMass = 0.93827;
constexpr double kNeutronMass = 0.93957;

inline double nucleon_mass(Nucleon n) {
    return n == Nucleon::Proton ? kProtonMass : kNeutronMass;
}

inline const char* nucleon_name(Nucleon n) {
    return n == Nucleon::Proton ? "proton" : "neutron";
}

// Daughter-particle rest masses [GeV]. Neutrino and photon are massless.
constexpr double kMassNu = 0.0;
constexpr double kMassElectron = 0.000511;
constexpr double kMassMuon = 0.10566;
constexpr double kMassPi0 = 0.13498;
constexpr double kMassPiCharged = 0.13957;
constexpr double kMassEta = 0.54786;
constexpr double kMassKCharged = 0.49368;
constexpr double kMassK0 = 0.49761;

// Kept for backward references to the original p -> K+ nu mode.
constexpr double kKaonMass = kMassKCharged;

// A two-body nucleon-decay mode N -> (lepton-side) + (hadron-side).
struct Channel {
    const char* key;     // CLI key, e.g. "pToKnu"
    const char* pretty;  // human-readable, e.g. "p -> K+ nubar"
    Nucleon parent;      // decaying nucleon
    double m1;           // lepton-side daughter mass [GeV]
    double m2;           // hadron-side daughter mass [GeV]
    const char* lab1;    // lepton-side label, e.g. "nubar"
    const char* lab2;    // hadron-side label, e.g. "K+"
    bool vis1;           // lepton-side visible? (nu-bar is not)
    bool vis2;           // hadron-side visible?
};

// Registry of supported channels (SUSY/GUT "classic" set, proton + neutron).
inline const Channel* channel_registry(std::size_t& n) {
    static const Channel kChannels[] = {
        // proton modes
        {"pToKnu",   "p -> K+ nubar",  Nucleon::Proton, kMassNu,       kMassKCharged, "nubar", "K+",  false, true},
        {"pToEPi0",  "p -> e+ pi0",    Nucleon::Proton, kMassElectron, kMassPi0,      "e+",    "pi0", true,  true},
        {"pToMuPi0", "p -> mu+ pi0",   Nucleon::Proton, kMassMuon,     kMassPi0,      "mu+",   "pi0", true,  true},
        {"pToNuPip", "p -> nubar pi+", Nucleon::Proton, kMassNu,       kMassPiCharged,"nubar", "pi+", false, true},
        {"pToEEta",  "p -> e+ eta",    Nucleon::Proton, kMassElectron, kMassEta,      "e+",    "eta", true,  true},
        {"pToMuEta", "p -> mu+ eta",   Nucleon::Proton, kMassMuon,     kMassEta,      "mu+",   "eta", true,  true},
        {"pToEK0",   "p -> e+ K0",     Nucleon::Proton, kMassElectron, kMassK0,       "e+",    "K0",  true,  true},
        {"pToMuK0",  "p -> mu+ K0",    Nucleon::Proton, kMassMuon,     kMassK0,       "mu+",   "K0",  true,  true},
        // neutron modes
        {"nToEPim",  "n -> e+ pi-",    Nucleon::Neutron, kMassElectron, kMassPiCharged,"e+",    "pi-", true,  true},
        {"nToMuPim", "n -> mu+ pi-",   Nucleon::Neutron, kMassMuon,     kMassPiCharged,"mu+",   "pi-", true,  true},
        {"nToNuPi0", "n -> nubar pi0", Nucleon::Neutron, kMassNu,       kMassPi0,      "nubar", "pi0", false, true},
        {"nToNuEta", "n -> nubar eta", Nucleon::Neutron, kMassNu,       kMassEta,      "nubar", "eta", false, true},
        {"nToNuK0",  "n -> nubar K0",  Nucleon::Neutron, kMassNu,       kMassK0,       "nubar", "K0",  false, true},
        {"nToEKm",   "n -> e+ K-",     Nucleon::Neutron, kMassElectron, kMassKCharged, "e+",    "K-",  true,  true},
    };
    n = sizeof(kChannels) / sizeof(kChannels[0]);
    return kChannels;
}

// Parse a channel key (case-insensitive). Returns false on no match.
inline bool parse_channel(std::string s, Channel& out) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    std::size_t n = 0;
    const Channel* reg = channel_registry(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::string key = reg[i].key;
        std::transform(key.begin(), key.end(), key.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (s == key) {
            out = reg[i];
            return true;
        }
    }
    return false;
}

}  // namespace pdk

#endif  // PDK_CHANNELS_H
