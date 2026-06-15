// Generates bound-nucleon decay events N -> d1 + d2: samples the initial nucleon
// momentum from the chosen nuclear model, performs the off-shell two-body decay
// for the selected channel, and boosts the daughters to the lab.
//
// Output columns (whitespace-separated): event  nucleon_p  d1_p  d2_p  e_rem
//   event     : 1-based event index
//   nucleon_p : parent nucleon momentum magnitude [GeV/c]
//   d1_p      : lepton-side daughter momentum, lab frame [GeV/c]
//   d2_p      : hadron-side daughter momentum, lab frame [GeV/c]
//   e_rem     : removal (separation) energy of the bound nucleon [GeV]
//
// Usage: LunarPDKGenerator [options]
//   --events N      number of events (default 10000)
//   --channel KEY   decay channel (default pToKnu); --help lists all
//   --model NAME    nucleon momentum model (default polynomial)
//   --binding NAME  mean-field removal-energy model (default potential)
//   --seed S        RNG seed (default: non-deterministic)
//   --config PATH   params file for the polynomial model (default config/params.dat)
//   -h, --help      print this help (with the supported channels/models/bindings)

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <optional>
#include <string>

#include "PDKCascade.h"
#include "PDKChannels.h"
#include "PDKConfig.h"
#include "PDKDecay.h"
#include "PDKKinematics.h"
#include "PDKMomentum.h"

namespace {

// Map a channel's hadron-side label to a PDG code for the FSI cascade.
int hadron_pdg(const std::string& lab) {
    if (lab == "K+") return pdk::fsi::kPdgKPlus;
    if (lab == "K0") return pdk::fsi::kPdgK0;
    if (lab == "K-") return pdk::fsi::kPdgKMinus;
    if (lab == "pi+") return pdk::fsi::kPdgPiPlus;
    if (lab == "pi0") return pdk::fsi::kPdgPi0;
    if (lab == "pi-") return pdk::fsi::kPdgPiMinus;
    if (lab == "eta") return pdk::fsi::kPdgEta;
    return 0;
}

// Map a channel's lepton-side label to a PDG code (nu-bar is invisible -> 0).
int lepton_pdg(const std::string& lab) {
    if (lab == "e+") return -11;
    if (lab == "mu+") return -13;
    return 0;
}

void print_help(const char* prog) {
    std::cout
        << "Usage: " << prog << " [options]\n"
        << "  --events N      number of events (default 10000)\n"
        << "  --channel KEY   decay channel (default pToKnu)\n"
        << "  --model NAME    momentum model: polynomial|gfg|lfg|src|sf|\n"
        << "                  hosm|br|gauss|cfg|benhar|ankowski\n"
        << "  --binding NAME  removal-energy model: potential|constant|shell\n"
        << "                  (ignored for benhar/ankowski; the table sets E)\n"
        << "  --sf-file PATH  spectral-function grid for benhar/ankowski\n"
        << "                  (default config/sf/gsf_Ar40{P,N}.grid by nucleon)\n"
        << "  --fsi on|off    final-state interactions of the hadron (default on)\n"
        << "  --fsi-pot       apply a nucleon exit potential in the FSI cascade\n"
        << "  --decay-mesons  decay escaped pi0/eta/K0 (-> gamma, K_S/K_L, pi); FSI-on only\n"
        << "  --seed S        RNG seed (default: non-deterministic)\n"
        << "  --config PATH   polynomial-model params file (default config/params.dat)\n"
        << "  -h, --help      show this help\n\n"
        << "Supported channels:\n";
    std::size_t n = 0;
    const pdk::Channel* reg = pdk::channel_registry(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::printf("  %-10s  %s\n", reg[i].key, reg[i].pretty);
    }
}

// Pop the value following a "--flag" option, or report a missing argument.
bool take_value(int argc, char* argv[], int& i, std::string& out) {
    if (i + 1 >= argc) {
        std::cerr << "Error: missing value for " << argv[i] << "\n";
        return false;
    }
    out = argv[++i];
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    int n_events = 10000;
    pdk::Channel channel;
    pdk::parse_channel("pToKnu", channel);  // default
    pdk::MomentumModel model = pdk::MomentumModel::Polynomial;
    std::string model_key = "polynomial";  // CLI key, for the forbidden summary
    pdk::BindingModel binding = pdk::BindingModel::Potential;
    std::optional<std::uint32_t> seed;
    std::string config_path = "config/params.dat";
    std::string sf_path;     // spectral-function grid (benhar/ankowski)
    bool fsi_on = true;      // run the FSI cascade on the hadron daughter
    bool fsi_pot = false;    // apply the nucleon exit potential in the cascade
    bool decay_mesons = false;  // decay escaped pi0/eta/K0 secondaries

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i], val;
        try {
            if (arg == "-h" || arg == "--help") {
                print_help(argv[0]);
                return 0;
            } else if (arg == "--events") {
                if (!take_value(argc, argv, i, val)) return 1;
                n_events = std::stoi(val);
                if (n_events <= 0) {
                    std::cerr << "Error: --events must be positive.\n";
                    return 1;
                }
            } else if (arg == "--channel") {
                if (!take_value(argc, argv, i, val)) return 1;
                if (!pdk::parse_channel(val, channel)) {
                    std::cerr << "Error: unknown channel '" << val
                              << "'. Use --help to list channels.\n";
                    return 1;
                }
            } else if (arg == "--model") {
                if (!take_value(argc, argv, i, val)) return 1;
                if (!pdk::parse_model(val, model)) {
                    std::cerr << "Error: unknown model '" << val
                              << "'. Use polynomial|gfg|lfg|src|sf|hosm|br|gauss|cfg|benhar|ankowski.\n";
                    return 1;
                }
                model_key = val;
            } else if (arg == "--binding") {
                if (!take_value(argc, argv, i, val)) return 1;
                if (!pdk::parse_binding(val, binding)) {
                    std::cerr << "Error: unknown binding '" << val
                              << "'. Use potential|constant|shell.\n";
                    return 1;
                }
            } else if (arg == "--seed") {
                if (!take_value(argc, argv, i, val)) return 1;
                seed = static_cast<std::uint32_t>(std::stoul(val));
            } else if (arg == "--config") {
                if (!take_value(argc, argv, i, val)) return 1;
                config_path = val;
            } else if (arg == "--sf-file") {
                if (!take_value(argc, argv, i, val)) return 1;
                sf_path = val;
            } else if (arg == "--fsi") {
                if (!take_value(argc, argv, i, val)) return 1;
                if (val == "on" || val == "1" || val == "true") fsi_on = true;
                else if (val == "off" || val == "0" || val == "false") fsi_on = false;
                else {
                    std::cerr << "Error: --fsi expects on|off.\n";
                    return 1;
                }
            } else if (arg == "--fsi-pot") {
                fsi_pot = true;
            } else if (arg == "--decay-mesons") {
                decay_mesons = true;
            } else {
                std::cerr << "Error: unknown option '" << arg
                          << "'. Use --help.\n";
                return 1;
            }
        } catch (const std::exception&) {
            std::cerr << "Error: bad value for " << arg << ".\n";
            return 1;
        }
    }

    const double m_nucleon = pdk::nucleon_mass(channel.parent);
    const bool tabulated = model == pdk::MomentumModel::Benhar ||
                           model == pdk::MomentumModel::Ankowski;

    // Default the spectral-function grid by parent nucleon for tabulated models.
    if (tabulated && sf_path.empty()) {
        sf_path = channel.parent == pdk::Nucleon::Proton
                      ? "config/sf/gsf_Ar40P.grid"
                      : "config/sf/gsf_Ar40N.grid";
    }

    pdk::Config cfg = pdk::load_config(config_path);
    pdk::NuclearParams np;
    pdk::NucleonMomentumSampler sample_nucleon(model, cfg, np, binding,
                                               channel.parent, sf_path);
    std::mt19937 gen = pdk::make_rng(seed);

    std::cerr << "Channel: " << channel.pretty << " | parent: "
              << pdk::nucleon_name(channel.parent) << " | model: "
              << pdk::model_name(model);
    if (tabulated)
        std::cerr << " | sf-file: " << sf_path << " (binding ignored)";
    else
        std::cerr << " | binding: " << pdk::binding_name(binding);
    std::cerr << " | FSI: " << (fsi_on ? "on" : "off");
    if (fsi_on && fsi_pot) std::cerr << " (+exit potential)";
    if (fsi_on && decay_mesons) std::cerr << " (+meson decays)";
    std::cerr << "\n";

    // FSI cascade and the channel's hadron / lepton PDG codes.
    pdk::Cascade cascade(np, gen, fsi_pot);
    const int hpdg = hadron_pdg(channel.lab2);
    const int lpdg = lepton_pdg(channel.lab1);
    if (fsi_on) {
        std::printf("# post-FSI final state; columns: event  pdg  px  py  pz  E  outcome\n");
        std::printf("# per-event header: # event <i>: nucleon_p=<GeV> e_rem=<GeV> outcome=<...>\n");
    }

    long forbidden = 0;
    for (int i = 0; i < n_events; ++i) {
        // Draw a bound nucleon; resample if it is too deeply bound / too fast for
        // the off-shell decay N -> d1 + d2 to be kinematically allowed.
        pdk::NucleonState ns;
        do {
            ns = sample_nucleon.sample(gen);
            if (!pdk::nucleon_can_decay(m_nucleon, ns.p, ns.e_rem, channel.m1,
                                        channel.m2))
                ++forbidden;
        } while (!pdk::nucleon_can_decay(m_nucleon, ns.p, ns.e_rem, channel.m1,
                                         channel.m2));

        double ux, uy, uz;
        pdk::random_direction(gen, ux, uy, uz);

        pdk::LorentzVector d1, d2;
        pdk::decay_two_body(m_nucleon, ns.p, ux, uy, uz, ns.e_rem, channel.m1,
                            channel.m2, gen, d1, d2);

        if (!fsi_on) {
            std::printf("%5d  %14.4f  %14.4f  %14.4f  %14.4f\n", i + 1, ns.p,
                        d1.p_mag(), d2.p_mag(), ns.e_rem);
            continue;
        }

        // Run the hadron daughter through the residual-nucleus cascade. The
        // lepton daughter does not interact strongly and is reported as-is.
        double vx, vy, vz;
        cascade.sample_vertex(vx, vy, vz);
        pdk::FsiResult fr =
            cascade.propagate(pdk::Particle{hpdg, d2}, vx, vy, vz);
        const char* oc = pdk::outcome_name(fr.outcome);

        std::printf("# event %d: nucleon_p=%.4f e_rem=%.4f outcome=%s\n", i + 1,
                    ns.p, ns.e_rem, oc);
        if (channel.vis1 && lpdg != 0) {
            std::printf("%5d  %6d  %12.5f  %12.5f  %12.5f  %12.5f  %s\n", i + 1,
                        lpdg, d1.px, d1.py, d1.pz, d1.E, oc);
        }
        // Optionally decay the escaped unstable mesons (pi0/eta/K0) into their
        // observable daughters (photons, K_S/K_L, pions).
        const std::vector<pdk::Particle>& finals =
            decay_mesons ? (fr.out = pdk::decay_final_mesons(fr.out, gen))
                         : fr.out;
        for (const pdk::Particle& q : finals) {
            std::printf("%5d  %6d  %12.5f  %12.5f  %12.5f  %12.5f  %s\n", i + 1,
                        q.pdg, q.p4.px, q.p4.py, q.p4.pz, q.p4.E, oc);
        }
    }

    // Structured run summary on stderr (stdout/data is untouched): the number of
    // kinematically forbidden draws that had to be resampled, as a fraction of
    // all attempts. Parsed by macros/offshell_W.C via report/forbidden_log.txt.
    const double forbidden_frac =
        static_cast<double>(forbidden) / (forbidden + n_events);
    std::cerr << "# forbidden channel=" << channel.key
              << " model=" << model_key << " generated=" << n_events
              << " forbidden=" << forbidden << " frac=" << forbidden_frac
              << "\n";
    return 0;
}
