#ifndef PDK_CASCADE_H
#define PDK_CASCADE_H

// Semi-classical intranuclear cascade (final-state interactions) for the PDK
// generator.
//
// The hadron produced in a bound-nucleon decay (K+, K0, pi, eta, ...) is
// transported out through the residual argon nucleus before it is observed. The
// model follows the same philosophy as NuWro's `kaskada`:
//
//   1. Sample the decay vertex from the Woods-Saxon nucleon density rho(r)
//      (shared with the initial-state sampler, see PDKMomentum.h).
//   2. Step the hadron in small straight segments. At each step the local mean
//      free path is lambda = 1 / (rho * sigma), with sigma the hadron-nucleon
//      total cross section (PDKFsiXsec.h) split into proton/neutron pieces.
//   3. On interaction, draw a target nucleon with a local-Fermi-gas momentum,
//      pick a sub-channel (elastic / charge-exchange / absorption / production),
//      generate the kinematics by boosting to the projectile+target rest frame,
//      and apply Pauli blocking to every outgoing nucleon.
//   4. Cascade all secondaries (knocked-out nucleons, produced/converted pions,
//      hyperons) with the same loop until everything escapes or is absorbed.
//
// Leptons (e+, mu+, nu) do not interact strongly and never enter the cascade.
//
// This is a compact, header-only model meant for spectrum / efficiency studies,
// not a precision transport code: the angular distributions are isotropic in the
// centre of mass and the multi-body final states use a sequential (approximate)
// phase space. Four-momentum and charge are conserved exactly in every channel.

#include <cmath>
#include <random>
#include <vector>

#include "PDKFsiXsec.h"
#include "PDKKinematics.h"  // LorentzVector, boost, kallen, random_direction
#include "PDKMomentum.h"    // NuclearParams, nucleon_density, local_fermi_momentum

namespace pdk {

// A cascade particle: PDG id + lab-frame four-momentum.
struct Particle {
    int pdg;
    LorentzVector p4;
};

// Per-event classification of the PRIMARY meson's fate.
enum class FsiOutcome {
    None,      // left the nucleus without interacting
    Elastic,   // interacted but kept its identity, no extra mesons made
    Cex,       // left as a different meson (single charge exchange)
    Produced,  // survived but created extra meson(s) (inelastic production)
    Absorbed,  // did not leave the nucleus (absorbed or converted)
};

inline const char* outcome_name(FsiOutcome o) {
    switch (o) {
        case FsiOutcome::None: return "none";
        case FsiOutcome::Elastic: return "elastic";
        case FsiOutcome::Cex: return "cex";
        case FsiOutcome::Produced: return "produced";
        case FsiOutcome::Absorbed: return "absorbed";
    }
    return "unknown";
}

struct FsiResult {
    std::vector<Particle> out;       // all hadrons that escaped the nucleus
    FsiOutcome outcome = FsiOutcome::None;
};

// ---- PDG helpers -----------------------------------------------------------

inline double pdg_mass(int pdg) {
    switch (std::abs(pdg)) {
        case fsi::kPdgPiPlus: return kMassPiCharged;
        case fsi::kPdgPi0: return kMassPi0;
        case fsi::kPdgKPlus: return kMassKCharged;
        case fsi::kPdgK0: return kMassK0;
        case fsi::kPdgEta: return kMassEta;
        case fsi::kPdgProton: return kProtonMass;
        case fsi::kPdgNeutron: return kNeutronMass;
        case 3122: return 1.115683;  // Lambda
    }
    return 0.0;
}

inline int pdg_charge(int pdg) {
    switch (pdg) {
        case fsi::kPdgPiPlus: return +1;
        case fsi::kPdgPiMinus: return -1;
        case fsi::kPdgKPlus: return +1;
        case fsi::kPdgKMinus: return -1;
        case fsi::kPdgProton: return +1;
        default: return 0;  // pi0, K0, K0bar, eta, neutron, Lambda
    }
}

inline bool is_pion(int pdg) {
    return pdg == fsi::kPdgPiPlus || pdg == fsi::kPdgPi0 || pdg == fsi::kPdgPiMinus;
}
inline bool is_kaon(int pdg) {
    return pdg == fsi::kPdgKPlus || pdg == fsi::kPdgKMinus ||
           pdg == fsi::kPdgK0 || pdg == fsi::kPdgK0bar;
}
inline bool is_eta(int pdg) { return pdg == fsi::kPdgEta; }
inline bool is_nucleon(int pdg) {
    return pdg == fsi::kPdgProton || pdg == fsi::kPdgNeutron;
}

// ---- kinematic primitives --------------------------------------------------

// Decay/split a parent four-vector into two daughters of mass ma, mb, isotropic
// in the parent rest frame, returned in the lab frame. False if the parent
// invariant mass cannot accommodate the two masses.
inline bool decay_into(const LorentzVector& parent, double ma, double mb,
                       std::mt19937& gen, LorentzVector& a, LorentzVector& b) {
    double W = parent.mass();
    if (W <= ma + mb) return false;
    double lam = kallen(W * W, ma * ma, mb * mb);
    if (lam <= 0.0) return false;
    double pcm = std::sqrt(lam) / (2.0 * W);
    double ux, uy, uz;
    random_direction(gen, ux, uy, uz);
    LorentzVector acm{std::sqrt(pcm * pcm + ma * ma), pcm * ux, pcm * uy, pcm * uz};
    LorentzVector bcm{std::sqrt(pcm * pcm + mb * mb), -pcm * ux, -pcm * uy, -pcm * uz};
    double bx = parent.px / parent.E, by = parent.py / parent.E,
           bz = parent.pz / parent.E;
    a = boost(acm, bx, by, bz);
    b = boost(bcm, bx, by, bz);
    return true;
}

// Sample cos(theta*) in the CM with a forward-peaked weight exp(c*(u-1)),
// c >= 0 (c = 0 -> isotropic). Larger c favours the forward direction (u->+1).
// Used for diffractive-like quasi-elastic / charge-exchange scattering, whose
// real angular distribution is forward-peaked rather than isotropic.
inline double sample_cos_forward(double c, std::mt19937& gen) {
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    if (c < 1e-6) return 2.0 * u01(gen) - 1.0;  // isotropic limit
    // pdf ~ exp(c*u) on [-1,1]; inverse-CDF sampling.
    double r = u01(gen);
    double em = std::exp(-c), ep = std::exp(c);
    double u = std::log(em + r * (ep - em)) / c;
    return std::max(-1.0, std::min(1.0, u));
}

// Two-body scatter parent -> a (mass ma) + b (mass mb), forward-peaked in the CM
// about the incoming projectile direction. `p_proj` is the lab projectile
// four-vector and `slope` B [GeV^-2] sets the diffractive peaking through the
// momentum transfer t ~ -2 p*^2 (1 - cos theta*); slope = 0 reduces to the
// isotropic decay_into(). Returns the lab-frame daughters a, b; false below
// threshold.
inline bool decay_into_aniso(const LorentzVector& parent, double ma, double mb,
                             const LorentzVector& p_proj, double slope,
                             std::mt19937& gen, LorentzVector& a,
                             LorentzVector& b) {
    double W = parent.mass();
    if (W <= ma + mb) return false;
    double lam = kallen(W * W, ma * ma, mb * mb);
    if (lam <= 0.0) return false;
    double pcm = std::sqrt(lam) / (2.0 * W);

    // Incoming projectile direction in the parent (CM) rest frame = the axis the
    // scattering angle is measured from.
    double bx = parent.px / parent.E, by = parent.py / parent.E,
           bz = parent.pz / parent.E;
    LorentzVector proj_cm = boost(p_proj, -bx, -by, -bz);
    double pin = proj_cm.p_mag();
    double ax, ay, az;
    if (pin > 1e-12) {
        ax = proj_cm.px / pin; ay = proj_cm.py / pin; az = proj_cm.pz / pin;
    } else {
        random_direction(gen, ax, ay, az);
    }

    const double c = 2.0 * slope * pcm * pcm;          // peaking strength
    const double cth = sample_cos_forward(c, gen);
    const double sth = std::sqrt(std::max(0.0, 1.0 - cth * cth));
    std::uniform_real_distribution<double> dphi(0.0, 2.0 * M_PI);
    const double phi = dphi(gen);

    // Orthonormal frame (e1, e2) perpendicular to the axis.
    double e1x, e1y, e1z;
    if (std::fabs(ax) < 0.9) { e1x = 0.0; e1y = -az; e1z = ay; }
    else { e1x = -az; e1y = 0.0; e1z = ax; }
    double n1 = std::sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
    e1x /= n1; e1y /= n1; e1z /= n1;
    double e2x = ay * e1z - az * e1y;
    double e2y = az * e1x - ax * e1z;
    double e2z = ax * e1y - ay * e1x;

    double dx = cth * ax + sth * (std::cos(phi) * e1x + std::sin(phi) * e2x);
    double dy = cth * ay + sth * (std::cos(phi) * e1y + std::sin(phi) * e2y);
    double dz = cth * az + sth * (std::cos(phi) * e1z + std::sin(phi) * e2z);

    LorentzVector acm{std::sqrt(pcm * pcm + ma * ma), pcm * dx, pcm * dy, pcm * dz};
    LorentzVector bcm{std::sqrt(pcm * pcm + mb * mb), -pcm * dx, -pcm * dy, -pcm * dz};
    a = boost(acm, bx, by, bz);
    b = boost(bcm, bx, by, bz);
    return true;
}

// Sequential (approximate) n-body phase space. Conserves four-momentum exactly;
// the invariant-mass chain is sampled uniformly rather than with the exact
// phase-space weight, which is adequate for an FSI cascade. False if below
// threshold.
inline bool phase_space(const LorentzVector& parent,
                        const std::vector<double>& m, std::mt19937& gen,
                        std::vector<LorentzVector>& out) {
    const int n = static_cast<int>(m.size());
    out.assign(n, LorentzVector{});
    double sum = 0.0;
    for (double x : m) sum += x;
    if (parent.mass() <= sum) return false;
    if (n == 2) return decay_into(parent, m[0], m[1], gen, out[0], out[1]);

    LorentzVector cur = parent;
    double rest = sum;
    for (int i = 0; i < n - 1; ++i) {
        rest -= m[i];
        double M;
        if (i == n - 2) {
            M = m[n - 1];
        } else {
            double Mmin = rest, Mmax = cur.mass() - m[i];
            if (Mmax <= Mmin) return false;
            std::uniform_real_distribution<double> u(Mmin, Mmax);
            M = u(gen);
        }
        LorentzVector piece, remainder;
        if (!decay_into(cur, m[i], M, gen, piece, remainder)) return false;
        out[i] = piece;
        cur = remainder;
    }
    out[n - 1] = cur;
    return true;
}

// ---- charge-exchange / channel identity maps -------------------------------

inline void pion_cex(int pi, int N, int& out_pi, int& out_N) {
    if (pi == fsi::kPdgPiPlus && N == fsi::kPdgNeutron) {
        out_pi = fsi::kPdgPi0; out_N = fsi::kPdgProton;
    } else if (pi == fsi::kPdgPiMinus && N == fsi::kPdgProton) {
        out_pi = fsi::kPdgPi0; out_N = fsi::kPdgNeutron;
    } else if (pi == fsi::kPdgPi0 && N == fsi::kPdgProton) {
        out_pi = fsi::kPdgPiPlus; out_N = fsi::kPdgNeutron;
    } else if (pi == fsi::kPdgPi0 && N == fsi::kPdgNeutron) {
        out_pi = fsi::kPdgPiMinus; out_N = fsi::kPdgProton;
    } else {
        out_pi = pi; out_N = N;  // no CEX channel; falls back to elastic
    }
}

inline void kaon_cex(int k, int N, int& out_k, int& out_N) {
    if (k == fsi::kPdgKPlus && N == fsi::kPdgNeutron) {
        out_k = fsi::kPdgK0; out_N = fsi::kPdgProton;
    } else if (k == fsi::kPdgK0 && N == fsi::kPdgProton) {
        out_k = fsi::kPdgKPlus; out_N = fsi::kPdgNeutron;
    } else if (k == fsi::kPdgKMinus && N == fsi::kPdgProton) {
        out_k = fsi::kPdgK0bar; out_N = fsi::kPdgNeutron;
    } else if (k == fsi::kPdgK0bar && N == fsi::kPdgNeutron) {
        out_k = fsi::kPdgKMinus; out_N = fsi::kPdgProton;
    } else {
        out_k = k; out_N = N;  // no CEX channel; falls back to elastic
    }
}

// ---- the cascade -----------------------------------------------------------

class Cascade {
public:
    Cascade(const NuclearParams& np, std::mt19937& gen, bool apply_potential = false)
        : np_(np), gen_(gen), apply_pot_(apply_potential) {
        const int steps = 2000;
        for (int i = 0; i <= steps; ++i) {
            double r = np_.r_max * i / steps;
            r2rho_max_ = std::max(r2rho_max_, r * r * nucleon_density(r, np_));
        }
    }

    // Sample a decay vertex (x,y,z) [fm] from the r^2 rho(r) radial weight.
    void sample_vertex(double& x, double& y, double& z) {
        std::uniform_real_distribution<double> dr(0.0, np_.r_max);
        std::uniform_real_distribution<double> dy(0.0, r2rho_max_);
        double r;
        while (true) {
            r = dr(gen_);
            if (dy(gen_) < r * r * nucleon_density(r, np_)) break;
        }
        double ux, uy, uz;
        random_direction(gen_, ux, uy, uz);
        x = r * ux; y = r * uy; z = r * uz;
    }

    // Transport the primary meson from the vertex and return the full post-FSI
    // hadronic final state plus the primary's outcome classification.
    FsiResult propagate(const Particle& meson, double vx, double vy, double vz) {
        FsiResult res;
        std::vector<Track> queue;
        Fate fate;
        fate.init_pdg = meson.pdg;

        LineResult primary = transport(meson, vx, vy, vz, queue, &fate);
        if (primary.escaped) {
            res.out.push_back(primary.particle);
            if (primary.particle.pdg != fate.init_pdg)
                res.outcome = FsiOutcome::Cex;
            else if (fate.produced)
                res.outcome = FsiOutcome::Produced;
            else if (fate.interacted)
                res.outcome = FsiOutcome::Elastic;
            else
                res.outcome = FsiOutcome::None;
        } else {
            res.outcome = FsiOutcome::Absorbed;
        }

        // Drain the secondary queue (their individual fates do not relabel the
        // event, but they cascade fully and may spawn further secondaries).
        std::size_t guard = 0;
        while (!queue.empty() && guard++ < kMaxParticles) {
            Track t = queue.back();
            queue.pop_back();
            LineResult lr = transport(t.p, t.x, t.y, t.z, queue, nullptr);
            if (lr.escaped) res.out.push_back(lr.particle);
        }
        return res;
    }

private:
    struct Track {
        Particle p;
        double x, y, z;
    };
    struct LineResult {
        bool escaped;
        Particle particle;
    };
    struct Fate {
        bool interacted = false;
        bool produced = false;
        int init_pdg = 0;
    };
    enum class IAct { NoInteraction, Continued, Removed };

    static constexpr double kStep = 0.05;          // propagation step [fm]
    static constexpr int kMaxSteps = 200000;       // per-particle safety cap
    static constexpr std::size_t kMaxParticles = 400;  // total cascade safety cap

    // Diffractive forward-peaking slopes B [GeV^-2] for quasi-elastic / CEX
    // two-body scattering (decay_into_aniso); low-energy scattering stays nearly
    // isotropic since the peaking ~ exp(2 B p*^2 (cos-1)) vanishes as p* -> 0.
    static constexpr double kSlopePi = 4.0;        // pi  N
    static constexpr double kSlopeK = 3.0;         // K   N
    static constexpr double kSlopeN = 5.0;         // N   N

    double rng() { return uni_(gen_); }

    // Total hadron-nucleon cross section [mb] of `p` on a nucleon `Npdg`, target
    // at rest (Fermi motion enters only the interaction kinematics).
    double total_xsec(const Particle& p, int Npdg) const {
        if (is_pion(p.pdg)) {
            double T = p.p4.E - pdg_mass(p.pdg);
            return fsi::pion_nucleon(T, p.pdg, Npdg).total;
        }
        if (is_kaon(p.pdg))
            return fsi::kaon_nucleon(p.p4.p_mag(), p.pdg, Npdg).total;
        if (is_eta(p.pdg)) return fsi::eta_nucleon(p.p4.p_mag()).total;
        if (is_nucleon(p.pdg)) {
            double T = p.p4.E - pdg_mass(p.pdg);
            return fsi::nucleon_nucleon(T, p.pdg, Npdg).total;
        }
        return 0.0;  // hyperons, leptons: no FSI
    }

    // Sample a target nucleon four-vector from the local Fermi sphere at radius r.
    LorentzVector sample_target(int Npdg, double r) {
        double frac = (Npdg == fsi::kPdgProton) ? np_.frac_p : np_.frac_n;
        double kF = local_fermi_momentum(r, np_, frac);
        double pmag = kF * std::cbrt(rng());
        double ux, uy, uz;
        random_direction(gen_, ux, uy, uz);
        double m = pdg_mass(Npdg);
        return LorentzVector{std::sqrt(pmag * pmag + m * m), pmag * ux, pmag * uy,
                             pmag * uz};
    }

    bool pauli_blocked(const LorentzVector& v, double kF) const {
        return v.p_mag() < kF;
    }

    // Apply an optional constant nuclear potential step to an escaping nucleon:
    // it loses |V| of kinetic energy at the surface. Returns false if it is
    // trapped (kinetic energy would go negative).
    bool apply_exit_potential(Particle& p) const {
        if (!apply_pot_ || !is_nucleon(p.pdg)) return true;
        double m = pdg_mass(p.pdg);
        double T = p.p4.E - m - kPotDepth;
        if (T <= 0.0) return false;
        double pnew = std::sqrt(T * (T + 2.0 * m));
        double pmag = p.p4.p_mag();
        if (pmag <= 0.0) return false;
        double s = pnew / pmag;
        p.p4 = LorentzVector{T + m, p.p4.px * s, p.p4.py * s, p.p4.pz * s};
        return true;
    }
    static constexpr double kPotDepth = 0.040;  // nucleon exit potential [GeV]

    // Transport a single particle line until it escapes (returned) or is removed
    // (escaped=false). Secondaries are pushed onto `queue`; `fate` (if non-null,
    // i.e. the primary line) records whether it interacted / produced mesons.
    LineResult transport(Particle p, double x, double y, double z,
                         std::vector<Track>& queue, Fate* fate) {
        for (int step = 0; step < kMaxSteps; ++step) {
            double pm = p.p4.p_mag();
            if (pm <= 0.0) return {true, p};  // at rest: cannot propagate
            double ux = p.p4.px / pm, uy = p.p4.py / pm, uz = p.p4.pz / pm;
            x += ux * kStep; y += uy * kStep; z += uz * kStep;
            double r = std::sqrt(x * x + y * y + z * z);
            if (r >= np_.r_max) {  // reached the surface
                if (apply_exit_potential(p)) return {true, p};
                return {false, p};  // trapped by the potential
            }

            double rho = nucleon_density(r, np_);
            if (rho <= 0.0) continue;
            double xs_p = total_xsec(p, fsi::kPdgProton);
            double xs_n = total_xsec(p, fsi::kPdgNeutron);
            double wp = rho * np_.frac_p * xs_p;
            double wn = rho * np_.frac_n * xs_n;
            double macro = (wp + wn) * fsi::kMbToFm2;  // fm^-1
            if (macro <= 0.0) continue;
            if (rng() >= 1.0 - std::exp(-kStep * macro)) continue;

            int Npdg = (rng() < wp / (wp + wn)) ? fsi::kPdgProton : fsi::kPdgNeutron;
            IAct act = interact(p, Npdg, r, x, y, z, queue, fate);
            if (act == IAct::Removed) return {false, p};
            // Continued or NoInteraction: keep propagating with (possibly) updated p.
        }
        return {true, p};  // safety: emit if the step cap is hit
    }

    // Resolve one interaction of projectile `p` with a nucleon `Npdg` at radius r.
    IAct interact(Particle& p, int Npdg, double r, double x, double y, double z,
                  std::vector<Track>& queue, Fate* fate) {
        double frac = (Npdg == fsi::kPdgProton) ? np_.frac_p : np_.frac_n;
        double kF = local_fermi_momentum(r, np_, frac);
        LorentzVector tgt = sample_target(Npdg, r);

        if (is_pion(p.pdg)) return interact_pion(p, Npdg, tgt, kF, x, y, z, queue, fate);
        if (is_kaon(p.pdg)) return interact_kaon(p, Npdg, tgt, kF, x, y, z, queue, fate);
        if (is_eta(p.pdg)) return interact_eta(p, Npdg, tgt, kF, x, y, z, queue, fate);
        if (is_nucleon(p.pdg)) return interact_nucleon(p, Npdg, tgt, kF, x, y, z, queue, fate);
        return IAct::NoInteraction;
    }

    void push(std::vector<Track>& q, int pdg, const LorentzVector& p4, double x,
              double y, double z) {
        if (q.size() < kMaxParticles) q.push_back({Particle{pdg, p4}, x, y, z});
    }

    LorentzVector total4(const LorentzVector& a, const LorentzVector& b) const {
        return LorentzVector{a.E + b.E, a.px + b.px, a.py + b.py, a.pz + b.pz};
    }

    int pick_channel(double f0, double f1, double f2) {
        // Returns 0,1,2,3 by cumulative fractions (f0+f1+f2+rest == 1).
        double u = rng();
        if (u < f0) return 0;
        if (u < f0 + f1) return 1;
        if (u < f0 + f1 + f2) return 2;
        return 3;
    }

    IAct interact_pion(Particle& p, int Npdg, const LorentzVector& tgt, double kF,
                       double x, double y, double z, std::vector<Track>& queue,
                       Fate* fate) {
        double T = p.p4.E - pdg_mass(p.pdg);
        fsi::PionXsec xs = fsi::pion_nucleon(T, p.pdg, Npdg);
        int ch = pick_channel(xs.f_elastic, xs.f_cex, xs.f_abs);

        if (ch == 2) {  // absorption: pi N N -> N N
            return absorb_pion(p, Npdg, tgt, kF, x, y, z, queue, fate);
        }
        if (ch == 3) {  // production: pi N -> pi pi N
            LorentzVector W = total4(p.p4, tgt);
            std::vector<LorentzVector> o;
            if (phase_space(W, {pdg_mass(p.pdg), kMassPi0, pdg_mass(Npdg)}, gen_, o)) {
                if (pauli_blocked(o[2], kF)) return IAct::NoInteraction;
                p.p4 = o[0];
                push(queue, fsi::kPdgPi0, o[1], x, y, z);
                push(queue, Npdg, o[2], x, y, z);
                if (fate) { fate->interacted = true; fate->produced = true; }
                return IAct::Continued;
            }
            ch = 0;  // below threshold: fall back to elastic
        }

        int out_pi = p.pdg, out_N = Npdg;
        if (ch == 1) pion_cex(p.pdg, Npdg, out_pi, out_N);  // charge exchange

        LorentzVector W = total4(p.p4, tgt), meson, nucleon;
        if (!decay_into_aniso(W, pdg_mass(out_pi), pdg_mass(out_N), p.p4,
                              kSlopePi, gen_, meson, nucleon))
            return IAct::NoInteraction;
        if (pauli_blocked(nucleon, kF)) return IAct::NoInteraction;
        p.pdg = out_pi;
        p.p4 = meson;
        push(queue, out_N, nucleon, x, y, z);
        if (fate) fate->interacted = true;
        return IAct::Continued;
    }

    IAct absorb_pion(Particle& p, int Npdg, const LorentzVector& tgt, double kF,
                     double x, double y, double z, std::vector<Track>& queue,
                     Fate* fate) {
        int qpi = pdg_charge(p.pdg), q1 = pdg_charge(Npdg);
        // Choose the second nucleon so the two outgoing nucleons can carry the
        // total charge (must be 0, 1 or 2).
        int cand[2] = {fsi::kPdgProton, fsi::kPdgNeutron};
        double w[2] = {np_.frac_p, np_.frac_n};
        int N2 = 0;
        double wsum = 0.0;
        for (int i = 0; i < 2; ++i) {
            int Q = qpi + q1 + pdg_charge(cand[i]);
            if (Q >= 0 && Q <= 2) wsum += w[i];
        }
        if (wsum <= 0.0) return IAct::NoInteraction;  // no charge-allowed pairing
        double u = rng() * wsum, acc = 0.0;
        for (int i = 0; i < 2; ++i) {
            int Q = qpi + q1 + pdg_charge(cand[i]);
            if (Q >= 0 && Q <= 2) { acc += w[i]; if (u <= acc) { N2 = cand[i]; break; } }
        }
        int Qtot = qpi + q1 + pdg_charge(N2);
        int Na = (Qtot >= 1) ? fsi::kPdgProton : fsi::kPdgNeutron;
        int Nb = (Qtot == 2) ? fsi::kPdgProton : fsi::kPdgNeutron;

        LorentzVector tgt2 = sample_target(N2, std::sqrt(x * x + y * y + z * z));
        LorentzVector W = total4(total4(p.p4, tgt), tgt2);
        LorentzVector a, b;
        if (!decay_into(W, pdg_mass(Na), pdg_mass(Nb), gen_, a, b))
            return IAct::NoInteraction;
        if (pauli_blocked(a, kF) || pauli_blocked(b, kF)) return IAct::NoInteraction;
        push(queue, Na, a, x, y, z);
        push(queue, Nb, b, x, y, z);
        if (fate) fate->interacted = true;
        return IAct::Removed;  // pion absorbed
    }

    IAct interact_kaon(Particle& p, int Npdg, const LorentzVector& tgt, double kF,
                       double x, double y, double z, std::vector<Track>& queue,
                       Fate* fate) {
        fsi::KaonXsec xs = fsi::kaon_nucleon(p.p4.p_mag(), p.pdg, Npdg);
        int ch = pick_channel(xs.f_elastic, xs.f_cex, xs.f_abs);

        if (ch == 2) {  // Kbar N -> Lambda + pi (strangeness-exchange absorption)
            int Q = pdg_charge(p.pdg) + pdg_charge(Npdg);
            int out_pi = (Q == 1) ? fsi::kPdgPiPlus
                                  : (Q == -1) ? fsi::kPdgPiMinus : fsi::kPdgPi0;
            LorentzVector W = total4(p.p4, tgt), lam, pion;
            if (!decay_into(W, pdg_mass(3122), pdg_mass(out_pi), gen_, lam, pion))
                return IAct::NoInteraction;
            push(queue, 3122, lam, x, y, z);
            push(queue, out_pi, pion, x, y, z);
            if (fate) fate->interacted = true;
            return IAct::Removed;  // kaon absorbed
        }

        int out_k = p.pdg, out_N = Npdg;
        if (ch == 1) kaon_cex(p.pdg, Npdg, out_k, out_N);  // charge exchange

        LorentzVector W = total4(p.p4, tgt), meson, nucleon;
        if (!decay_into_aniso(W, pdg_mass(out_k), pdg_mass(out_N), p.p4, kSlopeK,
                              gen_, meson, nucleon))
            return IAct::NoInteraction;
        if (pauli_blocked(nucleon, kF)) return IAct::NoInteraction;
        p.pdg = out_k;
        p.p4 = meson;
        push(queue, out_N, nucleon, x, y, z);
        if (fate) fate->interacted = true;
        return IAct::Continued;
    }

    IAct interact_eta(Particle& p, int Npdg, const LorentzVector& tgt, double kF,
                      double x, double y, double z, std::vector<Track>& queue,
                      Fate* fate) {
        fsi::EtaXsec xs = fsi::eta_nucleon(p.p4.p_mag());
        bool convert = rng() < xs.f_conv;
        // eta N -> pi0 N (conversion) keeps the nucleon species and total charge.
        double m_meson = convert ? kMassPi0 : pdg_mass(p.pdg);
        LorentzVector W = total4(p.p4, tgt), meson, nucleon;
        if (!decay_into(W, m_meson, pdg_mass(Npdg), gen_, meson, nucleon))
            return IAct::NoInteraction;
        if (pauli_blocked(nucleon, kF)) return IAct::NoInteraction;
        push(queue, Npdg, nucleon, x, y, z);
        if (fate) fate->interacted = true;
        if (convert) {  // primary eta is gone, replaced by a pi0 secondary
            push(queue, fsi::kPdgPi0, meson, x, y, z);
            return IAct::Removed;
        }
        p.p4 = meson;  // elastic
        return IAct::Continued;
    }

    IAct interact_nucleon(Particle& p, int Npdg, const LorentzVector& tgt,
                          double kF, double x, double y, double z,
                          std::vector<Track>& queue, Fate* fate) {
        double T = p.p4.E - pdg_mass(p.pdg);
        fsi::NNXsec xs = fsi::nucleon_nucleon(T, p.pdg, Npdg);
        bool produce = rng() < xs.f_prod;

        if (produce) {  // N N -> N N pi0
            LorentzVector W = total4(p.p4, tgt);
            std::vector<LorentzVector> o;
            if (phase_space(W, {pdg_mass(p.pdg), pdg_mass(Npdg), kMassPi0}, gen_, o)) {
                if (pauli_blocked(o[0], kF) || pauli_blocked(o[1], kF))
                    return IAct::NoInteraction;
                p.p4 = o[0];
                push(queue, Npdg, o[1], x, y, z);
                push(queue, fsi::kPdgPi0, o[2], x, y, z);
                if (fate) { fate->interacted = true; fate->produced = true; }
                return IAct::Continued;
            }
        }
        // elastic N N
        LorentzVector W = total4(p.p4, tgt), n1, n2;
        if (!decay_into_aniso(W, pdg_mass(p.pdg), pdg_mass(Npdg), p.p4, kSlopeN,
                              gen_, n1, n2))
            return IAct::NoInteraction;
        if (pauli_blocked(n1, kF) || pauli_blocked(n2, kF)) return IAct::NoInteraction;
        p.p4 = n1;
        push(queue, Npdg, n2, x, y, z);
        if (fate) fate->interacted = true;
        return IAct::Continued;
    }

    NuclearParams np_;
    std::mt19937& gen_;
    bool apply_pot_;
    double r2rho_max_ = 0.0;
    std::uniform_real_distribution<double> uni_{0.0, 1.0};
};

}  // namespace pdk

#endif  // PDK_CASCADE_H
