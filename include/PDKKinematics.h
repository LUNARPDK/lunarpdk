#ifndef PDK_KINEMATICS_H
#define PDK_KINEMATICS_H

// Relativistic kinematics for the two-body decay N -> d1 + d2 of a bound nucleon
// carrying Fermi momentum inside the nucleus. The nucleon and daughter rest
// masses live in PDKChannels.h; the bound nucleon's nuclear binding
// (momentum-dependent optical potential) is applied in PDKMomentum.h and enters
// here through the removal energy and the off-shell invariant mass W below.

#include <cmath>
#include <random>

#include "PDKChannels.h"  // Nucleon, particle masses

namespace pdk {

struct LorentzVector {
    double E, px, py, pz;

    double p_mag() const { return std::sqrt(px * px + py * py + pz * pz); }
    double mass() const {
        double m2 = E * E - (px * px + py * py + pz * pz);
        return m2 > 0.0 ? std::sqrt(m2) : 0.0;
    }
};

// Sample an isotropic unit direction.
inline void random_direction(std::mt19937& gen, double& ux, double& uy,
                             double& uz) {
    std::uniform_real_distribution<double> distCos(-1.0, 1.0);
    std::uniform_real_distribution<double> distPhi(0.0, 2.0 * M_PI);
    double cosT = distCos(gen);
    double sinT = std::sqrt(1.0 - cosT * cosT);
    double phi = distPhi(gen);
    ux = sinT * std::cos(phi);
    uy = sinT * std::sin(phi);
    uz = cosT;
}

// Active Lorentz boost of a 4-vector by velocity (bx, by, bz).
inline LorentzVector boost(const LorentzVector& v, double bx, double by,
                           double bz) {
    double b2 = bx * bx + by * by + bz * bz;
    if (b2 <= 0.0) return v;
    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double bdotp = bx * v.px + by * v.py + bz * v.pz;
    double factor = (gamma - 1.0) * bdotp / b2 + gamma * v.E;
    return LorentzVector{gamma * (v.E + bdotp), v.px + factor * bx,
                         v.py + factor * by, v.pz + factor * bz};
}

// Kallen (triangle) function lambda(a,b,c) = a^2 + b^2 + c^2 - 2(ab+bc+ca).
inline double kallen(double a, double b, double c) {
    return a * a + b * b + c * c - 2.0 * (a * b + b * c + c * a);
}

// Effective (off-shell) invariant mass squared of a bound nucleon of rest mass
// m_nucleon carrying Fermi momentum p_mag and removed with removal energy
// e_removal. In the nuclear rest frame the bound nucleon four-momentum is
// (m_nucleon - e_removal, p), so
//   W^2 = (m_nucleon - e_removal)^2 - p^2.
// The energy e_removal is carried away by the residual nucleus. W^2 can fall
// below (m1+m2)^2 (or even go negative) for deeply-bound, high-momentum
// nucleons, in which case the decay N -> d1 + d2 is kinematically forbidden.
inline double offshell_mass2(double m_nucleon, double p_mag, double e_removal) {
    double E = m_nucleon - e_removal;
    return E * E - p_mag * p_mag;
}

// Two-body decay N -> d1 + d2 is allowed only if the off-shell mass W exceeds
// the sum of the daughter masses.
inline bool nucleon_can_decay(double m_nucleon, double p_mag, double e_removal,
                              double m1, double m2) {
    double W2 = offshell_mass2(m_nucleon, p_mag, e_removal);
    double thr = m1 + m2;
    return W2 > thr * thr;
}

// Given a bound nucleon (rest mass + Fermi momentum magnitude + direction +
// removal energy), produce the lab-frame 4-vectors of the two daughters d1, d2
// (rest masses m1, m2) from the off-shell two-body decay N -> d1 + d2. The
// nucleon is off-shell with invariant mass W; its energy in the nuclear frame is
// reduced by the removal energy, which sets both the decay momentum and the
// boost. Caller must ensure nucleon_can_decay() is true.
inline void decay_two_body(double m_nucleon, double p_mag, double ux, double uy,
                           double uz, double e_removal, double m1, double m2,
                           std::mt19937& gen, LorentzVector& d1,
                           LorentzVector& d2) {
    double E = m_nucleon - e_removal;                       // off-shell energy, nuclear frame
    double W = std::sqrt(offshell_mass2(m_nucleon, p_mag, e_removal));

    // Daughter momentum in the off-shell nucleon rest frame (back-to-back).
    double p_rest = std::sqrt(kallen(W * W, m1 * m1, m2 * m2)) / (2.0 * W);
    double E1 = std::sqrt(p_rest * p_rest + m1 * m1);
    double E2 = std::sqrt(p_rest * p_rest + m2 * m2);

    double dx, dy, dz;
    random_direction(gen, dx, dy, dz);
    LorentzVector d1_rest{E1, p_rest * dx, p_rest * dy, p_rest * dz};
    LorentzVector d2_rest{E2, -p_rest * dx, -p_rest * dy, -p_rest * dz};

    // Boost along the nucleon's velocity beta = p / E into the lab frame.
    double bx = p_mag * ux / E, by = p_mag * uy / E, bz = p_mag * uz / E;
    d1 = boost(d1_rest, bx, by, bz);
    d2 = boost(d2_rest, bx, by, bz);
}

}  // namespace pdk

#endif  // PDK_KINEMATICS_H
