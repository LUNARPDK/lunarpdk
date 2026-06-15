#ifndef PDK_FSI_XSEC_H
#define PDK_FSI_XSEC_H

// Hadron-nucleon cross sections for the PDK final-state-interaction cascade
// (PDKCascade.h).
//
// The pion-nucleon and nucleon-nucleon cross sections are the Metropolis tables
// used by NuWro's intranuclear cascade (kaskada):
//   pi N : N. Metropolis et al., Phys. Rev. 110 (1958) 204-219, as coded in
//          NuWro src/Interaction.cc PiData::setMetropolis() (xsec model 0).
//   N N  : the same reference, NuWro data/input/kaskada_NN_xsec_0.dat and
//          kaskada_NN_inel_0.dat (xsec model 0).
// The numeric tables are embedded verbatim below (the raw NuWro files are kept
// for provenance under config/fsi/) so the PDK generator stays header-only and
// needs no external data files at run time.
//
// The K+/K0 (S = +1, "transparent"), Kbar (S = -1, strongly absorbed) and eta
// (through the N*(1535)) cross sections have no Metropolis table. The kaon
// cross sections are short interpolation tables of the measured K-nucleon totals
// (Dover & Walker, Phys. Rept. 89 (1982) 1; Friedman & Gal, Phys. Rept. 452
// (2007) 89; PDG), in the same tabulated spirit as the Metropolis data above;
// the eta cross section is an N*(1535) Breit-Wigner. They are documented and
// tunable at each call site.
//
// Conventions:
//   * cross sections returned in millibarn [mb]; 1 mb = 0.1 fm^2 (kMbToFm2).
//   * the kinetic-energy / momentum argument is that of the projectile in the
//     REST FRAME of the struck nucleon (the cascade boosts before calling).
//   * channel fractions returned alongside the total sum to 1.

#include <cmath>

#include "PDKChannels.h"  // PDG-free; only used for shared constants if needed

namespace pdk {
namespace fsi {

constexpr double kMbToFm2 = 0.1;  // 1 mb = 0.1 fm^2

// PDG codes used across the cascade.
enum Pdg {
    kPdgPiPlus = 211,
    kPdgPi0 = 111,
    kPdgPiMinus = -211,
    kPdgKPlus = 321,
    kPdgK0 = 311,
    kPdgK0bar = -311,
    kPdgKMinus = -321,
    kPdgEta = 221,
    kPdgProton = 2212,
    kPdgNeutron = 2112,
};

// Linear interpolation on a monotonically increasing grid xs[0..n-1] with flat
// extrapolation past either end (the projectile energy can fall below the first
// tabulated point, where the table is least reliable; clamping is the standard
// cascade choice and avoids a spurious zero cross section).
inline double interp(const double* xs, const double* ys, int n, double x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    int hi = 1;
    while (hi < n && xs[hi] < x) ++hi;
    double t = (x - xs[hi - 1]) / (xs[hi] - xs[hi - 1]);
    return ys[hi - 1] + t * (ys[hi] - ys[hi - 1]);
}

// ----------------------------- pion - nucleon ------------------------------

struct PionXsec {
    double total = 0.0;     // mb
    double f_elastic = 1.0; // quasi-elastic (no charge change)
    double f_cex = 0.0;     // single charge exchange
    double f_abs = 0.0;     // true absorption (pi N N -> N N), pion removed
    double f_prod = 0.0;    // inelastic pion production (extra pion)
};

// pi N cross section from the Metropolis table. T_pi is the pion kinetic energy
// in the nucleon rest frame [GeV]. The isospin assignment follows Metropolis:
//   "ii" (resonant I=3/2): pi+ p  or  pi- n
//   "ij"                  : pi+ n  or  pi- p
//   pi0 sees the average (sii+sij)/2 with its own production/CEX fractions.
inline PionXsec pion_nucleon(double T_pi, int pion_pdg, int nucleon_pdg) {
    // Metropolis et al. (1958); pion kinetic energy grid [MeV].
    static const double E[]    = {0, 49, 85, 128, 184, 250, 350, 540, 1300, 1e10};
    static const double sii[]  = {16, 16, 50, 114, 200, 110, 51, 20, 30, 30};
    static const double sij[]  = {15, 15, 21, 43, 66, 44, 23, 22, 30, 30};
    static const double sabs[] = {20, 20, 32, 45, 36, 18, 0, 0, 0, 0};
    static const double fxii[] = {0, 0, 0, 0, 0.03, 0.06, 0.16, 0.30, 0.88, 0.88};
    static const double fxij[] = {0.45, 0.45, 0.57, 0.62, 0.64, 0.62, 0.56, 0.58, 0.94, 0.94};
    static const double fx0[]  = {0.42, 0.42, 0.36, 0.36, 0.37, 0.40, 0.50, 0.59, 0.94, 0.94};
    static const double fceij[]= {0.80, 1.00, 1.00, 1.00, 0.95, 0.89, 0.72, 0.51, 0.06, 0.06};
    static const double fce0[] = {0.80, 1.00, 1.00, 1.00, 0.90, 0.84, 0.67, 0.50, 0.05, 0.05};
    const int n = sizeof(E) / sizeof(E[0]);

    const double Tmev = T_pi * 1.0e3;
    const double s_ii = interp(E, sii, n, Tmev);
    const double s_ij = interp(E, sij, n, Tmev);
    const double s_abs = interp(E, sabs, n, Tmev);

    double s_scat, f_x, f_ce;
    if (pion_pdg == kPdgPi0) {
        s_scat = 0.5 * (s_ii + s_ij);
        f_x = interp(E, fx0, n, Tmev);
        f_ce = interp(E, fce0, n, Tmev);
    } else {
        const bool resonant =
            (pion_pdg == kPdgPiPlus && nucleon_pdg == kPdgProton) ||
            (pion_pdg == kPdgPiMinus && nucleon_pdg == kPdgNeutron);
        if (resonant) {
            s_scat = s_ii;
            f_x = interp(E, fxii, n, Tmev);
            f_ce = 0.0;  // fceii == 0 in the table (CEX forbidden for pi+ p / pi- n)
        } else {
            s_scat = s_ij;
            f_x = interp(E, fxij, n, Tmev);
            f_ce = interp(E, fceij, n, Tmev);
        }
    }

    PionXsec out;
    out.total = s_scat + s_abs;
    if (out.total <= 0.0) {
        out.f_elastic = 1.0;
        return out;
    }
    // Within the scattering part: a fraction f_x is inelastic (pion production);
    // the rest splits into single charge exchange (f_ce) and quasi-elastic.
    const double w_prod = s_scat * f_x;
    const double w_rest = s_scat * (1.0 - f_x);
    const double w_cex = w_rest * f_ce;
    const double w_el = w_rest * (1.0 - f_ce);
    out.f_prod = w_prod / out.total;
    out.f_cex = w_cex / out.total;
    out.f_elastic = w_el / out.total;
    out.f_abs = s_abs / out.total;
    return out;
}

// --------------------------- nucleon - nucleon -----------------------------

struct NNXsec {
    double total = 0.0;      // mb
    double f_elastic = 1.0;
    double f_prod = 0.0;     // inelastic (single-pion production), extra pion
};

// N N cross section, NuWro kaskada_NN_{xsec,inel}_0.dat (Metropolis). T_N is the
// projectile nucleon kinetic energy in the target rest frame [GeV]. "ii" is the
// like-pair (pp or nn), "ij" the unlike pair (pn). Below the first tabulated
// point (335 MeV) the values are clamped (flat) rather than extrapolated to the
// table's (0,0) sentinel, which would zero out low-energy rescattering.
inline NNXsec nucleon_nucleon(double T_N, int n1_pdg, int n2_pdg) {
    static const double E[]    = {335, 410, 510, 660, 840, 1160, 1780, 3900, 1e10};
    static const double xs_ii[]= {24.5, 26.4, 30.4, 41.2, 47.2, 48.0, 44.2, 41.0, 41.0};
    static const double xs_ij[]= {33.0, 34.0, 35.1, 36.5, 37.9, 40.2, 42.7, 42.0, 42.0};
    static const double in_ii[]= {0.07, 0.20, 0.31, 0.43, 0.58, 0.65, 0.69, 0.69, 0.69};
    static const double in_ij[]= {0.04, 0.07, 0.15, 0.27, 0.37, 0.36, 0.35, 0.35, 0.35};
    const int n = sizeof(E) / sizeof(E[0]);

    const double Tmev = T_N * 1.0e3;
    const bool like = (n1_pdg == n2_pdg);
    NNXsec out;
    out.total = interp(E, like ? xs_ii : xs_ij, n, Tmev);
    out.f_prod = interp(E, like ? in_ii : in_ij, n, Tmev);
    out.f_elastic = 1.0 - out.f_prod;
    return out;
}

// ------------------------------ kaon - nucleon -----------------------------

struct KaonXsec {
    double total = 0.0;     // mb
    double f_elastic = 1.0;
    double f_cex = 0.0;     // K+ n <-> K0 p (or Kbar charge exchange)
    double f_abs = 0.0;     // Kbar N -> Y pi (strangeness-exchange absorption)
};

// K-meson - nucleon cross section. p_lab is the kaon momentum in the nucleon
// rest frame [GeV/c]. The totals are short interpolation tables of the measured
// K-nucleon cross sections (Dover & Walker 1982; Friedman & Gal 2007; PDG):
//   * K+/K0 component (S=+1): nearly transparent, sigma ~ 10-18 mb with no
//     absorption (no S=+1 baryon resonance); the only inelasticity is charge
//     exchange K+ n -> K0 p (and K0 p -> K+ n). K+ n (I=0+I=1) is larger than
//     the pure-I=1 K+ p.
//   * Kbar component (S=-1): strong, sigma ~ 30-100 mb rising toward threshold,
//     with dominant strangeness-exchange absorption Kbar N -> Lambda/Sigma + pi
//     whose fraction also rises toward threshold.
// A neutral kaon (K0/K0bar, |pdg| = 311) is treated as a 50/50 incoherent mix
// of the S=+1 and S=-1 behaviours, the standard cascade approximation.
inline KaonXsec kaon_nucleon(double p_lab, int kaon_pdg, int nucleon_pdg) {
    // S=+1 totals [mb] vs kaon momentum [GeV/c] (K+ p pure I=1; K+ n I=0+I=1).
    // K+ p is the textbook "small and flat" cross section, ~10-12 mb essentially
    // from threshold (s-wave; Dover & Walker 1982; PDG), so the table plateaus
    // immediately rather than turning on from zero.
    static const double pK[]   = {0.0, 0.1, 0.2, 0.4, 0.6, 1.0, 1.5, 1e9};
    static const double sKpp[] = {10.0, 10.5, 11.0, 11.5, 12.0, 12.0, 13.0, 13.0};
    static const double sKpn[] = {6.0,  7.0,  9.0,  13.0, 16.0, 18.0, 18.0, 18.0};
    // S=-1 totals [mb] vs momentum (K- p, K- n): strong, falling from threshold.
    static const double pKb[]  = {0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 1e9};
    static const double sKmp[] = {100.0, 95.0, 70.0, 48.0, 42.0, 40.0, 38.0, 33.0, 33.0};
    static const double sKmn[] = {60.0, 55.0, 45.0, 38.0, 35.0, 33.0, 32.0, 30.0, 30.0};
    // Kbar absorption fraction (Kbar N -> Y pi) vs momentum: large near threshold.
    static const double fap[]  = {0.0, 0.2, 0.4, 0.6, 1.0, 1e9};
    static const double fav[]  = {0.65, 0.62, 0.58, 0.52, 0.45, 0.45};
    const int nK = sizeof(pK) / sizeof(pK[0]);
    const int nKb = sizeof(pKb) / sizeof(pKb[0]);
    const int nfa = sizeof(fap) / sizeof(fap[0]);

    // S = +1 piece (K+ - like), on a proton or a neutron target.
    auto kplus_like = [&](KaonXsec& k, bool on_proton) {
        k.total = interp(pK, on_proton ? sKpp : sKpn, nK, p_lab);
        // Charge exchange only off the unlike nucleon: K+ n -> K0 p, K0 p -> K+ n.
        const bool cex_allowed =
            (kaon_pdg == kPdgKPlus && nucleon_pdg == kPdgNeutron) ||
            (kaon_pdg == kPdgK0 && nucleon_pdg == kPdgProton);
        k.f_cex = cex_allowed ? 0.10 : 0.0;
        k.f_abs = 0.0;
        k.f_elastic = 1.0 - k.f_cex;
    };
    // S = -1 piece (Kbar - like): strong, absorption-dominated, falling with p.
    auto kbar_like = [&](KaonXsec& k, bool on_proton) {
        k.total = interp(pKb, on_proton ? sKmp : sKmn, nKb, p_lab);
        k.f_abs = interp(fap, fav, nfa, p_lab);  // Kbar N -> Y pi
        k.f_cex = 0.10;                          // K- p <-> K0bar n
        k.f_elastic = 1.0 - k.f_abs - k.f_cex;
        if (k.f_elastic < 0.0) k.f_elastic = 0.0;
    };

    KaonXsec out;
    const bool on_proton = (nucleon_pdg == kPdgProton);
    if (kaon_pdg == kPdgKPlus) {
        kplus_like(out, on_proton);
    } else if (kaon_pdg == kPdgKMinus) {
        kbar_like(out, on_proton);
    } else {  // neutral kaon: 50/50 incoherent K0 (S=+1) / K0bar (S=-1) mix
        KaonXsec a, b;
        kplus_like(a, on_proton);
        kbar_like(b, on_proton);
        out.total = 0.5 * (a.total + b.total);
        if (out.total <= 0.0) { out.f_elastic = 1.0; return out; }
        out.f_elastic =
            0.5 * (a.total * a.f_elastic + b.total * b.f_elastic) / out.total;
        out.f_cex = 0.5 * (a.total * a.f_cex + b.total * b.f_cex) / out.total;
        out.f_abs = 0.5 * (a.total * a.f_abs + b.total * b.f_abs) / out.total;
    }
    return out;
}

// ------------------------------- eta - nucleon -----------------------------

struct EtaXsec {
    double total = 0.0;     // mb
    double f_elastic = 1.0;
    double f_conv = 0.0;    // eta N -> pi N (eta converted to a pion), eta removed
};

// eta - nucleon cross section, dominated by the s-wave N*(1535) resonance, which
// couples strongly to both eta N and pi N. p_lab is the eta momentum in the
// nucleon rest frame [GeV/c]. Modelled as a Breit-Wigner in the eta-N invariant
// mass W (M_R = 1535 MeV, Gamma = 150 MeV) on a small non-resonant background;
// the eta production threshold sits right under the resonance, so the cross
// section is already large near threshold. The inelastic eta N -> pi N
// conversion (the eta is absorbed and replaced by a single pion) takes most of
// the cross section.
inline EtaXsec eta_nucleon(double p_lab) {
    constexpr double mR = 1.535, gR = 0.150;  // N*(1535) mass / width [GeV]
    constexpr double sigma_peak = 25.0;       // on-resonance eta-N total [mb]
    constexpr double sigma_bg = 4.0;          // non-resonant background [mb]
    const double mEta = kMassEta, mN = kProtonMass;
    const double Eeta = std::sqrt(p_lab * p_lab + mEta * mEta);
    const double W = std::sqrt(mEta * mEta + mN * mN + 2.0 * mN * Eeta);
    const double bw =
        (gR * gR / 4.0) / ((W - mR) * (W - mR) + gR * gR / 4.0);
    EtaXsec out;
    out.total = sigma_peak * bw + sigma_bg;
    out.f_conv = 0.65;  // eta N -> pi N (via the N*(1535))
    out.f_elastic = 1.0 - out.f_conv;
    return out;
}

}  // namespace fsi
}  // namespace pdk

#endif  // PDK_FSI_XSEC_H
