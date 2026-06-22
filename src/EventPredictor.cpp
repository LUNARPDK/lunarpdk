// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Expected number of nucleon-decay events in DUNE, for every channel in the
// generator, folding in the Monte-Carlo final-state-interaction (FSI) survival
// and the nuclear-model spread of the signal-containment efficiency.
//
//   N = N_nuc * (1 - exp(-T/tau)) * eps
//
// where N_nuc is the number of target protons (proton modes) or neutrons
// (neutron modes) in the fiducial mass, T the exposure, tau the assumed
// lifetime, and eps the detection efficiency. We take tau at each channel's
// CURRENT Super-Kamiokande 90% C.L. lower limit, so N is the MAXIMUM number of
// events consistent with present data ("discovery-reach" framing).
//
// Efficiency, avoiding a double-count of FSI: the DUNE TDR signal efficiencies
// for the two documented modes (p->K+ nu, p->e+ pi0) are evaluated on a GENIE
// sample that ALREADY includes Fermi motion and FSI -- e.g. the p->e+ pi0 number
// is "limited by inelastic intra-nuclear scattering". For those modes we
// therefore use eps = eps_det directly (fsi_inclusive = true) and do NOT re-apply
// our cascade survival eps_FSI. For every other mode, eps_det is an FSI-exclusive
// reconstruction assumption, so eps = eps_det * eps_FSI applies the cascade
// survival once. The nuclear model enters as a relative band on N from the spread
// of the signal-window (containment) efficiency across the ten momentum models
// (window_summary.txt). N_off is the FSI-unfolded count (= N_cen / eps_FSI): for
// the factorized modes it is rate*eps_det, and for the FSI-inclusive modes it
// unfolds eps_FSI from eps_det so the cascade suppression stays visible (rather
// than being hidden inside eps_det). The gap N_off - N_cen is the FSI effect.
//
// Detector assumptions (2 x 10 kt DUNE LAr modules, TDR-nominal exposure):
//   fiducial = 20 kt, exposure = 20 yr  -> 400 kt.yr.
//   Ar-40: 2.71e32 protons/kt (Z=18), 3.31e32 neutrons/kt (N=22).
//
// Reads report/{fsi_summary,window_summary}.txt (written by the FSI and window
// macros). Writes report/event_predictions.txt and prints a summary table.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "PDKChannels.h"

namespace {

// --- detector / exposure constants ----------------------------------------
constexpr double kProtonsPerKt = 2.71e32;   // Ar-40, Z = 18
constexpr double kNeutronsPerKt = 3.31e32;  // Ar-40, N = 22
constexpr double kFiducialKt = 20.0;        // 2 x 10 kt modules
constexpr double kExposureYr = 20.0;        // -> 400 kt.yr (TDR nominal)

// Per-channel current 90% C.L. lower limit on tau/B [years] and the detector
// efficiency eps_det. Limits are the PDG-2024-compiled Super-Kamiokande results.
// For the two DUNE-documented modes (p->K+ nu ~30%, p->e+ pi0 ~40%) eps_det is
// the DUNE TDR signal efficiency, which is FSI-inclusive (fsi_inclusive = true,
// so eps_FSI is not re-applied). Elsewhere eps_det is a flat 0.30 FSI-exclusive
// reconstruction assumption and the cascade survival eps_FSI is applied.
// A limit of 0 means "no dedicated experimental limit".
struct ChannelInput {
   double tau_limit_yr;
   double eps_det;
   bool fsi_inclusive = false;  // eps_det already includes FSI (DUNE TDR modes)
};
const std::map<std::string, ChannelInput> kInputs = {
    // proton modes
    {"pToKnu", {5.9e33, 0.30, true}},   // SuperK2014   (DUNE TDR eff ~30%, FSI-incl.)
    {"pToEPi0", {2.4e34, 0.40, true}},  // SuperK2020   (DUNE TDR eff ~40%, FSI-incl.)
    {"pToMuPi0", {1.6e34, 0.35}}, // SuperK2020
    {"pToNuPip", {3.9e32, 0.30}}, // SuperK2014nupi
    {"pToEEta", {1.4e34, 0.30}},  // SuperK2024eta
    {"pToMuEta", {7.3e33, 0.30}}, // SuperK2024eta
    {"pToEK0", {1.0e33, 0.30}},   // SuperK2005susy
    {"pToMuK0", {3.6e33, 0.30}},  // SuperK2022mK0
    // neutron modes
    {"nToEPim", {5.3e33, 0.30}},  // SuperK2017lmeson
    {"nToMuPim", {3.5e33, 0.30}}, // SuperK2017lmeson
    {"nToNuPi0", {1.1e33, 0.30}}, // SuperK2014nupi
    {"nToNuEta", {1.6e32, 0.30}}, // SuperK2017lmeson
    {"nToNuK0", {1.3e32, 0.30}},  // SuperK2005susy
    {"nToEKm", {0.0, 0.30}},      // no dedicated limit
};

// eps_FSI per channel = primary-meson survival (none + elastic + produced),
// read from report/fsi_summary.txt.
std::map<std::string, double> read_fsi_survival(const std::string& path) {
   std::map<std::string, double> out;
   std::ifstream in(path);
   std::string line;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string key, meson;
      double none, elastic, cex, produced, absorbed;
      if (ss >> key >> meson >> none >> elastic >> cex >> produced >> absorbed)
         out[key] = none + elastic + produced;
   }
   return out;
}

// Per-channel min/mean/max signal-window fraction across nuclear models, read
// from report/window_summary.txt (the per-channel summary lines, which have a
// numeric 2nd field; the indented per-model detail lines do not).
struct WindowBand {
   double min, mean, max;
};
std::map<std::string, WindowBand> read_window_band(const std::string& path) {
   std::map<std::string, WindowBand> out;
   std::ifstream in(path);
   std::string line;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string key;
      double pf, mn, mean, mx;
      if (ss >> key >> pf >> mn >> mean >> mx) out[key] = {mn, mean, mx};
   }
   return out;
}

}  // namespace

int main() {
   const auto fsi = read_fsi_survival("report/fsi_summary.txt");
   const auto win = read_window_band("report/window_summary.txt");

   std::ofstream f("report/event_predictions.txt");
   f << "# DUNE nucleon-decay event predictions\n";
   f << "# fiducial=" << kFiducialKt << " kt, exposure=" << kExposureYr
     << " yr (" << kFiducialKt * kExposureYr << " kt.yr), tau = current SuperK "
        "90% C.L. limit per channel.\n";
   f << "# eps = eps_det * eps_FSI; N_lo/N_hi = nuclear-model band (window "
        "spread); N_off = FSI off (eps_FSI=1).\n";
   f << "# columns: key parent tau_limit eps_det eps_FSI w_lo w_hi N_cen N_lo "
        "N_hi N_off\n";

   printf("\n  DUNE event predictions  (%.0f kt x %.0f yr = %.0f kt.yr, "
          "tau = current SuperK limit)\n",
          kFiducialKt, kExposureYr, kFiducialKt * kExposureYr);
   printf("  %-9s %-15s %9s %6s %6s   %8s  %-16s %8s\n", "key", "mode",
          "tau[yr]", "e_det", "e_FSI", "N_cen", "[model band]", "N(FSIoff)");
   printf("  %s\n", std::string(92, '-').c_str());

   std::size_t nch = 0;
   const pdk::Channel* reg = pdk::channel_registry(nch);
   for (std::size_t i = 0; i < nch; ++i) {
      const pdk::Channel& c = reg[i];
      const auto in_it = kInputs.find(c.key);
      if (in_it == kInputs.end()) continue;
      const ChannelInput& ci = in_it->second;

      double n_nuc = kFiducialKt * (c.parent == pdk::Nucleon::Proton
                                        ? kProtonsPerKt
                                        : kNeutronsPerKt);
      double eps_fsi = fsi.count(c.key) ? fsi.at(c.key) : 1.0;
      double w_lo = 1.0, w_hi = 1.0;
      if (win.count(c.key) && win.at(c.key).mean > 0) {
         const WindowBand& wb = win.at(c.key);
         w_lo = wb.min / wb.mean;  // relative nuclear-model band
         w_hi = wb.max / wb.mean;
      }

      // Write machine-readable row (N's = -1 when no lifetime limit exists).
      double n_cen = -1, n_lo = -1, n_hi = -1, n_off = -1;
      if (ci.tau_limit_yr > 0) {
         // -expm1(-T/tau) = 1 - exp(-T/tau), accurate for the tiny T/tau here
         // (the naive 1 - exp(...) cancels to 0 in double precision).
         double rate = n_nuc * (-std::expm1(-kExposureYr / ci.tau_limit_yr));
         // DUNE-documented modes: eps_det already includes FSI, so do not
         // re-apply eps_FSI. Other modes: apply the cascade survival once.
         n_cen = rate * ci.eps_det * (ci.fsi_inclusive ? 1.0 : eps_fsi);
         n_lo = n_cen * w_lo;
         n_hi = n_cen * w_hi;
         // FSI-off = FSI-unfolded count, uniformly N_cen / eps_FSI. For the
         // factorized modes this is rate*eps_det; for the FSI-inclusive modes it
         // unfolds eps_FSI from the DUNE efficiency so the cascade suppression
         // (which is otherwise embedded in eps_det) remains visible.
         n_off = (eps_fsi > 0.0) ? n_cen / eps_fsi : n_cen;
      }
      f << c.key << ' ' << pdk::nucleon_name(c.parent) << ' '
        << ci.tau_limit_yr << ' ' << ci.eps_det << ' ' << eps_fsi << ' ' << w_lo
        << ' ' << w_hi << ' ' << n_cen << ' ' << n_lo << ' ' << n_hi << ' '
        << n_off << '\n';

      if (ci.tau_limit_yr > 0)
         printf("  %-9s %-15s %9.1e %6.2f %6.2f   %8.2f  [%6.2f, %6.2f]  %8.2f\n",
                c.key, c.pretty, ci.tau_limit_yr, ci.eps_det, eps_fsi, n_cen,
                n_lo, n_hi, n_off);
      else
         printf("  %-9s %-15s %9s %6.2f %6.2f   %8s  %-16s %8s\n", c.key,
                c.pretty, "(none)", ci.eps_det, eps_fsi, "-", "", "-");
   }
   printf("  %s\n", std::string(92, '-').c_str());
   printf("  Wrote report/event_predictions.txt\n\n");
   return 0;
}
