// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Analytic effective Ankowski-Sobczyk spectral function S(p, E) of argon-40
// (proton and neutron).
//
// Companion to plot_spectral.C (which shows the tabulated NuWro grid behind the
// `benhar` model). This visualises the S(p,E) that the `ankowski` model builds
// analytically -- there is no grid file to read, so the map is evaluated from
// the same formulas the sampler uses: a shell mean field (each argon shell's
// harmonic-oscillator momentum profile |phi_{nl}(p)|^2 paired with its Gaussian-
// smeared separation energy) plus a correlated 1/p^4 tail carrying the two-
// nucleon removal energy E = E_offset + p^2/2M. The shell band appears as
// horizontal stripes at the shell separation energies; the correlated part as a
// thin ridge curving up in E with p.
//
// The physics MIRRORS the Ankowski case of NucleonMomentumSampler::sample() and
// the shell tables in include/PDKMomentum.h -- keep the two in sync.
//
// Run from the project root:  root -l -b -q macros/plot_ankowski_sf.C
// Writes plots/ankowski_sf.png.
#include <cmath>
#include <vector>

#include "pdk_style.h"

namespace {

// Model constants, mirroring NuclearParams / PDKMomentum.h (GeV, fm units).
constexpr double kHbarC = 0.197327;       // GeV*fm
constexpr double kHoB = 1.9;              // HO oscillator length b [fm]
constexpr double kKmax = 0.65;            // SRC tail upper edge [GeV/c]
constexpr double kSrcOffset = 0.020;      // two-nucleon removal-energy offset [GeV]
constexpr double kCorrFraction = 0.20;    // as_corr_fraction
constexpr double kProtonMass = 0.938272;  // [GeV]
constexpr double kNeutronMass = 0.939565; // [GeV]
constexpr double kKfProton = 0.217;       // proton global Fermi momentum [GeV/c]
constexpr double kKfNeutron = 0.230;      // neutron global Fermi momentum [GeV/c]

struct Shell {
   int occ;
   double e_sep;    // [GeV]
   double e_width;  // [GeV]
   int n, l;
};

// argon_proton_shells() / argon_neutron_shells() from include/PDKMomentum.h.
const std::vector<Shell>& proton_shells() {
   static const std::vector<Shell> s = {
       {2, 0.0520, 0.0090, 0, 0}, {4, 0.0360, 0.0070, 0, 1},
       {2, 0.0320, 0.0060, 0, 1}, {6, 0.0180, 0.0040, 0, 2},
       {2, 0.0130, 0.0030, 1, 0}, {2, 0.0125, 0.0030, 0, 2},
   };
   return s;
}
const std::vector<Shell>& neutron_shells() {
   static const std::vector<Shell> s = {
       {2, 0.0480, 0.0090, 0, 0}, {4, 0.0320, 0.0070, 0, 1},
       {2, 0.0280, 0.0060, 0, 1}, {6, 0.0150, 0.0040, 0, 2},
       {2, 0.0110, 0.0030, 1, 0}, {4, 0.0100, 0.0030, 0, 2},
       {2, 0.0099, 0.0030, 0, 3},
   };
   return s;
}

// Generalized Laguerre L_n^alpha(x), mirroring laguerre() in PDKMomentum.h.
double laguerre(int n, double alpha, double x) {
   double lkm1 = 0.0, lk = 1.0;
   for (int k = 0; k < n; ++k) {
      double lkp1 = ((2 * k + 1 + alpha - x) * lk - (k + alpha) * lkm1) / (k + 1);
      lkm1 = lk;
      lk = lkp1;
   }
   return lk;
}

// p^2 |phi_{nl}(p)|^2, mirroring ho_weight() in PDKMomentum.h (p in GeV/c).
double ho_weight(int n, int l, double p) {
   double x = p * kHoB / kHbarC;
   x = x * x;
   double lag = laguerre(n, l + 0.5, x);
   return std::pow(x, l + 1) * lag * lag * std::exp(-x);
}

double gauss(double e, double mu, double sigma) {
   double z = (e - mu) / sigma;
   return std::exp(-0.5 * z * z) / sigma;
}

// Build the analytic S(p,E) on a (p[MeV], E[MeV]) grid for one nucleon. The
// mean-field and correlated parts are each normalised in sum_ij p^2 S so the
// correlated strength is exactly kCorrFraction, matching the sampler.
TH2D* build_sf(bool proton, const char* name, const char* title) {
   const int nP = 200, nE = 200;
   const double pMaxMeV = 800.0, eMaxMeV = 200.0;
   TH2D* h = new TH2D(name, title, nP, 0.0, pMaxMeV, nE, 0.0, eMaxMeV);

   const auto& shells = proton ? proton_shells() : neutron_shells();
   const double m = proton ? kProtonMass : kNeutronMass;
   const double kF = proton ? kKfProton : kKfNeutron;

   std::vector<double> mf(static_cast<std::size_t>(nP) * nE, 0.0);
   std::vector<double> tail(static_cast<std::size_t>(nP) * nE, 0.0);
   double sumMF = 0.0, sumTail = 0.0;

   for (int i = 0; i < nP; ++i) {
      double p = (h->GetXaxis()->GetBinCenter(i + 1)) * 1.0e-3;  // GeV/c
      if (p <= 0.0) continue;
      double p2 = p * p;

      // Mean field: sum over shells of occ * |phi_nl(p)|^2 * Gauss(E; e_sep).
      for (int j = 0; j < nE; ++j) {
         double e = (h->GetYaxis()->GetBinCenter(j + 1)) * 1.0e-3;  // GeV
         double s = 0.0;
         for (const Shell& sh : shells)
            s += sh.occ * (ho_weight(sh.n, sh.l, p) / p2) *
                 gauss(e, sh.e_sep, sh.e_width);
         mf[static_cast<std::size_t>(i) * nE + j] = s;
         sumMF += p2 * s;
      }

      // Correlated 1/p^4 tail: a delta in E at E = offset + p^2/2M, placed in the
      // bin containing it (no invented spread).
      if (p > kF && p < kKmax) {
         double e_tail = kSrcOffset + p2 / (2.0 * m);  // GeV
         int j = h->GetYaxis()->FindBin(e_tail * 1.0e3) - 1;
         if (j >= 0 && j < nE) {
            double s = 1.0 / (p2 * p2);
            tail[static_cast<std::size_t>(i) * nE + j] += s;
            sumTail += p2 * s;
         }
      }
   }

   double fMF = (sumMF > 0.0) ? (1.0 - kCorrFraction) / sumMF : 0.0;
   double fTail = (sumTail > 0.0) ? kCorrFraction / sumTail : 0.0;
   for (int i = 0; i < nP; ++i)
      for (int j = 0; j < nE; ++j) {
         std::size_t k = static_cast<std::size_t>(i) * nE + j;
         h->SetBinContent(i + 1, j + 1, fMF * mf[k] + fTail * tail[k]);
      }

   h->GetXaxis()->SetTitle("p [MeV/c]");
   h->GetYaxis()->SetTitle("E_{rem} [MeV]");
   h->GetZaxis()->SetTitle("S(p,E)  [arb.]");
   return h;
}

}  // namespace

void plot_ankowski_sf() {
   pdk_set_style();
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   TH2D* hp = build_sf(true, "ankP", "Argon-40 proton  S(p,E)  (Ankowski)");
   TH2D* hn = build_sf(false, "ankN", "Argon-40 neutron  S(p,E)  (Ankowski)");

   TCanvas* cv = new TCanvas("c", "ankowski spectral function", 1400, 560);
   cv->Divide(2, 1);
   for (TH2D* h : {hp, hn}) {
      h->GetXaxis()->SetRangeUser(0, 500);
      h->GetYaxis()->SetRangeUser(0, 150);
      h->GetZaxis()->SetTitleOffset(1.1);
      // Floor the log scale ~5 decades below the peak: the Gaussian shell tails
      // and empty tail bins otherwise drag the range down past 1e-30 and wash
      // out the shell band.
      h->SetMinimum(1.0e-5 * h->GetMaximum());
   }

   cv->cd(1);
   gPad->SetRightMargin(0.16);
   gPad->SetLogz();
   hp->Draw("COLZ");

   cv->cd(2);
   gPad->SetRightMargin(0.16);
   gPad->SetLogz();
   hn->Draw("COLZ");

   cv->cd(0);
   pdk_label("Ankowski-Sobczyk effective SF (analytic)");
   cv->SaveAs("plots/ankowski_sf.png");
   printf("Wrote plots/ankowski_sf.png\n");
}
