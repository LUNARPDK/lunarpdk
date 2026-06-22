// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Signal-window (containment) efficiency per channel and nuclear model.
//
// For each decay channel the observable hadron-meson momentum peaks at the
// free-decay value p_free = (M_N^2 - m_had^2)/(2 M_N) and is smeared by Fermi
// motion and binding. As a toy reconstruction-efficiency proxy we count the
// fraction of decays whose hadron-daughter momentum falls inside a fixed window
// [p_free - dP, p_free + dP], dP = 0.04 GeV/c (the same half-width as the kaon
// signal window [0.30,0.38] used in plot_kaon_models.C). The spread of this
// fraction across the ten nuclear models is the nuclear-model variation band
// that feeds the DUNE event-rate prediction (EventPredictor).
//
// Reads data/win_<channel>_<model>.txt (legacy 5-column FSI-off table,
// event nucleon_p d1_p d2_p e_rem; d2 = hadron-daughter momentum) for all
// channels x models written by make_plots.sh.
// Run from the project root:  root -l -b -q macros/window_eff.C
// Writes report/window_summary.txt and plots/window_eff.png.
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "pdk_style.h"

static const double kHalfWidth = 0.04;  // GeV/c, window half-width about p_free

struct Chan {
   const char* key;
   const char* label;
   double m_nucleon;  // parent nucleon mass [GeV]
   double m_hadron;   // hadron-daughter mass [GeV]
};

// Masses from include/PDKChannels.h.
static const double Mp = 0.93827, Mn = 0.93957;
static const double mPi0 = 0.13498, mPiC = 0.13957, mEta = 0.54786,
                    mKch = 0.49368, mK0 = 0.49761;

static double p_free(const Chan& c) {
   return (c.m_nucleon * c.m_nucleon - c.m_hadron * c.m_hadron) /
          (2.0 * c.m_nucleon);
}

// In-window fraction of the hadron-daughter momentum for one channel/model file.
static double window_fraction(const std::string& fname, double plo, double phi) {
   std::ifstream in(fname);
   if (!in.is_open()) {
      printf("Warning: could not open %s (run make_plots.sh first)\n",
             fname.c_str());
      return -1.0;
   }
   int ev;
   float np, d1, d2, er;
   long n = 0, win = 0;
   while (in >> ev >> np >> d1 >> d2 >> er) {
      ++n;
      if (d2 >= plo && d2 <= phi) ++win;
   }
   return (n > 0) ? double(win) / n : -1.0;
}

void window_eff() {
   pdk_set_style();

   const Chan chans[] = {
       {"pToKnu", "p#rightarrowK^{+}#bar{#nu}", Mp, mKch},
       {"pToMuK0", "p#rightarrow#mu^{+}K^{0}", Mp, mK0},
       {"pToEK0", "p#rightarrowe^{+}K^{0}", Mp, mK0},
       {"pToEEta", "p#rightarrowe^{+}#eta", Mp, mEta},
       {"pToMuEta", "p#rightarrow#mu^{+}#eta", Mp, mEta},
       {"pToNuPip", "p#rightarrow#bar{#nu}#pi^{+}", Mp, mPiC},
       {"pToEPi0", "p#rightarrowe^{+}#pi^{0}", Mp, mPi0},
       {"pToMuPi0", "p#rightarrow#mu^{+}#pi^{0}", Mp, mPi0},
       {"nToNuK0", "n#rightarrow#bar{#nu}K^{0}", Mn, mK0},
       {"nToEKm", "n#rightarrowe^{+}K^{-}", Mn, mKch},
       {"nToNuEta", "n#rightarrow#bar{#nu}#eta", Mn, mEta},
       {"nToEPim", "n#rightarrowe^{+}#pi^{-}", Mn, mPiC},
       {"nToMuPim", "n#rightarrow#mu^{+}#pi^{-}", Mn, mPiC},
       {"nToNuPi0", "n#rightarrow#bar{#nu}#pi^{0}", Mn, mPi0},
   };
   const int nchan = sizeof(chans) / sizeof(chans[0]);
   const char* models[] = {"gfg",  "lfg", "src",    "sf",      "hosm",
                           "br",   "gauss", "cfg",  "benhar",  "ankowski"};
   const int nmod = sizeof(models) / sizeof(models[0]);

   FILE* sum = fopen("report/window_summary.txt", "w");
   if (sum) {
      fprintf(sum, "# Signal-window (containment) efficiency: fraction of the "
                   "hadron-daughter\n");
      fprintf(sum, "# momentum within [p_free-%.2f, p_free+%.2f] GeV/c, per "
                   "channel and model.\n", kHalfWidth, kHalfWidth);
      fprintf(sum, "# %-9s %8s %8s %8s %8s   %-8s %-8s %-8s\n", "channel",
              "p_free", "min", "mean", "max", "min_mod", "max_mod", "");
   }

   // For the plot: per-channel min / mean / max across models.
   std::vector<double> vmin(nchan), vmax(nchan), vmean(nchan), vfree(nchan);

   for (int c = 0; c < nchan; ++c) {
      double pf = p_free(chans[c]);
      double plo = pf - kHalfWidth, phi = pf + kHalfWidth;
      vfree[c] = pf;
      double lo = 1e9, hi = -1e9, sumf = 0.0;
      int nok = 0;
      std::string lo_mod = "-", hi_mod = "-";
      // Per-model detail lines accumulate into the summary too.
      std::string detail;
      for (int m = 0; m < nmod; ++m) {
         double f = window_fraction(
             Form("data/win_%s_%s.txt", chans[c].key, models[m]), plo, phi);
         if (f < 0) continue;
         detail += Form("    %-9s %-8s %8.4f\n", chans[c].key, models[m], f);
         if (f < lo) { lo = f; lo_mod = models[m]; }
         if (f > hi) { hi = f; hi_mod = models[m]; }
         sumf += f;
         ++nok;
      }
      double mean = (nok > 0) ? sumf / nok : 0.0;
      vmin[c] = (nok > 0) ? lo : 0.0;
      vmax[c] = (nok > 0) ? hi : 0.0;
      vmean[c] = mean;
      if (sum) {
         fprintf(sum, "  %-9s %8.4f %8.4f %8.4f %8.4f   %-8s %-8s\n",
                 chans[c].key, pf, vmin[c], mean, vmax[c], lo_mod.c_str(),
                 hi_mod.c_str());
         fputs(detail.c_str(), sum);
      }
   }
   if (sum) fclose(sum);

   // Plot: per-channel mean with a min->max band (asymmetric error bars).
   TH1F* frame = new TH1F("win_frame",
                          "Signal-window efficiency vs nuclear model;;"
                          "in-window fraction  (p_{free} #pm 0.04 GeV/c)",
                          nchan, 0.0, nchan);
   frame->SetMinimum(0.0);
   frame->SetMaximum(0.75);
   for (int c = 0; c < nchan; ++c)
      frame->GetXaxis()->SetBinLabel(c + 1, chans[c].label);
   frame->GetXaxis()->SetLabelSize(0.048);
   frame->GetXaxis()->LabelsOption("v");  // rotate the 14 channel names vertical
   frame->GetYaxis()->SetTitleOffset(0.90);

   TGraphAsymmErrors* g = new TGraphAsymmErrors(nchan);
   for (int c = 0; c < nchan; ++c) {
      g->SetPoint(c, c + 0.5, vmean[c]);
      g->SetPointError(c, 0.0, 0.0, vmean[c] - vmin[c], vmax[c] - vmean[c]);
   }
   g->SetMarkerStyle(20);
   g->SetMarkerSize(1.1);
   g->SetMarkerColor(kAzure + 2);
   g->SetLineColor(kAzure + 2);
   g->SetLineWidth(2);

   TCanvas* cv = new TCanvas("c", "window efficiency", 1300, 600);
   gPad->SetBottomMargin(0.26);
   frame->Draw("AXIS");
   g->Draw("P SAME");
   pdk_label("FSI off  #bullet  band = 10 nuclear models");
   cv->SaveAs("plots/window_eff.png");
   printf("Wrote plots/window_eff.png and report/window_summary.txt\n");
}
