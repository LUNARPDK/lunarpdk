// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Per-channel signal-meson yield (transparency) after FSI.
//
// For each benchmark channel, the yield is the number of signal mesons of the
// original species that escape the nucleus, per generated decay:
//   yield = (# escaped mesons with the signal PDG) / (# decays).
// A yield near 1 means the nucleus is transparent to that meson; a low yield
// means FSI absorbs or charge-exchanges away the signal. This compact bar chart
// orders the channels from transparent (K+) to opaque (K-, pi0).
//
// Reads data/fsi_<channel>_on.txt written by make_plots.sh.
// Run from the project root:  root -l -b -q macros/fsi_yield.C
// Writes plots/fsi_yield.png.
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "pdk_style.h"

// Count escaped mesons matching any PDG in `pdgs`; returns yield per decay.
static double meson_yield(const std::string& fname,
                          const std::vector<int>& pdgs) {
   std::ifstream in(fname);
   if (!in.is_open()) {
      printf("Warning: could not open %s (run make_plots.sh first)\n",
             fname.c_str());
      return 0.0;
   }
   long nev = 0, filled = 0;
   std::string line;
   while (std::getline(in, line)) {
      if (line.empty()) continue;
      if (line[0] == '#') {
         if (line.rfind("# event", 0) == 0) ++nev;
         continue;
      }
      int e, p;
      float px, py, pz, en;
      char oc[40];
      if (std::sscanf(line.c_str(), "%d %d %f %f %f %f %39s", &e, &p, &px, &py,
                      &pz, &en, oc) >= 6) {
         for (int want : pdgs)
            if (p == want) {
               ++filled;
               break;
            }
      }
   }
   return (nev > 0) ? double(filled) / nev : 0.0;
}

struct Chan {
   const char* key;
   const char* label;
   std::vector<int> pdgs;
};

void fsi_yield() {
   pdk_set_style();

   // Ordered transparent -> opaque.
   std::vector<Chan> chans = {
       {"pToKnu", "K^{+}", {321}},      {"pToMuK0", "K^{0}", {311, -311}},
       {"pToEEta", "#eta", {221}},      {"pToNuPip", "#pi^{+}", {211}},
       {"pToEPi0", "#pi^{0}", {111}},   {"nToEKm", "K^{-}", {-321}},
   };
   const int nchan = chans.size();

   TH1F* h = new TH1F("h_yield",
                      "Signal-meson yield after FSI in argon;decay channel;"
                      "escaped mesons per decay",
                      nchan, 0.0, nchan);
   for (int c = 0; c < nchan; ++c) {
      double y = meson_yield(Form("data/fsi_%s_on.txt", chans[c].key),
                             chans[c].pdgs);
      h->SetBinContent(c + 1, y);
      h->GetXaxis()->SetBinLabel(c + 1, chans[c].label);
   }
   h->SetMaximum(1.15);
   h->SetMinimum(0.0);
   h->SetBarWidth(0.7);
   h->SetBarOffset(0.15);
   h->SetFillColorAlpha(kAzure + 1, 0.55);
   h->SetLineColor(kAzure + 2);
   h->SetLineWidth(2);
   h->GetXaxis()->SetLabelSize(0.052);

   TCanvas* cv = new TCanvas("c", "FSI yield", 900, 600);
   h->Draw("bar2");

   // Annotate each bar with its percentage.
   TLatex tx;
   tx.SetTextFont(42);
   tx.SetTextSize(0.034);
   tx.SetTextAlign(21);  // centered, bottom
   for (int c = 0; c < nchan; ++c) {
      double y = h->GetBinContent(c + 1);
      tx.DrawLatex(h->GetBinCenter(c + 1), y + 0.02,
                   Form("%.0f%%", 100.0 * y));
   }

   pdk_label("argon, FSI on");
   cv->SaveAs("plots/fsi_yield.png");
   printf("Wrote plots/fsi_yield.png\n");
}
