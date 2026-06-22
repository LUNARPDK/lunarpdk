// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Final-state-interaction outcome fractions across decay channels.
//
// For each benchmark channel the FSI cascade labels every event by the fate of
// the PRIMARY meson, written on the per-event header line
//   # event <i>: nucleon_p=<GeV> e_rem=<GeV> outcome=<none|elastic|cex|produced|absorbed>
// This macro tallies those outcomes per channel and draws a stacked bar chart
// (one bar per channel, split into the five outcome fractions), so the contrast
// between the transparent K+ and the strongly-reworked pi0 / antikaon is visible
// at a glance. It also writes report/fsi_summary.txt for the report prose.
//
// Reads data/fsi_<channel>_on.txt written by make_plots.sh.
// Run from the project root:  root -l -b -q macros/fsi_outcomes.C
// Writes plots/fsi_outcomes.png and report/fsi_summary.txt.
#include <cstdio>
#include <fstream>
#include <string>

#include "pdk_style.h"

// Outcome order (bottom -> top of each stacked bar).
static const int kNoutc = 5;
static const char* kOutcName[kNoutc] = {"none", "elastic", "cex", "produced",
                                        "absorbed"};
static const char* kOutcLabel[kNoutc] = {"none (escaped)", "elastic",
                                         "charge exch.", "produced extra",
                                         "absorbed"};
static Color_t outc_color(int i) {
   const Color_t c[kNoutc] = {kGreen + 2, kAzure + 1, kOrange + 7, kViolet + 1,
                              kRed + 1};
   return c[i];
}

// One channel: file key, the analysed channel name, and a display label.
struct Chan {
   const char* key;    // data/fsi_<key>_on.txt
   const char* label;  // axis label
};

// Tally the five outcome fractions for one file; returns the number of events.
static long read_outcomes(const std::string& fname, double frac[kNoutc]) {
   for (int i = 0; i < kNoutc; ++i) frac[i] = 0.0;
   std::ifstream in(fname);
   if (!in.is_open()) {
      printf("Warning: could not open %s (run make_plots.sh first)\n",
             fname.c_str());
      return 0;
   }
   long counts[kNoutc] = {0, 0, 0, 0, 0};
   long n = 0;
   std::string line;
   while (std::getline(in, line)) {
      if (line.rfind("# event", 0) != 0) continue;
      auto pos = line.find("outcome=");
      if (pos == std::string::npos) continue;
      std::string oc = line.substr(pos + 8);
      // strip trailing whitespace
      while (!oc.empty() && (oc.back() == ' ' || oc.back() == '\n' ||
                             oc.back() == '\r' || oc.back() == '\t'))
         oc.pop_back();
      for (int i = 0; i < kNoutc; ++i) {
         if (oc == kOutcName[i]) {
            ++counts[i];
            ++n;
            break;
         }
      }
   }
   if (n > 0)
      for (int i = 0; i < kNoutc; ++i) frac[i] = double(counts[i]) / n;
   return n;
}

void fsi_outcomes() {
   pdk_set_style();

   // All 14 channels, grouped by meson species (the outcome depends on the
   // meson, not the lepton), proton modes then neutron modes.
   const Chan chans[] = {
       {"pToKnu", "K^{+}#bar{#nu}"},   {"pToMuK0", "#mu^{+}K^{0}"},
       {"pToEK0", "e^{+}K^{0}"},       {"pToEEta", "e^{+}#eta"},
       {"pToMuEta", "#mu^{+}#eta"},    {"pToNuPip", "#bar{#nu}#pi^{+}"},
       {"pToEPi0", "e^{+}#pi^{0}"},    {"pToMuPi0", "#mu^{+}#pi^{0}"},
       {"nToNuK0", "#bar{#nu}K^{0}"},  {"nToEKm", "e^{+}K^{-}"},
       {"nToNuEta", "#bar{#nu}#eta"},  {"nToEPim", "e^{+}#pi^{-}"},
       {"nToMuPim", "#mu^{+}#pi^{-}"}, {"nToNuPi0", "#bar{#nu}#pi^{0}"},
   };
   const int nchan = sizeof(chans) / sizeof(chans[0]);

   // One histogram per outcome, nchan category bins; stacked.
   TH1F* h[kNoutc];
   for (int o = 0; o < kNoutc; ++o) {
      h[o] = new TH1F(Form("h_outc_%d", o), "", nchan, 0.0, nchan);
      h[o]->SetFillColor(outc_color(o));
      h[o]->SetLineColor(kBlack);
      h[o]->SetLineWidth(1);
   }

   FILE* sum = fopen("report/fsi_summary.txt", "w");
   if (sum) {
      fprintf(sum, "# FSI primary-meson outcome fractions (per generated decay)\n");
      fprintf(sum, "# %-9s %8s %8s %8s %8s %8s %8s %10s\n", "channel", "meson",
              "none", "elastic", "cex", "produced", "absorbed", "n_events");
   }

   for (int c = 0; c < nchan; ++c) {
      double frac[kNoutc];
      long n = read_outcomes(Form("data/fsi_%s_on.txt", chans[c].key), frac);
      for (int o = 0; o < kNoutc; ++o) {
         h[o]->SetBinContent(c + 1, frac[o]);
         h[o]->GetXaxis()->SetBinLabel(c + 1, chans[c].label);
      }
      if (sum)
         fprintf(sum, "  %-9s %8s %8.4f %8.4f %8.4f %8.4f %8.4f %10ld\n",
                 chans[c].key, chans[c].label, frac[0], frac[1], frac[2],
                 frac[3], frac[4], n);
   }
   if (sum) fclose(sum);

   THStack* st = new THStack(
       "st",
       "FSI fate of the primary meson by channel;decay channel;fraction of decays");
   for (int o = 0; o < kNoutc; ++o) st->Add(h[o]);

   TCanvas* cv = new TCanvas("c", "FSI outcomes", 1300, 620);
   st->Draw("hist bar2");
   st->SetMaximum(1.28);
   st->GetXaxis()->SetLabelSize(0.050);
   st->GetYaxis()->SetTitleOffset(0.85);
   st->GetYaxis()->SetRangeUser(0.0, 1.28);
   gPad->SetTopMargin(0.08);

   // Legend in the empty band above the (unit-height) bars.
   TLegend* leg = new TLegend(0.16, 0.82, 0.93, 0.90);
   leg->SetNColumns(5);
   leg->SetFillColorAlpha(kWhite, 0.0);
   leg->SetTextSize(0.028);
   for (int o = 0; o < kNoutc; ++o) leg->AddEntry(h[o], kOutcLabel[o], "f");
   leg->Draw();

   pdk_label("argon, FSI on");
   cv->SaveAs("plots/fsi_outcomes.png");
   printf("Wrote plots/fsi_outcomes.png and report/fsi_summary.txt\n");
}
