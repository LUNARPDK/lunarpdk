// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// FSI on-vs-off meson momentum spectra for the channels NOT shown in
// fsi_compare.png (which covers K+ and pi0). Three panels:
//   left   : p -> nubar pi+   (pi+, strongly interacting)
//   middle : p -> e+ eta      (eta, converts via the N*(1535))
//   right  : p -> mu+ K0      (neutral kaon, mildly interacting; K0 <-> K0bar)
//
// Same construction as macros/compare_fsi.C: each curve is normalized PER
// GENERATED DECAY, so the FSI-on histogram shows both the depletion (smaller
// area = mesons lost to absorption / charge exchange) and the softening (a
// low-momentum tail from quasi-elastic energy loss).
//
// Reads the off/on files written by make_plots.sh:
//   data/fsi_pToNuPip_{off,on}.txt  data/fsi_pToEEta_{off,on}.txt
//   data/fsi_pToMuK0_{off,on}.txt
// Run from the project root:  root -l -b -q macros/fsi_spectra_channels.C
// Writes plots/fsi_spectra_channels.png.
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "pdk_style.h"

// Legacy 5-column file: hadron-daughter momentum (column d2). Returns # events.
static long fill_off(const char* fname, TH1F* h) {
   std::ifstream in(fname);
   if (!in.is_open()) {
      printf("Warning: could not open %s (run make_plots.sh first)\n", fname);
      return 0;
   }
   int ev;
   float np, d1, d2, er;
   long n = 0;
   while (in >> ev >> np >> d1 >> d2 >> er) {
      h->Fill(d2);
      ++n;
   }
   return n;
}

// Post-FSI per-particle file: momentum of every particle whose PDG is in `pdgs`
// (a list, so the neutral kaon K0/K0bar can be collected together). Returns the
// number of generated decays (# event headers); `filled` = # target mesons.
static long fill_on(const char* fname, const std::vector<int>& pdgs, TH1F* h,
                    long& filled) {
   std::ifstream in(fname);
   filled = 0;
   if (!in.is_open()) {
      printf("Warning: could not open %s (run make_plots.sh first)\n", fname);
      return 0;
   }
   std::string line;
   long nev = 0;
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
         for (int want : pdgs) {
            if (p == want) {
               h->Fill(std::sqrt(px * px + py * py + pz * pz));
               ++filled;
               break;
            }
         }
      }
   }
   return nev;
}

static void draw_panel(const char* off_file, const char* on_file,
                       const std::vector<int>& pdgs, const char* title,
                       double p_free, const char* corner) {
   TH1F* hoff = new TH1F(Form("hoff_%s", corner), title, 110, 0.0, 0.7);
   TH1F* hon = new TH1F(Form("hon_%s", corner), title, 110, 0.0, 0.7);

   long noff = fill_off(off_file, hoff);
   long filled = 0;
   long non = fill_on(on_file, pdgs, hon, filled);

   if (noff > 0) hoff->Scale(1.0 / noff, "width");
   double yield = (non > 0) ? double(filled) / non : 0.0;
   if (non > 0) hon->Scale(1.0 / non, "width");

   double ymax = TMath::Max(hoff->GetMaximum(), hon->GetMaximum());
   if (ymax <= 0) ymax = 1.0;
   hoff->SetMaximum(ymax * 1.30);
   hoff->SetMinimum(0.0);
   hoff->SetLineColor(kAzure + 1);
   hoff->SetLineWidth(3);
   hon->SetLineColor(kRed + 1);
   hon->SetLineWidth(3);
   hon->SetFillColorAlpha(kRed + 1, 0.18);

   hoff->Draw("HIST");
   hon->Draw("HIST SAME");

   TLine* ref = new TLine(p_free, 0.0, p_free, ymax * 1.30);
   ref->SetLineStyle(2);
   ref->SetLineColor(kGray + 2);
   ref->Draw();
   TLatex tx;
   tx.SetTextFont(42);
   tx.SetTextColor(kGray + 2);
   tx.SetTextSize(0.034);
   tx.SetTextAngle(90);
   tx.DrawLatex(p_free - 0.030, ymax * 0.28, "free decay");

   TLegend* leg = new TLegend(0.14, 0.74, 0.62, 0.90);
   leg->AddEntry(hoff, Form("FSI off  (#LTp#GT = %.3f)", hoff->GetMean()), "l");
   leg->AddEntry(hon,
                 Form("FSI on  (yield %.0f%%, #LTp#GT = %.3f)", 100.0 * yield,
                      hon->GetMean()),
                 "f");
   leg->Draw();
   pdk_label(corner);
}

void fsi_spectra_channels() {
   pdk_set_style();

   const double Mp = 0.93827, mPiC = 0.13957, mEta = 0.54786, mK0 = 0.49761;
   auto pfree = [&](double m) { return (Mp * Mp - m * m) / (2.0 * Mp); };

   TCanvas* c = new TCanvas("c", "FSI spectra by channel", 1500, 480);
   c->Divide(3, 1);

   c->cd(1);
   draw_panel("data/fsi_pToNuPip_off.txt", "data/fsi_pToNuPip_on.txt", {211},
              "p #rightarrow #bar{#nu}#pi^{+};"
              "p_{#pi^{+}} [GeV/c];(1/N_{dec}) dN/dp  [(GeV/c)^{-1}]",
              pfree(mPiC), "#pi^{+} FSI in argon");

   c->cd(2);
   draw_panel("data/fsi_pToEEta_off.txt", "data/fsi_pToEEta_on.txt", {221},
              "p #rightarrow e^{+}#eta;"
              "p_{#eta} [GeV/c];(1/N_{dec}) dN/dp  [(GeV/c)^{-1}]",
              pfree(mEta), "#eta FSI in argon");

   c->cd(3);
   draw_panel("data/fsi_pToMuK0_off.txt", "data/fsi_pToMuK0_on.txt", {311, -311},
              "p #rightarrow #mu^{+}K^{0};"
              "p_{K^{0}} [GeV/c];(1/N_{dec}) dN/dp  [(GeV/c)^{-1}]",
              pfree(mK0), "K^{0} FSI in argon");

   c->SaveAs("plots/fsi_spectra_channels.png");
   printf("Wrote plots/fsi_spectra_channels.png\n");
}
