// Compare the lab-frame meson momentum spectrum with and without final-state
// interactions (FSI), for the two benchmark channels:
//   left  : p -> K+ nu   (kaon, nearly transparent)
//   right : p -> e+ pi0  (pion, strongly absorbed near the Delta)
//
// Reads the four files written by make_plots.sh:
//   data/fsi_pToKnu_off.txt   data/fsi_pToKnu_on.txt
//   data/fsi_pToEPi0_off.txt  data/fsi_pToEPi0_on.txt
// The "off" files are the legacy 5-column table (event nucleon_p d1_p d2_p e_rem);
// the "on" files are the post-FSI per-particle format (event pdg px py pz E
// outcome) with a "# event <i>: ..." header line per event.
//
// Both curves are normalized PER GENERATED DECAY (not to unit area), so the
// integral of each curve is the mean number of signal mesons (K+ or pi0) that
// leave the nucleus per decay. FSI-off integrates to 1 (one meson per event);
// FSI-on sits below it (mesons lost to absorption / charge exchange) and is
// shifted softer (elastic/quasi-elastic energy loss).
//
// Run in batch from the project root:  root -l -b -q macros/compare_fsi.C
// Writes plots/fsi_compare.png.
#include <cstdio>
#include <fstream>
#include <string>

#include "pdk_style.h"

// Fill h with the hadron-daughter momentum from a legacy 5-column file; returns
// the number of events (= rows).
static long fill_off(const char* fname, TH1F* h) {
   std::ifstream in(fname);
   if (!in.is_open()) {
      printf("Warning: could not open %s (run compare_fsi.sh first)\n", fname);
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

// Fill h with the momentum of every particle of PDG `pdg` from a post-FSI
// per-particle file; returns the number of generated decays (# event headers)
// and reports how many target mesons were found via `filled`.
static long fill_on(const char* fname, int pdg, TH1F* h, long& filled) {
   std::ifstream in(fname);
   filled = 0;
   if (!in.is_open()) {
      printf("Warning: could not open %s (run compare_fsi.sh first)\n", fname);
      return 0;
   }
   std::string line;
   long nev = 0;
   while (std::getline(in, line)) {
      if (line.empty()) continue;
      if (line[0] == '#') {
         if (line.rfind("# event", 0) == 0) ++nev;  // one header per decay
         continue;
      }
      int e, p;
      float px, py, pz, en;
      char oc[40];
      if (std::sscanf(line.c_str(), "%d %d %f %f %f %f %39s", &e, &p, &px, &py,
                      &pz, &en, oc) >= 6 &&
          p == pdg) {
         h->Fill(std::sqrt(px * px + py * py + pz * pz));
         ++filled;
      }
   }
   return nev;
}

static void draw_panel(const char* off_file, const char* on_file, int pdg,
                       const char* title, double p_free, const char* free_tag) {
   TH1F* hoff = new TH1F(Form("hoff_%d", pdg), title, 110, 0.0, 0.7);
   TH1F* hon = new TH1F(Form("hon_%d", pdg), title, 110, 0.0, 0.7);

   long noff = fill_off(off_file, hoff);
   long filled = 0;
   long non = fill_on(on_file, pdg, hon, filled);

   if (noff > 0) hoff->Scale(1.0 / noff, "width");  // -> per-decay density
   double yield = (non > 0) ? double(filled) / non : 0.0;
   if (non > 0) hon->Scale(1.0 / non, "width");

   double ymax = TMath::Max(hoff->GetMaximum(), hon->GetMaximum());
   hoff->SetMaximum(ymax * 1.30);
   hoff->SetMinimum(0.0);
   hoff->SetLineColor(kAzure + 1);
   hoff->SetLineWidth(3);
   hon->SetLineColor(kRed + 1);
   hon->SetLineWidth(3);
   hon->SetFillColorAlpha(kRed + 1, 0.18);

   hoff->Draw("HIST");
   hon->Draw("HIST SAME");

   // Free-decay reference (proton at rest): p = (M_p^2 - m^2)/(2 M_p).
   TLine* ref = new TLine(p_free, 0.0, p_free, ymax * 1.30);
   ref->SetLineStyle(2);
   ref->SetLineColor(kGray + 2);
   ref->Draw();
   TLatex tx;
   tx.SetTextFont(42);
   tx.SetTextColor(kGray + 2);
   tx.SetTextSize(0.030);
   tx.SetTextAngle(90);
   tx.DrawLatex(p_free - 0.028, ymax * 0.30, free_tag);

   TLegend* leg = new TLegend(0.16, 0.74, 0.62, 0.90);
   leg->AddEntry(hoff, Form("FSI off  (#LTp#GT = %.3f)", hoff->GetMean()), "l");
   leg->AddEntry(hon,
                 Form("FSI on  (yield %.0f%%, #LTp#GT = %.3f)", 100.0 * yield,
                      hon->GetMean()),
                 "f");
   leg->Draw();
}

void compare_fsi() {
   pdk_set_style();

   const double Mp = 0.93827, mK = 0.49368, mPi0 = 0.13498;

   TCanvas* c = new TCanvas("c", "FSI on vs off", 1200, 520);
   c->Divide(2, 1);

   c->cd(1);
   draw_panel("data/fsi_pToKnu_off.txt", "data/fsi_pToKnu_on.txt", 321,
              "p #rightarrow K^{+}#bar{#nu}: kaon spectrum;"
              "p_{K^{+}} [GeV/c];(1/N_{dec}) dN/dp  [(GeV/c)^{-1}]",
              (Mp * Mp - mK * mK) / (2.0 * Mp), "free decay");
   pdk_label("K^{+} FSI in argon");

   c->cd(2);
   draw_panel("data/fsi_pToEPi0_off.txt", "data/fsi_pToEPi0_on.txt", 111,
              "p #rightarrow e^{+}#pi^{0}: pion spectrum;"
              "p_{#pi^{0}} [GeV/c];(1/N_{dec}) dN/dp  [(GeV/c)^{-1}]",
              (Mp * Mp - mPi0 * mPi0) / (2.0 * Mp), "free decay");
   pdk_label("#pi^{0} FSI in argon");

   c->SaveAs("plots/fsi_compare.png");
   printf("Wrote plots/fsi_compare.png\n");
}
