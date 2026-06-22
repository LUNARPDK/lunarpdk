// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Plot the proton removal (separation) energy distribution.
// Run in batch from the project root:  root -l -b -q macros/plot_removal.C
//
// Reads data/output.root (TTree "T") and writes plots/removal_E.png.
#include "pdk_style.h"

void plot_removal() {
   pdk_set_style();
   gStyle->SetOptStat("emr");  // entries, mean, rms (single-distribution plot)

   TFile *f = TFile::Open("data/output.root", "read");
   if (!f || f->IsZombie()) {
      printf("Error: could not open data/output.root\n");
      return;
   }
   TTree *T = nullptr;
   f->GetObject("T", T);
   if (!T) {
      printf("Error: TTree 'T' not found in data/output.root\n");
      return;
   }

   TCanvas *c = new TCanvas("c", "PDK removal energy", 800, 600);

   TH1F *h = new TH1F("h_e_rem",
                      "Proton removal energy;Removal energy E_{rem} [GeV];Events",
                      80, 0.0, 0.25);
   T->Draw("e_rem >> h_e_rem", "", "goff");
   pdk_style_hist(h, pdk_model_color("benhar"));  // orange signature
   h->Draw("HIST");
   pdk_label("removal energy");

   c->SaveAs("plots/removal_E.png");
   printf("Wrote plots/removal_E.png  (%lld events, mean = %.4f GeV)\n",
          (long long)h->GetEntries(), h->GetMean());
   f->Close();
}
