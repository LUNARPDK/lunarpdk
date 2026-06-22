// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Plot the kaon momentum distribution from the generated events.
// Run in batch from the project root:  root -l -b -q macros/plot_kaon.C
//
// Reads data/output.root (TTree "T") and writes plots/kaon_p.png.
#include "pdk_style.h"

void plot_kaon() {
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

   TCanvas *c = new TCanvas("c", "PDK kaon momentum", 800, 600);

   TH1F *h = new TH1F("h_d2_p",
                      "Hadron-side daughter momentum;"
                      "Daughter momentum p_{d2} [GeV/c];Events",
                      80, 0.0, 0.6);
   T->Draw("d2_p >> h_d2_p", "", "goff");
   pdk_style_hist(h, pdk_model_color("gfg"));  // blue signature
   h->Draw("HIST");
   pdk_label("lab-frame hadron daughter");

   c->SaveAs("plots/kaon_p.png");
   printf("Wrote plots/kaon_p.png  (%lld events, mean = %.4f GeV/c)\n",
          (long long)h->GetEntries(), h->GetMean());
   f->Close();
}
