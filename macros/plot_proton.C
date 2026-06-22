// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Plot the initial proton (Fermi) momentum distribution from the generated events.
// Run in batch from the project root:  root -l -b -q macros/plot_proton.C
//
// Reads data/output.root (TTree "T") and writes plots/proton_p.png.
#include "pdk_style.h"

void plot_proton() {
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

   TCanvas *c = new TCanvas("c", "PDK proton momentum", 800, 600);

   TH1F *h = new TH1F("h_nucleon_p",
                      "Initial nucleon Fermi momentum;"
                      "Nucleon momentum p_{N} [GeV/c];Events",
                      50, 0.0, 0.3);
   T->Draw("nucleon_p >> h_nucleon_p", "", "goff");
   pdk_style_hist(h, pdk_model_color("src"));  // red signature
   h->Draw("HIST");
   pdk_label("initial nucleon");

   c->SaveAs("plots/proton_p.png");
   printf("Wrote plots/proton_p.png  (%lld events, mean = %.4f GeV/c)\n",
          (long long)h->GetEntries(), h->GetMean());
   f->Close();
}
