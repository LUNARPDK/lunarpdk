// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Compare the removal-energy distribution of the tabulated argon spectral
// function for the proton vs the neutron grid (gsf_Ar40P.grid / gsf_Ar40N.grid).
// Reads data/erem_sf_proton.txt and data/erem_sf_neutron.txt
// (columns: event nucleon_p d1_p d2_p e_rem) produced by make_plots.sh.
//
// Run in batch from the project root:  root -l -b -q macros/compare_erem.C
// Writes plots/erem_proton_neutron.png.
#include "pdk_style.h"

void compare_erem() {
   pdk_set_style();

   const int nf = 2;
   const char* files[nf] = {"data/erem_sf_proton.txt", "data/erem_sf_neutron.txt"};
   const char* labels[nf] = {"argon proton SF", "argon neutron SF"};
   Color_t colors[nf] = {pdk_model_color("src"), pdk_model_color("gfg")};

   TCanvas* c = new TCanvas("c", "tabulated SF removal energy: p vs n", 800, 600);
   TLegend* leg = new TLegend(0.52, 0.74, 0.93, 0.90);

   double ymax = 0.0;
   TH1F* hists[nf] = {nullptr};
   for (int i = 0; i < nf; ++i) {
      std::ifstream in(files[i]);
      if (!in.is_open()) { printf("Warning: could not open %s\n", files[i]); continue; }
      TH1F* h = new TH1F(TString::Format("h_%d", i),
                         "Tabulated argon spectral function;"
                         "Removal energy E_{rem} [GeV];Probability density",
                         120, 0.0, 0.12);
      int ev; float np, d1, d2, er;
      while (in >> ev >> np >> d1 >> d2 >> er) h->Fill(er);
      in.close();
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
      h->SetLineColor(colors[i]);
      h->SetLineWidth(2);
      hists[i] = h;
      ymax = TMath::Max(ymax, h->GetMaximum());
      leg->AddEntry(h, labels[i], "l");
   }
   bool first = true;
   for (int i = 0; i < nf; ++i) {
      if (!hists[i]) continue;
      if (first) { hists[i]->SetMaximum(ymax * 1.25); hists[i]->Draw("HIST"); first = false; }
      else hists[i]->Draw("HIST SAME");
   }
   leg->Draw();
   pdk_label("tabulated S(p,E)");

   c->SaveAs("plots/erem_proton_neutron.png");
   printf("Wrote plots/erem_proton_neutron.png\n");
}
