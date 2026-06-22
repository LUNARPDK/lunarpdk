// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Tabulated spectral function S(p, E) of argon-40 (proton and neutron grids).
//
// Directly visualises the input that drives the `benhar` / `ankowski` momentum
// models: the JLab E12-14-012 NuWro grids in config/sf/. The map shows the
// mean-field shell band at low removal energy and low momentum plus the
// correlated tail extending to high p and high E -- the structure the analytic
// models only approximate. Parsing mirrors load_nuwro_grid in
// include/PDKSpectral.h.
//
// Run from the project root:  root -l -b -q macros/plot_spectral.C
// Writes plots/spectral_pe.png.
#include <fstream>
#include <string>

#include "pdk_style.h"

// Load a NuWro grid file into a TH2D(p[MeV], E[MeV]) of S(p,E). Returns nullptr
// on failure.
static TH2D* load_grid(const std::string& path, const char* name,
                       const char* title) {
   std::ifstream f(path);
   if (!f.is_open()) {
      printf("Warning: could not open %s\n", path.c_str());
      return nullptr;
   }
   int eRes = 0, pRes = 0;
   double eMin, pMin, eMax, pMax;
   if (!(f >> eRes >> pRes >> eMin >> pMin >> eMax >> pMax) || eRes <= 0 ||
       pRes <= 0) {
      printf("Warning: malformed header in %s\n", path.c_str());
      return nullptr;
   }
   // Bin edges span the axis range; grid points are bin centres.
   double dp = (pMax - pMin) / (pRes - 1), de = (eMax - eMin) / (eRes - 1);
   TH2D* h = new TH2D(name, title, pRes, pMin - 0.5 * dp, pMax + 0.5 * dp, eRes,
                      eMin - 0.5 * de, eMax + 0.5 * de);
   for (int i = 0; i < pRes; ++i) {
      double p_mev;
      if (!(f >> p_mev)) { printf("Truncated %s\n", path.c_str()); break; }
      for (int j = 0; j < eRes; ++j) {
         double e_mev, val;
         if (!(f >> e_mev >> val)) { printf("Truncated %s\n", path.c_str()); break; }
         if (val < 0) val = 0;
         h->SetBinContent(i + 1, j + 1, val);
      }
   }
   h->GetXaxis()->SetTitle("p [MeV/c]");
   h->GetYaxis()->SetTitle("E_{rem} [MeV]");
   h->GetZaxis()->SetTitle("S(p,E)  [arb.]");
   return h;
}

void plot_spectral() {
   pdk_set_style();
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   TH2D* hp = load_grid("config/sf/gsf_Ar40P.grid", "sfP",
                        "Argon-40 proton  S(p,E)");
   TH2D* hn = load_grid("config/sf/gsf_Ar40N.grid", "sfN",
                        "Argon-40 neutron  S(p,E)");
   if (!hp || !hn) {
      printf("Could not load spectral grids; aborting.\n");
      return;
   }

   TCanvas* cv = new TCanvas("c", "spectral functions", 1400, 560);
   cv->Divide(2, 1);
   // Zoom to the physically populated region (the grids run to 800/400 MeV but
   // the strength is concentrated below ~500/150 MeV).
   for (TH2D* h : {hp, hn}) {
      h->GetXaxis()->SetRangeUser(0, 500);
      h->GetYaxis()->SetRangeUser(0, 150);
      h->GetZaxis()->SetTitleOffset(1.1);
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
   pdk_label("JLab argon SF (NuWro grid)");
   cv->SaveAs("plots/spectral_pe.png");
   printf("Wrote plots/spectral_pe.png\n");
}
