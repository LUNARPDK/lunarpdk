// Off-shell invariant mass W of the decaying nucleon and the kinematically
// forbidden fraction, across nuclear models.
//
// A bound nucleon of momentum p and removal energy E_rem has nuclear-frame
// energy M_N - E_rem, hence an off-shell invariant mass
//   W = sqrt((M_N - E_rem)^2 - p^2),
// which Fermi motion and binding push below the free mass M_N. When W falls
// below m1 + m2 the two-body decay is forbidden and the draw is resampled. This
// macro (left) overlays the W distribution per momentum model, computed from the
// existing legacy (FSI-off) columns (event nucleon_p d1_p d2_p e_rem), and
// (right) shows the per-model forbidden fraction read from the generator's
// structured stderr summary (report/forbidden_log.txt, lines beginning
// "# forbidden ... model=<key> ... frac=<f>").
//
// Run from the project root:  root -l -b -q macros/offshell_W.C
// Writes plots/offshell_W.png.
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "pdk_style.h"

static const double Mp = 0.93827;        // proton mass [GeV]
static const double mKch = 0.49368;      // K+ mass [GeV] (pToKnu threshold)

// Fill a W histogram from a legacy FSI-off file (W from nucleon_p, e_rem).
static long fill_W(const std::string& fname, TH1* h) {
   std::ifstream in(fname);
   if (!in.is_open()) return 0;
   int ev;
   double np, d1, d2, er;
   long n = 0;
   while (in >> ev >> np >> d1 >> d2 >> er) {
      double E = Mp - er;
      double W2 = E * E - np * np;
      if (W2 > 0) h->Fill(std::sqrt(W2));
      ++n;
   }
   return n;
}

// Parse report/forbidden_log.txt -> map model key -> forbidden fraction (for the
// pToKnu sweep). Lines: "# forbidden channel=pToKnu model=gfg ... frac=0.0123".
static std::map<std::string, double> read_forbidden(const std::string& path,
                                                    const std::string& chan) {
   std::map<std::string, double> out;
   std::ifstream in(path);
   std::string line;
   while (std::getline(in, line)) {
      if (line.rfind("# forbidden", 0) != 0) continue;
      std::string ch, mod;
      double frac = -1;
      std::istringstream ss(line);
      std::string tok;
      while (ss >> tok) {
         if (tok.rfind("channel=", 0) == 0) ch = tok.substr(8);
         else if (tok.rfind("model=", 0) == 0) mod = tok.substr(6);
         else if (tok.rfind("frac=", 0) == 0) frac = atof(tok.c_str() + 5);
      }
      if (ch == chan && frac >= 0) out[mod] = frac;  // last wins
   }
   return out;
}

void offshell_W() {
   pdk_set_style();

   const char* models[] = {"gfg", "lfg", "src",    "sf",     "hosm",
                           "br",  "gauss", "cfg",  "benhar", "ankowski"};
   const int nmod = sizeof(models) / sizeof(models[0]);

   TCanvas* cv = new TCanvas("c", "off-shell W", 1500, 600);
   cv->Divide(2, 1);

   // ---- left: W distribution per model (pToKnu) ---------------------------
   cv->cd(1);
   gPad->SetLogy();
   TLegend* leg = new TLegend(0.15, 0.55, 0.45, 0.90);
   leg->SetTextSize(0.028);
   double ymax = 0;
   TH1F* hs[16];
   int nh = 0;
   for (int m = 0; m < nmod; ++m) {
      TH1F* h = new TH1F(Form("hW_%s", models[m]), "", 132, 0.45, 1.0);
      long n = fill_W(Form("data/win_pToKnu_%s.txt", models[m]), h);
      if (n == 0) { delete h; continue; }
      h->Scale(1.0 / h->Integral());
      h->SetLineColor(pdk_model_color(models[m]));
      h->SetLineWidth(2);
      ymax = TMath::Max(ymax, h->GetMaximum());
      hs[nh++] = h;
      leg->AddEntry(h, pdk_model_label(models[m]), "l");
   }
   if (nh == 0) {
      printf("No data/win_pToKnu_*.txt files; run make_plots.sh first.\n");
      return;
   }
   hs[0]->SetTitle("Off-shell nucleon mass W (p#rightarrowK^{+}#bar{#nu});"
                   "W [GeV];area-normalised");
   hs[0]->SetMaximum(2.0 * ymax);
   hs[0]->SetMinimum(1e-5);
   hs[0]->Draw("hist");
   for (int i = 1; i < nh; ++i) hs[i]->Draw("hist same");
   // free proton mass and the K+ threshold (W must exceed m_K).
   TLine* lMp = new TLine(Mp, 1e-5, Mp, 2.0 * ymax);
   lMp->SetLineStyle(2); lMp->SetLineColor(kGray + 2); lMp->Draw();
   TLine* lthr = new TLine(mKch, 1e-5, mKch, 2.0 * ymax);
   lthr->SetLineStyle(3); lthr->SetLineColor(kRed + 1); lthr->Draw();
   leg->Draw();
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.026); tx.SetTextColor(kGray + 2);
   tx.DrawLatex(0.70, 0.45, "M_{p}");
   tx.SetTextColor(kRed + 1);
   tx.DrawLatex(0.17, 0.45, "m_{K} (threshold)");

   // ---- right: forbidden fraction per model -------------------------------
   cv->cd(2);
   gPad->SetBottomMargin(0.16);
   auto fb = read_forbidden("report/forbidden_log.txt", "pToKnu");
   TH1F* hf = new TH1F("hf", "Kinematically forbidden fraction (p#rightarrowK^{+}#bar{#nu});;"
                             "resampled fraction", nmod, 0, nmod);
   for (int m = 0; m < nmod; ++m) {
      double v = fb.count(models[m]) ? fb[models[m]] : 0.0;
      hf->SetBinContent(m + 1, v);
      hf->GetXaxis()->SetBinLabel(m + 1, pdk_model_label(models[m]));
   }
   hf->SetFillColorAlpha(kAzure + 2, 0.55);
   hf->SetLineColor(kAzure + 2);
   hf->SetBarWidth(0.7);
   hf->SetBarOffset(0.15);
   hf->GetXaxis()->LabelsOption("v");
   hf->GetXaxis()->SetLabelSize(0.034);
   hf->SetMinimum(0);
   hf->Draw("bar");
   if (fb.empty())
      printf("Note: report/forbidden_log.txt empty/missing; right panel is blank.\n");

   cv->cd(0);
   pdk_label("FSI off  #bullet  binding = optical potential");
   cv->SaveAs("plots/offshell_W.png");
   printf("Wrote plots/offshell_W.png\n");
}
