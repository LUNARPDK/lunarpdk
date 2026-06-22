// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Overlay the nucleon momentum distribution for the different nuclear models.
// Expects one ascii file per model at data/proton_<model>.txt (produced by
// compare_models.sh), each with columns: event nucleon_p d1_p d2_p e_rem.
// Run in batch from the project root:  root -l -b -q macros/plot_models.C
//
// Writes plots/proton_models.png (area-normalized, log-y to show the SRC tail)
// and report/momentum_summary.txt (per-model mean/RMS/peak and tail fraction).
#include "pdk_style.h"

void plot_models() {
   pdk_set_style();

   const double kF = 0.217;  // argon proton Fermi momentum [GeV/c]

   const int nm = 10;
   const char* models[nm] = {"gfg", "lfg", "src",  "sf",     "hosm",
                             "br",  "gauss", "cfg", "benhar", "ankowski"};

   TCanvas *c = new TCanvas("c", "proton momentum models", 800, 600);
   c->SetLogy();

//   TLegend *leg = new TLegend(0.58, 0.50, 0.93, 0.90);
   TLegend *leg = new TLegend(0.68, 0.50, 0.95, 0.90);

   FILE* rep = fopen("report/momentum_summary.txt", "w");
   auto emit = [&](TString s) {
      printf("%s", s.Data());
      if (rep) fprintf(rep, "%s", s.Data());
   };
   emit("Initial nucleon momentum distribution by nuclear model\n");
   emit(TString::Format("Argon proton Fermi momentum: k_F = %.3f GeV/c\n\n", kF));
   emit(TString::Format("%-18s  %8s  %8s  %8s  %9s\n", "model", "mean", "RMS",
                        "peak", "above-kF%"));
   emit(TString::Format("%-18s  %8s  %8s  %8s  %9s\n", "------------------",
                        "--------", "--------", "--------", "---------"));

   double ymax = 0.0;
   TH1F *hists[nm] = {nullptr};
   for (int i = 0; i < nm; ++i) {
      TString fname = TString::Format("data/proton_%s.txt", models[i]);
      std::ifstream in(fname.Data());
      if (!in.is_open()) {
         printf("Warning: could not open %s (run compare_models.sh first)\n",
                fname.Data());
         continue;
      }
      TH1F *h = new TH1F(TString::Format("h_%s", models[i]),
                         "Nucleon momentum by nuclear model;"
                         "Nucleon momentum p_{N} [GeV/c];Probability density",
                         100, 0.0, 0.85);
      int event; float np, d1, d2, er;
      long n = 0, above = 0;
      while (in >> event >> np >> d1 >> d2 >> er) {
         h->Fill(np);
         ++n;
         if (np > kF) ++above;
      }
      in.close();

      emit(TString::Format("%-18s  %8.4f  %8.4f  %8.4f  %8.2f%%\n",
                           pdk_model_label(models[i]), h->GetMean(), h->GetRMS(),
                           h->GetBinCenter(h->GetMaximumBin()),
                           n ? 100.0 * above / n : 0.0));

      if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
      h->SetLineColor(pdk_model_color(models[i]));
      h->SetLineWidth(2);
      hists[i] = h;
      ymax = TMath::Max(ymax, h->GetMaximum());
      leg->AddEntry(h, pdk_model_label(models[i]), "l");
   }
   if (rep) fclose(rep);

   bool first = true;
   for (int i = 0; i < nm; ++i) {
      if (!hists[i]) continue;
      if (first) {
         hists[i]->SetMaximum(ymax * 3.0);
         hists[i]->SetMinimum(1e-3);
         hists[i]->Draw("HIST");
         first = false;
      } else {
         hists[i]->Draw("HIST SAME");
      }
   }
   leg->Draw();
   pdk_label("proton momentum");

   c->SaveAs("plots/proton_models.png");
   printf("\nWrote plots/proton_models.png and report/momentum_summary.txt\n");
}
