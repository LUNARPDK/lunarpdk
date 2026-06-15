// Compare the lab-frame kaon momentum spectrum across the three binding
// (removal-energy) models, for two representative nucleon-momentum models.
// Reads data/kbind_<model>_<binding>.txt (columns: event nucleon_p d1_p d2_p e_rem).
//
// The binding model sets the proton removal energy, hence the off-shell mass
// W < M_p, hence the kaon momentum: deeper binding -> smaller W -> softer kaon.
//
// Run in batch from the project root:  root -l -b -q macros/compare_kaon_binding.C
// Writes plots/kaon_binding.png and report/kaon_binding_summary.txt.
#include "pdk_style.h"

FILE* g_rep = nullptr;
void emit(TString s) {
   printf("%s", s.Data());
   if (g_rep) fprintf(g_rep, "%s", s.Data());
}

void draw_panel(const char* model, const char* title, double pK_free) {
   const int nb = 3;
   const char* bindings[nb] = {"potential", "constant", "shell"};

   TLegend* leg = new TLegend(0.15, 0.70, 0.55, 0.90);
   leg->SetHeader(title);

   emit(TString::Format("\n[%s]\n", title));
   emit(TString::Format("  %-18s  %8s  %8s  %8s  %9s\n", "binding", "mean",
                        "RMS", "peak", "in-win%"));

   double ymax = 0.0;
   TH1F* hists[nb] = {nullptr};
   for (int i = 0; i < nb; ++i) {
      TString fname = TString::Format("data/kbind_%s_%s.txt", model, bindings[i]);
      std::ifstream in(fname.Data());
      if (!in.is_open()) { printf("Warning: missing %s\n", fname.Data()); continue; }
      TH1F* h = new TH1F(TString::Format("h_%s_%s", model, bindings[i]),
                         ";Kaon momentum p_{K} [GeV/c];Probability density",
                         120, 0.0, 0.6);
      int ev; float np, d1, d2, er; long n = 0, in_win = 0;
      while (in >> ev >> np >> d1 >> d2 >> er) {
         h->Fill(d2); ++n;
         if (d2 >= 0.30 && d2 <= 0.38) ++in_win;
      }
      in.close();
      emit(TString::Format("  %-18s  %8.4f  %8.4f  %8.4f  %8.2f%%\n",
                           pdk_binding_label(bindings[i]),
                           h->GetMean(), h->GetRMS(),
                           h->GetBinCenter(h->GetMaximumBin()),
                           100.0 * in_win / n));
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral(), "width");
      h->SetLineColor(pdk_binding_color(bindings[i]));
      h->SetLineWidth(2);
      hists[i] = h;
      ymax = TMath::Max(ymax, h->GetMaximum());
      leg->AddEntry(h, pdk_binding_label(bindings[i]), "l");
   }
   bool first = true;
   for (int i = 0; i < nb; ++i) {
      if (!hists[i]) continue;
      if (first) { hists[i]->SetMaximum(ymax * 1.25); hists[i]->Draw("HIST"); first = false; }
      else hists[i]->Draw("HIST SAME");
   }
   TLine* ref = new TLine(pK_free, 0.0, pK_free, ymax * 1.25);
   ref->SetLineStyle(2); ref->SetLineColor(kGray + 2); ref->Draw();
   leg->Draw();
   pdk_label();
}

void compare_kaon_binding() {
   pdk_set_style();
   const double Mp = 0.93827, mK = 0.49368;
   const double pK_free = (Mp * Mp - mK * mK) / (2.0 * Mp);

   g_rep = fopen("report/kaon_binding_summary.txt", "w");
   emit("Kaon momentum vs binding model (free-decay p_K = 0.339 GeV/c, "
        "window [0.30,0.38])\n");

   TCanvas* c = new TCanvas("c", "kaon spectra vs binding", 1200, 500);
   c->Divide(2, 1);
   c->cd(1); draw_panel("lfg", "Local Fermi gas", pK_free);
   c->cd(2); draw_panel("sf", "Spectral function (toy)", pK_free);

   if (g_rep) fclose(g_rep);
   c->SaveAs("plots/kaon_binding.png");
   printf("\nWrote plots/kaon_binding.png and report/kaon_binding_summary.txt\n");
}
