// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Analyze the lab-frame hadron-daughter momentum distribution (the kaon, for the
// default p -> K+ nu channel) for every nucleon-momentum model. Reads
// data/kaon_<model>.txt (columns: event nucleon_p d1_p d2_p e_rem) produced by
// analyze_kaon.sh.
//
// Run in batch from the project root:  root -l -b -q macros/plot_kaon_models.C
// Writes plots/kaon_models.png (overlay) and report/kaon_summary.txt (table).
//
// Physics reference: a free proton at rest decaying p -> K+ nu emits a
// monochromatic kaon at p_K = (M_p^2 - m_K^2)/(2 M_p). Fermi motion broadens it
// and the binding (off-shell W < M_p) shifts the peak downward; high-momentum
// (SRC) tails feed a low-p shoulder. The signal-window fraction below quantifies
// how much of the rate stays near the free-decay line.
#include "pdk_style.h"

void plot_kaon_models() {
   pdk_set_style();

   const double Mp = 0.93827, mK = 0.49368;
   const double pK_free = (Mp * Mp - mK * mK) / (2.0 * Mp);  // ~0.339 GeV/c
   const double win_lo = 0.30, win_hi = 0.38;                // nominal signal window

   const int nm = 10;
   const char* models[nm] = {"gfg", "lfg", "src",  "sf",     "hosm",
                             "br",  "gauss", "cfg", "benhar", "ankowski"};

   TCanvas* c = new TCanvas("c", "kaon momentum models", 800, 600);

   TLegend* leg = new TLegend(0.15, 0.50, 0.45, 0.90);

   FILE* rep = fopen("report/kaon_summary.txt", "w");
   auto emit = [&](TString s) {
      printf("%s", s.Data());
      if (rep) fprintf(rep, "%s", s.Data());
   };
   emit("Lab-frame kaon momentum from p -> K+ nu, by nuclear model\n");
   emit(TString::Format("Free-decay reference: p_K = %.4f GeV/c\n", pK_free));
   emit(TString::Format("Signal window: [%.2f, %.2f] GeV/c\n\n", win_lo, win_hi));
   emit(TString::Format("%-18s  %8s  %8s  %8s  %9s  %9s\n", "model", "mean",
                        "RMS", "peak", "in-win%", "below.25%"));
   emit(TString::Format("%-18s  %8s  %8s  %8s  %9s  %9s\n", "------------------",
                        "--------", "--------", "--------", "---------",
                        "---------"));

   double ymax = 0.0;
   TH1F* hists[nm] = {nullptr};
   for (int i = 0; i < nm; ++i) {
      TString fname = TString::Format("data/kaon_%s.txt", models[i]);
      std::ifstream in(fname.Data());
      if (!in.is_open()) {
         printf("Warning: could not open %s (run analyze_kaon.sh first)\n",
                fname.Data());
         continue;
      }
      TH1F* h = new TH1F(TString::Format("hk_%s", models[i]),
                         "Lab-frame kaon momentum by nuclear model;"
                         "Kaon momentum p_{K} [GeV/c];Probability density",
                         120, 0.0, 0.6);
      int ev; float np, d1, d2, er;
      long n = 0, in_win = 0, below = 0;
      while (in >> ev >> np >> d1 >> d2 >> er) {
         h->Fill(d2);  // hadron-side daughter (kaon for the default channel)
         ++n;
         if (d2 >= win_lo && d2 <= win_hi) ++in_win;
         if (d2 < 0.25) ++below;
      }
      in.close();

      emit(TString::Format("%-18s  %8.4f  %8.4f  %8.4f  %8.2f%%  %8.2f%%\n",
                           pdk_model_label(models[i]), h->GetMean(), h->GetRMS(),
                           h->GetBinCenter(h->GetMaximumBin()),
                           100.0 * in_win / n, 100.0 * below / n));

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
         hists[i]->SetMaximum(ymax * 1.25);
         hists[i]->Draw("HIST");
         first = false;
      } else {
         hists[i]->Draw("HIST SAME");
      }
   }

   // Free-decay reference line at p_K ~ 0.339 GeV/c.
   TLine* ref = new TLine(pK_free, 0.0, pK_free, ymax * 1.25);
   ref->SetLineStyle(2);
   ref->SetLineColor(kGray + 2);
   ref->Draw();
   TLatex tx;
   tx.SetTextFont(42);
   tx.SetTextColor(kGray + 2);
   tx.SetTextSize(0.028);
   tx.SetTextAngle(90);
   tx.DrawLatex(pK_free + 0.006, ymax * 0.35, "free decay  p_{K} = 0.339");

   leg->Draw();
   pdk_label("kaon momentum");
   c->SaveAs("plots/kaon_models.png");
   printf("\nWrote plots/kaon_models.png and report/kaon_summary.txt\n");
}
