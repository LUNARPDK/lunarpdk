// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Hadron-nucleon cross sections used by the FSI cascade.
//
// This macro includes the generator's own cross-section header
// (include/PDKFsiXsec.h) and plots the very functions the cascade calls, so the
// figure is a faithful picture of the FSI input (and a validation that the
// data-driven K / eta forms reproduce the expected magnitudes): the Delta(1232)
// bump in pi N, the rise of N N inelasticity, the transparent K+ vs strongly
// absorbed K-, and the N*(1535) peak in eta N.
//
// Run from the project root:  root -l -b -q macros/plot_xsec.C
// Writes plots/fsi_xsec.png.
#include "../include/PDKFsiXsec.h"
#include "pdk_style.h"

using namespace pdk::fsi;

static TGraph* mkgraph(int n, Color_t col, int style = 1) {
   TGraph* g = new TGraph(n);
   g->SetLineColor(col);
   g->SetLineWidth(3);
   g->SetLineStyle(style);
   return g;
}

// Enlarge the axis fonts on a pad of the 2x2 grid: each pad is only half the
// canvas height, so the global label/title sizes render small once the figure is
// scaled to the column/text width. Call on the lead graph after Draw("AL").
static void bump_axes(TGraph* g) {
   for (TAxis* a : {g->GetXaxis(), g->GetYaxis()}) {
      a->SetLabelSize(0.055);
      a->SetTitleSize(0.060);
   }
   g->GetXaxis()->SetTitleOffset(1.00);
   g->GetYaxis()->SetTitleOffset(1.05);
}

void plot_xsec() {
   pdk_set_style();
   const int N = 240;

   TCanvas* cv = new TCanvas("c", "FSI cross sections", 1300, 1000);
   cv->Divide(2, 2);

   // ---- pad 1: pi-N total cross section vs pion kinetic energy -------------
   cv->cd(1);
   TGraph* g_pip_p = mkgraph(N, kRed + 1);     // pi+ p (resonant I=3/2)
   TGraph* g_pip_n = mkgraph(N, kAzure + 2);   // pi+ n
   TGraph* g_pi0_p = mkgraph(N, kGreen + 2);   // pi0 p
   for (int i = 0; i < N; ++i) {
      double T = 0.001 + 1.0 * i / (N - 1);    // GeV
      g_pip_p->SetPoint(i, T, pion_nucleon(T, kPdgPiPlus, kPdgProton).total);
      g_pip_n->SetPoint(i, T, pion_nucleon(T, kPdgPiPlus, kPdgNeutron).total);
      g_pi0_p->SetPoint(i, T, pion_nucleon(T, kPdgPi0, kPdgProton).total);
   }
   g_pip_p->SetTitle(";T_{#pi} [GeV];#sigma_{#pi N} [mb]");
   g_pip_p->SetMaximum(220);
   g_pip_p->SetMinimum(0);
   g_pip_p->Draw("AL");
   bump_axes(g_pip_p);
   g_pip_n->Draw("L SAME");
   g_pi0_p->Draw("L SAME");
   TLegend* l1 = new TLegend(0.55, 0.68, 0.92, 0.90);
   l1->AddEntry(g_pip_p, "#pi^{+}p (#Delta resonance)", "l");
   l1->AddEntry(g_pip_n, "#pi^{+}n", "l");
   l1->AddEntry(g_pi0_p, "#pi^{0}p", "l");
   l1->Draw();

   // ---- pad 2: N-N total + inelastic vs nucleon kinetic energy -------------
   cv->cd(2);
   // The two inelastic curves are dashed, and a dash pattern only survives if the
   // polyline segments are longer than the dashes themselves: sampled on the same
   // dense grid as the totals they render as solid lines, which is what made the
   // legend disagree with the plot. They are piecewise-linear interpolations of a
   // short table, so a coarse grid loses no information.
   const int Nc = N / 6;
   TGraph* g_pp = mkgraph(N, kRed + 1);           // pp/nn total
   TGraph* g_pn = mkgraph(N, kAzure + 2);         // pn total
   TGraph* g_pp_in = mkgraph(Nc, kRed + 1, 2);    // pp/nn inelastic (dashed)
   TGraph* g_pn_in = mkgraph(Nc, kAzure + 2, 2);  // pn inelastic (dashed)
   for (int i = 0; i < N; ++i) {
      double T = 0.001 + 2.0 * i / (N - 1);  // GeV
      NNXsec pp = nucleon_nucleon(T, kPdgProton, kPdgProton);
      NNXsec pn = nucleon_nucleon(T, kPdgProton, kPdgNeutron);
      g_pp->SetPoint(i, T, pp.total);
      g_pn->SetPoint(i, T, pn.total);
   }
   for (int i = 0; i < Nc; ++i) {
      double T = 0.001 + 2.0 * i / (Nc - 1);  // GeV
      NNXsec pp = nucleon_nucleon(T, kPdgProton, kPdgProton);
      NNXsec pn = nucleon_nucleon(T, kPdgProton, kPdgNeutron);
      g_pp_in->SetPoint(i, T, pp.total * pp.f_prod);
      g_pn_in->SetPoint(i, T, pn.total * pn.f_prod);
   }
   g_pp->SetTitle(";T_{N} [GeV];#sigma_{NN} [mb]");
   g_pp->SetMaximum(55);
   g_pp->SetMinimum(0);
   g_pp->Draw("AL");
   bump_axes(g_pp);
   g_pn->Draw("L SAME");
   g_pp_in->Draw("L SAME");
   g_pn_in->Draw("L SAME");
   TLegend* l2 = new TLegend(0.52, 0.18, 0.92, 0.42);
   l2->AddEntry(g_pp, "pp/nn total", "l");
   l2->AddEntry(g_pn, "pn total", "l");
   l2->AddEntry(g_pp_in, "pp/nn inelastic", "l");
   l2->AddEntry(g_pn_in, "pn inelastic", "l");
   l2->Draw();

   // ---- pad 3: K-N total vs kaon momentum ---------------------------------
   cv->cd(3);
   TGraph* g_Kpp = mkgraph(N, kGreen + 2);    // K+ p
   TGraph* g_Kpn = mkgraph(N, kSpring + 4);   // K+ n
   TGraph* g_Kmp = mkgraph(N, kRed + 1);      // K- p
   TGraph* g_Kmn = mkgraph(N, kOrange + 7);   // K- n
   for (int i = 0; i < N; ++i) {
      double p = 0.001 + 1.5 * i / (N - 1);  // GeV/c
      g_Kpp->SetPoint(i, p, kaon_nucleon(p, kPdgKPlus, kPdgProton).total);
      g_Kpn->SetPoint(i, p, kaon_nucleon(p, kPdgKPlus, kPdgNeutron).total);
      g_Kmp->SetPoint(i, p, kaon_nucleon(p, kPdgKMinus, kPdgProton).total);
      g_Kmn->SetPoint(i, p, kaon_nucleon(p, kPdgKMinus, kPdgNeutron).total);
   }
   g_Kmp->SetTitle(";p_{K} [GeV/c];#sigma_{KN} [mb]");
   g_Kmp->SetMaximum(105);
   g_Kmp->SetMinimum(0);
   g_Kmp->Draw("AL");
   bump_axes(g_Kmp);
   g_Kmn->Draw("L SAME");
   g_Kpp->Draw("L SAME");
   g_Kpn->Draw("L SAME");
   TLegend* l3 = new TLegend(0.52, 0.60, 0.92, 0.90);
   l3->AddEntry(g_Kmp, "K^{-}p (S=#minus1)", "l");
   l3->AddEntry(g_Kmn, "K^{-}n", "l");
   l3->AddEntry(g_Kpp, "K^{+}p (S=+1)", "l");
   l3->AddEntry(g_Kpn, "K^{+}n", "l");
   l3->Draw();

   // ---- pad 4: eta-N total vs eta momentum --------------------------------
   cv->cd(4);
   TGraph* g_eta = mkgraph(N, kViolet + 1);
   for (int i = 0; i < N; ++i) {
      double p = 0.001 + 1.0 * i / (N - 1);  // GeV/c
      g_eta->SetPoint(i, p, eta_nucleon(p).total);
   }
   g_eta->SetTitle(";p_{#eta} [GeV/c];#sigma_{#eta N} [mb]");
   g_eta->SetMinimum(0);
   g_eta->Draw("AL");
   bump_axes(g_eta);
   TLegend* l4 = new TLegend(0.45, 0.74, 0.92, 0.90);
   l4->AddEntry(g_eta, "#eta N  (N*(1535) Breit-Wigner)", "l");
   l4->Draw();

   cv->cd(0);
   pdk_label("FSI cross-section inputs");
   cv->SaveAs("plots/fsi_xsec.png");
   printf("Wrote plots/fsi_xsec.png\n");
}
