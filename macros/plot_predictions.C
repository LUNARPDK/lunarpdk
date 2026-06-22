// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Predicted number of nucleon-decay events in DUNE per channel.
//
// Reads report/event_predictions.txt (written by build/EventPredictor) and draws,
// for every channel with a current experimental lifetime limit, the expected
// count in 400 kt.yr at tau = the SuperK 90% C.L. limit:
//   * filled marker + horizontal error bar = central count with the nuclear-model
//     band (spread of the signal-window efficiency across the ten models);
//   * open marker = the FSI-unfolded count (N_cen / eps_FSI), so the gap between
//     the two shows the cascade's suppression. For the two DUNE-documented modes
//     (p->K+ nu, p->e+ pi0) the published DUNE efficiency already includes FSI;
//     unfolding eps_FSI keeps that suppression visible (e.g. the p->e+ pi0 gap),
//     rather than hiding it inside eps_det.
// A dashed line at N = 1 marks the one-event threshold. Log-x, since the counts
// span several orders of magnitude.
//
// Run from the project root:  root -l -b -q macros/plot_predictions.C
// Writes plots/event_predictions.png.
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "pdk_style.h"

// Pretty TLatex label per channel key.
static const char* chan_label(const std::string& k) {
   if (k == "pToKnu") return "p#rightarrowK^{+}#bar{#nu}";
   if (k == "pToEPi0") return "p#rightarrowe^{+}#pi^{0}";
   if (k == "pToMuPi0") return "p#rightarrow#mu^{+}#pi^{0}";
   if (k == "pToNuPip") return "p#rightarrow#bar{#nu}#pi^{+}";
   if (k == "pToEEta") return "p#rightarrowe^{+}#eta";
   if (k == "pToMuEta") return "p#rightarrow#mu^{+}#eta";
   if (k == "pToEK0") return "p#rightarrowe^{+}K^{0}";
   if (k == "pToMuK0") return "p#rightarrow#mu^{+}K^{0}";
   if (k == "nToEPim") return "n#rightarrowe^{+}#pi^{-}";
   if (k == "nToMuPim") return "n#rightarrow#mu^{+}#pi^{-}";
   if (k == "nToNuPi0") return "n#rightarrow#bar{#nu}#pi^{0}";
   if (k == "nToNuEta") return "n#rightarrow#bar{#nu}#eta";
   if (k == "nToNuK0") return "n#rightarrow#bar{#nu}K^{0}";
   if (k == "nToEKm") return "n#rightarrowe^{+}K^{-}";
   return k.c_str();
}

void plot_predictions() {
   pdk_set_style();

   std::ifstream in("report/event_predictions.txt");
   if (!in.is_open()) {
      printf("Error: report/event_predictions.txt not found (run EventPredictor).\n");
      return;
   }

   std::vector<std::string> keys;
   std::vector<double> ncen, nlo, nhi, noff;
   std::string line;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string key, parent;
      double tau, edet, efsi, wlo, whi, nc, nl, nh, no;
      if (!(ss >> key >> parent >> tau >> edet >> efsi >> wlo >> whi >> nc >>
            nl >> nh >> no))
         continue;
      if (nc <= 0) continue;  // no limit -> skip
      keys.push_back(key);
      ncen.push_back(nc);
      nlo.push_back(nl);
      nhi.push_back(nh);
      noff.push_back(no);
   }
   const int n = keys.size();
   if (n == 0) {
      printf("Error: no channels with predictions found.\n");
      return;
   }

   // x-range across all positive values.
   double xmin = 1e30, xmax = -1e30;
   for (int i = 0; i < n; ++i) {
      for (double v : {nlo[i], ncen[i], nhi[i], noff[i]})
         if (v > 0) { xmin = TMath::Min(xmin, v); xmax = TMath::Max(xmax, v); }
   }
   xmin *= 0.4;
   xmax *= 2.5;

   // Frame: log-x, channel modes as y labels (top = first row).
   TH2D* fr = new TH2D("fr",
                       "DUNE nucleon-decay predictions (400 kt#upoint yr, "
                       "#tau = SuperK limit);expected events;",
                       100, xmin, xmax, n, 0, n);
   for (int i = 0; i < n; ++i)
      fr->GetYaxis()->SetBinLabel(n - i, chan_label(keys[i]));  // first row on top
   fr->GetYaxis()->SetLabelSize(0.044);
   fr->GetXaxis()->SetTitleOffset(1.10);
   fr->SetStats(0);

   TGraphAsymmErrors* g_on = new TGraphAsymmErrors(n);  // FSI on, with model band
   TGraphAsymmErrors* g_off = new TGraphAsymmErrors(n);  // FSI off
   for (int i = 0; i < n; ++i) {
      double y = n - i - 0.5;  // first row near top
      g_on->SetPoint(i, ncen[i], y);
      g_on->SetPointError(i, ncen[i] - nlo[i], nhi[i] - ncen[i], 0.0, 0.0);
      g_off->SetPoint(i, noff[i], y);
      g_off->SetPointError(i, 0.0, 0.0, 0.0, 0.0);
   }
   g_on->SetMarkerStyle(20);
   g_on->SetMarkerSize(1.2);
   g_on->SetMarkerColor(kAzure + 2);
   g_on->SetLineColor(kAzure + 2);
   g_on->SetLineWidth(2);
   g_off->SetMarkerStyle(24);
   g_off->SetMarkerSize(1.2);
   g_off->SetMarkerColor(kRed + 1);
   g_off->SetLineColor(kRed + 1);

   TCanvas* cv = new TCanvas("c", "predictions", 1000, 680);
   cv->SetLogx();
   cv->SetLeftMargin(0.20);
   cv->SetGridx();
   fr->Draw("AXIS");
   g_on->Draw("P SAME");
   g_off->Draw("P SAME");

   // One-event threshold.
   TLine* one = new TLine(1.0, 0.0, 1.0, n);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->SetLineWidth(2);
   one->Draw();

   TLegend* leg = new TLegend(0.58, 0.78, 0.93, 0.90);
   leg->AddEntry(g_on, "FSI on (band = nucl. models)", "pl");
   leg->AddEntry(g_off, "cascade off", "p");
   leg->AddEntry(one, "1 event", "l");
   leg->Draw();

   pdk_label("2#times10 kt LAr  #bullet  #tau at current limit");
   cv->SaveAs("plots/event_predictions.png");
   printf("Wrote plots/event_predictions.png\n");
}
