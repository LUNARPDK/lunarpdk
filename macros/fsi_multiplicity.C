// Final-state multiplicity after the FSI cascade.
//
// The cascade knocks nucleons out of the residual nucleus and can produce or
// charge-exchange extra pions, so the observable final state is richer than the
// single decay meson. For each channel this macro tallies, per generated decay,
// the mean number of escaping protons, neutrons, charged/neutral pions and
// kaons, and draws (left) a stacked bar of the mean multiplicity per channel and
// (right) the distribution of the knocked-out nucleon multiplicity for two
// benchmark channels (transparent K+ vs strongly-interacting pi0).
//
// Reads the FSI-on per-particle files data/fsi_<channel>_on.txt written by
// make_plots.sh. Leptons (PDG -11/-13) are excluded.
// Run from the project root:  root -l -b -q macros/fsi_multiplicity.C
// Writes plots/fsi_multiplicity.png and report/fsi_multiplicity_summary.txt.
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "pdk_style.h"

struct Chan {
   const char* key;
   const char* label;
};

// Species bins for the stacked bar (order = bottom -> top).
enum { kProt, kNeut, kPiC, kPi0, kKaon, kNspec };
static const char* kSpecLabel[kNspec] = {"protons", "neutrons", "#pi^{#pm}",
                                         "#pi^{0}", "kaons"};
static Color_t spec_color(int i) {
   const Color_t c[kNspec] = {kRed + 1, kAzure + 2, kOrange + 7, kViolet + 1,
                              kGreen + 2};
   return c[i];
}

static int species_of(int pdg) {
   switch (pdg) {
      case 2212: return kProt;
      case 2112: return kNeut;
      case 211:
      case -211: return kPiC;
      case 111: return kPi0;
      case 321:
      case -321:
      case 311:
      case -311:
      case 310:
      case 130: return kKaon;
      default: return -1;  // eta, Lambda, photon, leptons: not stacked
   }
}

// Accumulate per-decay mean species multiplicity; optionally fill a nucleon
// (p+n) knockout-multiplicity histogram. Returns the number of events.
static long read_mult(const std::string& fname, double mean[kNspec],
                      TH1* h_nucl) {
   for (int i = 0; i < kNspec; ++i) mean[i] = 0.0;
   std::ifstream in(fname);
   if (!in.is_open()) {
      printf("Warning: could not open %s (run make_plots.sh first)\n",
             fname.c_str());
      return 0;
   }
   long counts[kNspec] = {0, 0, 0, 0, 0};
   long n = 0;
   int cur_nucl = 0;
   bool have_event = false;
   std::string line;
   auto flush = [&]() {
      if (have_event && h_nucl) h_nucl->Fill(cur_nucl);
   };
   while (std::getline(in, line)) {
      if (line.rfind("# event", 0) == 0) {
         flush();
         ++n;
         have_event = true;
         cur_nucl = 0;
         continue;
      }
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      int ev, pdg;
      double px, py, pz, E;
      if (!(ss >> ev >> pdg >> px >> py >> pz >> E)) continue;
      int s = species_of(pdg);
      if (s >= 0) ++counts[s];
      if (pdg == 2212 || pdg == 2112) ++cur_nucl;
   }
   flush();  // last event
   if (n > 0)
      for (int i = 0; i < kNspec; ++i) mean[i] = double(counts[i]) / n;
   return n;
}

void fsi_multiplicity() {
   pdk_set_style();

   const Chan chans[] = {
       {"pToKnu", "K^{+}#bar{#nu}"},   {"pToMuK0", "#mu^{+}K^{0}"},
       {"pToEK0", "e^{+}K^{0}"},       {"pToEEta", "e^{+}#eta"},
       {"pToMuEta", "#mu^{+}#eta"},    {"pToNuPip", "#bar{#nu}#pi^{+}"},
       {"pToEPi0", "e^{+}#pi^{0}"},    {"pToMuPi0", "#mu^{+}#pi^{0}"},
       {"nToNuK0", "#bar{#nu}K^{0}"},  {"nToEKm", "e^{+}K^{-}"},
       {"nToNuEta", "#bar{#nu}#eta"},  {"nToEPim", "e^{+}#pi^{-}"},
       {"nToMuPim", "#mu^{+}#pi^{-}"}, {"nToNuPi0", "#bar{#nu}#pi^{0}"},
   };
   const int nchan = sizeof(chans) / sizeof(chans[0]);

   TH1F* h[kNspec];
   for (int s = 0; s < kNspec; ++s) {
      h[s] = new TH1F(Form("h_spec_%d", s), "", nchan, 0.0, nchan);
      h[s]->SetFillColor(spec_color(s));
      h[s]->SetLineColor(kBlack);
      h[s]->SetLineWidth(1);
   }
   // Nucleon-knockout multiplicity for two benchmark channels.
   TH1F* hk_kp = new TH1F("hk_kp", "", 8, -0.5, 7.5);
   TH1F* hk_pi = new TH1F("hk_pi", "", 8, -0.5, 7.5);

   FILE* sum = fopen("report/fsi_multiplicity_summary.txt", "w");
   if (sum) {
      fprintf(sum, "# Mean post-FSI final-state multiplicity per decay (FSI on)\n");
      fprintf(sum, "# %-9s %8s %8s %8s %8s %8s %10s\n", "channel", "protons",
              "neutrons", "pi_pm", "pi0", "kaons", "n_events");
   }

   for (int c = 0; c < nchan; ++c) {
      double mean[kNspec];
      TH1* target = (std::string(chans[c].key) == "pToKnu")   ? hk_kp
                    : (std::string(chans[c].key) == "pToEPi0") ? hk_pi
                                                               : nullptr;
      long n = read_mult(Form("data/fsi_%s_on.txt", chans[c].key), mean, target);
      for (int s = 0; s < kNspec; ++s) {
         h[s]->SetBinContent(c + 1, mean[s]);
         h[s]->GetXaxis()->SetBinLabel(c + 1, chans[c].label);
      }
      if (sum)
         fprintf(sum, "  %-9s %8.4f %8.4f %8.4f %8.4f %8.4f %10ld\n",
                 chans[c].key, mean[kProt], mean[kNeut], mean[kPiC], mean[kPi0],
                 mean[kKaon], n);
   }
   if (sum) fclose(sum);

   TCanvas* cv = new TCanvas("c", "FSI multiplicity", 1500, 600);
   cv->Divide(2, 1);

   cv->cd(1);
   gPad->SetBottomMargin(0.14);
   THStack* st = new THStack("st",
                             "Mean post-FSI final-state multiplicity by channel;;"
                             "mean particles per decay");
   for (int s = 0; s < kNspec; ++s) st->Add(h[s]);
   st->Draw("hist bar2");
   st->GetXaxis()->SetLabelSize(0.034);
   st->GetYaxis()->SetTitleOffset(0.95);
   TLegend* leg = new TLegend(0.15, 0.70, 0.55, 0.90);
   leg->SetNColumns(2);
   leg->SetTextSize(0.030);
   for (int s = 0; s < kNspec; ++s) leg->AddEntry(h[s], kSpecLabel[s], "f");
   leg->Draw();

   cv->cd(2);
   hk_pi->SetTitle("Knocked-out nucleon multiplicity;N_{p}+N_{n} escaping;"
                   "fraction of decays");
   if (hk_kp->Integral() > 0) hk_kp->Scale(1.0 / hk_kp->Integral());
   if (hk_pi->Integral() > 0) hk_pi->Scale(1.0 / hk_pi->Integral());
   hk_pi->SetLineColor(kOrange + 7);
   hk_pi->SetFillColorAlpha(kOrange + 7, 0.30);
   hk_kp->SetLineColor(kGreen + 2);
   hk_kp->SetFillColorAlpha(kGreen + 2, 0.30);
   hk_pi->SetMaximum(1.15 * TMath::Max(hk_pi->GetMaximum(), hk_kp->GetMaximum()));
   hk_pi->Draw("hist");
   hk_kp->Draw("hist same");
   TLegend* l2 = new TLegend(0.45, 0.74, 0.92, 0.90);
   l2->AddEntry(hk_kp, "p#rightarrowK^{+}#bar{#nu}", "f");
   l2->AddEntry(hk_pi, "p#rightarrowe^{+}#pi^{0}", "f");
   l2->Draw();

   cv->cd(0);
   pdk_label("argon, FSI on");
   cv->SaveAs("plots/fsi_multiplicity.png");
   printf("Wrote plots/fsi_multiplicity.png and report/fsi_multiplicity_summary.txt\n");
}
