#ifndef PDK_STYLE_H
#define PDK_STYLE_H
//
// Shared plotting style for all PDK macros. Including this header and calling
// pdk_set_style() once at the top of a macro gives every figure the same fonts,
// margins, ticks, line widths and palette. The colour/label helpers guarantee a
// given nuclear model (or binding model) keeps the same colour and name across
// every plot it appears in.
//
#include <map>
#include <string>

#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

// ---- global look -----------------------------------------------------------
inline void pdk_set_style() {
   TStyle* s = new TStyle("PDK", "PDK shared style");
   const int font = 42;  // Helvetica

   s->SetCanvasColor(0);
   s->SetPadColor(0);
   s->SetFrameFillColor(0);
   s->SetFrameBorderMode(0);
   s->SetCanvasBorderMode(0);
   s->SetPadBorderMode(0);

   s->SetPadTopMargin(0.08);
   s->SetPadRightMargin(0.05);
   s->SetPadBottomMargin(0.13);
   s->SetPadLeftMargin(0.13);

   s->SetTextFont(font);
   s->SetLabelFont(font, "xyz");
   s->SetTitleFont(font, "xyz");
   s->SetTitleFont(font, "");  // pad title
   s->SetLegendFont(font);
   s->SetStatFont(font);

   s->SetTitleSize(0.045, "xyz");
   s->SetLabelSize(0.040, "xyz");
   s->SetTitleSize(0.050, "");  // pad title
   s->SetTitleOffset(1.20, "x");
   s->SetTitleOffset(1.35, "y");

   s->SetHistLineWidth(2);
   s->SetFrameLineWidth(1);
   s->SetPadTickX(1);
   s->SetPadTickY(1);
   s->SetNdivisions(508, "xyz");

   s->SetOptTitle(1);
   s->SetTitleBorderSize(0);
   s->SetTitleFillColor(0);
   s->SetTitleAlign(23);
   s->SetTitleX(0.5);

   s->SetOptStat(0);  // overlays use a legend; single plots re-enable below
   s->SetStatBorderSize(1);
   s->SetStatColor(0);
   s->SetStatFont(font);
   s->SetStatX(0.95);
   s->SetStatY(0.92);
   s->SetStatW(0.18);
   s->SetStatH(0.14);

   s->SetLegendBorderSize(0);
   s->SetLegendFillColor(0);
   s->SetPalette(kBird);

   gROOT->SetStyle("PDK");
   gROOT->ForceStyle();
}

// ---- consistent per-model colour and label ---------------------------------
inline Color_t pdk_model_color(const std::string& key) {
   static const std::map<std::string, Color_t> m = {
       {"gfg", kBlue + 1},      {"lfg", kGreen + 2},  {"src", kRed + 1},
       {"sf", kMagenta + 1},    {"benhar", kOrange + 7}, {"ankowski", kCyan + 2},
       {"hosm", kViolet + 1},   {"br", kSpring + 4},  {"gauss", kGray + 2},
       {"cfg", kPink + 7},
   };
   auto it = m.find(key);
   return it == m.end() ? kBlack : it->second;
}

inline const char* pdk_model_label(const std::string& key) {
   static const std::map<std::string, std::string> m = {
       {"gfg", "Global Fermi gas"}, {"lfg", "Local Fermi gas"},
       {"src", "Short-range corr."}, {"sf", "Spectral function"},
       {"benhar", "Benhar SF"},     {"ankowski", "Ankowski SF"},
       {"hosm", "HO shell model"},  {"br", "Bodek-Ritchie"},
       {"gauss", "Gaussian"},       {"cfg", "Correlated FG"},
   };
   auto it = m.find(key);
   return it == m.end() ? "unknown" : it->second.c_str();
}

// ---- consistent per-binding colour and label -------------------------------
inline Color_t pdk_binding_color(const std::string& key) {
   static const std::map<std::string, Color_t> m = {
       {"potential", kAzure + 1}, {"constant", kOrange + 7}, {"shell", kGreen + 2},
   };
   auto it = m.find(key);
   return it == m.end() ? kBlack : it->second;
}

inline const char* pdk_binding_label(const std::string& key) {
   static const std::map<std::string, std::string> m = {
       {"potential", "Optical potential"}, {"constant", "Constant 30 MeV"},
       {"shell", "Shell levels"},
   };
   auto it = m.find(key);
   return it == m.end() ? "unknown" : it->second.c_str();
}

// ---- per-model fill colour for the single-distribution histograms -----------
inline void pdk_style_hist(TH1* h, Color_t col) {
   h->SetLineColor(col);
   h->SetLineWidth(2);
   h->SetFillColorAlpha(col, 0.30);
}

// ---- a consistent corner tag drawn on every figure -------------------------
inline void pdk_label(const char* extra = "") {
   TLatex* t = new TLatex();
   t->SetNDC();
   t->SetTextFont(42);
   t->SetTextSize(0.030);
   t->SetTextColor(kGray + 2);
   t->SetTextAlign(31);  // bottom-right
   TString s = "#bf{LUNAR} toy MC";
   if (extra && extra[0]) s += TString::Format("  #bullet  %s", extra);
   t->DrawLatex(0.95, 0.16, s);
}

#endif  // PDK_STYLE_H
