#ifndef PDK_STYLE_H
#define PDK_STYLE_H
// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
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

   s->SetTitleSize(0.050, "xyz");
   s->SetLabelSize(0.045, "xyz");
   s->SetTitleSize(0.052, "");  // pad title
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

// Line style per model, so the ten overlaid curves stay distinguishable when the
// figure is printed in greyscale or read by a colour-blind reader: colour alone
// cannot separate ten curves that partly overlap.
//
// Only five styles are used --- solid, dashed, dotted, long-dash and
// dash-dot-dot --- because ROOT's remaining patterns are too fine to survive the
// dense polyline of a 100-bin histogram outline and end up looking alike. Each
// style is therefore shared by exactly two models, which are always given very
// different colours, so no two curves share both attributes.
inline Style_t pdk_model_style(const std::string& key) {
   static const std::map<std::string, Style_t> m = {
       {"gfg", 1},  {"br", 1},          // solid:          blue     / olive
       {"lfg", 2},  {"cfg", 2},         // dashed:         green    / pink
       {"src", 3},  {"ankowski", 3},    // dotted:         red      / cyan
       {"sf", 7},   {"gauss", 7},       // long dash:      magenta  / grey
       {"hosm", 9}, {"benhar", 9},      // dash-dot-dot:   violet   / orange
   };
   auto it = m.find(key);
   return it == m.end() ? 1 : it->second;
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

// Line style per binding model (same rationale as pdk_model_style).
inline Style_t pdk_binding_style(const std::string& key) {
   static const std::map<std::string, Style_t> m = {
       {"potential", 1}, {"constant", 2}, {"shell", 7},
   };
   auto it = m.find(key);
   return it == m.end() ? 1 : it->second;
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
// The short brand tag sits in the bottom-right margin, below the frame and clear
// of the (centred) axis title. The optional `extra` context is intentionally not
// drawn: it made the tag long enough to collide with axis titles on narrow
// multi-panel pads, and it duplicates information already in the title/caption.
inline void pdk_label(const char* /*extra*/ = "") {
   TLatex* t = new TLatex();
   t->SetNDC();
   t->SetTextFont(42);
   t->SetTextSize(0.024);
   t->SetTextColor(kGray + 2);
   t->SetTextAlign(33);  // right-top: tucked into the extreme top-right corner,
                         // clear of the (centred) title and all in-frame content
   t->DrawLatex(0.985, 0.988, "#bf{LUNAR} PDK MC");
}

#endif  // PDK_STYLE_H
