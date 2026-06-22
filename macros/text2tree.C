// Author: Jarek Nowak <lunar_pdk@proton.me>, 2026
//
// Read the generator's ascii output into a ROOT TTree.
// Columns: event nucleon_p d1_p d2_p e_rem (d1 = lepton-side, d2 = hadron-side).
// Run from the project root:  root -l -q macros/text2tree.C
void text2tree() {
   TFile *f = TFile::Open("data/output.root", "recreate");
   TTree *T = new TTree("T", "nucleon decay events");
   T->ReadFile("data/output.txt",
               "event/I:nucleon_p/F:d1_p/F:d2_p/F:e_rem/F");
   T->Print();
   T->Write();
}
