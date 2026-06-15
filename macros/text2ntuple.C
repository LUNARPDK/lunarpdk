// Alternative reader: load the generator's ascii output into a TNtuple.
// Run from the project root:  root -l -q macros/text2ntuple.C
#include "Riostream.h"

void text2ntuple() {
   std::ifstream in;
   in.open("data/output.txt");

   Int_t event;
   Float_t nucleon_p, d1_p, d2_p, e_rem;
   Int_t nlines = 0;

   TFile *f = new TFile("data/output.root", "RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple", "nucleon decay events",
                                 "event:nucleon_p:d1_p:d2_p:e_rem");

   while (in >> event >> nucleon_p >> d1_p >> d2_p >> e_rem) {
      if (nlines < 5)
         printf("event=%d, nucleon_p=%8f, d1_p=%8f, d2_p=%8f, e_rem=%8f\n",
                event, nucleon_p, d1_p, d2_p, e_rem);
      ntuple->Fill(event, nucleon_p, d1_p, d2_p, e_rem);
      nlines++;
   }
   printf(" found %d points\n", nlines);

   in.close();
   f->Write();
}
