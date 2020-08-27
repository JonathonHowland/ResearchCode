#include <iostream>
#include "TH1D.h"
#include "TRandom.h"
#include "TCanvas.h"

void Gauss(){
  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
  TH1D* h1 = new TH1D("h1", "h1", 120, -6., 6.);
  TRandom* Rndm = new TRandom();
  for(int i=0; i<1000000; i++){
    double sum = 0;
    for(int j=0; j<12; j++){
      double k = Rndm->Rndm();
      sum = sum+k;
    }
    double gauss = sum - 6;
    h1->Fill(gauss);
  }
  c1->cd();
  h1->Draw();
}
