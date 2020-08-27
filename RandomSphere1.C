#include <iostream>
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
//plot phi, plot cos(phi), plot 2d histogram, play with plots
void RandomSphere1(){
  double rand;
  TRandom* Random = new TRandom();
  double theta;
  double phi;
  TH1D* hist = new TH1D("hist","hist",100, 0., 2.*acos(-1));
  TH2D* hist2 = new TH2D("hist2", "A Random Sphere Using Theta and Cos(Phi)", 100, 0., acos(-1), 100,  0., 2.*acos(-1));
  for (int i = 0; i<1000; i++){
     rand = Random->Rndm();
     theta = rand*acos(-1)*2;
     rand = Random->Rndm();
     phi = acos(2.*rand - 1.);
     // cout << theta << " " << phi << endl;
     hist2->Fill(cos(phi), theta);
     }
  gStyle->SetOptStat(0);
  hist2->Draw("SPH SURF5");
  cout << "OK" << endl;
}
