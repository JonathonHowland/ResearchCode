
#include <iostream>
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TStyle.h"
//Make histograms of photons and change parameters
void HiggsDecay(TLorentzVector H, TLorentzVector& P1, TLorentzVector& P2);

void RandomSphere(){
  double rand;
  double rand2;
  TRandom* Random = new TRandom();
  TH1D* hist = new TH1D("hist","hist",1000, 0., 200.);
  TH2D* hist2 = new TH2D("hist2", "Momentum of Photon 1 versus Momentum of Photon 2", 1000, 0., 250., 1000,  0., 250.);
  double PH = 100.;
  double MH = 125.;
  gStyle->SetOptStat(0);
  TLorentzVector H, P1, P2;
   for(int j = 0; j<1000; j++){
     double theta = gRandom->Rndm()*2.*acos(-1.);
     double phi = gRandom->Rndm()*acos(-1.);
     double Px = PH*cos(theta)*sin(phi);
     double Py = PH*sin(theta)*sin(phi);
     double Pz = PH*cos(phi);
     H.SetPxPyPzE(Px,Py,Pz,sqrt(MH*MH+Px*Px+Py*Py+Pz*Pz));
     HiggsDecay(H, P1, P2);
     hist2->Fill(P1.P(), P2.P());
}
  hist2->GetXaxis()->SetTitle("Momentum of Photon 1");
  hist2->GetYaxis()->SetTitle("Momentum of Photon 2");
  hist2->GetYaxis()->SetTitleOffset(1.0);
  // hist2D->GetZaxis()->SetTitle("Number of events");
  hist2->Draw();
  /* hist2->SetMarkerSize(2);
  hist2->SetMarkerStyle(22); 
  hist2->SetMarkerColor(kBlue+3);*/
  cout << "OK" << endl;
}
void HiggsDecay(TLorentzVector H, TLorentzVector& P1, TLorentzVector& P2){
  double MH = H.M();
  TVector3 BH = (1/H.E())*H.Vect();
  double Epho = .5*MH;
  TVector3 Ppho;
  double theta = gRandom->Rndm()*2.0*acos(-1);
  double phi = acos(2.0*gRandom->Rndm()-1);
  
  Ppho.SetXYZ(cos(theta)*sin(phi)*Epho, sin(theta)*sin(phi)*Epho, cos(phi)*Epho);
  P1.SetVectM(Ppho, 0.);
  P2.SetVectM(-Ppho, 0.);

  P1.Boost(BH);
  P2.Boost(BH);
}
