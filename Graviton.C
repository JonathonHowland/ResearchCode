
#include <iostream>
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
//Make histograms of photons and change parameters
void HiggsDecay1(TLorentzVector H1, TLorentzVector& P11, TLorentzVector& P12);
void HiggsDecay2(TLorentzVector H2, TLorentzVector& P21, TLorentzVector& P22);
void GravDecay(TLorentzVector G, TLorentzVector& H1, TLorentzVector& H2);

void Graviton(){
  double w = 1200;
  double h = 900;
  TCanvas* cGrav = new TCanvas("c","Old Graviton",w,h);
  double rand;
  double rand2;
  TRandom* Random = new TRandom();
  TH1D* hist = new TH1D("hist","hist",1000, 0., 200.);
  TH2D* hist2 = new TH2D("hist2", "hist2", 1000, 0., 250., 1000,  0., 250.);
  double MG, MH1, MH2, PG;
  MG = 250;
  PG = 200;
  // double PH = 100.;
  // double MH = 125.;
  TLorentzVector G, H1, H2, P11, P12, P21, P22;
   for(int j = 0; j<10000; j++){
     double theta = gRandom->Rndm()*2.*acos(-1.);
     double phi = acos(2.0*gRandom->Rndm()-1);
     double Px = PG*cos(theta)*sin(phi);
     double Py = PG*sin(theta)*sin(phi);
     double Pz = PG*cos(phi);
     G.SetPxPyPzE(Px,Py,Pz,sqrt(MG*MG+Px*Px+Py*Py+Pz*Pz));
     GravDecay(G, H1, H2);
     HiggsDecay1(H1, P11, P12);
     HiggsDecay2(H2, P21, P22);
     hist2->Fill(H1.P(), H2.P());
     // cout << "G mass " << G.M() << " " << (H1 +H2).M() << endl;
     //cout << "G pt " << G.Pt() << " " << (H1 + H2).Pt() << endl;
     //cout << "H1 mass " << H1.M() << " " << "H2 mass " << H2.M() << endl;
     // cout << P11.Px() << " " << P12.Px() << endl;
 }

   cGrav->cd();
  hist2->Draw();
  hist2->SetMarkerSize(3);
  hist2->SetMarkerStyle(22); 
  hist2->SetMarkerColor(kBlue+3);
  cout << "OK" << endl;
}

void GravDecay(TLorentzVector G, TLorentzVector& H1, TLorentzVector& H2){
  double MH1 = 100;
  double MH2 = 100;
  double MG = G.M();
  TVector3 BG = (1/G.E())*G.Vect();
  double EH1 = (MG*MG + MH1*MH1 - MH2*MH2)/2./MG;
  double EH2 = (MG*MG - MH1*MH1 + MH2*MH2)/2./MG;
  double PH = sqrt( (MG*MG - (MH1+MH2)*(MH1+MH2))*
		    (MG*MG - (MH1-MH2)*(MH1-MH2)) )/2./MG;
  TVector3 vPH;
  double theta = gRandom->Rndm()*2.0*acos(-1);
  double phi = acos(2.0*gRandom->Rndm()-1);
 
  double Px = PH*cos(theta)*sin(phi); 
  double Py = PH*sin(theta)*sin(phi);
  double Pz = PH*cos(phi);
  //cout << MG*MG  << " " << (MH1+MH2)*(MH1+MH2) << endl;
  H1.SetPxPyPzE(Px, Py, Pz, sqrt(MH1*MH1 + PH*PH));
  H2.SetPxPyPzE(-Px, -Py, -Pz, sqrt(MH2*MH2 + PH*PH));

  H1.Boost(BG);
  H2.Boost(BG);
   
  }
void HiggsDecay1(TLorentzVector H1, TLorentzVector& P11, TLorentzVector& P12){
  double MH = H1.M();
  TVector3 BH = (1/H1.E())*H1.Vect();
  double Epho = .5*MH;
  TVector3 Ppho;
  double theta = gRandom->Rndm()*2.0*acos(-1);
  double phi = acos(2.0*gRandom->Rndm()-1);
  
  Ppho.SetXYZ(cos(theta)*sin(phi)*Epho, sin(theta)*sin(phi)*Epho, cos(phi)*Epho);
  P11.SetVectM(Ppho, 0.);
  P12.SetVectM(-Ppho, 0.);

  P11.Boost(BH);
  P12.Boost(BH);
}
void HiggsDecay2(TLorentzVector H2, TLorentzVector& P21, TLorentzVector& P22){
  double MH = H2.M();
  TVector3 BH = (1/H2.E())*H2.Vect();
  double Epho = .5*MH;
  TVector3 Ppho;
  double theta = gRandom->Rndm()*2.0*acos(-1);
  double phi = acos(2.0*gRandom->Rndm()-1);
  
  Ppho.SetXYZ(cos(theta)*sin(phi)*Epho, sin(theta)*sin(phi)*Epho, cos(phi)*Epho);
  P21.SetVectM(Ppho, 0.);
  P22.SetVectM(-Ppho, 0.);

  P21.Boost(BH);
  P22.Boost(BH);
  }
