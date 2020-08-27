
#include <iostream>
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
//Analyze three body decay in three different ways. Find mc and cos(phi)c. Find mp and cos(phi)p.
void TwoBodyDecay(TLorentzVector, double, double, 
		  TLorentzVector&, TLorentzVector&);
void GravitonDecay(double MG, double MH1, double MH2, double PG, TLorentzVector& vG, TLorentzVector& vH1,
		   TLorentzVector& vH2, TLorentzVector& vGam11, TLorentzVector& vGam12,
		   TLorentzVector& vGam21, TLorentzVector& vGam22);

void DecayTwoBody(){
  double rand;
  double rand2;
  TRandom* Random = new TRandom();
  TH1D* hist = new TH1D("hist","hist",1000, 0., 200.);
  TH2D* hist2 = new TH2D("hist2", "hist2", 1000, 0., 500., 1000,  0., 500.);

  TH2D* hist2D = new TH2D("hist2D","hist2D",100,-1.,1,100,-1.,1.);
  double MG, MH1, MH2, PG;
  // double PH = 100.;
  // double MH = 125.;
  TLorentzVector vG, vH1, vH2, vGam11, vGam12, vGam21, vGam22;
  for (int j = 0; j<1000; j++){
  GravitonDecay(MG, MH1, MH2, PG, vG, vH1, vH2, vGam11, vGam12, vGam21, vGam22);
  hist2->Fill((vGam11+vGam21).M(), (vGam12+vGam22).M());
  //cout << vG.Px() <<  " " <<  (vGam11 + vGam12 +vGam21 +vGam22).Px() << endl;

  TLorentzVector temp = vGam11;
  vGam11 = vGam22;
  vGam22 = temp;

  TLorentzVector vGr = (vGam11+vGam12+vGam21+vGam22);
  double mG = vGr.M();
  TVector3 vBeta_G = vGr.BoostVector();
  // moving to G rest frame
  vGam11.Boost(-vBeta_G);
  vGam12.Boost(-vBeta_G);
  vGam21.Boost(-vBeta_G);
  vGam22.Boost(-vBeta_G);
  TLorentzVector vHiggs1 = vGam11+vGam12;
  TLorentzVector vHiggs2 = vGam21+vGam22;
  double cosG = vBeta_G.Unit().Dot(vHiggs1.Vect().Unit());

  double mHiggs1 = vHiggs1.M();
  double mHiggs2 = vHiggs2.M();
  TVector3 vBeta_H1 = vHiggs1.BoostVector();
  TVector3 vBeta_H2 = vHiggs2.BoostVector();
  // moving to H1 rest frame
  vGam11.Boost(-vBeta_H1);
  vGam12.Boost(-vBeta_H1);
  // moving to H2 rest frame
  vGam21.Boost(-vBeta_H2);
  vGam22.Boost(-vBeta_H2);

  double cosH1 = vBeta_H1.Unit().Dot(vGam11.Vect().Unit());
  double cosH2 = vBeta_H2.Unit().Dot(vGam21.Vect().Unit());

  hist2D->Fill(cosH1,cosH2);

}
  /* for(int j = 0; j<10000; j++){
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
     cout << "G mass " << G.M() << " " << (H1 +H2).M() << endl;
     cout << "G pt " << G.Pt() << " " << (H1 + H2).Pt() << endl;
     cout << "H1 mass " << H1.M() << " " << "H2 mass " << H2.M() << endl;
     }*/
  hist2D->Draw();
  /* hist2->SetMarkerSize(3);
  hist2->SetMarkerStyle(22); 
  hist2->SetMarkerColor(kBlue+3);*/
  cout << "OK" << endl;
}
void GravitonDecay(double MG, double MH1, double MH2, double PG, TLorentzVector& vG, TLorentzVector& vH1,
		   TLorentzVector& vH2, TLorentzVector& vGam11, TLorentzVector& vGam12,
		   TLorentzVector& vGam21, TLorentzVector& vGam22){
  MG = 500;
  MH1 = 125;
  MH2 = 125;
  PG = 200;
  double theta = gRandom->Rndm()*2.*acos(-1.);
  double phi = acos(2.0*gRandom->Rndm()-1);
  double Px = PG*cos(theta)*sin(phi);
  double Py = PG*sin(theta)*sin(phi);
  double Pz = PG*cos(phi);
  vG.SetPxPyPzE(Px,Py,Pz,sqrt(MG*MG+Px*Px+Py*Py+Pz*Pz));
  TwoBodyDecay(vG, MH1, MH2, vH1, vH2);
  TwoBodyDecay(vH1, 0., 0., vGam11, vGam12);
  TwoBodyDecay(vH2, 0., 0., vGam21, vGam22);
}
void TwoBodyDecay(TLorentzVector vP, double mC1, double mC2, TLorentzVector& vC1, TLorentzVector& vC2){
  double mvP = vP.M();
  TVector3 BvP = (1/vP.E())*vP.Vect();
  double EvC1 = (mvP*mvP + mC1*mC1 - mC2*mC2)/2./mvP;
  double EvC2 = (mvP*mvP - mC1*mC1 + mC2*mC2)/2./mvP;
  double PvC = sqrt((mvP*mvP - (mC1+mC2)*(mC1+mC2))*(mvP*mvP - (mC1-mC2)*(mC1-mC2)))/2./mvP;
  double theta = gRandom->Rndm()*2.0*acos(-1);
  double phi = acos(2.0*gRandom->Rndm()-1);
 
  double Px = PvC*cos(theta)*sin(phi); 
  double Py = PvC*sin(theta)*sin(phi);
  double Pz = PvC*cos(phi);
  vC1.SetPxPyPzE(Px, Py, Pz, sqrt(mC1*mC1 + PvC*PvC));
  vC2.SetPxPyPzE(-Px, -Py, -Pz, sqrt(mC2*mC2 + PvC*PvC));
  vC1.Boost(BvP);
  vC2.Boost(BvP);
}
