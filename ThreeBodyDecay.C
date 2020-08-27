#include <iostream>
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TCanvas.h"
//Simulate intrinsic width and experimental width. Write a function that takes a photon and smears it. Make sure E2-P2 = 0. Make sure each reconstructed photon corresponds to one event. Plot quantities like Preco/Ptrue, MGreco/MGtrue, sin(0true-0reco), etc.
void TwoBodyDecay(TLorentzVector, double, double, 
		  TLorentzVector&, TLorentzVector&);
void GravitonDecay(double MG, double MP, double PG, TLorentzVector& vG,
		   TLorentzVector&, TLorentzVector&, TLorentzVector&,
		   TLorentzVector&);
void ThreeBodyDecay(){
  gStyle->SetPalette(51);
  double rand;
  double rand2;
  TRandom* Random = new TRandom();
  TH1D* hist = new TH1D("Statistics","Cos(P)",100, -1., 1.);
  TH2D* hist2 = new TH2D("hist2", "hist2", 1000, 0., 500., 1000,  0., 500.);

  TH2D* hist2D = new TH2D("Statistics","MP versus cos(P)",30,-1., 1.,30,0.,250.);
  TH2D* Wrong1 = new TH2D("Statistics","MP versus cos(P) with incorrect decay",30,-1., 1.,30,0.,250.);
  TH2D* Wrong2 = new TH2D("Statistics","MP versus cos(P) with incorrect decay",30,-1., 1.,30,0.,250.);
  TH1D* wrong1 = new TH1D("Statistics","Cos(P) with incorrect decay",100, -1., 1.);
  TH1D* wrong2 = new TH1D("Statistics","Cos(P) with incorrect decay",100, -1., 1.);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.23);
   Double_t w = 1800;
   Double_t h = 600;
   TCanvas* c1 = new TCanvas("c", "c", w, h);
   c1->Divide(3,2, 0.01, 0.01, 0);
  double MG = 250, MP = 125, PG = 300;
  // double PH = 100.;
  // double MH = 125.;
  TLorentzVector vG, VP, VC1, VC2, VC3;
  for (int j = 0; j<100000; j++){
    double JMG = -1;
    double JMP = -.5;
   while(JMP>JMG || JMP<=0){
    JMG =  Random->Gaus(MG,MG*0.1);    
    JMP = Random->Gaus(MP,MP*0.1);
    }
    GravitonDecay(JMG, JMP, PG, vG, VP, VC1, VC2, VC3);
    hist2->Fill((VP).P(), (VC3).P());

  TLorentzVector vGr = (VC1+VC2+VC3);
  double mG = vGr.M();
  TVector3 vBeta_G = vGr.BoostVector();
  // moving to G rest frame
  VC1.Boost(-vBeta_G);
  VC2.Boost(-vBeta_G);
  VC3.Boost(-vBeta_G);
  TLorentzVector vParent = VC1 +VC2;
  double cosG = vBeta_G.Unit().Dot(vParent.Vect().Unit());

  double mParent = vParent.M();
  double mChild3 = VC3.M();
  TVector3 vBeta_P = vParent.BoostVector();
  TVector3 vBeta_C3 = VC3.BoostVector();
  // moving to H1 rest frame
  VC1.Boost(-vBeta_P);
  VC2.Boost(-vBeta_P);
  double cosP = vBeta_P.Unit().Dot(VC1.Vect().Unit());
  double cosC3 = vBeta_C3.Unit().Dot(VC3.Vect().Unit());
 
  hist2D->Fill(cosP, mParent);
  hist ->Fill(cosP);
 
}
  c1->cd(1);
  hist2D->Draw("colz");
  hist2D->GetXaxis()->SetTitle("Cos(#theta_{P})");
  hist2D->GetYaxis()->SetTitle("MP");
  hist2D->GetYaxis()->SetTitleOffset(1.5);
  hist2D->GetZaxis()->SetTitle("Number of events");
  c1->cd(4);
  hist->Draw();
  hist->GetXaxis()->SetTitle("Cos(#theta_{P})");
  hist->GetYaxis()->SetTitle("Number of events");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetRangeUser(0,hist->GetMaximum()*1.1);
  for (int j = 0; j<100000; j++){
    GravitonDecay(MG, MP, PG, vG, VP, VC1, VC2, VC3);
    hist2->Fill((VP).P(), (VC3).P());

  TLorentzVector temp = VC2;
  VC2 = VC3;
  VC3 = temp;

  TLorentzVector vGr = (VC1+VC2+VC3);
  double mG = vGr.M();
  TVector3 vBeta_G = vGr.BoostVector();
  // moving to G rest frame
  VC1.Boost(-vBeta_G);
  VC2.Boost(-vBeta_G);
  VC3.Boost(-vBeta_G);
  TLorentzVector vParent = VC1 +VC2;
  double cosG = vBeta_G.Unit().Dot(vParent.Vect().Unit());

  double mParent = vParent.M();
  double mChild3 = VC3.M();
  TVector3 vBeta_P = vParent.BoostVector();
  TVector3 vBeta_C3 = VC3.BoostVector();
  // moving to H1 rest frame
  VC1.Boost(-vBeta_P);
  VC2.Boost(-vBeta_P);

  double cosP = vBeta_P.Unit().Dot(VC1.Vect().Unit());
  double cosC3 = vBeta_C3.Unit().Dot(VC3.Vect().Unit());
 
  Wrong1->Fill(cosP, mParent);
  wrong1->Fill(cosP);
 
}
  c1->cd(2);
  Wrong1->Draw("colz");
  Wrong1->GetXaxis()->SetTitle("Cos(#theta_{P})");
  Wrong1->GetYaxis()->SetTitle("MP");
  Wrong1->GetYaxis()->SetTitleOffset(1.5);
  c1->cd(5);
  wrong1->Draw();
  wrong1->GetXaxis()->SetTitle("Cos(#theta_{P})");
  wrong1->GetYaxis()->SetTitle("Number of events");
  wrong1->GetYaxis()->SetTitleOffset(1.5);
  for (int j = 0; j<100000; j++){
    GravitonDecay(MG, MP, PG, vG, VP, VC1, VC2, VC3);
    hist2->Fill((VP).P(), (VC3).P());

  TLorentzVector temp = VC1;
  VC1 = VC3;
  VC3 = temp;

  TLorentzVector vGr = (VC1+VC2+VC3);
  double mG = vGr.M();
  TVector3 vBeta_G = vGr.BoostVector();
  // moving to G rest frame
  VC1.Boost(-vBeta_G);
  VC2.Boost(-vBeta_G);
  VC3.Boost(-vBeta_G);
  TLorentzVector vParent = VC1 +VC2;
  double cosG = vBeta_G.Unit().Dot(vParent.Vect().Unit());

  double mParent = vParent.M();
  double mChild3 = VC3.M();
  TVector3 vBeta_P = vParent.BoostVector();
  TVector3 vBeta_C3 = VC3.BoostVector();
  // moving to H1 rest frame
  VC1.Boost(-vBeta_P);
  VC2.Boost(-vBeta_P);

  double cosP = vBeta_P.Unit().Dot(VC1.Vect().Unit());
  double cosC3 = vBeta_C3.Unit().Dot(VC3.Vect().Unit());
 
  Wrong2->Fill(cosP, mParent);
  wrong2->Fill(cosP);
 
}
  c1->cd(3);
  Wrong2->Draw("colz");
  Wrong2->GetXaxis()->SetTitle("Cos(#theta_{P})");
  Wrong2->GetYaxis()->SetTitle("MP");
  Wrong2->GetYaxis()->SetTitleOffset(1.5);
  c1->cd(6);
  wrong2->Draw();
  wrong2->GetXaxis()->SetTitle("Cos(#theta_{P})");
  wrong2->GetYaxis()->SetTitle("Number of events");
  wrong2->GetYaxis()->SetTitleOffset(1.5);

  c1->SaveAs("Decays.jpg");
 /* hist2->SetMarkerSize(3);
  hist2->SetMarkerStyle(22); 
  hist2->SetMarkerColor(kBlue+3);*/
  cout << "OK" << endl;
}
void GravitonDecay(double MG, double MP, double PG, TLorentzVector& vG, TLorentzVector& VP,
		   TLorentzVector& VC1, TLorentzVector& VC2,
		   TLorentzVector& VC3){
  double theta = gRandom->Rndm()*2.*acos(-1.);
  double phi = acos(2.0*gRandom->Rndm()-1);
  double Px = PG*cos(theta)*sin(phi);
  double Py = PG*sin(theta)*sin(phi);
  double Pz = PG*cos(phi);
  vG.SetPxPyPzE(Px,Py,Pz,sqrt(MG*MG+Px*Px+Py*Py+Pz*Pz));
  TwoBodyDecay(vG, MP, 0., VP, VC3);
  TwoBodyDecay(VP, 0., 0., VC1, VC2);
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
