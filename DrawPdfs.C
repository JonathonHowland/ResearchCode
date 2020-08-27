
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <vector>
#include <cmath>

using namespace RooFit;

void DrawPdfs(){

  //////////////////////////////////////////////////////////////////////////////////////////
  //Generate TTree from file
  TFile *f = new TFile("RooFitData.root");
  TTree *t1 = (TTree*)f->Get("t1");
  double cos1, cos2, cos3, MHa1, MHa2, MHa3, cosSum, massSum;
  t1->SetBranchAddress("cos1", &cos1);
  t1->SetBranchAddress("cos2", &cos2);
  t1->SetBranchAddress("cos3", &cos3);
  t1->SetBranchAddress("MHa", &MHa1);
  t1->SetBranchAddress("MHa2", &MHa2);
  t1->SetBranchAddress("MHa3", &MHa3);
  t1->SetBranchAddress("CosSum", &cosSum);
  t1->SetBranchAddress("MassSum", &massSum);

  cout << "TTree imported successfully" << endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Make 1D histograms
  
  Int_t palette[5];
  palette[0] = (kMagenta+2);
  palette[1] = (kRed);
  palette[2] = (kYellow);
  palette[3] = (kAzure+6);
  palette[4] = kBlue;
  gStyle->SetPalette(5,palette);
  double w = 1500;
  double h = 900;
  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("New Method","c1",w,h);
  c1->Divide(3,2,0.005,0.01,0);

  //Mass histograms
  TH1D* origm = new TH1D ("hist","Reconstructed Mass of H_{a}",100,0.,500.);
  TH1D* origm2 = new TH1D ("hist","Reconstructed Mass of H_{a}",100,0.,500.);
  TH1D* origm3 = new TH1D ("hist","Reconstructed Mass of H_{a}",100,0.,500.);
  
  //Decay angle histograms
  TH1D* origcos = new TH1D ("hist","Decay Angle of H_{a} to #it{l}_{a} and #it{l}_{b}",100,-1.,1.);
  TH1D* origcos2 = new TH1D ("hist","Decay Angle of H_{a} to #it{l}_{a} and #it{l}_{c}",100,-1.,1.);
  TH1D* origcos3 = new TH1D ("hist","Decay Angle of H_{a} to #it{l}_{b} and #it{l}_{c}",100,-1.,1.);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Fill 1D histograms
  int Nentries = t1->GetEntries();

  for(int i=0; i<Nentries; i++){
    t1->GetEntry(i);
    origm->Fill(MHa1);
    origm2->Fill(MHa2);
    origm3->Fill(MHa3);
    origcos->Fill(cos1);
    origcos2->Fill(cos2);
    origcos3->Fill(cos3);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Draw 1D histograms
  c1->cd(1);

  origm->GetXaxis()->SetTitle("Reconstructed Mass of Ha");
  origm->GetYaxis()->SetTitle("Counts");
  origm->GetYaxis()->SetRangeUser(0,origm->GetMaximum()*1.1);
  origm->SetFillColor(palette[4]);
  origm->SetFillStyle(3001);
  origm->SetLineColor(kBlack);
  origm->Draw();

  c1->cd(2);

  origm2->GetXaxis()->SetTitle("Reconstructed Mass of Ha");
  origm2->GetYaxis()->SetTitle("Counts");
  origm2->SetFillColor(palette[1]);
  origm2->SetFillStyle(3001);
  origm2->SetLineColor(kBlack);
  origm2->Draw();

  c1->cd(3);

  origm3->GetXaxis()->SetTitle("Reconstructed Mass of Ha");
  origm3->GetYaxis()->SetTitle("Counts");
  origm3->SetFillColor(palette[0]);
  origm3->SetFillStyle(3001);
  origm3->SetLineColor(kBlack);
  origm3->Draw();
  
  c1->cd(4);

  origcos->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos->GetYaxis()->SetTitle("Counts");
  origcos->GetYaxis()->SetRangeUser(0,origcos->GetMaximum()*1.1);
  origcos->SetFillColor(palette[4]);
  origcos->SetFillStyle(3001);
  origcos->SetLineColor(kBlack);
  origcos->Draw();

  c1->cd(5);

  origcos2->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos2->GetYaxis()->SetTitle("Counts");
  origcos2->SetFillColor(palette[1]);
  origcos2->SetFillStyle(3001);
  origcos2->SetLineColor(kBlack);
  origcos2->Draw();

  c1->cd(6);

  origcos3->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos3->GetYaxis()->SetTitle("Counts");
  origcos3->SetFillColor(palette[0]);
  origcos3->SetFillStyle(3001);
  origcos3->SetLineColor(kBlack);
  origcos3->Draw();

  
  c1->SaveAs("1DHists.jpg");
}
