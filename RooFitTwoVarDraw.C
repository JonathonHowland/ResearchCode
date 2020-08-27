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
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <vector>
#include <cmath>

using namespace std;

void RooFitTwoVarDraw(){  
  gStyle->SetPalette(51);
  gStyle->SetOptStat(0);
  
  //Import histograms

  TFile *fin2 = new TFile("2dHist.root","READ");
  TH2D *cosMassHist1 = (TH2D*) fin2->Get("cosMassHist1");
  TH2D *cosMassHist2 = (TH2D*) fin2->Get("cosMassHist2");
  TH2D *cosMassHist3 = (TH2D*) fin2->Get("cosMassHist3");
  TH2D *cosMassHistBack = new TH2D("cosMassHistBack", "", 50, -1., 1., 50, 0., 500.);

  std::cout << "Histograms imported successfully" << std::endl;

  cosMassHist1->SetTitle("2D histogram of signal events' decay angle and mass");
  cosMassHist2->SetTitle("2D histogram of background events' decay angle and mass");
  cosMassHist3->SetTitle("2D histogram of background events' decay angle and mass");

  //Create a histogram of all background events
  cosMassHistBack->SetTitle("2D histogram of background events' decay angle and mass");
  cosMassHistBack->Add(cosMassHist2);
  cosMassHistBack->Add(cosMassHist3);
  
  TCanvas* c1 = new TCanvas("c1", "2D histogram of signal events decay angle and mass", 1200, 600);
  TCanvas* c2 = new TCanvas("c2", "2D histogram of background events decay angle and mass", 1200, 600);
  TCanvas* c3 = new TCanvas("c3", "2D histogram of background decay events angle and mass", 1200, 600);
  TCanvas* c4 = new TCanvas("c4", "2D histogram of background decay events angle and mass", 1200, 600);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //Draw histograms
  
  c1->cd();
  c1->SetLogz();
  cosMassHist1->GetXaxis()->SetTitle("Decay angle");
  cosMassHist1->GetYaxis()->SetTitle("Mass");
  cosMassHist1->Draw("colz");

  c2->cd();
  c2->SetLogz();
  cosMassHist2->GetXaxis()->SetTitle("Decay angle");
  cosMassHist2->GetYaxis()->SetTitle("Mass");
  cosMassHist2->Draw("colz");

  c3->cd();
  c3->SetLogz();
  cosMassHist3->GetXaxis()->SetTitle("Decay angle");
  cosMassHist3->GetYaxis()->SetTitle("Mass");
  cosMassHist3->Draw("colz");

  c4->cd();
  c4->SetLogz();
  cosMassHistBack->GetXaxis()->SetTitle("Decay angle");
  cosMassHistBack->GetYaxis()->SetTitle("Mass");
  cosMassHistBack->Draw("colz");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //Save histograms

  c1->SaveAs("SignalHist.jpg");
  c2->SaveAs("BackHist1.jpg");
  c3->SaveAs("BackHist2.jpg");
  c4->SaveAs("BackHist.jpg");
}
