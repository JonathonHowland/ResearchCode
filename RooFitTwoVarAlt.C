#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TRandom.h"
#include "TGraph.h"
#include <vector>
#include <cmath>

using namespace RooFit;

void RooFitTwoVarAlt(){
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


  //Observables
  RooRealVar cosR("cosR", "cosR", -1.,1.);
  
  RooRealVar MHaR("MHa", "MHa", 0., 500.);

  //Datasets
  //Treat cos1D and cos2D as backgroundd
  RooDataSet cos1D("name","title", RooArgSet(cosR));
  RooDataSet cos2D("name","title", RooArgSet(cosR));
  RooDataSet cos3D("name","title", RooArgSet(cosR));
  //RooDataSet TotalCosBkg("bkg_cos","bkg_cos",RooArgSet(cosR));
  RooDataSet TotalCos("name","title", RooArgSet(cosR));
  RooDataSet MHa1D("name","title", RooArgSet(MHaR));
  RooDataSet MHa2D("name","title", RooArgSet(MHaR));
  RooDataSet MHa3D("name","title", RooArgSet(MHaR));
  RooDataSet TotalMass("totalm", "totalm", RooArgSet(MHaR));

  //RooDataSet for both
  RooDataSet cosMass1D("cosMass1D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet cosMass2D("cosMass2D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet cosMass3D("cosMass3D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet TotalCosMass("Total both", "Total both", RooArgSet(cosR, MHaR));

  double prob1 = .1;
  double prob2 = .2;
  double prob3 = 1. - prob1 - prob2;
  
  //Build the data sets for the pdf
  int Nentries = t1->GetEntries();

  cout << Nentries << endl;
  for(int i = 0; i < Nentries/100; i++){
    t1->GetEntry(i);
    //Decay angle
    cosR.setVal(cos1);
    cos1D.add(cosR);
    cosMass1D.add(cosR);
    cosR.setVal(cos2);
    cos2D.add(cosR);
    cosMass2D.add(cosR);
    cosR.setVal(cos3);
    cos3D.add(cosR);
    cosMass3D.add(cosR);

    //Mass
    MHaR.setVal(MHa1);
    MHa1D.add(MHaR);
    cosMass1D.add(MHaR);
    MHaR.setVal(MHa2);
    MHa2D.add(MHaR);
    cosMass2D.add(MHaR);
    MHaR.setVal(MHa3);
    MHa3D.add(MHaR);
    cosMass3D.add(MHaR);
  }
  std::cout << "Data sets built successfully" << std::endl;

  //Build a 2D RooKeysPdf

  RooNDKeysPdf cosMass1Pdf("cosMass1", "cosMass1", RooArgList(cosR, MHaR), cosMass1D);
  RooNDKeysPdf cosMass2Pdf("cosMass2", "cosMass2", RooArgList(cosR, MHaR), cosMass2D);
  RooNDKeysPdf cosMass3Pdf("cosMass3", "cosMass3", RooArgList(cosR, MHaR), cosMass3D);

  std::cout << "Pdfs built successfully" << std::endl;
  
  TH2D* cosMassHist1 = new TH2D("cosMassHist1","cosMassHist1",200,-1,1.,200,0.,500.);
  TH2D* cosMassHist2 = new TH2D("cosMassHist2","cosMassHist2",200,-1,1.,200,0.,500.);
  TH2D* cosMassHist3 = new TH2D("cosMassHist3","cosMassHist3",200,-1,1.,200,0.,500.);

  int NbinX = cosMassHist1->GetNbinsX();
  int NbinY = cosMassHist1->GetNbinsY();

  for(int ix = 0; ix < NbinX; ix++){
    cosR.setVal( cosMassHist1->GetXaxis()->GetBinCenter(ix+1) );
    for(int iy = 0; iy < NbinY; iy++){
      MHaR.setVal( cosMassHist1->GetYaxis()->GetBinCenter(iy+1) );

      cosMassHist1->SetBinContent(ix+1,iy+1, cosMass1Pdf.getVal());
      cosMassHist2->SetBinContent(ix+1,iy+1, cosMass2Pdf.getVal());
      cosMassHist3->SetBinContent(ix+1,iy+1, cosMass3Pdf.getVal());
    }
  }

  TFile* fout = new TFile("2dHist2.root","RECREATE");
  fout->cd();
  cosMassHist1->Write("",TObject::kOverwrite);
  cosMassHist2->Write("",TObject::kOverwrite);
  cosMassHist3->Write("",TObject::kOverwrite);
  fout->Close();
  
  std::cout << "Histogram successfully saved to file." << std::endl;

}
