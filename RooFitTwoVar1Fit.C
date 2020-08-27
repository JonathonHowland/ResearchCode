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

using namespace RooFit;

void RooFitTwoVar1Fit(){
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

  /////////////////////////////////////////////////////////////////////////////////////////////
  RooMsgService::instance().setSilentMode(1);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  //Declare RooFit variables and datasets
  //Observables
  RooRealVar cosR("cosR", "cosR", -1.,1.);
  
  RooRealVar MHaR("MHaR", "MHa", 0., 500.);

  //Datasets
  //Treat cos1D and cos2D as backgroundd
  RooDataSet cos1D("name","title", RooArgSet(cosR));
  RooDataSet cos2D("name","title", RooArgSet(cosR));
  RooDataSet cos3D("name","title", RooArgSet(cosR));
  RooDataSet MHa1D("name","title", RooArgSet(MHaR));
  RooDataSet MHa2D("name","title", RooArgSet(MHaR));
  RooDataSet MHa3D("name","title", RooArgSet(MHaR));

  //RooDataSet for both
  RooDataSet cosMass1D("cosMass1D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet cosMass2D("cosMass2D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet cosMass3D("cosMass3D", "cosMassD", RooArgSet(cosR, MHaR));

  //Decay angle
  
  //Parameters
  RooRealVar ncos1("ncos1","signal events", 5000., 0., 300000);
  RooRealVar ncos2("ncos2", "signal events", 5000., 0., 300000);
  RooRealVar ncos3("ncos3","signal events", 5000., 0., 300000);
  
  /////////////////////////////////////////////////////////////////////////////
  //Read 1D pdfs from file

  TFile *fin1 = new TFile("1DPdf.root", "READ");

  RooWorkspace* w = (RooWorkspace*) fin1->Get("w");
  
  RooAbsPdf* cos1Pdf = w->pdf("cos1Pdf");
  RooAbsPdf* cos2Pdf = w->pdf("cos2Pdf");
  RooAbsPdf* cos3Pdf = w->pdf("cos3Pdf");
  RooAbsPdf* MHa1Pdf = w->pdf("MHa1Pdf");
  RooAbsPdf* MHa2Pdf = w->pdf("MHa2Pdf");
  RooAbsPdf* MHa3Pdf = w->pdf("MHa3Pdf");
  
  cout <<"1D Pdfs imported successfully" << endl;
  
  //Add the pdfs together

  RooAddPdf cosModel("modelCos", "modelCos", RooArgList(*cos1Pdf, *cos2Pdf, *cos3Pdf), RooArgList(ncos1, ncos2, ncos3));
  
  RooAddPdf massModel("massModel", "massModel", RooArgList(*MHa1Pdf, *MHa2Pdf, *MHa3Pdf), RooArgList(ncos1, ncos2, ncos3));

  cout << "1D models built successfully" << endl;

  //////////////////////////////////////////////////////////////////////////////////////
  //Build the 2D pdf
  //Import histograms

  TFile *fin2 = new TFile("2dHist.root","READ");
  TH2D *cosMassHist1 = (TH2D*) fin2->Get("cosMassHist1");
  TH2D *cosMassHist2 = (TH2D*) fin2->Get("cosMassHist2");
  TH2D *cosMassHist3 = (TH2D*) fin2->Get("cosMassHist3");

  std::cout << "Histograms imported successfully" << std::endl;

  //Change to RooDataHist
  RooDataHist *cosMassDataH1 = new RooDataHist("DataHist1", "DataHist1", RooArgList(cosR, MHaR), cosMassHist1);
  RooDataHist *cosMassDataH2 = new RooDataHist("DataHist2", "DataHist2", RooArgList(cosR, MHaR), cosMassHist2);
  RooDataHist *cosMassDataH3 = new RooDataHist("DataHist3", "DataHist3", RooArgList(cosR, MHaR), cosMassHist3);

  std::cout << "RooDataHists created successfully" << std::endl;

  //Create RooHistPdfs
  RooHistPdf *cosMassPdf1 = new RooHistPdf("cosMassPdf1", "cosMassPdf1", RooArgList(cosR, MHaR), *cosMassDataH1);
  RooHistPdf *cosMassPdf2 = new RooHistPdf("cosMassPdf2", "cosMassPdf2", RooArgList(cosR, MHaR), *cosMassDataH2);
  RooHistPdf *cosMassPdf3 = new RooHistPdf("cosMassPdf3", "cosMassPdf3", RooArgList(cosR, MHaR), *cosMassDataH3);

  cout << "RooHistPdfs created successfully" << endl;

  //Add the pdfs together to create the model

  RooAddPdf *cosMassModel = new RooAddPdf("cosMassModel", "cosMassModel", RooArgList(*cosMassPdf1, *cosMassPdf2, *cosMassPdf3), RooArgList(ncos1, ncos2, ncos3));

  cout << "Pdfs added together. Model built successfully. " << endl;

  ///////////////////////////////////////////////////////////////////////////////////
  //Do the fit

  

  //Build data sets
  
  int N = 1000;
  double prob1 = 0.01;
  double prob2 = (1.0-prob1)/2.0;
  double prob3 = 1.0 - prob1 - prob2;

  int N1 = prob1*N;
  int N2 = prob2*N;
  int N3 = prob3*N;
  
  RooDataSet TotalCos("name","title", RooArgSet(cosR));
  RooDataSet TotalMass("totalm", "totalm", RooArgSet(MHaR));
  RooDataSet TotalCosMass("Total both", "Total both", RooArgSet(cosR, MHaR));

  //Fill the datasets from the histograms
  for(int k=0; k<N1; k++){
    double cos = 0.0;
    double mass = 0.0;
    cosMassHist1->GetRandom2(cos, mass);
    cosR.setVal(cos);
    MHaR.setVal(mass);
    TotalCosMass.add(RooArgSet(cosR, MHaR));    
  }
  for(int k=0; k<(N2); k++){
    double cos = 0.0;
    double mass = 0.0;
    cosMassHist2->GetRandom2(cos, mass);
    cosR.setVal(cos);
    MHaR.setVal(mass);
    TotalCosMass.add(RooArgSet(cosR, MHaR));    
  }
  for(int k=0; k<(N3); k++){
    double cos = 0.0;
    double mass = 0.0;
    cosMassHist3->GetRandom2(cos, mass);
    cosR.setVal(cos);
    MHaR.setVal(mass);
    TotalCosMass.add(RooArgSet(cosR, MHaR));    
  }

  //Fit with cos only
  cosModel.fitTo(TotalCosMass);
  
  //Plotting
  TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
  c1->cd();
  RooPlot* xframe = MHaR.frame();
  TotalCosMass.plotOn(xframe);
  cosMassModel->plotOn(xframe);
  xframe->Draw();

  //Fit with mass only
  massModel.fitTo(TotalCosMass);
  
  //Plotting
  TCanvas* c1 = new TCanvas("c2", "c2", 1200, 600);
  c2->cd();
  RooPlot* yframe = cosR.frame();
  TotalCosMass.plotOn(yframe);
  cosMassModel->plotOn(yframe);
  yframe->Draw();
    
  //Fit with both cos and mass
  cosMassModel->fitTo(TotalCosMass);
  
  //Plotting
  TCanvas* c3 = new TCanvas("c3", "c3", 1200, 600);
  c3->cd();
  RooPlot* xframe = MHaR.frame();
  TotalCosMass.plotOn(zframe);
  cosMassModel->plotOn(zframe);
  zframe->Draw();

}
