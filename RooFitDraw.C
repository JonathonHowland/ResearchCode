//Code to draw fits

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
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMsgService.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TError.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace RooFit;

void RooFitDraw(){
  //////////////////////////////////////////////////////////////////////////////////////////
  // gStyle->SetOptStat(0);

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
  //RooMsgService::instance().setSilentMode(1);
  //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  //Declare RooFit variables and datasets
  //Observables
  RooRealVar cosR("cosR", "cosR", -1.,1.);
  
  RooRealVar MHaR("MHaR", "MHa", 0., 500.);

  //Datasets
  //Treat cos1D and cos2D as backgroundd
  RooDataSet cos1D("cos1D","title", RooArgSet(cosR));
  RooDataSet cos2D("cos2D","title", RooArgSet(cosR));
  RooDataSet cos3D("cos3D","title", RooArgSet(cosR));
  RooDataSet TotalCos("TotalCos","title", RooArgSet(cosR));
  RooDataSet MHa1D("MHa1D","title", RooArgSet(MHaR));
  RooDataSet MHa2D("MHa2D","title", RooArgSet(MHaR));
  RooDataSet MHa3D("MHa3D","title", RooArgSet(MHaR));
  RooDataSet TotalMass("TotalMass", "totalm", RooArgSet(MHaR));

  //RooDataSet for both
  RooDataSet cosMass1D("cosMass1D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet cosMass2D("cosMass2D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet cosMass3D("cosMass3D", "cosMassD", RooArgSet(cosR, MHaR));
  RooDataSet TotalCosMass("Total both", "Total both", RooArgSet(cosR, MHaR));
  RooDataSet TotalCosMassBack("TotalCosMassBack", "Total background", RooArgSet(cosR, MHaR));

  //Decay angle
  
  //Parameters
  RooRealVar ncos1("ncos1","signal events", 5000., 0., 1000000);
  RooRealVar ncos2("ncos2", "background events", 5000., 0., 1000000);
  RooRealVar ncos3("ncos3","background events", 5000., 0., 1000000);

  //Parameters for pdfs without signal
  RooRealVar ncos22("ncos22", "background events", 5000., 0., 1000000);
  RooRealVar ncos23("ncos23","background events", 5000., 0., 1000000);

  
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

  RooAddPdf cosModel("cosModel", "modelCos", RooArgList(*cos1Pdf, *cos2Pdf, *cos3Pdf), RooArgList(ncos1, ncos2, ncos3));
  
  RooAddPdf massModel("massModel", "massModel", RooArgList(*MHa1Pdf, *MHa2Pdf, *MHa3Pdf), RooArgList(ncos1, ncos2, ncos3));

  cout << "1D models built successfully" << endl;

  //Add just the background Pdfs together

  RooAddPdf cosBackModel("cosBackModel", "cosBackModel", RooArgList(*cos2Pdf, *cos3Pdf), RooArgList(ncos22, ncos23));

  RooAddPdf massBackModel("massBackModel", "cosBackModel", RooArgList(*MHa2Pdf, *MHa2Pdf), RooArgList(ncos22, ncos23));

  cout << "1D models (without signal) built successfully." << endl;

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

  ////////////////////////////////////////////////////////////////////////////////////////////
  //Fill datasets
  int N1 = 10;
  int N2 = 495;
  int N3 = 495;
  
    //Fill the datasets from the histograms (all 3 types of events)
    for(int k=0; k<N1; k++){
      double cos = 0.0;
      double mass = 0.0;
      cosMassHist1->GetRandom2(cos, mass);
      cosR.setVal(cos);
      MHaR.setVal(mass);
      TotalCos.add(cosR);
      TotalMass.add(MHaR);
      TotalCosMass.add(RooArgSet(cosR, MHaR));    
    }
    for(int k=0; k<N2; k++){
      double cos = 0.0;
      double mass = 0.0;
      cosMassHist2->GetRandom2(cos, mass);
      cosR.setVal(cos);
      MHaR.setVal(mass);
      TotalCos.add(cosR);
      TotalMass.add(MHaR);
      TotalCosMass.add(RooArgSet(cosR, MHaR));    
    }
    for(int k=0; k<N3; k++){
      double cos = 0.0;
      double mass = 0.0;
      cosMassHist3->GetRandom2(cos, mass);
      cosR.setVal(cos);
      MHaR.setVal(mass);
      TotalCos.add(cosR);
      TotalMass.add(MHaR);
      TotalCosMass.add(RooArgSet(cosR, MHaR));    
    }
    massModel.fitTo(TotalCosMass);
    cosModel.fitTo(TotalCosMass);
    cosMassModel->fitTo(TotalCosMass);

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //Draw the fits

    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
    RooPlot* xframe = MHaR.frame();
    TotalCosMass.plotOn(xframe);
    massModel.plotOn(xframe);
    xframe->Draw();
    c1->SaveAs("MassFit.jpg");

    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 600);
    RooPlot* x2frame = MHaR.frame();
    TotalCosMass.plotOn(xframe);
    cosMassModel->plotOn(xframe);
    x2frame->Draw();
    c2->SaveAs("MassFit2.jpg");

    TCanvas* c3 = new TCanvas("c3", "c3", 1200, 600);
    RooPlot* yframe = cosR.frame();
    TotalCosMass.plotOn(yframe);
    cosModel.plotOn(yframe);
    yframe->Draw();
    c3->SaveAs("CosFit.jpg");

    TCanvas* c4 = new TCanvas("c4", "c4", 1200, 600);
    RooPlot* y2frame = cosR.frame();
    TotalCosMass.plotOn(y2frame);
    cosMassModel->plotOn(y2frame);
    y2frame->Draw();
    c4->SaveAs("CosFit2.jpg");
    
    
}
