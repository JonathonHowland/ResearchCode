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

void RooFitTwoVar2(){
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

  //////////////////////////////////////////////////////////////////////////////////////////////////////
 
  //Fit multiple times

  //Vectors for storing data for graphing
  vector<double> vNcos;
  vector<double> vNcosErr;
  vector<double> vNcosRatio;
  vector<double> vNcosRealErr;
  vector<double> vNmass;
  vector<double> vNmassErr;
  vector<double> vNmassRatio;
  vector<double> vNmassRealErr;
  vector<double> vNcosMass;
  vector<double> vNcosMassErr;
  vector<double> vNcosMassRatio;
  vector<double> vNcosMassRealErr;

  //Counter
  int i = 0;

  //Build data sets
  
  int N = 10000;
  double prob1 = 0.1;
  double prob2 = (1.0-prob1)/2.0;
  double prob3 = 1.0 - prob1 - prob2;

  int N1 = prob1*N;
  int N2 = prob2*N;
  int N3 = prob3*N;
  
  while(N1<10000){
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
    
    cout << "Do fit number " << i+1 << "..." << endl;

    //Fit with cos only
    cosModel.fitTo(TotalCosMass);
    vNcos.push_back(ncos1.getVal());
    vNcosErr.push_back(ncos1.getError());
    vNcosRatio.push_back(ncos1.getError()/ncos1.getVal());
    vNcosRealErr.push_back(N1 - ncos1.getVal());

    //Fit with mass only
    massModel.fitTo(TotalCosMass);
    vNmass.push_back(ncos1.getVal());
    vNmassErr.push_back(ncos1.getError());
    vNmassRatio.push_back(ncos1.getError()/ncos1.getVal());
    vNmassRealErr.push_back(N1 - ncos1.getVal());

    //Fit with both cos and mass
    cosMassModel->fitTo(TotalCosMass);
    vNcosMass.push_back(ncos1.getVal());
    vNcosMassErr.push_back(ncos1.getError());
    vNcosMassRatio.push_back(ncos1.getError()/ncos1.getVal());
    vNcosMassRealErr.push_back(N1 - ncos1.getVal());
    std::cout << "Fit number " << i+1 << " completed successfully" << std::endl;
    TotalCosMass.reset();
    i++;
    prob1 = 1.5*prob1;
    prob2 = (1.0-prob1)/2.0;
    prob3 = 1.0 - prob1 - prob2;

    N1 = prob1*N;
    N2 = prob2*N;
    N3 = prob3*N;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////
  //Graphing
  Int_t palette[5];
  palette[0] = (kMagenta+2);
  palette[1] = (kRed);
  palette[2] = (kYellow);
  palette[3] = (kAzure+6);
  palette[4] = kBlue;
  gStyle->SetPalette(5,palette);
    

  //Make TGraphs
  TGraph *cosErr = new TGraph(i, &vNcos[0], &vNcosErr[0]);
  TGraph *cosRatio = new TGraph(i, &vNcos[0], &vNcosRatio[0]);
  TGraph *cosRealErr = new TGraph(i, &vNcos[0], &vNcosRealErr[0]);
  TGraph *massErr = new TGraph(i, &vNmass[0], &vNmassErr[0]);
  TGraph *massRatio = new TGraph(i, &vNmass[0], &vNmassRatio[0]);
  TGraph *massRealErr = new TGraph(i, &vNmass[0], &vNmassRealErr[0]);
  TGraph *cosMassErr = new TGraph(i, &vNcosMass[0], &vNcosMassErr[0]);
  TGraph *cosMassRatio = new TGraph(i, &vNcosMass[0], &vNcosMassRatio[0]);
  TGraph *cosMassRealErr = new TGraph(i, &vNcosMass[0], &vNcosMassRealErr[0]);

  //Make 'em purdy
  cosErr->SetLineColor(palette[0]);
  cosErr->SetLineWidth(3);
  cosRatio->SetLineColor(palette[0]);
  cosRatio->SetLineWidth(3);
  cosRealErr->SetLineColor(palette[0]);
  cosRealErr->SetLineWidth(3);
  
  massErr->SetLineColor(palette[1]);
  massErr->SetLineWidth(3);
  massRatio->SetLineColor(palette[1]);
  massRatio->SetLineWidth(3);
  massRealErr->SetLineColor(palette[1]);
  massRealErr->SetLineWidth(3);
  
  cosMassErr->SetLineColor(palette[3]);
  cosMassErr->SetLineWidth(3);
  cosMassRatio->SetLineColor(palette[3]);
  cosMassRatio->SetLineWidth(3);
  cosMassRealErr->SetLineColor(palette[3]);
  cosMassRealErr->SetLineWidth(3);

  cosRatio->SetTitle("cos");
  massRatio->SetTitle("mass");
  cosMassRatio->SetTitle("cosMass");
			
  //Build multi-graphs
  TMultiGraph *Errs = new TMultiGraph();
  Errs->Add(cosErr);
  Errs->Add(massErr);
  Errs->Add(cosMassErr);
  Errs->SetTitle("Errors on number of signal events; Number of signal events (true); Error");

  TMultiGraph *Ratio = new TMultiGraph();
  Ratio->Add(cosRatio);
  Ratio->Add(massRatio);
  Ratio->Add(cosMassRatio);
  Ratio->SetTitle("Ratio of errors to number of signal events; Number of signal events (true); Ratio");

  TMultiGraph *RealErr = new TMultiGraph();
  RealErr->Add(cosRealErr);
  RealErr->Add(massRealErr);
  RealErr->Add(cosMassRealErr);

  //Draw the graphs
  TCanvas *c1 = new TCanvas("c1", "Errors", 600, 600);
  TCanvas *c2 = new TCanvas("c2", "Error Ratios", 600, 600);
  TCanvas *c3 = new TCanvas("c3", "Real Errors", 600, 600);
  TCanvas *c4 = new TCanvas("c4", "Legend", 600, 600);

  //Make a legend
  TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  legend->SetHeader("Legend", "C");
  legend->AddEntry(cosErr, "Using only cos", "l");
  legend->AddEntry(massErr, "Using only mass");
  legend->AddEntry(cosMassErr, "Using both cos and mass");
  
  c1->cd();
  Errs->Draw("APC");
  //legend->Draw();

  c2->cd();
  Ratio->Draw("APC");
  //legend->Draw();

  c3->cd();
  RealErr->Draw("APC");
  //legend->Draw();

  c4->cd();
  legend->Draw();
  
  c1->SaveAs("ErrsSpread.jpg");
  c2->SaveAs("RatioSpread.jpg");
  c3->SaveAs("RealErr.jpg");
  c4->SaveAs("Legend.jpg");
}

  /*
  //Generate data from 2D Pdf
  
  TotalCosMass.append( *cosMassPdf1->generate(RooArgSet(cosR, MHaR),N1) );
  TotalCosMass.append( *cosMassPdf2->generate(RooArgSet(cosR, MHaR),N2) );
  TotalCosMass.append( *cosMassPdf3->generate(RooArgSet(cosR, MHaR),N3) );
  */

  /* 
  //Plotting
  RooPlot* xframe = MHaR.frame();
  TotalCosMass.plotOn(xframe);
  cosMassModel->plotOn(xframe);
  //model.plotOn(xframe, Components(cos1Pdf), LineColor(kRed));
  //model.plotOn(xframe, Components(cos2Pdf), LineColor(kBlue-8));
  //model.plotOn(xframe, Components(cos3Pdf), LineColor(kBlack));
  //gauss.plotOn(mframe, Components(background), LineStyle(ELineStyle::kDashed));
  //gauss.plotOn(xframe, Components(gauss), LineColor(kRed));
  
  xframe->Draw();
  */

/* for (int k = 0; k<10000; k++){
    t1->GetEntry(k);
    double r = gRandom->Rndm();
    if(r < prob1){
      cosR.setVal(cos1);
      MHaR.setVal(MHa1);
    }
    else if(r < prob1 + prob2){
      cosR.setVal(cos2);
      MHaR.setVal(MHa2);
    }
    else{
      cosR.setVal(cos3);
      MHaR.setVal(MHa3);
    }

    TotalCos.add(cosR);
    TotalMass.add(MHaR);
    TotalCosMass.add(RooArgSet(cosR, MHaR));
    }*/

  /* 
  //Fill the data sets from TTree

  double prob1 = 0.01;
  double prob2 = 0.3;
  double prob3 = 1.0 - prob1 - prob2;
  while(prob1<0.3){
    int j = i;
    int k = j;
    while((j+10000*i)<(10000+10000*k)){
      t1->GetEntry(j+10000*i);
	
      double r = gRandom->Rndm();
      if(r < prob1){
	cosR.setVal(cos1);
	MHaR.setVal(MHa1);
      }
      else if(r < prob1 + prob2){
	cosR.setVal(cos2);
	MHaR.setVal(MHa2);
      }
      else{
	cosR.setVal(cos3);
	MHaR.setVal(MHa3);
      }

      TotalCos.add(cosR);
      TotalMass.add(MHaR);
      TotalCosMass.add(RooArgSet(cosR, MHaR));

      j++;
    }
  */
