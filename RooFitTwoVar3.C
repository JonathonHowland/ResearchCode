//Run the experiments twice (once with pdfs with signal, once without). Generate entirely new pdfs (and new coefficients) for both. Look at the likelihood ratios. Plot both likelihood ratios (with different pdfs and with different data) vs number of signal events. Keep the same number of total events in all trials.

//For a given Pdf(x) and true number of events for each experiment generate the same number of events (either background only or all three) from the pdfs and fit them (background only and all three). Record ln(L) for both fits. Calculate the difference in the ln(L)'s for the background only and all three.

//Check order on TTrees and TFiles. Find the region with fairly small sigma(sqrt(-2LLR)), and run many simulations in that area. 
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

void RooFitTwoVar3(){
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

  //Add the pdfs to create the model without signal

  RooAddPdf *cosMassBackModel = new RooAddPdf("cosMassBackModel", "cosMassBackModel", RooArgList(*cosMassPdf2, *cosMassPdf3), RooArgList(ncos22, ncos23));

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //Fit multiple times for different total numbers of events
  double nsig3, nsig2, lnsig3, lnback3, lnsig2, lnback2, diff3, diff2;
  
  //Create a TTree and File
  TFile fout("Likelihood.root", "recreate");
  TTree tout("tout","tout");
    tout.Branch("nsig3", &nsig3);
    tout.Branch("lnsig3", &lnsig3);
    tout.Branch("lnback3", &lnback3);
    tout.Branch("diff3", &diff3);
    tout.Branch("nsig2", &nsig2);
    tout.Branch("lnsig2", &lnsig2);
    tout.Branch("lnback2", &lnback2);
    tout.Branch("diff2", &diff2);
  
    ofstream results;
    results.open("3DataSets.dat");
    
    int TotalN = 1000; //Total number of events for a trial
    for(int i=0; i<20; i++){

      //Probabilities
      double prob1 = 0.01*i; //Signal probability
      
      double prob2 = 0.5*(1.0 - prob1); //Background probability
      double prob3 = prob2;

      //Numbers of each type of event
      int N1 = TotalN*prob1;
      int N2 = TotalN*prob2;
      int N3 = TotalN*prob3;

      int N22 = TotalN*0.5;
      int N23 = TotalN*0.5;
  
      //Fill the datasets from the histograms (all 3 types of events)
      for(int k=0; k<N1; k++){
	double cos = 0.0;
	double mass = 0.0;
	cosMassHist1->GetRandom2(cos, mass);
	cosR.setVal(cos);
	MHaR.setVal(mass);
	TotalCosMass.add(RooArgSet(cosR, MHaR));    
      }
      for(int k=0; k<N2; k++){
	double cos = 0.0;
	double mass = 0.0;
	cosMassHist2->GetRandom2(cos, mass);
	cosR.setVal(cos);
	MHaR.setVal(mass);
	TotalCosMass.add(RooArgSet(cosR, MHaR));    
      }
      for(int k=0; k<N3; k++){
	double cos = 0.0;
	double mass = 0.0;
	cosMassHist3->GetRandom2(cos, mass);
	cosR.setVal(cos);
	MHaR.setVal(mass);
	TotalCosMass.add(RooArgSet(cosR, MHaR));    
      }

      //Fit with both cos and mass
      RooFitResult* rcosMass = cosMassModel->fitTo(TotalCosMass, Save());
      RooFitResult* rcosMassBack = cosMassBackModel->fitTo(TotalCosMass, Save());
      ///////////////////////////////////////////////////////////////////////////////////////////
      //Check likelihood differences
      nsig3 = N1;
      lnsig3 = rcosMass->minNll();
      lnback3 = rcosMassBack->minNll();
      diff3 = lnsig3 - lnback3;
  
      //Fill the datasets from the histograms (only background events)
      for(int k=0; k<N22; k++){
	double cos = 0.0;
	double mass = 0.0;
	cosMassHist2->GetRandom2(cos, mass);
	cosR.setVal(cos);
	MHaR.setVal(mass);
	TotalCosMassBack.add(RooArgSet(cosR, MHaR));    
      }
      for(int k=0; k<N23; k++){
	double cos = 0.0;
	double mass = 0.0;
	cosMassHist3->GetRandom2(cos, mass);
	cosR.setVal(cos);
	MHaR.setVal(mass);
	TotalCosMassBack.add(RooArgSet(cosR, MHaR));    
      }

      //Fit with both cos and mass
      RooFitResult* rcosMass2 = cosMassModel->fitTo(TotalCosMassBack, Save());
      RooFitResult* rcosMassBack2 = cosMassBackModel->fitTo(TotalCosMassBack, Save());
      ///////////////////////////////////////////////////////////////////////////////////////////
      //Check likelihood differences
      nsig2 = N1;
      lnsig2 = rcosMass2->minNll();
      lnback2 = rcosMassBack2->minNll();
      diff2 = lnsig2 - lnback2;
      tout.Fill();
      results << N1 << endl;
      results << lnsig3 << endl;
      results << lnback3 << endl;
      results << diff3 << endl;
      results << N1 << endl;
      results << lnsig2 << endl;
      results << lnback2 << endl;
      results << diff2 << endl;
    }
    results.close();
    fout.Close();
}
