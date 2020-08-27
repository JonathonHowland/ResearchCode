//Treat 2 as backgrounds, one as signal
//Keep background numbers constant, vary number of signal events, check accuracy (graph known event number versus measured event number). Also plot (signals observed)/error versus known signal
//RooNDKeysPdf for two dimensions
//Create multiple data sets (from RooKeysPdf or RestFrames) and look at distributions

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TGraph.h"
#include <vector>
#include <cmath>

using namespace RooFit;

void RooFitAnalysis(){
  
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

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  //Declaring RooFit stuff
  //Observables
  RooRealVar cosR("cosR", "cosR", -1.,1.);
  RooRealVar MHaR("MHa", "MHa", 0., 500.);

  //Datasets
  //Treat cos1D and cos2D as backgroundd
  RooDataSet cos1D("name","title", RooArgSet(cosR));
  RooDataSet cos2D("name","title", RooArgSet(cosR));
  RooDataSet cos3D("name","title", RooArgSet(cosR));
  RooDataSet TotalCosBkg("bkg_cos","bkg_cos",RooArgSet(cosR));
  RooDataSet TotalCos("name","title", RooArgSet(cosR));
  RooDataSet MHa1D("name","title", RooArgSet(MHaR));
  RooDataSet MHa2D("name","title", RooArgSet(MHaR));
  RooDataSet MHa3D("name","title", RooArgSet(MHaR));
  RooDataSet TotalMass("totalm", "totalm", RooArgSet(MHaR));

  ///////////////////////////////////////////////////////////////////////////////////////////

  //Build the data sets for the pdf
  int Nentries = t1->GetEntries();
  for(int i = 0; i < Nentries/10; i++){
    t1->GetEntry(i);
    //cout << "cos1 " << cos1 << " cos2 " << cos2 << " cos3 " << cos3 << endl;
    //Decay angle
    cosR.setVal(cos1);
    cos1D.add(cosR);
    cosR.setVal(cos2);
    cos2D.add(cosR);
    cosR.setVal(cos3);
    cos3D.add(cosR);

    //Mass
    MHaR.setVal(MHa1);
    MHa1D.add(MHaR);
    MHaR.setVal(MHa2);
    MHa2D.add(MHaR);
    MHaR.setVal(MHa3);
    MHa3D.add(MHaR);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////////////

  //Declaring Pdfs and parameters

  //Decay angle
  
  //Parameters
  RooRealVar mean("mean", "mean", 0., -1., 1.);
  RooRealVar RMS("RMS", "RMS", 5, 0., 100.);
  RooRealVar ncos1("ncos1","signal events", 5000., 0., 300000);
  RooRealVar ncos2("ncos2", "signal events", 5000., 0., 300000);
  RooRealVar ncos3("ncos3","signal events", 5000., 0., 300000);
  
  //Build a numeric pdf
  RooKeysPdf cos1Pdf("cos1Pdf", "title", cosR, cos1D); 
  RooKeysPdf cos2Pdf("cos2Pdf", "title", cosR, cos2D); 
  RooKeysPdf cos3Pdf("cos3Pdf", "title", cosR, cos3D);
  //RooKeysPdf TotalCosPdf("name", "title", cos1R, TotalCos);

  RooAddPdf model("model", "model", RooArgList(cos1Pdf, cos2Pdf, cos3Pdf), RooArgList(ncos1, ncos2, ncos3));

  //Mass
  
  //Parameters
  RooRealVar nMHa1("nMHa1","signal events", 5000., 0., 300000);
  RooRealVar nMHa2("nMHa2", "signal events", 5000., 0., 300000);
  RooRealVar nMHa3("nMHa3","signal events", 5000., 0., 300000);
  
  //Build a numeric pdf
  RooKeysPdf MHa1Pdf("name", "title", MHaR, MHa1D); 
  RooKeysPdf MHa2Pdf("name2", "title", MHaR, MHa2D); 
  RooKeysPdf MHa3Pdf("name3", "title", MHaR, MHa3D);

  //Create a RooWorkspace and save pdf to a file
  RooWorkspace *w = new RooWorkspace("w", "workspace");
  w->import(cos1Pdf);
  w->import(cos2Pdf);
  w->import(cos3Pdf);
  w->import(MHa1Pdf);
  w->import(MHa2Pdf);
  w->import(MHa3Pdf);
  w->writeToFile("1DPdf2.root", "RECREATE");

  cout <<"Pdfs written to 1DPdf.root" << endl;

  RooAddPdf modelMass("mass model", "mass model", RooArgList(MHa1Pdf, MHa2Pdf, MHa3Pdf), RooArgList(nMHa1, nMHa2, nMHa3));

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int NTOT = 10000;
  
  //Filling the data set to fit to (3 methods)

  double prob1 = .1;
  double prob2 = .2;
  double prob3 = 1. - prob1 - prob2;

  //Choose 1 to use the TTree. Choose 2 to use RestFrames. Choose anything else to use the Pdfs.
  int choice = 0;

  if (choice == 1){
    //Fill the data sets from TTree

    for (int k = 0; k<10000; k++){
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
    
    }
  }
  else if(choice==2){
    //Build the data set from other RestFrames events
    for(int i=Nentries/10; i<Nentries/5; i++){
      t1->GetEntry(i);
      double r = gRandom->Rndm();
      if (r<prob1){
	cosR.setVal(cos1);
	TotalCos.add(cosR);
	MHaR.setVal(MHa1);
	TotalMass.add(MHaR);
      }
      else if(r<(prob1 + prob2)){
	  cosR.setVal(cos2);
	  TotalCos.add(cosR);
	  MHaR.setVal(MHa2);
	  TotalMass.add(MHaR);
	}
	else{
	  cosR.setVal(cos3);
	  TotalCos.add(cosR);
	  MHaR.setVal(MHa3);
	  TotalMass.add(MHaR);
	}
	}
    }
    else{
      //Build the data set from the pdfs
      int N1 = prob1*double(NTOT);
      int N2 = prob2*double(NTOT);
      int N3 = prob3*double(NTOT);

      TotalCos.append( *cos1Pdf.generate(RooArgSet(cosR),N1) );
      TotalCos.append( *cos2Pdf.generate(RooArgSet(cosR),N2) );
      TotalCos.append( *cos3Pdf.generate(RooArgSet(cosR),N3) );

      TotalMass.append( *MHa1Pdf.generate(RooArgSet(MHaR),N1) );
      TotalMass.append( *MHa2Pdf.generate(RooArgSet(MHaR),N2) );
      TotalMass.append( *MHa3Pdf.generate(RooArgSet(MHaR),N3) );

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*//Fitting multiple times
  
    //Vectors for storing data for graphing
    vector<double> vNcos3;
    vector<double> vNcos3Err;
    vector<double> vNcos3Ratio;
    vector<double> vNcos3RealErr;

    //Counter for for loop
    int i = 0;
  
    //Build the background
    double mycos;
  
    for (int i =0; i < Nentries; i++){
      t1->GetEntry(i);
      cosR.setVal(cos1);
      TotalCosBkg.add(cosR);
      cosR.setVal(cos2);
      TotalCosBkg.add(cosR);
    }
    for (int Nsig = 1; Nsig < 10000; Nsig *= 2){

      RooDataSet TotalCos("total_cos","total_cos",RooArgSet(cosR));
      TotalCos.append(TotalCosBkg);
    
      //Build the signal
      for (int j = 0; j < Nsig; j++){
	t1->GetEntry(j);
	cosR.setVal(cos3);
	TotalCos.add(cosR);
      }
  
      //Fit
      model.fitTo(TotalCos);

      vNcos3.push_back(ncos3.getVal());
      vNcos3Err.push_back(ncos3.getError());
      vNcos3Ratio.push_back(ncos3.getError()/ncos3.getVal());
      vNcos3RealErr.push_back(pow(2.0, i-1) - ncos3.getVal());
      i++;
      }*/

  modelMass.fitTo(TotalMass);

    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Plotting

    RooPlot* xframe = MHaR.frame();
    TotalMass.plotOn(xframe);
    modelMass.plotOn(xframe);
    //model.plotOn(xframe, Components(cos1Pdf), LineColor(kRed));
    //model.plotOn(xframe, Components(cos2Pdf), LineColor(kBlue-8));
    //model.plotOn(xframe, Components(cos3Pdf), LineColor(kBlack));
    //gauss.plotOn(mframe, Components(background), LineStyle(ELineStyle::kDashed));
    //gauss.plotOn(xframe, Components(gauss), LineColor(kRed));
  
    xframe->Draw();
  
    /*
    //Plotting graphs
    TGraph* Errs = new TGraph(i, &vNcos3[0], &vNcos3Err[0]);
    TGraph* ErrRatio = new TGraph(i, &vNcos3[0], &vNcos3Ratio[0]);
    TGraph* RealErr = new TGraph(i, &vNcos3[0], &vNcos3RealErr[0]);

    ErrRatio->Draw();
    */
  }
