//Code to create the 1D pdfs

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
#include <vector>
#include <cmath>

using namespace RooFit;

void RooFit1DPdf(){

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

  ///////////////////////////////////////////////////////////////////////////////////
  //Declare RooFit variables and datasets
  //Observables
  RooRealVar cosR("cosR", "cosR", -1.,1.);
  
  RooRealVar MHaR("MHaR", "MHa", 0., 500.);

  //Datasets
  RooDataSet cos1D("cos1D","title", RooArgSet(cosR));
  RooDataSet cos2D("cos2D","title", RooArgSet(cosR));
  RooDataSet cos3D("cos3D","title", RooArgSet(cosR));
  RooDataSet TotalCos("TotalCos","title", RooArgSet(cosR));
  RooDataSet MHa1D("MHa1D","title", RooArgSet(MHaR));
  RooDataSet MHa2D("MHa2D","title", RooArgSet(MHaR));
  RooDataSet MHa3D("MHa3D","title", RooArgSet(MHaR));
  RooDataSet TotalMass("TotalMass", "totalm", RooArgSet(MHaR));

  //////////////////////////////////////////////////////////////////////////////////////
  //Build the 1D pdfs

  //Build the datasets for the pdfs
  int Nentries = t1->GetEntries();
  cout << "TTree contains " <<  Nentries << " entries" << endl;
  cout << "Using " << Nentries/1000 << " entries" << endl;
    
  for(int i=0; i<Nentries/1000; i++){
    t1->GetEntry(i);
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

  cout << "Data sets built successfully" << endl;
  //Build the component pdfs
  RooKeysPdf cos1Pdf("cos1Pdf", "cos1Pdf", cosR, cos1D);
  RooKeysPdf cos2Pdf("cos2Pdf", "cos2Pdf", cosR, cos2D);
  RooKeysPdf cos3Pdf("cos3Pdf", "cos3Pdf", cosR, cos3D);
  RooKeysPdf MHa1Pdf("MHa1Pdf", "MHa1Pdf", MHaR, MHa1D);
  RooKeysPdf MHa2Pdf("MHa2Pdf", "MHa2Pdf", MHaR, MHa2D);
  RooKeysPdf MHa3Pdf("MHa3Pdf", "MHa3Pdf", MHaR, MHa3D);


  cout <<"1D Pdfs built successfully" << endl;

  //////////////////////////////////////////////////////////////////////////////////////

  //Create a RooWorkspace and save pdf to a file
  RooWorkspace *w = new RooWorkspace("w", "workspace");
  w->import(cos1Pdf);
  w->import(cos2Pdf);
  w->import(cos3Pdf);
  w->import(MHa1Pdf);
  w->import(MHa2Pdf);
  w->import(MHa3Pdf);
  w->writeToFile("1DPdf.root");

  cout <<"Pdfs written to 1DPdf.root" << endl;
  
 
}
