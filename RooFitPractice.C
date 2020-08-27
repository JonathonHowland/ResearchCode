#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <fstream>

using namespace RooFit;

void RooFitPractice(){

  TH1D* cosHist = new TH1D("hist", "hist", 100, -1., 1.);
  ifstream cosData;
  cosData.open("cos.dat");
  double cos = 0;
  for(int i=0; i<30000; i++){
    cosData >> cos;
    cosHist->Fill(cos);
    }
    cosData.close();
    cosHist->Draw();

  //Observable
  RooRealVar m("m","m", 0., 10.);

  //Parameters
  RooRealVar mean("mean", "mean", 0.5, 0., 1.);
  RooRealVar RMS("RMS", "RMS", 5, 0., 100.);

  //Gaussian distribution
  RooGaussian gauss("gauss", "gauss", m, mean, RMS);

  //Background
  RooRealVar argpar("argpar", "argpar", -20., -100., -1.);
  RooArgusBG background("background", "background", m, RooConst(5.291), argpar);

  //Signal + Background
  RooRealVar ngauss("ngauss","signal events", 200., 0., 10000);
  RooRealVar nbg("nbg", "background events", 800., 0., 10000);
  RooAddPdf model("model", "model", RooArgList(gauss,background), RooArgList(ngauss, nbg));

  //Generate data
  RooDataSet *data = model.generate(m, 2000);

  //Fit
  model.fitTo(*data);

  //Plotting
  RooPlot* mframe = m.frame();
  data->plotOn(mframe);
  model.plotOn(mframe);
  model.plotOn(mframe, Components(background), LineStyle(ELineStyle::kDashed));
  model.plotOn(mframe, Components(gauss), LineColor(kRed));
  
  //mframe->Draw();

//RooKeysPdf
  

}
