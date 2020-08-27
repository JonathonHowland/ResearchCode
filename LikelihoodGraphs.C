
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"

using namespace std;

double w = 1200;
double h = 600;

TGraph* Graph(int, vector<double>, vector<double>, int);
TGraphErrors* ErrorGraph(int, vector<double>, vector<double>, vector<double>);
void DrawM(TMultiGraph*, TLegend*, int);
void DrawE(TGraphErrors* Eg, int i);

void LikelihoodGraphs(){
  /*
  //Generate TTree from file
  TFile* f = new TFile("Likelihood.root", "READ");
  cout << "File opened successfully" << endl;
  TTree* tout = (TTree*)f->Get("tout");
  cout << "TTree obtained successfully" << endl;
  double nsig3, nsig2, lnsig3, lnback3, lnsig2, lnback2, diff3, diff2;
  tout->SetBranchAddress("nsig3", &nsig3);
  tout->SetBranchAddress("lnsig3", &lnsig3);
  tout->SetBranchAddress("lnback3", &lnback3);
  tout->SetBranchAddress("diff3", &diff3);
  tout->SetBranchAddress("nsig2", &nsig2);
  tout->SetBranchAddress("lnsig2", &lnsig2);
  tout->SetBranchAddress("lnback2", &lnback2);
  tout->SetBranchAddress("diff2", &diff2);

  cout << "TTree imported successfully" << endl;*/
  
  //Keep track of things
  //int nentries = tout->GetEntries();

  ifstream data;
  data.open("3DataSets.dat");
  int count = 0;
  vector<double> vN3;
  vector<double> vSig3;
  vector<double> vBack3;
  vector<double> vDiff3;
  vector<double> vSigma3;
  vector<double> vN2;
  vector<double> vSig2;
  vector<double> vBack2;
  vector<double> vDiff2;
  vector<double> vSigma2;
  for (int i=0; i<20; i++){
    double x;
    data >> x;
    vN3.push_back(x);
    data >> x;
    vSig3.push_back(x);
    data >> x;
    vBack3.push_back(x);
    data >> x;
    vDiff3.push_back(x);
    vSigma3.push_back(sqrt(2.0*abs(x)));
    data >> x;
    vN2.push_back(x);
    data >> x;
    vSig2.push_back(x);
    data >> x;
    vBack2.push_back(x);
    data >> x;
    vDiff2.push_back(x);
    vSigma2.push_back(sqrt(2.0*abs(x)));
    count++;
  }
 ifstream data2;
  data2.open("SameProportion.dat");
  int count2 = 0;
  vector<double> vTotalN;
  vector<double> vMedSigProb;
  vector<double> vMedBackProb;
  vector<double> vIQSigProb;
  vector<double> vIQBackProb;
  for (int i=0; i<9; i++){
    double x;
    data2 >> x;
    vTotalN.push_back(x);
    data2 >> x;
    vMedSigProb.push_back(x);
    data2 >> x;
    vMedBackProb.push_back(x);
    data2 >> x;
    vIQSigProb.push_back(0.5*x);
    data2 >> x;
    vIQBackProb.push_back(0.5*x);
    count2++;
  }
  
 ifstream data3;
  data3.open("SameN.dat");
  int count3 = 0;
  vector<double> vTotalN2;
  vector<double> vMedSigN;
  vector<double> vMedBackN;
  vector<double> vIQSigN;
  vector<double> vIQBackN;
  for (int i=0; i<9; i++){
    double x;
    data3 >> x;
    vTotalN2.push_back(x);
    data3 >> x;
    vMedSigN.push_back(x);
    data3 >> x;
    vMedBackN.push_back(x);
    data3 >> x;
    vIQSigN.push_back(0.5*x);
    data3 >> x;
    vIQBackN.push_back(0.5*x);
    count3++;
  }

  /////////////////////////////////////////////////////////
  //Graphing
  //Graphs for one trial

  TGraph* gSig3 = Graph(count, vN3, vSig3, 0);
  TGraph* gBack3 = Graph(count, vN3, vBack3, 1);
  TGraph* gDiff3 = Graph(count, vN3, vDiff3, 0);
  TGraph* gSigma3 = Graph(count, vN3, vSigma3, 0);
  TGraph* gSig2 = Graph(count, vN2, vSig2, 0);
  TGraph* gBack2 = Graph(count, vN2, vBack2, 1);
  TGraph* gDiff2 = Graph(count, vN2, vDiff2, 1);
  TGraph* gSigma2 = Graph(count, vN2, vSigma2, 1);

  TMultiGraph* g3Sets = new TMultiGraph();
  g3Sets->Add(gSig3);
  g3Sets->Add(gBack3);
  g3Sets->SetTitle("Negative log likelihoods with 1,000 total events (signal and background); Number of signal events measured; Likelihood ratio");

  TMultiGraph* g2Sets = new TMultiGraph();
  g2Sets->Add(gSig2);
  g2Sets->Add(gBack2);
  g2Sets->SetTitle("Negative log likelihoods with 1,000 total events (background only); Number of signal events measured; Negative log likelihood");

  TMultiGraph* gDiffs = new TMultiGraph();
  gDiffs->Add(gDiff3);
  gDiffs->Add(gDiff2);
  gDiffs->SetTitle("Likelihood ratios with 1,000 total events; Number of signal events measured; Likelihood ratio");

  TMultiGraph* gSigmas = new TMultiGraph();
  gSigmas->Add(gSigma3);
  gSigmas->Add(gSigma2);
  gSigmas->SetTitle("Number of sigmas (derived from likelihood ratios) with 1,000 total events; Number of signal events measured; Sigmas");
  
  //Make a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->SetHeader("Legend", "C");
  legend->AddEntry(gSig3, "Fitting with signal and background", "l");
  legend->AddEntry(gBack3, "Fitting with only background");
  
  //Make another legend
  TLegend *legend2 = new TLegend(0.1, 0.1, 0.3, 0.3);
  legend2->SetHeader("Legend", "C");
  legend2->AddEntry(gDiff3, "Data contains signal and background", "l");
  legend2->AddEntry(gDiff2, "Data contains only background");
  
  //Make a third legend
  TLegend *legend3 = new TLegend(0.1, 0.9, 0.3, 0.7);
  legend3->SetHeader("Legend", "C");
  legend3->AddEntry(gSigma3, "Data contains signal and background", "l");
  legend3->AddEntry(gSigma2, "Data contains only background");

  DrawM(g3Sets, legend, 1);

  
  DrawM(g2Sets, legend, 2);

  DrawM(gDiffs, legend2, 3);

  DrawM(gSigmas, legend3, 6);
  
  ////////////////////////////////////////////////////////////////////////
  //Graphs for many trials
  TGraphErrors* SigEG = ErrorGraph(count2, vTotalN, vMedSigProb, vIQSigProb);
  TGraphErrors* SigEG2 = ErrorGraph(count3, vTotalN2, vMedSigN, vIQSigN);

  SigEG->SetTitle("Means of sigmas of distributions with different numbers of total events; Number of total events; Means of sigmas with standard deviations as error bars");
  SigEG2->SetTitle("Means of sigmas of distributions with different signal probabilities; Fraction of (signal events)/(total events); Means of sigmas with standard deviations as error bars");
  
  DrawE(SigEG, 4);
  DrawE(SigEG2, 5);
  
}

TGraph* Graph(int count, vector<double> vX, vector<double> vY, int choice){
  TGraph* g = new TGraph(count, &vX[0], &vY[0]);
  g->SetLineWidth(2);
  if(choice == 0){
  g->SetLineColor(kRed);
  }
  else if(choice == 1){
    g->SetLineColor(kBlue);
  }
  return g;
  
}

TGraphErrors* ErrorGraph(int count, vector<double> vX, vector<double> vY, vector<double> vE){
  TGraphErrors* eg = new TGraphErrors(count, &vX[0], &vY[0], 0, &vE[0]);
  eg->SetLineWidth(2);
  return eg;
  
}

void DrawM(TMultiGraph* mg, TLegend* leg, int i){
  string str = to_string(i);
  TCanvas* c1 = new TCanvas(&str[0], &str[0], w, h);
  c1->cd();
  mg->Draw("APC");
  leg->Draw();
  string name = "Graph";
  name += str;
  name += ".jpg";
  c1->SaveAs(&name[0]);
}

void DrawE(TGraphErrors* Eg, int i){
  string str = to_string(i);
  TCanvas* c1 = new TCanvas(&str[0], &str[0], w, h);
  c1->cd();
  Eg->SetMarkerColor(4);
  Eg->SetMarkerStyle(21);
  Eg->Draw("AP");
  string name = "Graph";
  name += str;
  name += ".jpg";
  c1->SaveAs(&name[0]);
}
