#include "RestFrames/RestFrames.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TColor.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include <fstream>

using namespace RestFrames;

void TreePlots(const std::string& output_name =
			   "output_H_to_WlnuWlnu.root"){
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 1200);
  //gStyle->SetOptStat(0);

    // set particle masses and widths
  double mW   = 80.385;  // GeV, PDG 2016
  double wW = 2.085;
  double mL   = 0.106;   // muons
  double mN   = 0.;
  std::vector<double> mG; // vary neutral Higgs mass
  std::vector<double> wG; // vary Graviton mass
  wG.push_back(0.1);
  wG.push_back(0.2);
  wG.push_back(0.3);
  wG.push_back(0.4);
  wG.push_back(0.5);
  mG.push_back(125.);
  mG.push_back(400.);
  mG.push_back(750.);
  mG.push_back(1000.);
  mG.push_back(1500.);
  int Nmass = 5;

  // Number of events to generate (per H mass)
  int Ngen = 10000000;

  /////////////////////////////////////////////////////////////////////////////////////////
  g_Log << LogInfo << "Initializing generator frames and tree..." << LogEnd;
  /////////////////////////////////////////////////////////////////////////////////////////
  ppLabGenFrame     LAB_Gen("LAB_Gen","LAB");
   ResonanceGenFrame     G_Gen("G","G");
  ResonanceGenFrame Ha_Gen("Ha_Gen","H_{a}");
  VisibleGenFrame   La_Gen("La_Gen","#it{l}_{a}");
  VisibleGenFrame Lb_Gen("Lb_Gen","#it{l}_{b}");
  VisibleGenFrame   Lc_Gen("Lc_Gen","#it{l}_{c}");

  ppLabGenFrame     LAB_Gen2("LAB_Gen","LAB");
   ResonanceGenFrame     G_Gen2("G","G");
  ResonanceGenFrame Ha_Gen2("Ha_Gen","H_{a}");
  VisibleGenFrame   La_Gen2("La_Gen","#it{l}_{a}");
  VisibleGenFrame Lb_Gen2("Lb_Gen","#it{l}_{b}");
  VisibleGenFrame   Lc_Gen2("Lc_Gen","#it{l}_{c}");

  ppLabGenFrame     LAB_Gen3("LAB_Gen","LAB");
   ResonanceGenFrame     G_Gen3("G","G");
  ResonanceGenFrame Ha_Gen3("Ha_Gen","H_{a}");
  VisibleGenFrame   La_Gen3("La_Gen","#it{l}_{a}");
  VisibleGenFrame Lb_Gen3("Lb_Gen","#it{l}_{b}");
  VisibleGenFrame   Lc_Gen3("Lc_Gen","#it{l}_{c}");

  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//
  //Tree 1
  
  LAB_Gen.SetChildFrame(G_Gen);
  G_Gen.AddChildFrame(Ha_Gen);
  G_Gen.AddChildFrame(Lc_Gen);
  Ha_Gen.AddChildFrame(La_Gen);
  Ha_Gen.AddChildFrame(Lb_Gen);

  if(LAB_Gen.InitializeTree())
    g_Log << LogInfo << "...Successfully initialized generator tree" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator tree" << LogEnd;
  
  //Tree 2
  
  LAB_Gen2.SetChildFrame(G_Gen2);
  G_Gen2.AddChildFrame(Ha_Gen2);
  G_Gen2.AddChildFrame(La_Gen2);
  Ha_Gen2.AddChildFrame(Lb_Gen2);
  Ha_Gen2.AddChildFrame(Lc_Gen2);

  if(LAB_Gen2.InitializeTree())
    g_Log << LogInfo << "...Successfully initialized generator tree" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator tree" << LogEnd;
    
  //Tree 3
  
  LAB_Gen3.SetChildFrame(G_Gen3);
  G_Gen3.AddChildFrame(Ha_Gen3);
  G_Gen3.AddChildFrame(Lb_Gen3);
  Ha_Gen3.AddChildFrame(Lc_Gen3);
  Ha_Gen3.AddChildFrame(La_Gen3);

  if(LAB_Gen3.InitializeTree())
    g_Log << LogInfo << "...Successfully initialized generator tree" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator tree" << LogEnd;
  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//

  if(mG.size() < 1) return;
  
  // set Higgs masses
  G_Gen.SetMass(mG[1]);
  G_Gen.SetWidth(wG[0]*mG[0]);
  // set W masses and widths
  Ha_Gen.SetMass(mW);                    Ha_Gen.SetWidth(wW);
  // set lepton and neutrino masses
  La_Gen.SetMass(0);                    Lb_Gen.SetMass(0);
  Lc_Gen.SetMass(0);                    


  if(LAB_Gen.InitializeAnalysis())
    g_Log << LogInfo << "...Successfully initialized generator analysis" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator analysis" << LogEnd;

    // set Higgs masses
  G_Gen2.SetMass(mG[1]);
  G_Gen2.SetWidth(wG[0]*mG[0]);
  // set W masses and widths
  Ha_Gen2.SetMass(mW);                    Ha_Gen2.SetWidth(wW);
  // set lepton and neutrino masses
  La_Gen2.SetMass(0);                    Lb_Gen2.SetMass(0);
  Lc_Gen2.SetMass(0);                    


  if(LAB_Gen2.InitializeAnalysis())
    g_Log << LogInfo << "...Successfully initialized generator analysis" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator analysis" << LogEnd;

    // set Higgs masses
  G_Gen3.SetMass(mG[1]);
  G_Gen3.SetWidth(wG[0]*mG[0]);
  // set W masses and widths
  Ha_Gen3.SetMass(mW);                    Ha_Gen3.SetWidth(wW);
  // set lepton and neutrino masses
  La_Gen3.SetMass(0);                    Lb_Gen3.SetMass(0);
  Lc_Gen3.SetMass(0);                    


  if(LAB_Gen3.InitializeAnalysis())
    g_Log << LogInfo << "...Successfully initialized generator analysis" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator analysis" << LogEnd;
  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////
  g_Log << LogInfo << "Initializing reconstruction frames and trees..." << LogEnd;
  //////////////////////////////////////////////////////////////////////////////////////////

  LabRecoFrame LAB_Reco("LAB_Reco","LAB");
  DecayRecoFrame G_Reco("G","G");
  DecayRecoFrame Ha_Reco("Ha","H_{a}");
  VisibleRecoFrame La_Reco("La","#it{l}_{a}");
  VisibleRecoFrame Lb_Reco("Lb","#it{l}_{b}");
  VisibleRecoFrame Lc_Reco("Lc","#it{l}_{c}");

  LAB_Reco.SetChildFrame(G_Reco);
  G_Reco.AddChildFrame(Ha_Reco);
  G_Reco.AddChildFrame(Lc_Reco);
  Ha_Reco.AddChildFrame(La_Reco);
  Ha_Reco.AddChildFrame(Lb_Reco);

  LAB_Reco.InitializeTree();
  LAB_Reco.InitializeAnalysis();

  /////////////////////////////////////////////////////////////////////////////////////////

  TreePlot* treePlotR = new TreePlot("TreePlotR","TreePlot");
  TreePlot* treePlotG1 = new TreePlot("TreePlotG1","TreePlot");
  TreePlot* treePlotG2 = new TreePlot("TreePlotG2","TreePlot");
  TreePlot* treePlotG3 = new TreePlot("TreePlotG3","TreePlot");
  c1->cd();
  
  treePlotR->SetTree(LAB_Reco);
  treePlotR->Draw();
  //c1->SaveAs("RecoTree.jpg");
  
  treePlotG1->SetTree(LAB_Gen);
  treePlotG1->Draw();
  // c1->SaveAs("GeneratorTree1.jpg");
  
  treePlotG2->SetTree(LAB_Gen2);
  treePlotG2->Draw();
  //c1->SaveAs("GeneratorTree2.jpg");
  
  treePlotG3->SetTree(LAB_Gen3);
  treePlotG3->Draw();
  //c1->SaveAs("GeneratorTree3.jpg");

}
