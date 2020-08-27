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

void RootFit(const std::string& output_name =
			   "output_H_to_WlnuWlnu.root"){

  Int_t palette[5];
  palette[0] = (kMagenta+2);
  palette[1] = (kRed);
  palette[2] = (kYellow);
  palette[3] = (kAzure+6);
  palette[4] = kBlue;
  gStyle->SetPalette(5,palette);
  double w = 1500;
  double h = 900;
  TCanvas* c1 = new TCanvas("New Method","New Method",w,h);
  TCanvas* c2 = new TCanvas("RooFit", "RooFit", w, h);
  c1->Divide(3,2,0.005,0.01,0);
  c2->Divide(2,1,0.005,0.01,0);
  //gStyle->SetOptStat(0);
  
  //Declare variables
  double cos1, cos2, cos3, MHa, MHa2, MHa3;
  double cosSum = cos1+cos2+cos3;
  double massSum = MHa+MHa2+MHa3;

  //Filling a TTree and File
  TFile f1("RooFitData.root", "recreate");
  TTree t1("t1","t1");
    t1.Branch("cos1", &cos1);
    t1.Branch("cos2", &cos2);
    t1.Branch("cos3", &cos3);
    t1.Branch("MHa", &MHa);
    t1.Branch("MHa2", &MHa2);
    t1.Branch("MHa3", &MHa3);
    t1.Branch("CosSum", &cosSum);
    t1.Branch("MassSum", &massSum);

  //Declaring histograms
   TH2D* newm2D = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Ha from La and Lb",32,-1.,1,32.,0.,250.);
   TH2D* newm2D2 = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Ha from Lb and Lc",32,-1.,1,32.,0.,250.);
   TH2D* newm2D3 = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Ha from Lc and Ld",32,-1.,1,32.,0.,250.);
   TH1D* origcos = new TH1D ("hist","Decay Angle of Ha to La and Lb",100,-1.,1.);
   TH1D* origcos2 = new TH1D ("hist","Decay Angle of Ha to Lb and Lc",100,-1.,1.);
   TH1D* origcos3 = new TH1D ("hist","Decay Angle of Ha to Lc and Ld",100,-1.,1.);
   TH2D* combinedm = new TH2D("hist", "Combined Mass versus Decay Angle", 100, -1., 1., 100, 0., 250.);
   TH1D* combinedcos = new TH1D ("hist", "Combined Decay Angle",100,-1.,1.);

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
  int Ngen = 100000000;

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
  G_Gen.SetWidth(wG[2]*mG[0]);
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
  G_Gen2.SetWidth(wG[2]*mG[0]);
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
  G_Gen3.SetWidth(wG[2]*mG[0]);
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
  DecayRecoFrame Ha_Reco("Ha","Ha");
  VisibleRecoFrame La_Reco("La","La");
  VisibleRecoFrame Lb_Reco("Lb","Lb");
  VisibleRecoFrame Lc_Reco("Lc","Lc");

  LAB_Reco.SetChildFrame(G_Reco);
  G_Reco.AddChildFrame(Ha_Reco);
  G_Reco.AddChildFrame(Lc_Reco);
  Ha_Reco.AddChildFrame(La_Reco);
  Ha_Reco.AddChildFrame(Lb_Reco);

  LAB_Reco.InitializeTree();
  LAB_Reco.InitializeAnalysis();

  /////////////////////////////////////////////////////////////////////////////////////////

  TreePlot* treePlot = new TreePlot("TreePlot","TreePlot");
 
  treePlot->SetTree(LAB_Gen);
  treePlot->Draw();


  /////////////////////////////////////////////////////////////////////////////////////////
   G_Gen.SetMass(mG[0]);
  G_Gen.SetWidth(wG[2]*mG[0]);
  LAB_Gen.InitializeAnalysis();

   for(int igen = 0; igen < Ngen; igen++){
      if(igen%((std::max(Ngen,10))/10) == 0)
	g_Log << LogInfo << "Generating event " << igen << " of " << Ngen << LogEnd;

      // generate event
      LAB_Gen.ClearEvent();                            // clear the gen tree
      
      LAB_Gen.AnalyzeEvent();                          // generate a new event

      // LAB Frame FourVectors
      TLorentzVector G = G_Gen.GetFourVector();
      TLorentzVector Ha = Ha_Gen.GetFourVector();
      TLorentzVector La = La_Gen.GetFourVector();
      TLorentzVector Lb = Lb_Gen.GetFourVector();
      TLorentzVector Lc = Lc_Gen.GetFourVector();

      LAB_Gen2.ClearEvent();                            // clear the gen tree
      
      LAB_Gen2.AnalyzeEvent();                          // generate a new event

      // LAB Frame FourVectors
      TLorentzVector G2 = G_Gen2.GetFourVector();
      TLorentzVector Ha2 = Ha_Gen2.GetFourVector();
      TLorentzVector La2 = La_Gen2.GetFourVector();
      TLorentzVector Lb2 = Lb_Gen2.GetFourVector();
      TLorentzVector Lc2 = Lc_Gen2.GetFourVector();

      LAB_Gen3.ClearEvent();                            // clear the gen tree
      
      LAB_Gen3.AnalyzeEvent();                          // generate a new event

      // LAB Frame FourVectors
      TLorentzVector G3 = G_Gen3.GetFourVector();
      TLorentzVector Ha3 = Ha_Gen3.GetFourVector();
      TLorentzVector La3 = La_Gen3.GetFourVector();
      TLorentzVector Lb3 = Lb_Gen3.GetFourVector();
      TLorentzVector Lc3 = Lc_Gen3.GetFourVector();


      LAB_Reco.ClearEvent();

      La_Reco.SetLabFrameFourVector(La);
      Lb_Reco.SetLabFrameFourVector(Lb);
      Lc_Reco.SetLabFrameFourVector(Lc);

      LAB_Reco.AnalyzeEvent();

      cos1 = Ha_Reco.GetCosDecayAngle();
      MHa = Ha_Reco.GetMass();


      LAB_Reco.ClearEvent();

      La_Reco.SetLabFrameFourVector(La2);
      Lb_Reco.SetLabFrameFourVector(Lb2);
      Lc_Reco.SetLabFrameFourVector(Lc2);

      LAB_Reco.AnalyzeEvent();

      cos2 = Ha_Reco.GetCosDecayAngle();
      MHa2= Ha_Reco.GetMass();


      LAB_Reco.ClearEvent();

      La_Reco.SetLabFrameFourVector(La3);
      Lb_Reco.SetLabFrameFourVector(Lb3);
      Lc_Reco.SetLabFrameFourVector(Lc3);

      LAB_Reco.AnalyzeEvent();

      cos3 = Ha_Reco.GetCosDecayAngle();
      MHa3 = Ha_Reco.GetMass();

      //Filling the original histograms
      origcos->Fill(cos1);
      origcos2->Fill(cos2);
      origcos3->Fill(cos3);
      newm2D->Fill(cos1, MHa);
      newm2D2->Fill(cos2, MHa2);
      newm2D3->Fill(cos3, MHa3);
      //Filling the combined histograms and data files
      combinedcos->Fill(cos1);
      combinedcos->Fill(cos2);
      combinedcos->Fill(cos3);
      combinedm->Fill(cos1, MHa);
      combinedm->Fill(cos2, MHa2);
      combinedm->Fill(cos3, MHa3);

      //Filling the TTree
      t1.Fill();


   }
   t1.Write();
  /////////////////////////////////////////////////////////////////////////////////////////
   //Histograms for known distributions
   c1->cd(1);

   newm2D->Draw("colz");

   c1->cd(2);

   newm2D2->Draw("colz");

   c1->cd(3);

   newm2D3->Draw("colz");

   c1->cd(4);

   origcos->Draw();

   c1->cd(5);

   origcos2->Draw();

   c1->cd(6);

   origcos3->Draw();

   //Histograms for combined distribution
   c2->cd(1);

   combinedm->Draw("colz");

   c2->cd(2);

   combinedcos->Draw();
}
