//Use different reco and gen tree. Look at plots of Hb mass versus Ha mass. Change Ha and Hb mass. Look at plots for Ha.  

#include "RestFrames/RestFrames.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TColor.h"
#include "TLegend.h"
#include "TVector3.h"

using namespace RestFrames;

void DecayTrees(const std::string& output_name =
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
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("New Method","c1",w,h);
  TCanvas* c2 = new TCanvas("Ha Plots", "c2", w, h);
  TCanvas* c3 = new TCanvas("Combined Ha and Hb Masses", "c3", w,0.5*h);
  c1->Divide(3,2,0.005,0.01,0);
  c2->Divide(3,2,0.005,0.01,0);
  c3->Divide(3,1,0.005,0.01,0);
    //Write histograms here
    THStack* stack1 = new THStack ("Stack", "Masses of Gravitons of Different Widths");
    THStack* stack2 = new THStack ("Stack", "Momentums of Gravitons of Different Widths");
    THStack* stack3 = new THStack ("Stack", "Graviton Masses (Width = 10%)");
    THStack* stack4 = new THStack ("Stack", "Graviton Momentums (Width = 10%)");
    THStack* stack5 = new THStack ("Stack", "Test");
    THStack* stack6 = new THStack ("Stack", "Test");
    TH1D* origm = new TH1D ("hist","Reconstructed Mass of Ha from La and Lb",100,0.,300.);
    TH1D* origm2 = new TH1D ("hist","Reconstructed Mass of Ha from Lb and Lc",100,0.,300.);
    TH1D* origm3 = new TH1D ("hist","Reconstructed Mass of Ha from Lc and La",100,0.,300.);
    TH1D* newm = new TH1D ("hist","hist",1000,0.,1600.);
    //Decay angle histograms
    TH1D* origcos = new TH1D ("hist","Decay Angle of Hb to La and Lb",100,-1.,1.);
    TH1D* origcos2 = new TH1D ("hist","Decay Angle of Hb to Lb and Lc",100,-1.,1.);
    TH1D* origcos3 = new TH1D ("hist","Decay Angle of Hb to Lc and Ld",100,-1.,1.);
    TH1D* origcosa = new TH1D ("hist","Decay Angle of Ha to Hb and Lc",100,-1.,1.);
    TH1D* origcos2a = new TH1D ("hist","Decay Angle of Ha to Hb and La",100,-1.,1.);
    TH1D* origcos3a = new TH1D ("hist","Decay Angle of Ha to Lc and Ld",100,-1.,1.);
    TH1D* newcos = new TH1D ("hist","Cosine of Theta with Spread",1000,-1.,1.);
    TH1D* origp = new TH1D ("hist","hist",1000,0.,1600.);
    TH1D* newp = new TH1D ("hist","hist",1000,0.,1600.);
    TH2D* hist2 = new TH2D ("hist2","hist2",1000,0.,200., 1000,0.,200.);
    // TH2D* origm2 = new TH2D ("hist","Ha and Hb Mass with Spread",1000,0.,1600.,1000.,0.,1600.);
    //Combined masses
    TH2D* newm1 = new TH2D ("hist","Combined Mass of Ha and Hb",32,0.,300,32.,0.,300.);
    TH2D* newm2 = new TH2D ("hist","Combined Mass of Ha and Hb",32,0.,300,32.,0.,300.);
    TH2D* newm3 = new TH2D ("hist","Combined Mass of Ha and Hb",32,0.,300,32.,0.,300.);
    //Mass versus decay angle histograms
    TH2D* newm2D = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Hb from La and Lb",32,-1.,1,32.,0.,250.);
    TH2D* newm2D2 = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Hb from Lb and Lc",32,-1.,1,32.,0.,250.);
    TH2D* newm2D3 = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Hb from Lc and Ld",32,-1.,1,32.,0.,250.);
    TH2D* newm2Da = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Ha from Hb and Lc",32,-1.,1,32.,0.,300.);
    TH2D* newm2D2a = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Ha from Hb and La",32,-1.,1,32.,0.,300.);
    TH2D* newm2D3a = new TH2D ("hist","Reconstructed Mass versus Decay Angle of Ha from Lc and Ld",32,-1.,1,32.,0.,300.);
    //TH2D* origcos2 = new TH2D ("hist","hist",1000,-1.,1.,1000.,-1.,1.);
    TH2D* newcos2 = new TH2D ("hist","hist",1000,-1.,1.,1000.,-1.,1.);
    TH2D* origp2 = new TH2D ("hist","Momentum of Ha and Hb without Spread",1000,0.,1600.,1000.,0.,1600.);
    TH2D* newp2 = new TH2D ("hist","Momentum of Ha + Hb versus Width",5,0.05,0.55,100.,0.,1600.);
    TH1D* width = new TH1D ("hist", "hist", 1000, 0., 1.);

    TH1D* newp_array[5];
    TH1D* newm_array[5];
    TH1D* newcos_array[5];
    TH1D* width_array[5];
    TH1D* origm2_array[5];
    TH1D* newpm_array[5];
    TH1D* newmm_array[5]; 

    for(int i = 0; i < 5; i++){
      newp_array[i] = new TH1D(Form("newp_array_%d",i), "Momentums of Gravitons of Different Widths",
			       64, 0., 2000);
      newm_array[i] = new TH1D(Form("newm_array_%d",i),Form("newm_array_%d",i),
			       64, 0., 1000);
      newcos_array[i] = new TH1D(Form("newcos_array_%d",i),Form("newcos_array_%d",i),
			       64, 0., 2000);
      width_array[i] = new TH1D(Form("width_array_%d",i),Form("width_array_%d",i),64, 0.,1.);

      origm2_array[i] = new TH1D(Form("origm2_array_%d",i),Form("origm2_array_%d",i),
			       64, 0., 1500);
      newpm_array[i] = new TH1D(Form("newpm_array_%d",i), "Masses of Gravitons of Different Widths",
      			       64, 0., 2000);
      newmm_array[i] = new TH1D(Form("newmm_array_%d",i),Form("newm_array_%d",i),
                               64, 0., 2000);
    }
    
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
  int Ngen = 10000;

  /////////////////////////////////////////////////////////////////////////////////////////
  g_Log << LogInfo << "Initializing generator frames and tree..." << LogEnd;
  /////////////////////////////////////////////////////////////////////////////////////////
  ppLabGenFrame     LAB_Gen("LAB_Gen","LAB");
   ResonanceGenFrame     G_Gen("G","G");
  ResonanceGenFrame Ha_Gen("Ha_Gen","H_{a}");
  ResonanceGenFrame Hb_Gen("Hb_Gen","H_{b}");
  VisibleGenFrame   La_Gen("La_Gen","#it{l}_{a}");
  VisibleGenFrame Lb_Gen("Lb_Gen","#it{l}_{b}");
  VisibleGenFrame   Lc_Gen("Lc_Gen","#it{l}_{c}");
  VisibleGenFrame Ld_Gen("Ld_Gen","#it{l}_{d}");

  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//

  LAB_Gen.SetChildFrame(G_Gen);
  G_Gen.AddChildFrame(Ha_Gen);
  G_Gen.AddChildFrame(Ld_Gen);
  Ha_Gen.AddChildFrame(Hb_Gen);
  Ha_Gen.AddChildFrame(Lc_Gen);
  Hb_Gen.AddChildFrame(La_Gen);
  Hb_Gen.AddChildFrame(Lb_Gen);

  if(LAB_Gen.InitializeTree())
    g_Log << LogInfo << "...Successfully initialized generator tree" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator tree" << LogEnd;								    
  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//

  if(mG.size() < 1) return;
  
  // set Higgs masses
  G_Gen.SetMass(mG[1]);
  G_Gen.SetWidth(wG[0]*mG[0]);
  // set W masses and widths
  Ha_Gen.SetMass(mW*3);                    Ha_Gen.SetWidth(wW);
  Hb_Gen.SetMass(mW);                    Hb_Gen.SetWidth(wW);
  // set lepton and neutrino masses
  La_Gen.SetMass(0);                    Lb_Gen.SetMass(0);
  Lc_Gen.SetMass(0);                    Ld_Gen.SetMass(0);


  if(LAB_Gen.InitializeAnalysis())
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
  DecayRecoFrame Hb_Reco("Hb", "Hb");
  VisibleRecoFrame La_Reco("La","La");
  VisibleRecoFrame Lb_Reco("Lb","Lb");
  VisibleRecoFrame Lc_Reco("Lc","Lc");
  VisibleRecoFrame Ld_Reco("Ld", "Ld");

  LAB_Reco.SetChildFrame(G_Reco);
  G_Reco.AddChildFrame(Ha_Reco);
  G_Reco.AddChildFrame(Ld_Reco);
  Ha_Reco.AddChildFrame(Lc_Reco);
  Ha_Reco.AddChildFrame(Hb_Reco);
  Hb_Reco.AddChildFrame(La_Reco);
  Hb_Reco.AddChildFrame(Lb_Reco);

  LabRecoFrame LAB_Reco2("LAB_Reco","LAB");
  DecayRecoFrame G_Reco2("G","G");
  DecayRecoFrame Ha_Reco2("Ha","Ha");
  DecayRecoFrame Hb_Reco2("Hb", "Hb");
  VisibleRecoFrame La_Reco2("La","La");
  VisibleRecoFrame Lb_Reco2("Lb","Lb");
  VisibleRecoFrame Lc_Reco2("Lc","Lc");
  VisibleRecoFrame Ld_Reco2("Ld", "Ld");

  LAB_Reco2.SetChildFrame(G_Reco2);
  G_Reco2.AddChildFrame(Ha_Reco2);
  G_Reco2.AddChildFrame(Hb_Reco2);
  Ha_Reco2.AddChildFrame(La_Reco2);
  Ha_Reco2.AddChildFrame(Lb_Reco2);
  Hb_Reco2.AddChildFrame(Lc_Reco2);
  Hb_Reco2.AddChildFrame(Ld_Reco2);

  LAB_Reco.InitializeTree();
  LAB_Reco.InitializeAnalysis();

  LAB_Reco2.InitializeTree();
  LAB_Reco2.InitializeAnalysis();

  /////////////////////////////////////////////////////////////////////////////////////////

  TreePlot* treePlot = new TreePlot("TreePlot","TreePlot");
 
  treePlot->SetTree(LAB_Gen);
  treePlot->Draw();

  treePlot->SetTree(LAB_Reco2);
  treePlot->Draw();

  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//


  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  //Keep mass and width the same
  G_Gen.SetMass(mG[0]);
  G_Gen.SetWidth(wG[0]*mG[0]);
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
      TLorentzVector Hb = Hb_Gen.GetFourVector();
      TLorentzVector La = La_Gen.GetFourVector();
      TLorentzVector Lb = Lb_Gen.GetFourVector();
      TLorentzVector Lc = Lc_Gen.GetFourVector();
      TLorentzVector Ld = Ld_Gen.GetFourVector();


      LAB_Reco.ClearEvent();

      La_Reco.SetLabFrameFourVector(La);
      Lb_Reco.SetLabFrameFourVector(Lb);
      Lc_Reco.SetLabFrameFourVector(Lc);
      Ld_Reco.SetLabFrameFourVector(Ld);

      LAB_Reco.AnalyzeEvent();
      double cos1 = Hb_Reco.GetCosDecayAngle();
      double cos1a = Ha_Reco.GetCosDecayAngle();

      
      //double Gmass = G.M();
      double Gmass = G_Reco.GetMass();
      double cosThetaG = G_Gen.GetCosDecayAngle();

      // Graviton RestFrame FourVectors
      TLorentzVector Hb_G = Hb_Gen.GetFourVector(G_Gen);

      // Lab frame
      double RecoHM = (La+Lb).M();
      double RecoHM2 = (Lb+Lc).M();
      double RecoHM3 = (Lc+Ld).M();
      double RecoHMa = (Hb+Lc).M();
      double RecoHM2a = (Hb+La).M();
      double RecoHM3a = (Lc+Ld).M();


      LAB_Reco.ClearEvent();
      
      La_Reco.SetLabFrameFourVector(Lb);
      Lb_Reco.SetLabFrameFourVector(Lc);
      Lc_Reco.SetLabFrameFourVector(La);
      Ld_Reco.SetLabFrameFourVector(Ld);

      LAB_Reco.AnalyzeEvent();


      double cos2 = Hb_Reco.GetCosDecayAngle();
      double cos2a = Ha_Reco.GetCosDecayAngle();

      LAB_Reco2.ClearEvent();

      La_Reco2.SetLabFrameFourVector(La);
      Lb_Reco2.SetLabFrameFourVector(Lb);
      Lc_Reco2.SetLabFrameFourVector(Lc);
      Ld_Reco2.SetLabFrameFourVector(Ld);

      LAB_Reco2.AnalyzeEvent();

      double cos3 = Hb_Reco2.GetCosDecayAngle();
      double cos3a = Ha_Reco2.GetCosDecayAngle();
      origcos->Fill(cos1);
      origcos2->Fill(cos2);
      origcos3->Fill(cos3);
      newm2D->Fill(cos1, RecoHM);
      newm2D2->Fill(cos2, RecoHM2);
      newm2D3->Fill(cos3, RecoHM3);
      origcosa->Fill(cos1a);
      origcos2a->Fill(cos2a);
      origcos3a->Fill(cos3a);
      newm2Da->Fill(cos1a, RecoHMa);
      newm2D2a->Fill(cos2a, RecoHM2a);
      newm2D3a->Fill(cos3a, RecoHM3a);
      newm1->Fill(RecoHMa, RecoHM);
      newm2->Fill(RecoHM2a, RecoHM2);
      newm3->Fill(RecoHM3a, RecoHM3);

      // lab frame
      // cout << "G mass = " << G.M() << " " << (Ha_G+Hb_G).M() << " momentum of Ha+Hb = " << (Ha_G+Hb_G).P() << endl;
      
      // analyze event
    }

    LAB_Gen.PrintGeneratorEfficiency();
  

		 
  c1->cd(1);
  
  newm2D->GetXaxis()->SetTitle("Decay Angle");
  newm2D->GetYaxis()->SetTitle("Reconstructed Mass of Hb");
  newm2D->Draw("colz");
		
  c1->cd(2);

  newm2D2->GetXaxis()->SetTitle("Decay Angle");
  newm2D2->GetYaxis()->SetTitle("Reconstructed Mass of Hb");
  newm2D2->Draw("colz");

  c1->cd(3);

  newm2D3->GetXaxis()->SetTitle("Decay Angle");
  newm2D3->GetYaxis()->SetTitle("Reconstructed Mass of Hb");
  newm2D3->Draw("colz");
  
  
  c1->cd(4);

  origcos->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos->GetYaxis()->SetTitle("Counts");
  origcos->GetYaxis()->SetRangeUser(0,origcos->GetMaximum()*1.1);
  origcos->SetFillColor(palette[4]);
  origcos->SetFillStyle(3001);
  origcos->SetLineColor(kBlack);
  origcos->Draw();

  c1->cd(5);

  origcos2->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos2->GetYaxis()->SetTitle("Counts");
  origcos2->SetFillColor(palette[1]);
  origcos2->SetFillStyle(3001);
  origcos2->SetLineColor(kBlack);
  origcos2->Draw();

  c1->cd(6);

  origcos3->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos3->GetYaxis()->SetTitle("Counts");
  origcos3->SetFillColor(palette[0]);
  origcos3->SetFillStyle(3001);
  origcos3->SetLineColor(kBlack);
  origcos3->Draw();
  
  //Plots for Ha
  c2->cd();
  c2->cd(1);
  
  newm2Da->GetXaxis()->SetTitle("Decay Angle");
  newm2Da->GetYaxis()->SetTitle("Reconstructed Mass of Hb");
  newm2Da->Draw("colz");
		
  c2->cd(2);

  newm2D2a->GetXaxis()->SetTitle("Decay Angle");
  newm2D2a->GetYaxis()->SetTitle("Reconstructed Mass of Hb");
  newm2D2a->Draw("colz");

  c2->cd(3);

  newm2D3a->GetXaxis()->SetTitle("Decay Angle");
  newm2D3a->GetYaxis()->SetTitle("Reconstructed Mass of Hb");
  newm2D3a->Draw("colz");
  
  c2->cd(4);

  origcosa->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcosa->GetYaxis()->SetTitle("Counts");
  origcosa->GetYaxis()->SetRangeUser(0,origcos->GetMaximum()*1.1);
  origcosa->SetFillColor(palette[4]);
  origcosa->SetFillStyle(3001);
  origcosa->SetLineColor(kBlack);
  origcosa->Draw();

  c2->cd(5);

  origcos2a->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos2a->GetYaxis()->SetTitle("Counts");
  origcos2a->SetFillColor(palette[1]);
  origcos2a->SetFillStyle(3001);
  origcos2a->SetLineColor(kBlack);
  origcos2a->Draw();

  c2->cd(6);

  origcos3a->GetXaxis()->SetTitle("Cosine of Decay Angle");
  origcos3a->GetYaxis()->SetTitle("Counts");
  origcos3a->SetFillColor(palette[0]);
  origcos3a->SetFillStyle(3001);
  origcos3a->SetLineColor(kBlack);
  origcos3a->Draw();

  c3->cd();
  c3->cd(1);
  newm1->GetXaxis()->SetTitle("Mass of Ha");
  newm1->GetYaxis()->SetTitle("Mass of Hb");
  newm1->Draw("colz");

  c3->cd(2);

  newm2->GetXaxis()->SetTitle("Mass of Ha");
  newm2->GetYaxis()->SetTitle("Mass of Hb");
  newm2->Draw("colz");

  c3->cd(3);

  newm3->GetXaxis()->SetTitle("Mass of Ha");
  newm3->GetYaxis()->SetTitle("Mass of Hb");
  newm3->Draw("colz");




  
   /* c1->cd(3);
  newm2->GetXaxis()->SetTitle("Width of Graviton's Mass (Proportion of Gmass)");
  newm2->GetYaxis()->SetTitle("Mass of the Graviton");
  newm2->Draw("colz");
  
  c1->cd(6);
  newp2->GetXaxis()->SetTitle("Width of Graviton's Mass (Proportion of Gmass)");

  newp2->GetYaxis()->SetTitle("Combined Momentum of Ha and Hb");
  newp2->GetZaxis()->SetNdivisions(20);
  gPad->SetLogz();
  newp2->Draw("colz");
  
  c1->cd(5);
  stack3->Draw("hist pfc");
  stack3->GetXaxis()->SetTitle("Mass of the Graviton");
  stack3->GetYaxis()->SetTitle("Counts (Normalized and Stacked)");
  
  c1->cd(2);
  stack4->Draw("hist pfc");
  stack4->GetXaxis()->SetTitle("Momentum of the Graviton");
  stack4->GetYaxis()->SetTitle("Counts (Normalized and Stacked)");
   */
  /* histPlot->Draw();

  TFile fout(output_name.c_str(),"RECREATE");
  fout.Close();
  histPlot->WriteOutput(output_name);
  histPlot->WriteHist(output_name);
  treePlot->WriteOutput(output_name);*/
  
  g_Log << LogInfo << "Finished" << LogEnd;
  }
