/////////////////////////////////////////////////////////////////////////
//   RestFrames: particle physics event analysis library
//   --------------------------------------------------------------------
//   Copyright (c) 2014-2016, Christopher Rogan
/////////////////////////////////////////////////////////////////////////
///
///  \file   example_H_to_WlnuWlnu.C
///
///  \author Christopher Rogan
///          (crogan@cern.ch)
///
///  \date   2015 June
//
//   This file is part of RestFrames.
//
//   RestFrames is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
// 
//   RestFrames is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
// 
//   You should have received a copy of the GNU General Public License
//   along with RestFrames. If not, see <http://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////////////////////////

//Write new code. Create 3 generator trees. Create a plot of the reconstructed intermediary mass for all 3 decay trees. Create plots of decay angle vs mass for each way. Change masses if time permits. 

#include "RestFrames/RestFrames.hh"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TColor.h"
#include "TLegend.h"

using namespace RestFrames;

void HDecay(const std::string& output_name =
			   "output_H_to_WlnuWlnu.root"){

  Int_t palette[5];
  palette[0] = (kMagenta+2);
  palette[1] = (kRed+3);
  palette[2] = (kYellow);
  palette[3] = (kAzure+6);
  palette[4] = kBlue;
  gStyle->SetPalette(5,palette);
  double w = 1500;
  double h = 900;
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c","c",w,h);
  c1->Divide(3,2,0.005,0.01,0);
    //Write histograms here
    THStack* stack1 = new THStack ("Stack", "Masses of Gravitons of Different Widths");
    THStack* stack2 = new THStack ("Stack", "Momentums of Gravitons of Different Widths");
    THStack* stack3 = new THStack ("Stack", "Graviton Masses (Width = 10%)");
    THStack* stack4 = new THStack ("Stack", "Graviton Momentums (Width = 10%)");
    THStack* stack5 = new THStack ("Stack", "Test");
    THStack* stack6 = new THStack ("Stack", "Test");
    TH1D* origm = new TH1D ("hist","hist",1000,0.,1600.);
    TH1D* newm = new TH1D ("hist","hist",1000,0.,1600.);
    TH1D* origcos = new TH1D ("hist","Cosine of Theta without Spread",1000,-1.,1.);
    TH1D* newcos = new TH1D ("hist","Cosine of Theta with Spread",1000,-1.,1.);
    TH1D* origp = new TH1D ("hist","hist",1000,0.,1600.);
    TH1D* newp = new TH1D ("hist","hist",1000,0.,1600.);
    TH2D* hist2 = new TH2D ("hist2","hist2",1000,0.,200., 1000,0.,200.);
    TH2D* origm2 = new TH2D ("hist","Ha and Hb Mass with Spread",1000,0.,1600.,1000.,0.,1600.);
    TH2D* newm2 = new TH2D ("hist","Graviton  Mass versus Width",1000,0.,0.6,1000.,0.,1600.);
    TH2D* origcos2 = new TH2D ("hist","hist",1000,-1.,1.,1000.,-1.,1.);
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
  G_Gen.AddChildFrame(Hb_Gen);
  Ha_Gen.AddChildFrame(La_Gen);
  Ha_Gen.AddChildFrame(Lb_Gen);
  Hb_Gen.AddChildFrame(Lc_Gen);
  Hb_Gen.AddChildFrame(Ld_Gen);

  if(LAB_Gen.InitializeTree())
    g_Log << LogInfo << "...Successfully initialized generator tree" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator tree" << LogEnd;								    
  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//

  if(mG.size() < 1) return;
  
  // set Higgs masses
  G_Gen.SetMass(mG[0]);
  G_Gen.SetWidth(wG[0]*mG[0]);
  // set W masses and widths
  Ha_Gen.SetMass(mW);                    Ha_Gen.SetWidth(wW);
  Hb_Gen.SetMass(mW);                    Hb_Gen.SetWidth(wW);
  // set lepton and neutrino masses
  La_Gen.SetMass(0);                    Lb_Gen.SetMass(0);
  Lc_Gen.SetMass(0);                    Ld_Gen.SetMass(0);

  /* // set lepton pT and eta cuts
  La_Gen.SetPtCut(10.);                  Lb_Gen.SetPtCut(10.);
  La_Gen.SetEtaCut(2.5);                 Lb_Gen.SetEtaCut(2.5);
  Lc_Gen.SetPtCut(10.);                  Ld_Gen.SetPtCut(10.);
  Lc_Gen.SetEtaCut(2.5);                 Ld_Gen.SetEtaCut(2.5);
  */
  if(LAB_Gen.InitializeAnalysis())
    g_Log << LogInfo << "...Successfully initialized generator analysis" << LogEnd;
  else
    g_Log << LogError << "...Failed initializing generator analysis" << LogEnd;
  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////
  g_Log << LogInfo << "Initializing reconstruction frames and trees..." << LogEnd;
  //////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  TreePlot* treePlot = new TreePlot("TreePlot","TreePlot");
 
  treePlot->SetTree(LAB_Gen);

  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//


  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////

  for(int m = 0; m < Nmass; m++){
    g_Log << LogInfo << "Generating events for H^{0} mass = " << mG[m] << LogEnd;

    G_Gen.SetMass(mG[0]);
    G_Gen.SetWidth(wG[m]*mG[m]);
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

      //double Gmass = G.M();
      double Gmass = G_Gen.GetMass();
      double cosThetaG = G_Gen.GetCosDecayAngle();

      // Graviton RestFrame FourVectors
      TLorentzVector Ha_G = Ha_Gen.GetFourVector(G_Gen);
      TLorentzVector Hb_G = Hb_Gen.GetFourVector(G_Gen);

      origm->Fill(Ha_G.M());
      newm->Fill(Gmass);
      newm_array[m]->Fill(Gmass);
      origcos->Fill(cosThetaG);
      newcos->Fill(cosThetaG);
      origp->Fill(G.P());
      newp->Fill(G.P());
      newp_array[m]->Fill(G.P());
      origm2->Fill(Ha_G.M(),Hb_G.M());
      newm2->Fill(wG[m], Gmass);
      origp2->Fill(Ha_G.P(),Hb_G.P());
      newp2->Fill(wG[m], Ha_G.P()+Hb_G.P());
      width->Fill(wG[m]); 

      // lab frame
      // cout << "G mass = " << G.M() << " " << (Ha_G+Hb_G).M() << " momentum of Ha+Hb = " << (Ha_G+Hb_G).P() << endl;
      
      // analyze event
    }

    LAB_Gen.PrintGeneratorEfficiency();
  }

  //-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//-//


  /////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  
  //Changing mass
  for(int m = 0; m < Nmass; m++){
    g_Log << LogInfo << "Generating events for H^{0} mass = " << mG[m] << LogEnd;

    G_Gen.SetMass(mG[m]);
    G_Gen.SetWidth(wG[0]*mG[m]);
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

      //double Gmass = G.M();
      double Gmass = G_Gen.GetMass();
      double cosThetaG = G_Gen.GetCosDecayAngle();

      // Graviton RestFrame FourVectors
      TLorentzVector Ha_G = Ha_Gen.GetFourVector(G_Gen);
      TLorentzVector Hb_G = Hb_Gen.GetFourVector(G_Gen);
      
      newmm_array[m]->Fill(Gmass);
      newpm_array[m]->Fill(G.P()); 

      // lab frame
      // cout << "G mass = " << G.M() << " " << (Ha_G+Hb_G).M() << " momentum of Ha+Hb = " << (Ha_G+Hb_G).P() << endl;
      
      // analyze event
    }

    LAB_Gen.PrintGeneratorEfficiency();
    }
  
  for(int i = 0; i < 5; i++){
    newm_array[i]->SetLineColor(kBlack);
    newm_array[i]->SetLineWidth(2);
    newm_array[i]->SetFillStyle(3001);
    newp_array[i]->SetLineColor(kBlack);
    newp_array[i]->SetLineWidth(2);
    newp_array[i]->SetFillStyle(3001);
    newm_array[i]->Scale(1./newm_array[i]->GetMaximum());
    newp_array[i]->Scale(1./newp_array[i]->GetMaximum());
    stack1->Add(newm_array[i]);
    stack2->Add(newp_array[i]);
    newmm_array[i]->SetLineColor(kBlack);
    newmm_array[i]->SetLineWidth(2);
    newmm_array[i]->SetFillStyle(3001);
    newpm_array[i]->SetLineColor(kBlack);
    newpm_array[i]->SetLineWidth(2);
    newpm_array[i]->SetFillStyle(3001);
    newmm_array[i]->Scale(1./newm_array[i]->GetMaximum());
    newpm_array[i]->Scale(1./newp_array[i]->GetMaximum());
    stack3->Add(newmm_array[i]);
    stack4->Add(newpm_array[i]);
    }
  c1->cd(1);
  

  //Drawing histograms with "stack"
  stack2->Draw("hist pfc");
  stack2->GetXaxis()->SetTitle("Momentum the Graviton");
  stack2->GetYaxis()->SetTitle("Counts (Normalized and Stacked)");
   
/* 
  //Drawing histograms with "same"

  newp_array[0]->GetXaxis()->SetTitle("Momentum of the Graviton");
  newp_array[0]->GetYaxis()->SetTitle("Normalized Counts");
  newp_array[0]->Draw();

  for(int i = 1; i < 5; i++)
   newp_array[i]->Draw("same hist");
*/
  c1->cd(4);

   stack1->Draw("nostack");
   stack1->GetXaxis()->SetTitle("Masses of the Gravitons");
   stack1->GetYaxis()->SetTitle("Counts (Normalized and Stacked)");
   TLegend* Legend4 = new TLegend(0.75,0.75,0.9,0.9);
   Legend4->AddEntry(newm_array[0], "Mg = 250");
   Legend4->AddEntry(newm_array[1], "Mg=300");
   Legend4->Draw();
  
  c1->cd(3);
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
  
  /* histPlot->Draw();

  TFile fout(output_name.c_str(),"RECREATE");
  fout.Close();
  histPlot->WriteOutput(output_name);
  histPlot->WriteHist(output_name);
  treePlot->WriteOutput(output_name);*/
  
  g_Log << LogInfo << "Finished" << LogEnd;
  }
