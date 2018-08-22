#include "interface/AnalysisUtils.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TDirectory.h"



int main(int argc, char** argv)
{
  // define samples
  std::map<std::string,TTree*> trees;
  
  TFile* inFile;
  
  inFile = TFile::Open("plots/plots_mc_ggH.root","READ");
  trees["sig_ggH"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_qqH.root","READ");
  trees["sig_qqH"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_WH.root","READ");
  trees["sig_WH"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_ZH.root","READ");
  trees["sig_ZH"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_DY.root","READ");
  trees["bkg_DY"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_tt.root","READ");
  trees["bkg_tt"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_st.root","READ");
  trees["bkg_st"] = (TTree*)( inFile->Get("outTree") );
  
  inFile = TFile::Open("plots/plots_mc_VV.root","READ");
  trees["bkg_VV"] = (TTree*)( inFile->Get("outTree") );
  
  
  // define variables
  float weight;
  float mu1_pt;
  float mu2_pt;
  float H_mass;
  int jets_all_bTagM_n;
  float met_pt;
  float jet_all_mass;
  float jet_all_Deta;
  
  
  // define histograms
  std::map<std::string,TH1F*> h1_sig;
  std::map<std::string,TH1F*> h1_bkg;
  
  std::vector<float> mjjcut;
  mjjcut.push_back(300);
  // mjjcut.push_back(400);
  // mjjcut.push_back(500);
  // mjjcut.push_back(600);
  // mjjcut.push_back(700);
  // mjjcut.push_back(800);
  std::vector<float> detacut;
  detacut.push_back(2);
  // detacut.push_back(3);
  // detacut.push_back(4);
  // detacut.push_back(5);
  // detacut.push_back(6);
  // detacut.push_back(7);
  
  TH2F* h2_totSignificance = new TH2F("h2_totSignificance","",6,250,850,6,1.5,7.5);
  
  
  // loop over samples
  for(auto mapIt : trees)
  {
    std::cout << ">>> sample: " << mapIt.first << std::endl;
    
    std::size_t found_sig = (mapIt.first).find("sig");
    bool isSig = (found_sig != std::string::npos);
    
    std::size_t found_bkg = (mapIt.first).find("bkg");
    bool isBkg = (found_bkg != std::string::npos);
            
    (mapIt.second) -> SetBranchStatus("*",0);
    (mapIt.second) -> SetBranchStatus("weight",1); (mapIt.second) -> SetBranchAddress("weight",&weight);
    (mapIt.second) -> SetBranchStatus("mu1_pt",1); (mapIt.second) -> SetBranchAddress("mu1_pt",&mu1_pt);
    (mapIt.second) -> SetBranchStatus("mu2_pt",1); (mapIt.second) -> SetBranchAddress("mu2_pt",&mu2_pt);
    (mapIt.second) -> SetBranchStatus("H_mass",1); (mapIt.second) -> SetBranchAddress("H_mass",&H_mass);
    (mapIt.second) -> SetBranchStatus("jets_all_bTagM_n",1); (mapIt.second) -> SetBranchAddress("jets_all_bTagM_n",&jets_all_bTagM_n);
    (mapIt.second) -> SetBranchStatus("met_pt",1); (mapIt.second) -> SetBranchAddress("met_pt",&met_pt);
    (mapIt.second) -> SetBranchStatus("jet_all_mass",1); (mapIt.second) -> SetBranchAddress("jet_all_mass",&jet_all_mass);
    (mapIt.second) -> SetBranchStatus("jet_all_Deta",1); (mapIt.second) -> SetBranchAddress("jet_all_Deta",&jet_all_Deta);
    
    // loop over entries
    int nEntries = (mapIt.second)->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%100000 == 0 ) std::cout << ">>>>>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      (mapIt.second)->GetEntry(entry);
      
      if( !(weight<10000.) ) weight = 0.;
      else weight *= 35.9;
      
      if( jets_all_bTagM_n > 0 ) continue;
      if( met_pt > 60. ) continue;
      if( std::max(mu1_pt,mu2_pt) < 26. ) continue;
      
      std::string label(Form("mjj%.1f_deta%.1f",mjjIt,detaIt));
            
      for(auto mjjIt : mjjcut)
        for(auto detaIt : detacut)
        {
          if( jet_all_mass >= mjjIt && jet_all_Deta >= detaIt)
          {            
            std::string label_cat = label + "_cat0";
            
            if( h1_sig[label_cat] == NULL ) h1_sig[label_cat] = new TH1F(Form("h1_sig_%s",label_cat.c_str()),"",60,110.,140.);
            if( h1_bkg[label_cat] == NULL ) h1_bkg[label_cat] = new TH1F(Form("h1_bkg_%s",label_cat.c_str()),"",60,110.,140.);
            
            if( isSig ) h1_sig[label_cat] -> Fill(H_mass,weight);
            if( isBkg ) h1_bkg[label_cat] -> Fill(H_mass,weight);
          }
          
          else
          {
            std::string label_cat = label + "_cat1";
            
            if( h1_sig[label_cat] == NULL ) h1_sig[label_cat] = new TH1F(Form("h1_sig_%s",label_cat.c_str()),"",60,110.,140.);
            if( h1_bkg[label_cat] == NULL ) h1_bkg[label_cat] = new TH1F(Form("h1_bkg_%s",label_cat.c_str()),"",60,110.,140.);

            if( isSig ) h1_sig[label_cat] -> Fill(H_mass,weight);
            if( isBkg ) h1_bkg[label_cat] -> Fill(H_mass,weight);            
          }
        }
    } // loop over entries
  } // loop over samples
  
}
