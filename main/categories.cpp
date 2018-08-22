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
#include "TF1.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraph2D.h"
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
  TFile* outFile = TFile::Open(Form("categories.root"),"RECREATE");
  outFile -> cd();
  
  std::map<std::pair<float,float>,std::map<std::string,TH1F*> > h1_sig;
  std::map<std::pair<float,float>,std::map<std::string,TH1F*> > h1_bkg;
  std::map<std::pair<float,float>,std::map<std::string,TH1F*> > h1_bkg_fit;
  
  std::vector<float> mjjcut;
  mjjcut.push_back(400.);
  mjjcut.push_back(500.);
  mjjcut.push_back(600.);
  mjjcut.push_back(700.);
  mjjcut.push_back(800.);
  mjjcut.push_back(900.);
  mjjcut.push_back(1000.);
  mjjcut.push_back(1100.);
  mjjcut.push_back(1200.);
  mjjcut.push_back(1300.);
  mjjcut.push_back(1500.);
  mjjcut.push_back(1700.);
  mjjcut.push_back(2000.);
  std::vector<float> detacut;
  detacut.push_back(2.);
  detacut.push_back(3.);
  detacut.push_back(4.);
  detacut.push_back(5.);
  detacut.push_back(6.);
  detacut.push_back(7.);
  
  std::vector<std::string> categories;
  categories.push_back("cat0");
  categories.push_back("cat1");
  
  for(auto mjjIt : mjjcut)
    for(auto detaIt : detacut)
      for(auto cat : categories)
      {            
        std::pair<float,float> dummy_key(mjjIt,detaIt);
        h1_sig[dummy_key][cat] = new TH1F(Form("h1_sig_mjj%.1f_deta%.1f_%s",mjjIt,detaIt,cat.c_str()),"",60,110.,140.);
        h1_bkg[dummy_key][cat] = new TH1F(Form("h1_bkg_mjj%.1f_deta%.1f_%s",mjjIt,detaIt,cat.c_str()),"",60,110.,140.);
        h1_bkg_fit[dummy_key][cat] = new TH1F(Form("h1_bkg_fit_mjj%.1f_deta%.1f_%s",mjjIt,detaIt,cat.c_str()),"",60,110.,140.);
      }
  
  
  
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
      if( entry%100000 == 0 ) std::cout << ">>>>>> 1st loop: reading entry "
                                        << std::fixed << std::setw(8) << entry
                                        << " / "
                                        << std::fixed << std::setw(8) << nEntries
                                        << " (" << std::fixed << std::setprecision(1) << std::setw(5) << float(entry)/float(nEntries)*100. << "%)"
                                        << "\r" << std::flush;
      (mapIt.second)->GetEntry(entry);
      
      if( !(weight<10000.) ) weight = 0.;
      else weight *= 35.9;
      
      if( jets_all_bTagM_n > 0 ) continue;
      if( met_pt > 60. ) continue;
      if( std::max(mu1_pt,mu2_pt) < 26. ) continue;
      
      for(auto mjjIt : mjjcut)
        for(auto detaIt : detacut)
        {
          std::pair<float,float> dummy_key(mjjIt,detaIt);
          
          if( jet_all_mass >= mjjIt && jet_all_Deta >= detaIt)
          {            
            if( isSig ) h1_sig[dummy_key]["cat0"] -> Fill(H_mass,weight);
            if( isBkg ) h1_bkg[dummy_key]["cat0"] -> Fill(H_mass,weight);
          }
          
          else
          {
            if( isSig ) h1_sig[dummy_key]["cat1"] -> Fill(H_mass,weight);
            if( isBkg ) h1_bkg[dummy_key]["cat1"] -> Fill(H_mass,weight);
          }
        }
    } // loop over entries
    std::cout << std::endl;
    
  } // loop over samples
  std::cout << std::endl;
  
  
  TGraph2D* g2_totSignificance_mode1 = new TGraph2D();
  TGraph2D* g2_totSignificance_mode2 = new TGraph2D();
  
  for(auto mjjIt : mjjcut)
    for(auto detaIt : detacut)
    {
      double totSignificance_mode1 = 0.;
      double totSignificance_mode2 = 0.;
      
      // loop over categories
      for(auto cat : categories )
      {
        std::pair<float,float> dummy_key(mjjIt,detaIt);
        
        std::string label(Form("mjj%.1f_deta%.1f_%s",mjjIt,detaIt,cat.c_str()));
        TF1* fitFunc = new TF1(Form("fitFunc_%s",label.c_str()),"[3]*(exp([0]*x+[1]*x*x))/((x-90)^[2]+(2.5/2.)^[2])",110.,140.);
        fitFunc -> SetParameters(5.88026e-02,-1.89398e-04,2.28365e+00,h1_bkg[dummy_key][cat]->Integral());
        h1_bkg[dummy_key][cat] -> Fit(fitFunc,"QNRLS+");
        fitFunc -> Write();
        
        for(int bin = 1; bin <= h1_bkg_fit[dummy_key][cat]->GetNbinsX(); ++bin)
          h1_bkg_fit[dummy_key][cat] -> SetBinContent(bin,fitFunc->Eval(h1_bkg_fit[dummy_key][cat]->GetBinCenter(bin)));
        
        double catSignificance_mode1 = ComputeSignificance(h1_sig[dummy_key][cat],h1_bkg_fit[dummy_key][cat],1);
        double catSignificance_mode2 = ComputeSignificance(h1_sig[dummy_key][cat],h1_bkg_fit[dummy_key][cat],2);
        totSignificance_mode1 += catSignificance_mode1*catSignificance_mode1;
        totSignificance_mode2 += catSignificance_mode2*catSignificance_mode2;
      } // loop over categories
      
      totSignificance_mode1 = sqrt(totSignificance_mode1);
      totSignificance_mode2 = sqrt(totSignificance_mode2);
      
      g2_totSignificance_mode1 -> SetPoint(g2_totSignificance_mode1->GetN(),mjjIt,detaIt,totSignificance_mode1);
      g2_totSignificance_mode2 -> SetPoint(g2_totSignificance_mode2->GetN(),mjjIt,detaIt,totSignificance_mode2);
    }
  
  
  
  g2_totSignificance_mode1 -> Write("g2_totSignificance_mode1");
  g2_totSignificance_mode2 -> Write("g2_totSignificance_mode2");
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
