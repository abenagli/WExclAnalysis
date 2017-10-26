#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
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



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> decayModes.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //--- open files and get trees
  std::vector<std::string> inputFiles = opts.GetOpt<std::vector<std::string> >("Input.inputFiles");
  
  TChain* chain_gen = new TChain("DumpGenParticles/gen_tree","DumpGenParticles/gen_tree");

  for(unsigned int fileIt = 0; fileIt < inputFiles.size(); ++fileIt)
    chain_gen -> Add(inputFiles.at(fileIt).c_str());
  
  TreeVars treeVars;
  InitTreeVars(NULL,chain_gen,treeVars);
  
  int nEntries_genTree = chain_gen -> GetEntries();
  std::cout << "Read " << nEntries_genTree  << " in tree gen_tree"  << std::endl;
  
  //--- get other variables
  int printGenEvent = opts.GetOpt<int>("Input.printGenEvent");
  
  
  //---------------------
  //--- define histograms
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  TFile* outFile = TFile::Open(Form("%s.root",outputFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  TH1F* h1_W_mass    = new TH1F("h1_W_mass",   "",1000,0.,200.);
  TH1F* h1_JPsi_mass = new TH1F("h1_JPsi_mass","",1000,0.,5.);
  TH1F* h1_mumu_mass = new TH1F("h1_mumu_mass","",1000,0.,5.);
  TH1F* h1_Ds_mass   = new TH1F("h1_Ds_mass",  "",1000,0.,5.);
  TH1F* h1_kepi_mass = new TH1F("h1_kepi_mass","",1000,0.,5.);
  
  
  //------------------------------------
  //--- container for decay mode counter
  std::vector<std::pair<std::string,int> > decayMode_WP;
  std::vector<std::pair<std::string,int> > decayMode_WM;
  
  
  //--------------------
  //--- loop over events
  int nEvents_WP = 0;
  int nEvents_WM = 0;
  
  for(int entry = 0; entry < chain_gen->GetEntries(); ++entry)
  {
    if( entry%100 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_genTree << "\r" << std::flush;
    chain_gen -> GetEntry(entry);
    
    
    //--------------
    //--- candidates
    
    particle gen_W;
    particle gen_JPsi;
    particle gen_Ds;
    std::vector<particle> gen_mu;
    std::vector<particle> gen_kepi;
    
    
    //---------------------
    //--- retrieve gen info
    
    gen_W.v.SetPtEtaPhiE(treeVars.reso_pt->at(0),treeVars.reso_eta->at(0),treeVars.reso_phi->at(0),treeVars.reso_energy->at(0));
    gen_W.charge = treeVars.reso_charge->at(0);
    gen_W.pdgId = treeVars.reso_pdgId->at(0);
    
    if( gen_W.charge > 0 ) ++nEvents_WP;
    else                   ++nEvents_WM;
    
    for(int ii = 0; ii < (*treeVars.resoDau1_n)[0]; ++ii)
    {
      if( (*treeVars.resoDau1_pdgId)[0][ii] == 443 )
      {
        gen_JPsi.v.SetPtEtaPhiE((*treeVars.resoDau1_pt)[0][ii],(*treeVars.resoDau1_eta)[0][ii],(*treeVars.resoDau1_phi)[0][ii],(*treeVars.resoDau1_energy)[0][ii]);
        gen_JPsi.charge = (*treeVars.resoDau1_charge)[0][ii];
        gen_JPsi.pdgId = (*treeVars.resoDau1_pdgId)[0][ii];
        
        for(int jj = 0; jj < (*treeVars.resoDau2_n)[0][ii]; ++jj)
        {
          if( abs((*treeVars.resoDau2_pdgId)[0][ii][jj]) == 13 )
          {
            particle temp;
            temp.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj],(*treeVars.resoDau2_eta)[0][ii][jj],(*treeVars.resoDau2_phi)[0][ii][jj],(*treeVars.resoDau2_energy)[0][ii][jj]);
            temp.charge = (*treeVars.resoDau2_charge)[0][ii][jj];
            temp.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj];
            
            gen_mu.push_back( temp );
          }
        }
      }
      
      if( abs((*treeVars.resoDau1_pdgId)[0][ii]) == 431 )
      {
        gen_Ds.v.SetPtEtaPhiE((*treeVars.resoDau1_pt)[0][ii],(*treeVars.resoDau1_eta)[0][ii],(*treeVars.resoDau1_phi)[0][ii],(*treeVars.resoDau1_energy)[0][ii]);
        gen_Ds.charge = (*treeVars.resoDau1_charge)[0][ii];
        gen_Ds.pdgId = (*treeVars.resoDau1_pdgId)[0][ii];
        
        for(int jj = 0; jj < (*treeVars.resoDau2_n)[0][ii]; ++jj)
        {
          particle temp;
          temp.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj],(*treeVars.resoDau2_eta)[0][ii][jj],(*treeVars.resoDau2_phi)[0][ii][jj],(*treeVars.resoDau2_energy)[0][ii][jj]);
          temp.charge = (*treeVars.resoDau2_charge)[0][ii][jj];
          temp.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj];
          
          if( (temp.pdgId == 22) && (jj < ((*treeVars.resoDau2_n)[0][ii]-1)) )
          {
            particle temp2;
            temp2.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj+1],(*treeVars.resoDau2_eta)[0][ii][jj+1],(*treeVars.resoDau2_phi)[0][ii][jj+1],(*treeVars.resoDau2_energy)[0][ii][jj+1]);
            temp2.charge = (*treeVars.resoDau2_charge)[0][ii][jj+1];
            temp2.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj+1];
            
            if( (temp2.pdgId == 22) && (fabs((temp.v+temp2.v).M()-MASS_PI0) < 0.001) )
            {
              particle pi0;
              pi0.v = temp.v + temp2.v;
              pi0.charge = 0;
              pi0.pdgId = 111;
              
              gen_kepi.push_back( pi0 );
              ++jj;
              continue;
            }
          }
          
          gen_kepi.push_back( temp );
        }
      }
    }
    
    std::sort(gen_mu.begin(),gen_mu.end(),NameSort());
    std::sort(gen_kepi.begin(),gen_kepi.end(),NameSort());
    
    particle mumu;
    for(unsigned int ii = 0; ii < gen_mu.size(); ++ii)
      mumu.v += gen_mu.at(ii).v;
    
    particle kepi;
    for(unsigned int ii = 0; ii < gen_kepi.size(); ++ii)
      kepi.v += gen_kepi.at(ii).v;
    
    h1_W_mass    -> Fill( gen_W.v.M() );
    h1_JPsi_mass -> Fill( gen_JPsi.v.M() );
    h1_mumu_mass -> Fill( mumu.v.M() );
    h1_Ds_mass   -> Fill( gen_Ds.v.M() );
    h1_kepi_mass -> Fill( kepi.v.M() );
    
    
    std::string decayMode = "";
    for(unsigned int ii = 0; ii < gen_kepi.size(); ++ii)
      decayMode += GetParticleName(gen_kepi.at(ii).pdgId) + " ";
    
    if( gen_W.charge > 0 )
    {
      std::pair<std::string,int> p(decayMode,1);
      std::vector<std::pair<std::string,int> >::iterator it = std::find_if(decayMode_WP.begin(),decayMode_WP.end(),FindPair(decayMode));
      if( it != decayMode_WP.end() ) it->second += 1;
      else decayMode_WP.push_back( p );
    }
    else
    {
      std::pair<std::string,int> p(decayMode,1);
      std::vector<std::pair<std::string,int> >::iterator it = std::find_if(decayMode_WM.begin(),decayMode_WM.end(),FindPair(decayMode));
      if( it != decayMode_WM.end() ) it->second += 1;
      else decayMode_WM.push_back( p );
    }
    
    
    //---------------
    //--- print event
    if( printGenEvent )
    {
      std::cout << "\n\n Event: " << entry << std::endl;
      {
        std::cout << "--------------------- GEN ---------------------" << std::endl;
        PrintEvent(gen_mu,gen_kepi);
      }
    }    
  }
  std::cout << std::endl;
  
  
  std::sort(decayMode_WP.begin(),decayMode_WP.end(),PairSort());
  std::sort(decayMode_WM.begin(),decayMode_WM.end(),PairSort());
  
  std::cout << std::fixed << std::endl;
  std::cout << "===================================================="  << std::endl;
  std::cout << "         Tot. events: "    << std::setw(6) << nEvents_WP + nEvents_WM << std::endl;
  
  std::cout << ">>> D_s+ decay modes: " << std::endl;
  for(unsigned int ii = 0; ii < decayMode_WP.size(); ++ii)
    std::cout << std::setw(75) << decayMode_WP.at(ii).first << ": " << std::setw(5) << decayMode_WP.at(ii).second << " (" << std::setprecision(2) << std::setw(5) << 100.*decayMode_WP.at(ii).second/nEvents_WP << "%)" << std::endl;
  
  std::cout << ">>> D_s- decay modes: " << std::endl;
  for(unsigned int ii = 0; ii < decayMode_WM.size(); ++ii)
    std::cout << std::setw(75) << decayMode_WM.at(ii).first << ": " << std::setw(5) << decayMode_WM.at(ii).second << " (" << std::setprecision(2) << std::setw(5) << 100.*decayMode_WM.at(ii).second/nEvents_WM << "%)" << std::endl;
  
  std::cout << std::endl;
  std::cout << "===================================================="  << std::endl;
  std::cout << std::endl;
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
