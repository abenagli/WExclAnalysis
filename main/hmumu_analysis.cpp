#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/PUReweighting.h"
#include "interface/OutTree.h"
#include "interface/StdHistoSet.h"
#include "src/RoccoR.cc"
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
  if( argc < 2 )
  {
    std::cerr << ">>>>> hmumu_analysis.cpp::usage:   " << argv[0] << " configFileName   [default=0/debug=1]" << std::endl;
    return -1;
  }
  TDirectory* baseDir = gDirectory;
  
  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  //--- open files and get trees
  int isData = opts.GetOpt<int>("Input.isData");
  int isSignal = opts.GetOpt<int>("Input.isSignal");
  std::string inputFileList = opts.GetOpt<std::string>("Input.inputFileList");
  std::string inputFileListPU = opts.GetOpt<std::string>("Input.inputFileListPU");
  
  TChain* chain_reco = new TChain("DumpReco/reco_tree","DumpReco/reco_tree");
  TChain* chain_gen = NULL; if( !isData ) chain_gen = new TChain("DumpGenParticles/gen_tree","DumpGenParticles/gen_tree");
  
  std::ifstream list(inputFileList.c_str(),std::ios::in);
  std::string fileName;
  while(1)
  {
    getline(list,fileName,'\n');
    if( !list.good() ) break;
    if( fileName.at(0) == '#' ) continue;
    std::cout << "Reading file " << fileName << std::endl;
    
    chain_reco -> Add(fileName.c_str());
    if( !isData ) chain_gen -> Add(fileName.c_str());
  }
  
  TreeVars treeVars;
  InitTreeVars(chain_reco,chain_gen,treeVars);
  
  int nEntries_recoTree = chain_reco -> GetEntries();
  int nEntries_genTree = 0; if( !isData ) nEntries_genTree = chain_gen -> GetEntries();
  std::cout << "Read " << nEntries_recoTree << " in tree reco_tree" << std::endl;
  std::cout << "Read " << nEntries_genTree  << " in tree gen_tree"  << std::endl;
  
  
  //--- get other variables
  int printGenEvent = opts.GetOpt<int>("Input.printGenEvent");
  int printRecoGenMatchEvent = opts.GetOpt<int>("Input.printRecoGenMatchEvent");
  int printRecoEvent = opts.GetOpt<int>("Input.printRecoGenMatchEvent");
  
  float gen_mu_ptMin = opts.GetOpt<float>("Cuts.gen_mu_ptMin");
  float gen_mu_etaMax = opts.GetOpt<float>("Cuts.gen_mu_etaMax");
  int genCut = opts.GetOpt<int>("Cuts.genCut");
  
  int genMatch = opts.GetOpt<int>("Cuts.genMatch");  
  int genMatchCut = opts.GetOpt<int>("Cuts.genMatchCut");
  
  int reco_HLTCut = opts.GetOpt<int>("Cuts.reco_HLTCut");
  std::vector<std::string> reco_HLT_paths = opts.GetOpt<std::vector<std::string> >("Cuts.reco_HLT_paths");
  int reco_HLT_OR = opts.GetOpt<int>("Cuts.reco_HLT_OR");
  
  float reco_ptMin1_mu = opts.GetOpt<float>("Cuts.reco_ptMin1_mu");
  float reco_ptMin2_mu = opts.GetOpt<float>("Cuts.reco_ptMin2_mu");
  
  
  //--- define histograms
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  TFile* outFile = TFile::Open(Form("%s.root",outputFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  OutTree outTree = OutTree("outTree","Final tree for H > mu mu studies");
  
  StdHistoSet* genHistos = isData? NULL : new StdHistoSet("gen",outFile);;
  StdHistoSet* genCutHistos = isData? NULL : new StdHistoSet("genCut",outFile);;
  StdHistoSet* recoGenMatchHistos = isData? NULL : new StdHistoSet("recoGenMatch",outFile);;
  StdHistoSet* recoHistos = new StdHistoSet("reco",outFile);
  
  TH1F* h1_nEvents = new TH1F("h1_nEvents",";;entries",11,0.,11.);
  
  TH2F* h2_mu_pt2_vs_pt1_total = new TH2F("h2_mu_pt2_vs_pt1_total","",100,-0.5,99.5,100,-0.5,99.5);
  TH2F* h2_mu_pt2_vs_pt1_pass  = new TH2F("h2_mu_pt2_vs_pt1_pass", "",100,-0.5,99.5,100,-0.5,99.5);
  TH2F* h2_mu_cuts_pt2_vs_pt1_total = new TH2F("h2_mu_cuts_pt2_vs_pt1_total","",100,-0.5,99.5,100,-0.5,99.5);
  TH2F* h2_mu_cuts_pt2_vs_pt1_pass  = new TH2F("h2_mu_cuts_pt2_vs_pt1_pass", "",100,-0.5,99.5,100,-0.5,99.5);
  
  TEfficiency* p1_eff_mu_pt = new TEfficiency("p1_eff_mu_pt",";muon p_{T} (GeV);#epsilon",200,0.,500.);
  TEfficiency* p1_eff_mu_eta = new TEfficiency("p1_eff_mu_eta",";muon #eta;#epsilon",200,0.,10.);
  
  TH1F* h1_reco_mu_gen_mu_DR = new TH1F("h1_reco_mu_gen_mu_DR","",10000,0.,10.);
  TH1F* h1_reco_mu_gen_mu_ptRatio = new TH1F("h1_reco_mu_gen_mu_ptRatio","",10000,0.,10.);
  
  std::vector<std::pair<std::string,int> > vec_triggerPass;
  
  int nTriggerEvents_noCuts = 0;
  int nTriggerEvents_genCuts = 0;
  int nTriggerEvents_recoMuSelection = 0;
  std::vector<std::pair<std::string,int> > vec_triggerPass_noCuts;
  std::vector<std::pair<std::string,int> > vec_triggerPass_genCuts;
  std::vector<std::pair<std::string,int> > vec_triggerPass_recoMuSelection;
  TH1F* h1_triggerEff_noCuts = new TH1F("h1_triggerEff_noCuts","",15,0.,15.);
  TH1F* h1_triggerEff_genCuts = new TH1F("h1_triggerEff_genCuts","",15,0.,15.);
  TH1F* h1_triggerEff_recoMuSelection = new TH1F("h1_triggerEff_recoMuSelection","",15,0.,15.);
  
  
  
  //----------------------
  //--- pileup reweighting
  PUReweighting* puReweighting = NULL;
  std::ifstream listPU(inputFileListPU.c_str(),std::ios::in);
  if( !isData )
  {
    std::string pileupFileName_data = opts.GetOpt<std::string>("Weights.pileupFileName_data");
    TFile* pileupFile_data = TFile::Open(pileupFileName_data.c_str(),"READ");
    TH1F* pileup_data = (TH1F*)( pileupFile_data->Get("pileup") );
    
    TH1F* pileup_mc = NULL;
    while(1)
    {
      getline(listPU,fileName,'\n');
      if( !listPU.good() ) break;
      
      TFile* tempInFile = TFile::Open(("root://eoscms/"+fileName).c_str(),"READ");
      TH1F* tempHisto = (TH1F*)( tempInFile->Get("DumpPU/pileup") );
      
      if( pileup_mc == NULL )
      {
        baseDir -> cd();
        pileup_mc = new TH1F("pileup","",tempHisto->GetNbinsX(),tempHisto->GetBinLowEdge(1),tempHisto->GetBinLowEdge(tempHisto->GetNbinsX()+tempHisto->GetBinWidth(tempHisto->GetNbinsX())));
      }
      pileup_mc -> Add(tempHisto);
      
      tempInFile -> Close();
    }
    puReweighting = new PUReweighting(pileup_data,pileup_mc);
    TH1F* pileup_weights = puReweighting -> GetWeightsHistogram();
    
    outFile -> cd();
    pileup_data -> Write();
    pileup_mc -> Write();
    pileup_weights -> Write();
  }
  
  
  
  //----------------
  //--- event weight
  float mcWeight = 1.;
  if( !isData )
  {
    float lumi = opts.GetOpt<float>("Weights.lumi");
    float xsec = opts.GetOpt<float>("Weights.xsec");
    float nevt = opts.GetOpt<int>("Weights.nevt");
    
    mcWeight = 1. * lumi * xsec / nevt;
    std::cout << "mcWeight:  "<< mcWeight << std::endl;
  }
  
  
  
  
  //------------------------------
  //--- muon Rochester corrections
  RoccoR rc("data/rcdata.2016.v3");
  
  
  
  //--------------------
  //--- loop over events
  int nEvents_tot = 0;
  int nEvents_genCut = 0;
  int nEvents_recoGenMatch = 0;
  int nEvents_selected = 0;
  int nEvents_selected_recoGenMatch = 0;
  
  // for(int entry = 0; entry < 1000; ++entry)
  for(int entry = 0; entry < chain_reco->GetEntries(); ++entry)
  {
    if( !debugMode && entry%100 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_recoTree << "\r" << std::flush;
    if(  debugMode && entry%1 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_recoTree << std::endl;
    chain_reco -> GetEntry(entry);
    if( !isData ) chain_gen -> GetEntry(entry);
    ++nEvents_tot;
    
    
    float weight = 1.;
    float puWeight = 1.;
    if( !isData )
    {
      puWeight = puReweighting -> GetPUWeight(treeVars.trueNumInteractions);
      weight = puWeight * mcWeight;
    }
    
    
    ++nTriggerEvents_noCuts;
    for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
    {
      std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
      std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_noCuts.begin(),vec_triggerPass_noCuts.end(),FindPair(treeVars.trgs_name->at(ii)));
      if( it != vec_triggerPass_noCuts.end() ) it->second += treeVars.trgs_pass->at(ii);
      else vec_triggerPass_noCuts.push_back( p );
    }
    
    
    //--------------
    //--- candidates
    particle gen_H;
    std::vector<particle> gen_mu;
    
    std::vector<particle> recoGenMatch_mu;
    
    particle reco_H;
    std::vector<particle> reco_mu;
    
    
    //-------------------------
    //--- fill particle vectors
    if( debugMode ) std::cout << ">>>>>> start filling particle vectors" << std::endl;
    std::vector<particle> temp_mu;
    for(unsigned int ii = 0; ii < treeVars.muons_pt->size(); ++ii)
    {
      particle temp;
      
      double SF = 1.;
      if( isData )
        SF = rc.kScaleDT(treeVars.muons_charge->at(ii),treeVars.muons_pt->at(ii),treeVars.muons_eta->at(ii),treeVars.muons_phi->at(ii));
      if( !isData && treeVars.muons_trackerLayersWithMeasurement->at(ii) >  0 )
        SF = rc.kScaleAndSmearMC(treeVars.muons_charge->at(ii),treeVars.muons_pt->at(ii),treeVars.muons_eta->at(ii),treeVars.muons_phi->at(ii),treeVars.muons_trackerLayersWithMeasurement->at(ii),gRandom->Rndm(),gRandom->Rndm());
      
      temp.v.SetPtEtaPhiE(SF*treeVars.muons_pt->at(ii),treeVars.muons_eta->at(ii),treeVars.muons_phi->at(ii),SF*treeVars.muons_energy->at(ii));
      temp.charge = treeVars.muons_charge->at(ii);
      temp.pdgId = -13 * temp.charge;
      temp.it = ii;
      temp_mu.push_back(temp);
    }
    std::sort(temp_mu.begin(),temp_mu.end(),PtSort());
    if( debugMode ) std::cout << ">>>>>> end filling particle vectors" << std::endl;
    
    
    //---------------------
    //--- retrieve gen info
    if( debugMode ) std::cout << ">>>>>> start retrieve gen info" << std::endl;
    if( isSignal )
    {
      gen_H.v.SetPtEtaPhiE(treeVars.reso_pt->at(0),treeVars.reso_eta->at(0),treeVars.reso_phi->at(0),treeVars.reso_energy->at(0));
      gen_H.charge = treeVars.reso_charge->at(0);
      gen_H.pdgId = treeVars.reso_pdgId->at(0);
      
      for(int ii = 0; ii < (*treeVars.resoDau1_n)[0]; ++ii)
      {
        if( fabs((*treeVars.resoDau1_pdgId)[0][ii]) == 13 )
        {
          particle temp;
          temp.v.SetPtEtaPhiE((*treeVars.resoDau1_pt)[0][ii],(*treeVars.resoDau1_eta)[0][ii],(*treeVars.resoDau1_phi)[0][ii],(*treeVars.resoDau1_energy)[0][ii]);
          temp.charge = (*treeVars.resoDau1_charge)[0][ii];
          temp.pdgId = (*treeVars.resoDau1_pdgId)[0][ii];
          
          gen_mu.push_back( temp );
          // for(int jj = 0; jj < (*treeVars.resoDau2_n)[0][ii]; ++jj)
          // {
          //   if( abs((*treeVars.resoDau2_pdgId)[0][ii][jj]) == 13 )
          //   {
          //     particle temp;
          //     temp.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj],(*treeVars.resoDau2_eta)[0][ii][jj],(*treeVars.resoDau2_phi)[0][ii][jj],(*treeVars.resoDau2_energy)[0][ii][jj]);
          //     temp.charge = (*treeVars.resoDau2_charge)[0][ii][jj];
          //     temp.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj];
              
          //     gen_mu.push_back( temp );
          //     gen_mu_etaSort.push_back( temp );
          //   }
          // }
        }
      }
      if( debugMode ) std::cout << ">>>>>> end retrieve gen info" << std::endl;
      
      std::sort(gen_mu.begin(),gen_mu.end(),PtSort());
      genHistos -> FillHistos(weight,gen_mu);
      
      
      bool accept_genCut = false;
      if( genCut )
      {
        //-----------------------------------
        //--- generator-level acceptance cuts
        if( gen_mu.at(0).v.Pt() > gen_mu_ptMin &&
            gen_mu.at(1).v.Pt() > gen_mu_ptMin &&
            fabs(gen_mu.at(0).v.Eta()) < gen_mu_etaMax &&
            fabs(gen_mu.at(1).v.Eta()) < gen_mu_etaMax )
        {
          ++nEvents_genCut;
          accept_genCut = true;
          
          
          ++nTriggerEvents_genCuts;
          for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
          {
            std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
            std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_genCuts.begin(),vec_triggerPass_genCuts.end(),FindPair(treeVars.trgs_name->at(ii)));
            if( it != vec_triggerPass_genCuts.end() ) it->second += treeVars.trgs_pass->at(ii);
            else vec_triggerPass_genCuts.push_back( p );
          }
        }
      }
      
      
      if( debugMode ) std::cout << ">>>>>> start gen match" << std::endl;
      if( genMatch )
      {
        bool recoGenMatch_mu_matching = false;
        
        //-----------------------
        //--- match reco with gen      
        std::vector<int> muIt_genMatch;
        muIt_genMatch.push_back( GetBestMatch(gen_mu.at(0),temp_mu,&muIt_genMatch) );
        muIt_genMatch.push_back( GetBestMatch(gen_mu.at(1),temp_mu,&muIt_genMatch) );
        if( muIt_genMatch.at(0) >= 0 ) recoGenMatch_mu.push_back( temp_mu.at(muIt_genMatch.at(0)) );
        if( muIt_genMatch.at(1) >= 0 ) recoGenMatch_mu.push_back( temp_mu.at(muIt_genMatch.at(1)) );
        
        if( recoGenMatch_mu.size() == 2 )
        {
          recoGenMatch_mu_matching = IsMatching(recoGenMatch_mu,gen_mu,0.01,0.1,false);
          if( (genCut && accept_genCut) || (!genCut) ) 
            if( recoGenMatch_mu_matching )
            {
              ++nEvents_recoGenMatch;
              
              for(int binx = 1; binx <= h2_mu_pt2_vs_pt1_pass->GetNbinsX(); ++binx)
                for(int biny = 1; biny <= h2_mu_pt2_vs_pt1_pass->GetNbinsY(); ++biny)
                {
                  float ptMin1 = h2_mu_pt2_vs_pt1_pass -> GetXaxis() -> GetBinCenter(binx);
                  float ptMin2 = h2_mu_pt2_vs_pt1_pass -> GetYaxis() -> GetBinCenter(biny);
                  
                  h2_mu_pt2_vs_pt1_total -> Fill(ptMin1,ptMin2);
                  
                  if( recoGenMatch_mu.at(0).v.Pt() > ptMin1 &&
                      recoGenMatch_mu.at(1).v.Pt() > ptMin2 )
                    h2_mu_pt2_vs_pt1_pass -> Fill(ptMin1,ptMin2);
                }
            }
        }
        
        
        // fill efficiency plots
        if( recoGenMatch_mu.size() > 0 )
        {
          if( IsMatching(gen_mu.at(0),recoGenMatch_mu.at(0),0.01,0.1) )
          {
            if( fabs(gen_mu.at(0).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(true,gen_mu.at(0).v.Pt());
            if( gen_mu.at(0).v.Pt() > 5. )         p1_eff_mu_eta -> Fill(true,fabs(gen_mu.at(0).v.Eta()));
          }
          else
          {
            if( fabs(gen_mu.at(0).v.Eta()) < 2.1 ) p1_eff_mu_pt ->  Fill(false,gen_mu.at(0).v.Pt());
            if( gen_mu.at(0).v.Pt() > 5.)          p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(0).v.Eta()));
          }
        }
        else
        {
          if( fabs(gen_mu.at(0).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(0).v.Pt());
          if( gen_mu.at(0).v.Pt() > 5. )         p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(0).v.Eta()));
          if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(1).v.Pt());
          if( gen_mu.at(1).v.Pt() > 5. )         p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(1).v.Eta()));
        }
        
        if( recoGenMatch_mu.size() > 1 )
        {
          if( IsMatching(gen_mu.at(1),recoGenMatch_mu.at(1),0.01,0.1) )
          {
            if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(true,gen_mu.at(1).v.Pt());
            if( gen_mu.at(1).v.Pt() > 5. )         p1_eff_mu_eta -> Fill(true,fabs(gen_mu.at(1).v.Eta()));
          }
          else
          {
            if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(1).v.Pt());
            if( gen_mu.at(1).v.Pt() > 5. )         p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(1).v.Eta()));
          }
        }
        else
        {
          if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(1).v.Pt());
          if( gen_mu.at(1).v.Pt() > 5. )         p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(1).v.Eta()));
        }
        
        
        // cut
        if( genMatchCut && !recoGenMatch_mu_matching )
          continue;
      }
    }
    if( debugMode ) std::cout << ">>>>>> end gen match" << std::endl;
    
    
    
    //---------------------
    //--- trigger selection
    if( debugMode ) std::cout << ">>>>>> start trigger selection" << std::endl;
    bool trgPassAll = true;
    bool trgPassOne = false;
    for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
    {
      // std::cout << "HLT name: " << treeVars.trgs_name->at(ii) << "   pass: " << treeVars.trgs_pass->at(ii) << "    prescale: " << treeVars.trgs_prescale->at(ii) << std::endl;
      for(unsigned int jj = 0; jj < reco_HLT_paths.size(); ++jj)
      {
        std::size_t found = treeVars.trgs_name->at(ii).find( reco_HLT_paths.at(jj) );
        if (found != std::string::npos )
        {
          // std::cout << "HLT name: " << treeVars.trgs_name->at(ii) << "   pass: " << treeVars.trgs_pass->at(ii) << "    prescale: " << treeVars.trgs_prescale->at(ii) << std::endl;
          if( treeVars.trgs_pass->at(ii) == 1 )
          {
            trgPassOne = true;
            std::pair<std::string,int> p(treeVars.trgs_name->at(ii),1);
            std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass.begin(),vec_triggerPass.end(),FindPair(treeVars.trgs_name->at(ii)));
            if( it != vec_triggerPass.end() ) it->second += 1;
            else vec_triggerPass.push_back( p );
          }
          else
          {
            trgPassAll = false;
          }
        }
      }
    }
    if( reco_HLTCut && (reco_HLT_OR  && !trgPassOne) ) continue;
    if( reco_HLTCut && (!reco_HLT_OR && !trgPassAll) ) continue;
    if( debugMode ) std::cout << ">>>>>> end trigger selection" << std::endl;
    
    
    
    //---------------------
    //--- select reco muons
    if( debugMode ) std::cout << ">>>>>> start reco mu selection" << std::endl;
    
    if( temp_mu.size() < 2 ) continue;
    
    float mu1_iso;
    for(unsigned int ii = 0; ii < temp_mu.size(); ++ii)
    {
      if( temp_mu.at(ii).v.Pt() < reco_ptMin1_mu ) continue;
      if( treeVars.muons_isLoose->at(temp_mu.at(ii).it) != 1 ) continue;
      if( fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)) > 0.05 ) continue;
      if( fabs(treeVars.muons_dz->at(temp_mu.at(ii).it))  > 0.10 ) continue;
      
      mu1_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
      std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
      if( mu1_iso/temp_mu.at(ii).v.Pt() > 0.25 ) continue;
      
      reco_mu.push_back( temp_mu.at(ii) );
      break;
    }
    if( reco_mu.size() < 1 ) continue;
    
    float mu2_iso;
    for(unsigned int ii = 0; ii < temp_mu.size(); ++ii)
    {
      if( reco_mu.at(0).charge*temp_mu.at(ii).charge == 1 ) continue;
      if( temp_mu.at(ii).v.Pt() < reco_ptMin2_mu ) continue;
      if( treeVars.muons_isLoose->at(temp_mu.at(ii).it) != 1 ) continue;
      if( fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)) > 0.05 ) continue;
      if( fabs(treeVars.muons_dz->at(temp_mu.at(ii).it))  > 0.10 ) continue;
      
      mu2_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
        std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
      if( mu2_iso/temp_mu.at(ii).v.Pt() > 0.25 ) continue;
      
      reco_mu.push_back( temp_mu.at(ii) );
      break;
    }
    if( reco_mu.size() < 2 ) continue;
    
    reco_H.v = reco_mu.at(0).v + reco_mu.at(1).v;
    reco_H.charge = reco_mu.at(0).charge + reco_mu.at(1).charge; 
    reco_H.pdgId = 25;
    if( debugMode ) std::cout << ">>>>>> end reco mu selection" << std::endl;
    
    
    
    if( debugMode ) std::cout << ">>>>>> start reco mu genMatch" << std::endl;
    bool reco_mu_genMatch = false;
    if( isSignal )
      reco_mu_genMatch = IsMatching(reco_mu,gen_mu,0.01,0.1,false);
    
    ++nEvents_selected;
    if( reco_mu_genMatch ) ++nEvents_selected_recoGenMatch;
    
    
    
    if( recoGenMatch_mu.size() == 2 )
    {
      bool recoGenMatch_mu_matching = IsMatching(recoGenMatch_mu,gen_mu,0.01,0.1,false);
      if( recoGenMatch_mu_matching )
      {
        for(int binx = 1; binx <= h2_mu_cuts_pt2_vs_pt1_pass->GetNbinsX(); ++binx)
          for(int biny = 1; biny <= h2_mu_cuts_pt2_vs_pt1_pass->GetNbinsY(); ++biny)
          {
            float ptMin1 = h2_mu_cuts_pt2_vs_pt1_pass -> GetXaxis() -> GetBinCenter(binx);
            float ptMin2 = h2_mu_cuts_pt2_vs_pt1_pass -> GetYaxis() -> GetBinCenter(biny);
            
            h2_mu_cuts_pt2_vs_pt1_total -> Fill(ptMin1,ptMin2);
            
            if( recoGenMatch_mu.at(0).v.Pt() > ptMin1 &&
                recoGenMatch_mu.at(1).v.Pt() > ptMin2 )
              h2_mu_cuts_pt2_vs_pt1_pass -> Fill(ptMin1,ptMin2);
          }
      }
    }
    if( debugMode ) std::cout << ">>>>>> end reco mu genMatch" << std::endl;
    
    
    
    //-----
    // jets
    if( debugMode ) std::cout << ">>>>>> start jets" << std::endl;
    
    int jets_all_n = 0;
    int jets_all_bTagL_n = 0;
    int jets_all_bTagM_n = 0;
    int jets_all_bTagT_n = 0;
    int jets_cen_n = 0;
    int jets_cen_bTagL_n = 0;
    int jets_cen_bTagM_n = 0;
    int jets_cen_bTagT_n = 0;
    int jets_fwd_n = 0;
    int jets_fwd_bTagL_n = 0;
    int jets_fwd_bTagM_n = 0;
    int jets_fwd_bTagT_n = 0;
    
    int jetId1_all = -1;
    int jetId2_all = -1;
    int jetId1_cen = -1;
    int jetId2_cen = -1;
    for(unsigned int ii = 0; ii < treeVars.jets_pt->size(); ++ii)
    {
      if( treeVars.jets_pt->at(ii) < 30. ) continue;
      if( fabs(treeVars.jets_eta->at(ii)) > 4.7 ) continue;
      
      bool skipJet = false;
      
      for(unsigned int jj = 0; jj < reco_mu.size(); ++jj)
        if( DeltaR(treeVars.jets_eta->at(ii),treeVars.jets_phi->at(ii),reco_mu.at(jj).v.Eta(),reco_mu.at(jj).v.Phi()) < 0.4 ) skipJet = true;
      
      if( fabs(treeVars.jets_eta->at(ii)) <= 2.7 )
      {
        if( treeVars.jets_NHF->at(ii)  >= 0.90 ) skipJet = true;
        if( treeVars.jets_NEMF->at(ii) >= 0.90 ) skipJet = true;
        if( treeVars.jets_MUF->at(ii)  >= 0.80 ) skipJet = true;
        if( (treeVars.jets_CM->at(ii)+treeVars.jets_NM->at(ii)) <= 1 ) skipJet = true;
      }
      
      if( fabs(treeVars.jets_eta->at(ii)) <= 2.4 )
      {
        if( treeVars.jets_CHF->at(ii)  <= 0.   ) skipJet = true;
        if( treeVars.jets_CM->at(ii)   <= 0    ) skipJet = true;
        if( treeVars.jets_CEMF->at(ii) >= 0.90 ) skipJet = true;
      }
      
      if( fabs(treeVars.jets_eta->at(ii)) >  2.7 &&
          fabs(treeVars.jets_eta->at(ii)) <= 3.0 )
      {
        if( treeVars.jets_NEMF->at(ii) <= 0.01 ) skipJet = true;
        if( treeVars.jets_NHF->at(ii)  >= 0.98 ) skipJet = true;
        if( treeVars.jets_NM->at(ii)   <= 1    ) skipJet = true;
      }
      
      if( fabs(treeVars.jets_eta->at(ii)) > 3.0 )
      {
        if( treeVars.jets_NEMF->at(ii) >= 0.90 ) skipJet = true;
        if( treeVars.jets_NM->at(ii)   <= 10   ) skipJet = true;
      }
      
      if( skipJet ) continue;
      
      
      ++jets_all_n;
      if( treeVars.jets_bTag->at(ii).at(0) > 0.460 ) ++jets_all_bTagL_n;
      if( treeVars.jets_bTag->at(ii).at(0) > 0.800 ) ++jets_all_bTagM_n;
      if( treeVars.jets_bTag->at(ii).at(0) > 0.935 ) ++jets_all_bTagT_n;
      
      if( fabs(treeVars.jets_eta->at(ii)) <= 2.4 )
      {
        ++jets_cen_n;
        if( treeVars.jets_bTag->at(ii).at(0) > 0.460 ) ++jets_cen_bTagL_n;
        if( treeVars.jets_bTag->at(ii).at(0) > 0.800 ) ++jets_cen_bTagM_n;
        if( treeVars.jets_bTag->at(ii).at(0) > 0.935 ) ++jets_cen_bTagT_n;
      }
      else
      {
        ++jets_fwd_n;
        if( treeVars.jets_bTag->at(ii).at(0) > 0.460 ) ++jets_fwd_bTagL_n;
        if( treeVars.jets_bTag->at(ii).at(0) > 0.800 ) ++jets_fwd_bTagM_n;
        if( treeVars.jets_bTag->at(ii).at(0) > 0.935 ) ++jets_fwd_bTagT_n;
      }
      
      if( jetId1_all != -1 && jetId2_all == -1 ) jetId2_all = ii;
      if( jetId1_all == -1 ) jetId1_all = ii;
      
      if( fabs(treeVars.jets_eta->at(ii)) <= 2.4 )
      {
        if( jetId1_cen != -1 && jetId2_cen == -1 ) jetId2_cen = ii;
        if( jetId1_cen == -1 ) jetId1_cen = ii;
      }
    }
    if( debugMode ) std::cout << ">>>>>> end jets" << std::endl;
    
    
    
    //---------------------------
    //--- fill out tree variables
    if( debugMode ) std::cout << ">>>>>> start filling out tree" << std::endl;
    outTree.vtxs_n = treeVars.vtxs_n;
    outTree.rho_all = treeVars.rho_all;
    
    outTree.mu1_pt = reco_mu.at(0).v.Pt();
    outTree.mu2_pt = reco_mu.at(1).v.Pt();
    outTree.mu1_eta = reco_mu.at(0).v.Eta();
    outTree.mu2_eta = reco_mu.at(1).v.Eta();
    outTree.mu1_phi = reco_mu.at(0).v.Phi();
    outTree.mu2_phi = reco_mu.at(1).v.Phi();
    outTree.mu1_dxy = treeVars.muons_dxy->at(reco_mu.at(0).it);
    outTree.mu2_dxy = treeVars.muons_dxy->at(reco_mu.at(1).it);
    outTree.mu1_dxyPull = fabs(treeVars.muons_dxy->at(reco_mu.at(0).it))/treeVars.muons_dxyErr->at(reco_mu.at(0).it);
    outTree.mu2_dxyPull = fabs(treeVars.muons_dxy->at(reco_mu.at(1).it))/treeVars.muons_dxyErr->at(reco_mu.at(1).it);
    outTree.mu1_dz = treeVars.muons_dz->at(reco_mu.at(0).it);
    outTree.mu2_dz = treeVars.muons_dz->at(reco_mu.at(1).it);
    outTree.mu1_dzPull = fabs(treeVars.muons_dz->at(reco_mu.at(0).it))/treeVars.muons_dzErr->at(reco_mu.at(0).it);
    outTree.mu2_dzPull = fabs(treeVars.muons_dz->at(reco_mu.at(1).it))/treeVars.muons_dzErr->at(reco_mu.at(1).it);
    outTree.mu1_relIso = mu1_iso/reco_mu.at(0).v.Pt();
    outTree.mu2_relIso = mu2_iso/reco_mu.at(1).v.Pt();
    outTree.mu1_isL = treeVars.muons_isLoose->at(reco_mu.at(0).it);
    outTree.mu2_isL = treeVars.muons_isLoose->at(reco_mu.at(1).it);
    outTree.mu1_isM = treeVars.muons_isMedium->at(reco_mu.at(0).it);
    outTree.mu2_isM = treeVars.muons_isMedium->at(reco_mu.at(1).it);
    outTree.mu1_isT = treeVars.muons_isTight->at(reco_mu.at(0).it);
    outTree.mu2_isT = treeVars.muons_isTight->at(reco_mu.at(1).it);
    outTree.mu_Deta = DeltaEta(reco_mu.at(0).v.Eta(),reco_mu.at(1).v.Eta());
    outTree.mu_Dphi = DeltaPhi(reco_mu.at(0).v.Phi(),reco_mu.at(1).v.Phi());
    outTree.mu_DR = DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi());
    
    outTree.H_pt = reco_H.v.Pt();
    outTree.H_eta = reco_H.v.Eta();
    outTree.H_phi = reco_H.v.Phi();
    outTree.H_mass = reco_H.v.M();
    
    outTree.jets_all_n = jets_all_n;
    outTree.jets_all_bTagL_n = jets_all_bTagL_n;
    outTree.jets_all_bTagM_n = jets_all_bTagM_n;
    outTree.jets_all_bTagT_n = jets_all_bTagT_n;
    outTree.jets_cen_n = jets_cen_n;
    outTree.jets_cen_bTagL_n = jets_cen_bTagL_n;
    outTree.jets_cen_bTagM_n = jets_cen_bTagM_n;
    outTree.jets_cen_bTagT_n = jets_cen_bTagT_n;
    outTree.jets_fwd_n = jets_fwd_n;
    outTree.jets_fwd_bTagL_n = jets_fwd_bTagL_n;
    outTree.jets_fwd_bTagM_n = jets_fwd_bTagM_n;
    outTree.jets_fwd_bTagT_n = jets_fwd_bTagT_n;
    
    TLorentzVector jet1_all;
    if( jetId1_all >= 0 )
    {
      jet1_all.SetPtEtaPhiE(treeVars.jets_pt->at(jetId1_all),
                            treeVars.jets_eta->at(jetId1_all),
                            treeVars.jets_phi->at(jetId1_all),
                            treeVars.jets_energy->at(jetId1_all));
      outTree.jet1_all_pt      = treeVars.jets_pt->at(jetId1_all);
      outTree.jet1_all_eta     = treeVars.jets_eta->at(jetId1_all);
      outTree.jet1_all_phi     = treeVars.jets_phi->at(jetId1_all);
      outTree.jet1_all_energy  = treeVars.jets_energy->at(jetId1_all);
    }
    else
    {
      outTree.jet1_all_pt      = -99.;
      outTree.jet1_all_eta     = -99.;
      outTree.jet1_all_phi     = -99.;
      outTree.jet1_all_energy  = -99.;
    }
    TLorentzVector jet2_all;
    if( jetId2_all >= 0 )
    {
      jet2_all.SetPtEtaPhiE(treeVars.jets_pt->at(jetId2_all),
                            treeVars.jets_eta->at(jetId2_all),
                            treeVars.jets_phi->at(jetId2_all),
                            treeVars.jets_energy->at(jetId2_all));
      outTree.jet2_all_pt      = treeVars.jets_pt->at(jetId2_all);
      outTree.jet2_all_eta     = treeVars.jets_eta->at(jetId2_all);
      outTree.jet2_all_phi     = treeVars.jets_phi->at(jetId2_all);
      outTree.jet2_all_energy  = treeVars.jets_energy->at(jetId2_all);
    }
    else
    {
      outTree.jet2_all_pt      = -99.;
      outTree.jet2_all_eta     = -99.;
      outTree.jet2_all_phi     = -99.;
      outTree.jet2_all_energy  = -99.;
    }
    if( jetId1_all >= 0 && jetId2_all >= 0 )
    {
      outTree.jet_all_Deta     = DeltaEta(treeVars.jets_eta->at(jetId1_all),treeVars.jets_eta->at(jetId2_all));
      outTree.jet_all_Dphi     = DeltaPhi(treeVars.jets_phi->at(jetId1_all),treeVars.jets_phi->at(jetId2_all));
      outTree.jet_all_mass     = (jet1_all+jet2_all).M();
    }
    else
    {
      outTree.jet_all_Deta     = -99.;
      outTree.jet_all_Dphi     = -99.;
      outTree.jet_all_mass     = -99.;
    }
    
    TLorentzVector jet1_cen;
    if( jetId1_cen >= 0 )
    {
      jet1_cen.SetPtEtaPhiE(treeVars.jets_pt->at(jetId1_cen),
                            treeVars.jets_eta->at(jetId1_cen),
                            treeVars.jets_phi->at(jetId1_cen),
                            treeVars.jets_energy->at(jetId1_cen));
      outTree.jet1_cen_pt      = treeVars.jets_pt->at(jetId1_cen);
      outTree.jet1_cen_eta     = treeVars.jets_eta->at(jetId1_cen);
      outTree.jet1_cen_phi     = treeVars.jets_phi->at(jetId1_cen);
      outTree.jet1_cen_energy  = treeVars.jets_energy->at(jetId1_cen);
    }
    else
    {
      outTree.jet1_cen_pt      = -99.;
      outTree.jet1_cen_eta     = -99.;
      outTree.jet1_cen_phi     = -99.;
      outTree.jet1_cen_energy  = -99.;
    }
    TLorentzVector jet2_cen;
    if( jetId2_cen >= 0 )
    {
      jet2_cen.SetPtEtaPhiE(treeVars.jets_pt->at(jetId2_cen),
                            treeVars.jets_eta->at(jetId2_cen),
                            treeVars.jets_phi->at(jetId2_cen),
                            treeVars.jets_energy->at(jetId2_cen));
      outTree.jet2_cen_pt      = treeVars.jets_pt->at(jetId2_cen);
      outTree.jet2_cen_eta     = treeVars.jets_eta->at(jetId2_cen);
      outTree.jet2_cen_phi     = treeVars.jets_phi->at(jetId2_cen);
      outTree.jet2_cen_energy  = treeVars.jets_energy->at(jetId2_cen);
    }
    else
    {
      outTree.jet2_cen_pt      = -99.;
      outTree.jet2_cen_eta     = -99.;
      outTree.jet2_cen_phi     = -99.;
      outTree.jet2_cen_energy  = -99.;
    }
    if( jetId1_cen >= 0 && jetId2_cen >= 0 )
    {
      outTree.jet_cen_Deta     = DeltaEta(treeVars.jets_eta->at(jetId1_cen),treeVars.jets_eta->at(jetId2_cen));
      outTree.jet_cen_Dphi     = DeltaPhi(treeVars.jets_phi->at(jetId1_cen),treeVars.jets_phi->at(jetId2_cen));
      outTree.jet_cen_mass     = (jet1_cen+jet2_cen).M();
    }
    else
    {
      outTree.jet_cen_Deta     = -99.;
      outTree.jet_cen_Dphi     = -99.;
      outTree.jet_cen_mass     = -99.;
    }
    
    outTree.met_pt = treeVars.met_pt;
    outTree.met_phi = treeVars.met_phi;
    outTree.met_sig = treeVars.met_sig;
    
    outTree.weight = weight;
    outTree.weight_MC = mcWeight;
    outTree.weight_PU = puWeight;
    
    outTree.GetTTreePtr()->Fill();
    if( debugMode ) std::cout << ">>>>>> end filling out tree" << std::endl;
    
    
    
    //---------------
    //--- print event
    if( debugMode ) std::cout << ">>>>>> start print event" << std::endl;
    if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
      std::cout << "\n\n Event: " << entry << std::endl;
    if( printGenEvent && isSignal )
    {
      std::cout << "--------------------- GEN ---------------------" << std::endl;
      PrintEvent(gen_mu);
    }
    if( printRecoGenMatchEvent && isSignal )
    {
      std::cout << "--------------------- GEN MATCH ---------------------" << std::endl;
      PrintEvent(recoGenMatch_mu);
    }
    if( printRecoEvent )
    {
      std::cout << "--------------------- RECO ---------------------" << std::endl;
      PrintEvent(reco_mu);
    }
    if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
      std::cout << std::endl;
    if( debugMode ) std::cout << ">>>>>> end print event" << std::endl;
    
    
    
    //-------------------
    //--- fill histograms
    if( isSignal )
      for(int ii = 0; ii < 2; ++ii)
      {
        h1_reco_mu_gen_mu_DR -> Fill( DeltaR(reco_mu.at(ii).v.Eta(),reco_mu.at(ii).v.Phi(),gen_mu.at(ii).v.Eta(),gen_mu.at(ii).v.Phi()),weight );
        h1_reco_mu_gen_mu_ptRatio -> Fill( reco_mu.at(ii).v.Pt()/gen_mu.at(ii).v.Pt(),weight );
      }
    
  } // loop over events
  std::cout << std::endl;
  
  
  
  // if( !isData )
  // {
  //   int N = 15;
  //   int v[15]={516,193,200,357,127,136,140,142,143,175,
  //              184,511,187,506, 26};
  //   for(int ii = 0; ii < N; ++ii)
  //   {
  //     std::pair<std::string,int> p = vec_triggerPass_noCuts.at(v[ii]); 
      
  //     std::string binLabel = p.first;
  //     binLabel.erase(0,4);
  //     binLabel.erase(binLabel.end()-3, binLabel.end());
      
  //     h1_triggerEff_noCuts -> SetBinContent(ii+1,1.*p.second/nTriggerEvents_noCuts);
  //     h1_triggerEff_noCuts -> GetXaxis() -> SetBinLabel(ii+1,binLabel.c_str());
      
  //     p = vec_triggerPass_genCuts.at(v[ii]); 
  //     h1_triggerEff_genCuts -> SetBinContent(ii+1,1.*p.second/nTriggerEvents_genCuts);
  //     h1_triggerEff_genCuts -> GetXaxis() -> SetBinLabel(ii+1,binLabel.c_str());
      
  //     p = vec_triggerPass_recoMuSelection.at(v[ii]); 
  //     h1_triggerEff_recoMuSelection -> SetBinContent(ii+1,1.*p.second/nTriggerEvents_recoMuSelection);
  //     h1_triggerEff_recoMuSelection -> GetXaxis() -> SetBinLabel(ii+1,binLabel.c_str());
      
  //     p = vec_triggerPass_recoJPsiMassSelection.at(v[ii]); 
  //     h1_triggerEff_recoJPsiMassSelection -> SetBinContent(ii+1,1.*p.second/nTriggerEvents_recoJPsiMassSelection);
  //     h1_triggerEff_recoJPsiMassSelection -> GetXaxis() -> SetBinLabel(ii+1,binLabel.c_str());
      
  //     p = vec_triggerPass_recoJPsiPtSelection.at(v[ii]); 
  //     h1_triggerEff_recoJPsiPtSelection -> SetBinContent(ii+1,1.*p.second/nTriggerEvents_recoJPsiPtSelection);
  //     h1_triggerEff_recoJPsiPtSelection -> GetXaxis() -> SetBinLabel(ii+1,binLabel.c_str());
  //   }
  
  
  // sort triggers
  std::sort(vec_triggerPass.begin(),vec_triggerPass.end(),PairSort());
  for(unsigned int ii = 0; ii < vec_triggerPass.size(); ++ii)
  {
    std::cout << std::fixed;
    std::cout << std::setw(100) << vec_triggerPass.at(ii).first << "   ";
    std::cout << std::setw(5)  << vec_triggerPass.at(ii).second << "   ";
    std::cout << std::setw(5)  << 1.*vec_triggerPass.at(ii).second/nEvents_selected << "   ";
    std::cout << std::endl;
  }
  
  
  
  h1_nEvents -> GetXaxis() -> SetBinLabel(1,"tot"); h1_nEvents -> SetBinContent(1,nEvents_tot);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(2,"gen. cut on muons"); h1_nEvents -> SetBinContent(2,nEvents_genCutMu);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(3,"2 gen.-matched muons"); h1_nEvents -> SetBinContent(3,nEvents_2RecoGenMatchMu);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(4,"gen. cut on k/#pi"); h1_nEvents -> SetBinContent(4,nEvents_genCutKepi);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(5,"3 gen.-matched k/#pi"); h1_nEvents -> SetBinContent(5,nEvents_3RecoGenMatchKepi);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(6,"full gen. cut"); h1_nEvents -> SetBinContent(6,nEvents_genCut);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(7,"full gen. match"); h1_nEvents -> SetBinContent(7,nEvents_recoGenMatch);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(8,"2 reco muons"); h1_nEvents -> SetBinContent(8,nEvents_2RecoMu);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(9,"2 gen.-matched reco muons"); h1_nEvents -> SetBinContent(9,nEvents_2RecoMu_genMatch);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(10,"3 reco k/#pi"); h1_nEvents -> SetBinContent(10,nEvents_3RecoKepi);
  // h1_nEvents -> GetXaxis() -> SetBinLabel(11,"3 gen.-matched reco k/#pi"); h1_nEvents -> SetBinContent(11,nEvents_3RecoKepi_genMatch);
  
  std::cout << std::fixed << std::endl;
  std::cout << "===================================================="  << std::endl;
  std::cout << "                              Tot. events: "    << std::setw(6) << nEvents_tot << std::endl;
  std::cout << std::endl;
  std::cout << "             Events with 2 gen mu in acceptance: " << std::setw(6) << nEvents_genCut          << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_genCut/nEvents_tot             << "%)" << std::endl;
  // std::cout << ">>>>>>       Events with 2 reco gen. matched mu: " << std::setw(6) << nEvents_2RecoGenMatchMu << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoGenMatchMu/nEvents_genCutMu << "%)" << std::endl;
  // std::cout << std::endl;
  // std::cout << "           Events with 2 reco. muons: " << std::setw(6) << nEvents_2RecoMu << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoMu/nEvents_tot << "%)"
  //           << "   + matching gen.: " << std::setw(6) << nEvents_2RecoMu_genMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoMu_genMatch/nEvents_2RecoMu << "%)" << std::endl;
  // std::cout << ">>>          Events with 3 reco. tracks: " << std::setw(6) << nEvents_3RecoKepi << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoKepi/nEvents_tot << "%)"
  //           << "   + matching gen.: " << std::setw(6) << nEvents_3RecoKepi_genMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoKepi_genMatch/nEvents_3RecoKepi << "%)" << std::endl;
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
