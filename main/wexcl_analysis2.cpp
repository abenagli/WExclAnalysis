#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "src/RoccoR.cc"
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
    std::cerr << ">>>>> wexcl_analysis2.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //--- open files and get trees
  int isData = opts.GetOpt<int>("Input.isData");
  std::vector<std::string> inputFiles = opts.GetOpt<std::vector<std::string> >("Input.inputFiles");
  
  TChain* chain_reco = new TChain("DumpReco/reco_tree","DumpReco/reco_tree");
  TChain* chain_gen = NULL; if( !isData ) chain_gen = new TChain("DumpGenParticles/gen_tree","DumpGenParticles/gen_tree");

  for(unsigned int fileIt = 0; fileIt < inputFiles.size(); ++fileIt)
  {
    chain_reco -> Add(inputFiles.at(fileIt).c_str());
    if( !isData ) chain_gen -> Add(inputFiles.at(fileIt).c_str());
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
  float gen_kepi_ptMin = opts.GetOpt<float>("Cuts.gen_kepi_ptMin");
  float gen_kepi_etaMax = opts.GetOpt<float>("Cuts.gen_kepi_etaMax");
  int genCut = opts.GetOpt<int>("Cuts.genCut");
  
  int DsDecayModeCut = opts.GetOpt<int>("Cuts.DsDecayModeCut");
  std::vector<std::string> DsDecayModes = opts.GetOpt<std::vector<std::string> >("Cuts.DsDecayModes");
  
  int genMatch = opts.GetOpt<int>("Cuts.genMatch");  
  int genMatchCut = opts.GetOpt<int>("Cuts.genMatchCut");
  
  int reco_HLTCut = opts.GetOpt<int>("Cuts.reco_HLTCut");
  std::vector<std::string> reco_HLT_paths = opts.GetOpt<std::vector<std::string> >("Cuts.reco_HLT_paths");
  int reco_HLT_OR = opts.GetOpt<int>("Cuts.reco_HLT_OR");
  
  float reco_ptMin1_mu = opts.GetOpt<float>("Cuts.reco_ptMin1_mu");
  float reco_ptMin2_mu = opts.GetOpt<float>("Cuts.reco_ptMin2_mu");
  float reco_ptMin_JPsi = opts.GetOpt<float>("Cuts.reco_ptMin_JPsi");
  
  float reco_ptMin1_kepi = opts.GetOpt<float>("Cuts.reco_ptMin1_kepi");
  float reco_ptMin2_kepi = opts.GetOpt<float>("Cuts.reco_ptMin2_kepi");
  float reco_ptMin3_kepi = opts.GetOpt<float>("Cuts.reco_ptMin3_kepi");
  float reco_DRMax_kepi = opts.GetOpt<float>("Cuts.reco_DRMax_kepi");
  float reco_DPhiMin_kepi_Jpsi = opts.GetOpt<float>("Cuts.reco_DPhiMin_kepi_Jpsi");
  float reco_ptMin_Ds = opts.GetOpt<float>("Cuts.reco_ptMin_Ds");
  
  int noMuonCuts = opts.GetOpt<int>("Cuts.noMuonCuts");
  int useResonances = opts.GetOpt<int>("Cuts.useResonances");
  int wrongDsCharge = opts.GetOpt<int>("Cuts.wrongDsCharge");
  int DsSidebands = opts.GetOpt<int>("Cuts.DsSidebands");
  int antiDsIso = opts.GetOpt<int>("Cuts.antiDsIso");
  
  
  //--- define histograms
  
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  std::string noMuonCutsLabel = noMuonCuts ? "__noMuonCuts" : "";
  std::string useResonancesLabel = useResonances ? "__resonances" : "";
  std::string wrongDsChargeLabel = wrongDsCharge ? "__wrongDsCharge" : "";
  std::string DsSidebandsLabel = DsSidebands ? "__DsSidebands" : "";
  std::string antiDsIsoLabel = antiDsIso ? "__antiDsIso" : "";
  TFile* outFile = TFile::Open(Form("%s__DRMin_Ds_JPsi_%.2f__DRMax_Ds_%.2f%s%s%s%s%s.root",outputFileName.c_str(),reco_DPhiMin_kepi_Jpsi,reco_DRMax_kepi,
                                    noMuonCutsLabel.c_str(),useResonancesLabel.c_str(),wrongDsChargeLabel.c_str(),DsSidebandsLabel.c_str(),antiDsIsoLabel.c_str()),"RECREATE");
  outFile -> cd();
  
  StdHistoSet* genHistos = isData? NULL : new StdHistoSet("gen",outFile);;
  StdHistoSet* genCutHistos = isData? NULL : new StdHistoSet("genCut",outFile);;
  StdHistoSet* recoGenMatchHistos = isData? NULL : new StdHistoSet("recoGenMatch",outFile);;
  StdHistoSet* recoGenMatchMassAssignHistos = isData? NULL : new StdHistoSet("recoGenMatchMassAssign",outFile);;
  StdHistoSet* recoHistos = new StdHistoSet("reco",outFile);
  StdHistoSet* recoMassAssignHistos = new StdHistoSet("recoMassAssign",outFile);
  StdHistoSet* recoMassAssignCutJPsiHistos = new StdHistoSet("recoMassAssignCutJPsi",outFile);
  StdHistoSet* recoMassAssignCutJPsiDsIso1Histos = new StdHistoSet("recoMassAssignCutJPsiDsIso1",outFile);
  StdHistoSet* recoMassAssignCutJPsiDsIso2Histos = new StdHistoSet("recoMassAssignCutJPsiDsIso2",outFile);
  StdHistoSet* recoMassAssignCutJPsiDsHistos = new StdHistoSet("recoMassAssignCutJPsiDs",outFile);
  
  TH1F* h1_nEvents = new TH1F("h1_nEvents",";;entries",11,0.,11.);
  
  TH2F* h2_mu_pt2_vs_pt1_total = new TH2F("h2_mu_pt2_vs_pt1_total","",100,-0.5,99.5,100,-0.5,99.5);
  TH2F* h2_mu_pt2_vs_pt1_pass  = new TH2F("h2_mu_pt2_vs_pt1_pass", "",100,-0.5,99.5,100,-0.5,99.5);
  TH2F* h2_mu_cuts_pt2_vs_pt1_total = new TH2F("h2_mu_cuts_pt2_vs_pt1_total","",100,-0.5,99.5,100,-0.5,99.5);
  TH2F* h2_mu_cuts_pt2_vs_pt1_pass  = new TH2F("h2_mu_cuts_pt2_vs_pt1_pass", "",100,-0.5,99.5,100,-0.5,99.5);
  
  TEfficiency* p1_eff_mu_pt = new TEfficiency("p1_eff_mu_pt",";muon p_{T} (GeV);#epsilon",200,0.,500.);
  TEfficiency* p1_eff_mu_eta = new TEfficiency("p1_eff_mu_eta",";muon #eta;#epsilon",200,0.,10.);
  
  TEfficiency* p1_eff_kepi_pt = new TEfficiency("p1_eff_kepi_pt",";K/#pi p_{T} (GeV);#epsilon",200,0.,500.);
  TEfficiency* p1_eff_kepi_eta = new TEfficiency("p1_eff_kepi_eta",";K/#pi #eta;#epsilon",200,0.,10.);
  
  TH1F* h1_reco_mu_gen_mu_DR = new TH1F("h1_reco_mu_gen_mu_DR","",10000,0.,10.);
  TH1F* h1_reco_mu_gen_mu_ptRatio = new TH1F("h1_reco_mu_gen_mu_ptRatio","",10000,0.,10.);
  
  TH1F* h1_kk_mass = new TH1F("h1_kk_mass","",5000,0.,5.);
  TH1F* h1_kpi_mass = new TH1F("h1_kpi_mass","",5000,0.,5.);
  
  std::vector<std::pair<std::string,int> > vec_triggerPass;
  
  
  
  //----------------
  //--- event weight

  float weight = 1.;
  if( !isData )
  {
    float lumi = opts.GetOpt<float>("Weights.lumi");
    float xsec = opts.GetOpt<float>("Weights.xsec");
    float nevt = opts.GetOpt<int>("Weights.nevt");
    
    float mcWeight = 1. * lumi * xsec / nevt;
    weight *= mcWeight;

  }
  
  
  
  // //----------------------
  // //--- pileup reweighting

  // if( !isData )
  // {
  //   std::string pileupFileName_data = opts.GetOpt<std::string>("Weights.pileupFileName_data");
  //   TFile* pileupFile_data = TFile::Open(pileupFileName_data.c_str(),"READ");
  //   TH1F* pileup_data = (TH1F*)( pileupFile_data->Get("pileup") );
    
  //   TH1F* pileup_mc;
  //   for(unsigned int ii = 0; ii < inputFiles.size(); ++ii)
  //   {
  //     TFile* tempInFle = TFile::Open(inputFiles.at(ii).c_str(),"READ");
  //     if( ii == 0 ) pileup_mc = (TH1F*)( tempInFile->Get("DumpPU/pileup") );
  //     else
  //     {
  //       TH1F* tempHisto = (TH1F*)( tempInFile->Get("DumpPU/pileup") );
  //       pileup_mc -> Add(tempHisto);
  //     }
  //   }
  // }
  
  
  
  //------------------------------
  //--- muon Rochester corrections
  RoccoR rc("data/rcdata.2016.v3");
  
  
  
  //--------------------
  //--- loop over events
  int nEvents_tot = 0;
  int nEvents_genCutMu = 0;
  int nEvents_genCutKepi = 0;
  int nEvents_genCut = 0;
  int nEvents_2RecoGenMatchMu = 0;
  int nEvents_3RecoGenMatchKepi = 0;
  int nEvents_recoGenMatch = 0;
  int nEvents_2RecoMu = 0; int nEvents_2RecoMu_genMatch = 0;
  int nEvents_3RecoKepi = 0; int nEvents_3RecoKepi_genMatch = 0;
  int nEvents_selected = 0;
  
  //for(int entry = 0; entry < 100; ++entry)
  for(int entry = 0; entry < chain_reco->GetEntries(); ++entry)
  {
    if( entry%100 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_recoTree << "\r" << std::flush;
    // if( entry%1 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_recoTree << std::endl;
    chain_reco -> GetEntry(entry);
    if( !isData ) chain_gen -> GetEntry(entry);
    ++nEvents_tot;
    
    
    //--------------
    //--- candidates
    
    particle gen_W;
    particle gen_JPsi;
    particle gen_Ds;
    std::vector<particle> gen_mu;
    std::vector<particle> gen_kepi;
    std::vector<particle> gen_mu_etaSort;
    std::vector<particle> gen_kepi_etaSort;
    
    std::vector<particle> recoGenMatch_mu;
    std::vector<particle> recoGenMatch_kepi;
    std::vector<particle> recoGenMatch_kepi_massAssign;

    particle reco_JPsi;
    particle reco_Ds;
    std::vector<particle> reco_mu;
    std::vector<particle> reco_kepi;
    std::vector<particle> reco_kepi_massAssign;
    
    
    
    //-----------------------------
    //--- fill particle collections
    std::vector<particle> temp_mu;
    for(int ii = 0; ii < treeVars.muons_n; ++ii)
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
    
    std::vector<particle> temp_chargedHadron;
    std::vector<particle> temp_neutralHadron;
    std::vector<particle> temp_photon;
    std::vector<particle> temp_PU;
    for(int ii = 0; ii < treeVars.pfCandidates_n; ++ii)
    {
      particle temp;
      temp.v.SetPtEtaPhiE(treeVars.pfCandidates_pt->at(ii),treeVars.pfCandidates_eta->at(ii),treeVars.pfCandidates_phi->at(ii),treeVars.pfCandidates_energy->at(ii));
      temp.charge = treeVars.pfCandidates_charge->at(ii);
      temp.pdgId = treeVars.pfCandidates_pdgId->at(ii);
      temp.it = ii;
      
      if( abs(treeVars.pfCandidates_pdgId->at(ii)) == 211 &&
          treeVars.pfCandidates_fromPV->at(ii) >= 2 )
        temp_chargedHadron.push_back(temp);
      
      if( abs(treeVars.pfCandidates_pdgId->at(ii)) == 130 &&
          treeVars.pfCandidates_pt->at(ii) > 0.5 )
        temp_neutralHadron.push_back(temp);
      
      if( abs(treeVars.pfCandidates_pdgId->at(ii)) == 22 &&
          treeVars.pfCandidates_pt->at(ii) > 0.5 )
        temp_photon.push_back(temp);
      
      if( abs(treeVars.pfCandidates_pdgId->at(ii)) == 211 &&
          treeVars.pfCandidates_fromPV->at(ii) < 2 &&
          treeVars.pfCandidates_pt->at(ii) > 0.5 )
        temp_PU.push_back(temp);        
    }
    std::sort(temp_chargedHadron.begin(),temp_chargedHadron.end(),PtSort());
    
    
    std::string DsDecayMode = "";
    if( !isData )
    {
      //---------------------
      //--- retrieve gen info
      
      gen_W.v.SetPtEtaPhiE(treeVars.reso_pt->at(0),treeVars.reso_eta->at(0),treeVars.reso_phi->at(0),treeVars.reso_energy->at(0));
      gen_W.charge = treeVars.reso_charge->at(0);
      gen_W.pdgId = treeVars.reso_pdgId->at(0);
      
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
              gen_mu_etaSort.push_back( temp );
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
            if( abs((*treeVars.resoDau2_pdgId)[0][ii][jj]) != 22 )
            {
              particle temp;
              temp.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj],(*treeVars.resoDau2_eta)[0][ii][jj],(*treeVars.resoDau2_phi)[0][ii][jj],(*treeVars.resoDau2_energy)[0][ii][jj]);
              temp.charge = (*treeVars.resoDau2_charge)[0][ii][jj];
              temp.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj];
              
              gen_kepi.push_back( temp );
              gen_kepi_etaSort.push_back( temp );
            }
          }
        }
      }
      
      
      std::sort(gen_kepi.begin(),gen_kepi.end(),NameSort());
      
      for(unsigned int ii = 0; ii < gen_kepi.size(); ++ii)
      {
        if( ii < gen_kepi.size()-1 )
          DsDecayMode += GetParticleName(gen_kepi.at(ii).pdgId) + "_";
        else
          DsDecayMode += GetParticleName(gen_kepi.at(ii).pdgId);
      }
      
      if( DsDecayModeCut )
      {
        bool skipEvent = true;
        for(unsigned int ii = 0; ii < DsDecayModes.size(); ++ii)
        {
          if( DsDecayMode == DsDecayModes.at(ii) )
          {
            skipEvent = false;
            break;
          }
        }
        // std::cout << "decayMode: " << DsDecayMode << "   skipEvent: " << skipEvent << std::endl;
        if( skipEvent ) continue;
      }
      
      
      std::sort(gen_mu.begin(),gen_mu.end(),PtSort());
      std::sort(gen_kepi.begin(),gen_kepi.end(),PtSort());
      genHistos -> FillHistos(weight,gen_mu,gen_kepi);
      
      
      bool accept_genCutMu = false;
      bool accept_genCutKepi = false;
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
          ++nEvents_genCutMu;
          accept_genCutMu = true;
        }
        
        if( gen_kepi.at(0).v.Pt() > gen_kepi_ptMin &&
            gen_kepi.at(1).v.Pt() > gen_kepi_ptMin &&
            gen_kepi.at(2).v.Pt() > gen_kepi_ptMin &&
            fabs(gen_kepi.at(0).v.Eta()) < gen_kepi_etaMax &&
            fabs(gen_kepi.at(1).v.Eta()) < gen_kepi_etaMax &&
            fabs(gen_kepi.at(2).v.Eta()) < gen_kepi_etaMax )
        {
          ++nEvents_genCutKepi;
          accept_genCutKepi = true;
        }
        
        if( gen_mu.at(0).v.Pt() > gen_mu_ptMin &&
            gen_mu.at(1).v.Pt() > gen_mu_ptMin &&
            fabs(gen_mu.at(0).v.Eta()) < gen_mu_etaMax &&
            fabs(gen_mu.at(1).v.Eta()) < gen_mu_etaMax &&
            gen_kepi.at(0).v.Pt() > gen_kepi_ptMin &&
            gen_kepi.at(1).v.Pt() > gen_kepi_ptMin &&
            gen_kepi.at(2).v.Pt() > gen_kepi_ptMin &&
            fabs(gen_kepi.at(0).v.Eta()) < gen_kepi_etaMax &&
            fabs(gen_kepi.at(1).v.Eta()) < gen_kepi_etaMax &&
            fabs(gen_kepi.at(2).v.Eta()) < gen_kepi_etaMax )
        {
          ++nEvents_genCut;
          accept_genCut = true;
          genCutHistos -> FillHistos(weight,gen_mu,gen_kepi);
        }
      }
      
      
      if( genMatch )
      {
        //-----------------------
        //--- match reco with gen      
        std::vector<int> muIt_genMatch;
        muIt_genMatch.push_back( GetBestMatch(gen_mu.at(0),temp_mu,&muIt_genMatch) );
        muIt_genMatch.push_back( GetBestMatch(gen_mu.at(1),temp_mu,&muIt_genMatch) );
        if( muIt_genMatch.at(0) >= 0 ) recoGenMatch_mu.push_back( temp_mu.at(muIt_genMatch.at(0)) );
        if( muIt_genMatch.at(1) >= 0 ) recoGenMatch_mu.push_back( temp_mu.at(muIt_genMatch.at(1)) );
        
        std::vector<int> kepiIt_genMatch;
        kepiIt_genMatch.push_back( GetBestMatch(gen_kepi.at(0),temp_chargedHadron,&kepiIt_genMatch) );
        kepiIt_genMatch.push_back( GetBestMatch(gen_kepi.at(1),temp_chargedHadron,&kepiIt_genMatch) );
        kepiIt_genMatch.push_back( GetBestMatch(gen_kepi.at(2),temp_chargedHadron,&kepiIt_genMatch) );
        if( kepiIt_genMatch.at(0) >= 0 ) recoGenMatch_kepi.push_back( temp_chargedHadron.at(kepiIt_genMatch.at(0)) );
        if( kepiIt_genMatch.at(1) >= 0 ) recoGenMatch_kepi.push_back( temp_chargedHadron.at(kepiIt_genMatch.at(1)) );
        if( kepiIt_genMatch.at(2) >= 0 ) recoGenMatch_kepi.push_back( temp_chargedHadron.at(kepiIt_genMatch.at(2)) );
        
        bool recoGenMatch_mu_matching = false;
        bool recoGenMatch_kepi_matching = false;
        
        if( recoGenMatch_mu.size() == 2 )
        {
          recoGenMatch_mu_matching = IsMatching(recoGenMatch_mu,gen_mu,0.01,0.1,false);
          if( (genCut && accept_genCutMu) || (!genCut) ) 
            if( recoGenMatch_mu_matching )
            {
              ++nEvents_2RecoGenMatchMu;
              
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
        
        if( recoGenMatch_kepi.size() == 3 )
        {
          recoGenMatch_kepi_matching = IsMatching(recoGenMatch_kepi,gen_kepi,0.01,0.1,false);
          if( (genCut && accept_genCutKepi) || (!genCut) )
            if( recoGenMatch_kepi_matching )
              ++nEvents_3RecoGenMatchKepi;
        }
        
        if( recoGenMatch_mu_matching && recoGenMatch_kepi_matching )
        {
          if( (genCut && accept_genCut) || (!genCut) )
            ++nEvents_recoGenMatch;
          
          recoGenMatchHistos -> FillHistos(weight,recoGenMatch_mu,recoGenMatch_kepi,&treeVars);
          
          recoGenMatch_kepi_massAssign = AssignMassDs(recoGenMatch_kepi);
          recoGenMatchMassAssignHistos -> FillHistos(weight,recoGenMatch_mu,recoGenMatch_kepi_massAssign,&treeVars);
        }
        
        
        if( recoGenMatch_mu.size() > 0 )
        {
          if( IsMatching(gen_mu.at(0),recoGenMatch_mu.at(0),0.01,0.1) )
          {
            if( fabs(gen_mu.at(0).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(true,gen_mu.at(0).v.Pt());
            if(gen_mu.at(0).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(true,fabs(gen_mu.at(0).v.Eta()));
          }
          else
          {
            if( fabs(gen_mu.at(0).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(0).v.Pt());
            if(gen_mu.at(0).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(0).v.Eta()));
          }
        }
        else
        {
          if( fabs(gen_mu.at(0).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(0).v.Pt());
          if(gen_mu.at(0).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(0).v.Eta()));
          if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(1).v.Pt());
          if(gen_mu.at(1).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(1).v.Eta()));
        }
        
        if( recoGenMatch_mu.size() > 1 )
        {
          if( IsMatching(gen_mu.at(1),recoGenMatch_mu.at(1),0.01,0.1) )
          {
            if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(true,gen_mu.at(1).v.Pt());
            if(gen_mu.at(1).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(true,fabs(gen_mu.at(1).v.Eta()));
          }
          else
          {
            if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(1).v.Pt());
            if(gen_mu.at(1).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(1).v.Eta()));
          }
        }
        else
        {
          if( fabs(gen_mu.at(1).v.Eta()) < 2.1 ) p1_eff_mu_pt  -> Fill(false,gen_mu.at(1).v.Pt());
          if(gen_mu.at(1).v.Pt() > 5. ) p1_eff_mu_eta -> Fill(false,fabs(gen_mu.at(1).v.Eta()));
        }
        
        
        if( recoGenMatch_kepi.size() > 0 )
        {
          if( IsMatching(gen_kepi.at(0),recoGenMatch_kepi.at(0),0.01,0.1) )
          {
            if(fabs(gen_kepi.at(0).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(true,gen_kepi.at(0).v.Pt());
            if( gen_kepi.at(0).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(true,fabs(gen_kepi.at(0).v.Eta()));
          }
          else
          {
            if(fabs(gen_kepi.at(0).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(0).v.Pt());
            if( gen_kepi.at(0).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(0).v.Eta()));
          }
        }
        else
        {
          if(fabs(gen_kepi.at(0).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(0).v.Pt());
          if( gen_kepi.at(0).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(0).v.Eta()));
          if(fabs(gen_kepi.at(1).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(1).v.Pt());
          if( gen_kepi.at(1).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(1).v.Eta()));
          if(fabs(gen_kepi.at(2).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(2).v.Pt());
          if( gen_kepi.at(2).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(2).v.Eta()));
        }
        
        if( recoGenMatch_kepi.size() > 1 )
        {
          if( IsMatching(gen_kepi.at(1),recoGenMatch_kepi.at(1),0.01,0.1) )
          {
            if(fabs(gen_kepi.at(1).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(true,gen_kepi.at(1).v.Pt());
            if( gen_kepi.at(1).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(true,fabs(gen_kepi.at(1).v.Eta()));
          }
          else
          {
            if(fabs(gen_kepi.at(1).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(1).v.Pt());
            if( gen_kepi.at(1).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(1).v.Eta()));
          }
        }
        else
        {
          if(fabs(gen_kepi.at(1).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(1).v.Pt());
          if( gen_kepi.at(1).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(1).v.Eta()));
          if(fabs(gen_kepi.at(2).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(2).v.Pt());
          if( gen_kepi.at(2).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(2).v.Eta()));
        }
        
        if( recoGenMatch_kepi.size() > 2 )
        {
          if( IsMatching(gen_kepi.at(2),recoGenMatch_kepi.at(2),0.01,0.1) )
          {
            if(fabs(gen_kepi.at(2).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(true,gen_kepi.at(2).v.Pt());
            if( gen_kepi.at(2).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(true,fabs(gen_kepi.at(2).v.Eta()));
          }
          else
          {
            if(fabs(gen_kepi.at(2).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(2).v.Pt());
            if( gen_kepi.at(2).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(2).v.Eta()));
          }
        }
        else
        {
          if(fabs(gen_kepi.at(2).v.Eta()) < 1.0) p1_eff_kepi_pt  -> Fill(false,gen_kepi.at(2).v.Pt());
          if( gen_kepi.at(2).v.Pt() > 5. ) p1_eff_kepi_eta -> Fill(false,fabs(gen_kepi.at(2).v.Eta()));
        }
        
        if( (genMatchCut) && (!recoGenMatch_mu_matching || !recoGenMatch_kepi_matching) )
          continue;
      }
    }
    
    
    
    //---------------------
    //--- trigger selection
    bool trgPassAll = true;
    bool trgPassOne = false;
    for(int ii = 0; ii < treeVars.trgs_n; ++ii)
    {
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
    
    
    
    //---------------------
    //--- select reco muons
    if( temp_mu.size() < 2 ) continue;
    
    for(unsigned int ii = 0; ii < temp_mu.size(); ++ii)
    {
      if( treeVars.muons_isTight->at(temp_mu.at(ii).it) != 1 ) continue;
      if( !noMuonCuts && fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)/treeVars.muons_dxyErr->at(temp_mu.at(ii).it)) > 5. ) continue;
      if( !noMuonCuts && fabs(treeVars.muons_dz->at(temp_mu.at(ii).it)/treeVars.muons_dzErr->at(temp_mu.at(ii).it)) > 5. ) continue;
      
      treeVars.mu1_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
      std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
      if( !noMuonCuts && treeVars.mu1_iso/temp_mu.at(ii).v.Pt() > 0.25 ) continue;
      
      reco_mu.push_back( temp_mu.at(ii) );
      break;
    }
    if( reco_mu.size() < 1 ) continue;
    
    for(unsigned int ii = 0; ii < temp_mu.size(); ++ii)
    {
      if( reco_mu.at(0).charge*temp_mu.at(ii).charge == 1 ) continue;
      if( treeVars.muons_isMedium->at(temp_mu.at(ii).it) != 1 ) continue;
      if( !noMuonCuts && fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)/treeVars.muons_dxyErr->at(temp_mu.at(ii).it)) > 5. ) continue;
      if( !noMuonCuts && fabs(treeVars.muons_dz->at(temp_mu.at(ii).it)/treeVars.muons_dzErr->at(temp_mu.at(ii).it)) > 5. ) continue;
      if( DeltaEta(temp_mu.at(ii).v.Eta(),reco_mu.at(0).v.Eta()) > 0.5 ) continue;
      
      treeVars.mu2_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
      std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
      if( !noMuonCuts && treeVars.mu2_iso/temp_mu.at(ii).v.Pt() > 0.25 ) continue;
      
      reco_mu.push_back( temp_mu.at(ii) );
      break;
    }
    if( reco_mu.size() < 2 ) continue;
    
    reco_JPsi.v = reco_mu.at(0).v + reco_mu.at(1).v;
    reco_JPsi.charge = reco_mu.at(0).charge + reco_mu.at(1).charge; 
    reco_JPsi.pdgId = 443;
    
    if( reco_mu.at(0).v.Pt() < reco_ptMin1_mu ) continue;
    if( reco_mu.at(1).v.Pt() < reco_ptMin2_mu ) continue;
    if( reco_JPsi.v.Pt() < reco_ptMin_JPsi ) continue;
    
    bool reco_mu_genMatch = false;
    if( !isData )
      reco_mu_genMatch = IsMatching(reco_mu,gen_mu,0.01,0.1,false);
    
    ++nEvents_2RecoMu;
    if( reco_mu_genMatch ) ++nEvents_2RecoMu_genMatch;
    
    
    
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
    
    
    
    //---------------------
    //--- select reco tracks    
    if( temp_chargedHadron.size() < 3 ) continue;
    
    int selId0 = -1;
    int charge0 = 0;
    for(unsigned int ii = 0; ii < temp_chargedHadron.size(); ++ii)
    {
      float DR_mu1 = DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi());
      float DR_mu2 = DeltaR(reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi(),temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi());
      if( DR_mu1 < 0.05 ) continue;
      if( DR_mu2 < 0.05 ) continue;
      if( temp_chargedHadron.at(ii).v.Pt() < reco_ptMin1_kepi ) continue;
      if( treeVars.pfCandidates_isHighPurity->at(temp_chargedHadron.at(ii).it) != 1 ) continue;
      //if( treeVars.pfCandidates_pvAssociationQuality->at(temp_chargedHadron.at(ii).it) < 4 ) continue;
      if( treeVars.pfCandidates_fromPV->at(temp_chargedHadron.at(ii).it) < 2 ) continue;
      
      float DPhiJPsi = DeltaPhi(temp_chargedHadron.at(ii).v.Phi(),reco_JPsi.v.Phi());
      if( DPhiJPsi < reco_DPhiMin_kepi_Jpsi ) continue;

      selId0 = int(ii);
      charge0 = temp_chargedHadron.at(ii).charge;
      break;
    }
    if( selId0 >= 0 ) reco_kepi.push_back( temp_chargedHadron.at(selId0) );
    else continue;
    
    float DRMin = 999999.;
    int selId1 = -1;
    for(unsigned int ii = 0; ii < temp_chargedHadron.size(); ++ii)
    {
      if( int(ii) == selId0 ) continue;
      if( wrongDsCharge && temp_chargedHadron.at(ii).charge != charge0 ) continue;
      
      float DR_mu1 = DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi());
      float DR_mu2 = DeltaR(reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi(),temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi());
      if( DR_mu1 < 0.05 ) continue;
      if( DR_mu2 < 0.05 ) continue;
      if( temp_chargedHadron.at(ii).v.Pt() < reco_ptMin2_kepi ) continue;
      if( treeVars.pfCandidates_isHighPurity->at(temp_chargedHadron.at(ii).it) != 1 ) continue;
      //if( treeVars.pfCandidates_pvAssociationQuality->at(temp_chargedHadron.at(ii).it) < 4 ) continue;
      if( treeVars.pfCandidates_fromPV->at(temp_chargedHadron.at(ii).it) < 2 ) continue;
      float DR = DeltaR(temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi(),reco_kepi.at(0).v.Eta(),reco_kepi.at(0).v.Phi());
      if( DR > reco_DRMax_kepi ) continue;
      
      
      if( DR < DRMin )
      {
        DRMin = DR;
        selId1 = ii;
      }
    }
    if( selId1 >= 0 ) reco_kepi.push_back( temp_chargedHadron.at(selId1) );
    else continue;
    
    float DRMin2 = 999999.;
    int selId2 = -1;
    for(unsigned int ii = 1; ii < temp_chargedHadron.size(); ++ii)
    {
      if( int(ii) == selId0 ) continue;
      if( int(ii) == selId1 ) continue;
      if( wrongDsCharge && temp_chargedHadron.at(ii).charge != charge0 ) continue;
      
      float DR_mu1 = DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi());
      float DR_mu2 = DeltaR(reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi(),temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi());
      if( DR_mu1 < 0.05 ) continue;
      if( DR_mu2 < 0.05 ) continue;
      if( temp_chargedHadron.at(ii).v.Pt() < reco_ptMin3_kepi ) continue;
      if( treeVars.pfCandidates_isHighPurity->at(temp_chargedHadron.at(ii).it) != 1 ) continue;
      //if( treeVars.pfCandidates_pvAssociationQuality->at(temp_chargedHadron.at(ii).it) < 4 ) continue;
      if( treeVars.pfCandidates_fromPV->at(temp_chargedHadron.at(ii).it) < 2 ) continue;
      if( !wrongDsCharge && ( abs(reco_kepi.at(0).charge + reco_kepi.at(1).charge + temp_chargedHadron.at(ii).charge) != 1 ) ) continue;
      float DR = DeltaR(temp_chargedHadron.at(ii).v.Eta(),temp_chargedHadron.at(ii).v.Phi(),reco_kepi.at(0).v.Eta(),reco_kepi.at(0).v.Phi());
      if( DR > reco_DRMax_kepi ) continue;
      
      if( DR < DRMin2 )
      {
        DRMin2 = DR;
        selId2 = ii;
      }
    }
    if( selId2 >= 0 ) reco_kepi.push_back( temp_chargedHadron.at(selId2) );
    else continue;
    
    if( reco_kepi.size() < 3 ) continue;
    
    
    bool phiFound = false;
    bool K0StarFound = false;
    std::vector<std::vector<particle> > reco_kepi_massCases = AssignMass(reco_kepi);
    for(int ii = 0; ii < 2; ++ii)
    {
      reco_kepi_massAssign = reco_kepi_massCases.at(ii);
      
      int charge_Ds = reco_kepi_massAssign.at(0).charge + reco_kepi_massAssign.at(1).charge + reco_kepi_massAssign.at(2).charge;
      particle temp_Km;
      particle temp_Kp;
      particle temp_pi;
      for(unsigned int jj = 0; jj < reco_kepi_massAssign.size(); ++jj)
      {
        if( reco_kepi_massAssign.at(jj).charge != charge_Ds )
          temp_Km = reco_kepi_massAssign.at(jj);
        else
        {
          if( fabs((reco_kepi_massAssign.at(jj).v.M()-MASS_K)/MASS_K) < 0.01 )
            temp_Kp = reco_kepi_massAssign.at(jj);
          else
            temp_pi = reco_kepi_massAssign.at(jj);
        }
      }
      
      h1_kk_mass -> Fill( (temp_Km.v+temp_Kp.v).M() );
      h1_kpi_mass -> Fill( (temp_Km.v+temp_pi.v).M() );
      
      if( (temp_Km.v+temp_Kp.v).M() > 1.010 && (temp_Km.v+temp_Kp.v).M() < 1.030 )
      {
        phiFound = true;
        break;
      }
    }
    
    if( phiFound == false )
    {
      std::vector<std::vector<particle> > reco_kepi_massCases = AssignMass(reco_kepi);
      for(int ii = 0; ii < 2; ++ii)
      {
        reco_kepi_massAssign = reco_kepi_massCases.at(ii);
        
        int charge_Ds = reco_kepi_massAssign.at(0).charge + reco_kepi_massAssign.at(1).charge + reco_kepi_massAssign.at(2).charge;
        particle temp_Km;
        particle temp_Kp;
        particle temp_pi;
        for(unsigned int jj = 0; jj < reco_kepi_massAssign.size(); ++jj)
        {
          if( reco_kepi_massAssign.at(jj).charge != charge_Ds )
            temp_Km = reco_kepi_massAssign.at(jj);
          else
          {
            if( fabs((reco_kepi_massAssign.at(jj).v.M()-MASS_K)/MASS_K) < 0.01 )
              temp_Kp = reco_kepi_massAssign.at(jj);
            else
              temp_pi = reco_kepi_massAssign.at(jj);
          }
        }
        
        if( (temp_Km.v+temp_pi.v).M() > 0.860 && (temp_Km.v+temp_pi.v).M() < 0.920 )
        {
          K0StarFound = true;
          break;
        }
      }
    }
    
    if( useResonances && (phiFound == false && K0StarFound == false) ) continue;
    
    if( !useResonances ) reco_kepi_massAssign = AssignMassDs(reco_kepi);
    reco_Ds.v = reco_kepi_massAssign.at(0).v + reco_kepi_massAssign.at(1).v + reco_kepi_massAssign.at(2).v;
    reco_Ds.charge = reco_kepi_massAssign.at(0).charge + reco_kepi_massAssign.at(1).charge + reco_kepi_massAssign.at(2).charge;
    reco_Ds.pdgId = reco_Ds.charge*431;
    
    if( reco_Ds.v.Pt() < reco_ptMin_Ds ) continue;
    
    
    recoHistos -> FillHistos(weight,reco_mu,reco_kepi,&treeVars);
    recoMassAssignHistos -> FillHistos(weight,reco_mu,reco_kepi_massAssign,&treeVars);
    
    bool reco_kepi_genMatch = false;
    if( !isData ) reco_kepi_genMatch = IsMatching(reco_kepi,gen_kepi,0.05,0.2,false);
    
    ++nEvents_3RecoKepi;
    if( reco_kepi_genMatch ) ++nEvents_3RecoKepi_genMatch;
    
    
    
    //----------------
    //--- Ds isolation
    
    float isoChargedHadron = 0.;
    for(unsigned int ii = 0; ii < temp_chargedHadron.size(); ++ii)
    {
      particle part = temp_chargedHadron.at(ii);
      
      if( int(ii) == selId0 ) continue;
      if( int(ii) == selId1 ) continue;
      if( int(ii) == selId2 ) continue;
      
      float DR = DeltaR(reco_Ds.v.Eta(),reco_Ds.v.Phi(),part.v.Eta(),part.v.Phi());
      if( DR < 0.3 && DR > 0.0001 ) isoChargedHadron += part.v.Pt();
    }
    float isoNeutralHadron = 0.;
    for(unsigned int ii = 0; ii < temp_neutralHadron.size(); ++ii)
    {
      particle part = temp_neutralHadron.at(ii);
      
      float DR = DeltaR(reco_Ds.v.Eta(),reco_Ds.v.Phi(),part.v.Eta(),part.v.Phi());
      if( DR < 0.3 && DR > 0.01 ) isoNeutralHadron += part.v.Pt();
    }
    float isoPhoton = 0.;
    for(unsigned int ii = 0; ii < temp_photon.size(); ++ii)
    {
      particle part = temp_photon.at(ii);
      
      float DR = DeltaR(reco_Ds.v.Eta(),reco_Ds.v.Phi(),part.v.Eta(),part.v.Phi());
      if( DR < 0.3 && DR > 0.01 ) isoPhoton += part.v.Pt();
    }
    float isoPU = 0.;
    for(unsigned int ii = 0; ii < temp_PU.size(); ++ii)
    {
      particle part = temp_PU.at(ii);
      
      float DR = DeltaR(reco_Ds.v.Eta(),reco_Ds.v.Phi(),part.v.Eta(),part.v.Phi());
      if( DR < 0.3 && DR > 0.01 ) isoPU += part.v.Pt();
    }
    
    treeVars.Ds_iso = isoChargedHadron + std::max(0.,isoNeutralHadron+isoPhoton-0.5*isoPU);
    
    
    
    //---------------------------------------
    //--- cuts on J/Psi and Ds invariant mass
    
    if( reco_JPsi.v.M() < 3.00 || reco_JPsi.v.M() > 3.20 ) continue;
    
    recoMassAssignCutJPsiHistos -> FillHistos(weight,reco_mu,reco_kepi_massAssign,&treeVars);
    
    if( !antiDsIso )
      if ( treeVars.Ds_iso/reco_Ds.v.Pt() > 0.50 ) continue;
      else ;
    else
      if ( treeVars.Ds_iso/reco_Ds.v.Pt() < 0.50 ) continue;
    
    recoMassAssignCutJPsiDsIso1Histos -> FillHistos(weight,reco_mu,reco_kepi_massAssign,&treeVars);
    
    if( !antiDsIso )
      if( treeVars.Ds_iso/reco_Ds.v.Pt() > 0.25 ) continue;
      else ;
    else
      if( treeVars.Ds_iso/reco_Ds.v.Pt() < 0.25 ) continue;
    
    recoMassAssignCutJPsiDsIso2Histos -> FillHistos(weight,reco_mu,reco_kepi_massAssign,&treeVars);
    
    if( !DsSidebands )
      if( reco_Ds.v.M() < 1.90 || reco_Ds.v.M() > 2.05 ) continue;
    if( DsSidebands )
    {
      if( reco_Ds.v.M() > 1.90 && reco_Ds.v.M() < 2.05 ) continue;
      if( reco_Ds.v.M() < 1.70 ) continue;
      if( reco_Ds.v.M() > 2.30 ) continue;
    }
    
    recoMassAssignCutJPsiDsHistos -> FillHistos(weight,reco_mu,reco_kepi_massAssign,&treeVars);
    
    
    
    //---------------
    //--- print event
    if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
      std::cout << "\n\n Event: " << entry << std::endl;
    if( printGenEvent && !isData )
    {
      std::cout << "--------------------- GEN ---------------------" << std::endl;
      PrintEvent(gen_mu,gen_kepi);
    }
    if( printRecoGenMatchEvent && !isData )
    {
      std::cout << "--------------------- GEN MATCH ---------------------" << std::endl;
      PrintEvent(recoGenMatch_mu,recoGenMatch_kepi_massAssign);
    }
    if( printRecoEvent )
    {
      std::cout << "--------------------- RECO ---------------------" << std::endl;
      PrintEvent(reco_mu,reco_kepi_massAssign);
    }
    if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
      std::cout << std::endl;
    
    
    
    //-------------------
    //--- fill histograms
    ++nEvents_selected;
    
    if( !isData )
      for(int ii = 0; ii < 2; ++ii)
      {
        h1_reco_mu_gen_mu_DR -> Fill( DeltaR(reco_mu.at(ii).v.Eta(),reco_mu.at(ii).v.Phi(),gen_mu.at(ii).v.Eta(),gen_mu.at(ii).v.Phi()),weight );
        h1_reco_mu_gen_mu_ptRatio -> Fill( reco_mu.at(ii).v.Pt()/gen_mu.at(ii).v.Pt(),weight );
      }
  }
  std::cout << std::endl;
  
  
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
  h1_nEvents -> GetXaxis() -> SetBinLabel(2,"gen. cut on muons"); h1_nEvents -> SetBinContent(2,nEvents_genCutMu);
  h1_nEvents -> GetXaxis() -> SetBinLabel(3,"2 gen.-matched muons"); h1_nEvents -> SetBinContent(3,nEvents_2RecoGenMatchMu);
  h1_nEvents -> GetXaxis() -> SetBinLabel(4,"gen. cut on k/#pi"); h1_nEvents -> SetBinContent(4,nEvents_genCutKepi);
  h1_nEvents -> GetXaxis() -> SetBinLabel(5,"3 gen.-matched k/#pi"); h1_nEvents -> SetBinContent(5,nEvents_3RecoGenMatchKepi);
  h1_nEvents -> GetXaxis() -> SetBinLabel(6,"full gen. cut"); h1_nEvents -> SetBinContent(6,nEvents_genCut);
  h1_nEvents -> GetXaxis() -> SetBinLabel(7,"full gen. match"); h1_nEvents -> SetBinContent(7,nEvents_recoGenMatch);
  h1_nEvents -> GetXaxis() -> SetBinLabel(8,"2 reco muons"); h1_nEvents -> SetBinContent(8,nEvents_2RecoMu);
  h1_nEvents -> GetXaxis() -> SetBinLabel(9,"2 gen.-matched reco muons"); h1_nEvents -> SetBinContent(9,nEvents_2RecoMu_genMatch);
  h1_nEvents -> GetXaxis() -> SetBinLabel(10,"3 reco k/#pi"); h1_nEvents -> SetBinContent(10,nEvents_3RecoKepi);
  h1_nEvents -> GetXaxis() -> SetBinLabel(11,"3 gen.-matched reco k/#pi"); h1_nEvents -> SetBinContent(11,nEvents_3RecoKepi_genMatch);
  
  std::cout << std::fixed << std::endl;
  std::cout << "===================================================="  << std::endl;
  std::cout << "                              Tot. events: "    << std::setw(6) << nEvents_tot << std::endl;
  std::cout << std::endl;
  std::cout << "             Events with 2 mu in acceptance: " << std::setw(6) << nEvents_genCutMu        << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_genCutMu/nEvents_tot             << "%)" << std::endl;
  std::cout << ">>>>>>   Events with 2 reco gen. matched mu: " << std::setw(6) << nEvents_2RecoGenMatchMu << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoGenMatchMu/nEvents_genCutMu << "%)" << std::endl;
  std::cout << std::endl;
  std::cout << "           Events with 3 kepi in acceptance: " << std::setw(6) << nEvents_genCutKepi        << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_genCutKepi/nEvents_tot               << "%)" << std::endl;
  std::cout << ">>>>>> Events with 3 reco gen. matched kepi: " << std::setw(6) << nEvents_3RecoGenMatchKepi << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoGenMatchKepi/nEvents_genCutKepi << "%)" << std::endl;
  std::cout << std::endl;
  std::cout << "  Events with 2 mu and 3 kepi in acceptance: " << std::setw(6) << nEvents_genCut       << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_genCut/nEvents_tot          << "%)" << std::endl;
  std::cout << ">>>>>>          Events with full gen. match: " << std::setw(6) << nEvents_recoGenMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_recoGenMatch/nEvents_genCut << "%)" << std::endl;
  std::cout << std::endl;
  std::cout << "           Events with 2 reco. muons: " << std::setw(6) << nEvents_2RecoMu << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoMu/nEvents_tot << "%)"
            << "   + matching gen.: " << std::setw(6) << nEvents_2RecoMu_genMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoMu_genMatch/nEvents_2RecoMu << "%)" << std::endl;
  std::cout << ">>>          Events with 3 reco. tracks: " << std::setw(6) << nEvents_3RecoKepi << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoKepi/nEvents_tot << "%)"
            << "   + matching gen.: " << std::setw(6) << nEvents_3RecoKepi_genMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoKepi_genMatch/nEvents_3RecoKepi << "%)" << std::endl;
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
