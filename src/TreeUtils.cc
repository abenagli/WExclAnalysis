#include "interface/TreeUtils.h"



void InitTreeVars(TChain* chain_reco, TChain* chain_gen, TreeVars& treeVars)
{
  treeVars.reso_pt = new std::vector<float>;
  treeVars.reso_eta = new std::vector<float>;
  treeVars.reso_phi = new std::vector<float>;
  treeVars.reso_energy = new std::vector<float>;
  treeVars.reso_charge = new std::vector<int>;
  treeVars.reso_pdgId = new std::vector<int>;
  
  treeVars.resoDau1_n = new std::vector<int>;
  treeVars.resoDau1_pt = new std::vector<std::vector<float> >;
  treeVars.resoDau1_eta = new std::vector<std::vector<float> >;
  treeVars.resoDau1_phi = new std::vector<std::vector<float> >;
  treeVars.resoDau1_energy = new std::vector<std::vector<float> >;
  treeVars.resoDau1_charge = new std::vector<std::vector<int> >;
  treeVars.resoDau1_pdgId = new std::vector<std::vector<int> >;
  
  treeVars.resoDau2_n = new std::vector<std::vector<int> >;
  treeVars.resoDau2_pt = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_eta = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_phi = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_energy = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_charge = new std::vector<std::vector<std::vector<int> > >;  
  treeVars.resoDau2_pdgId = new std::vector<std::vector<std::vector<int> > >;

  treeVars.trgs_name = new std::vector<std::string>;
  treeVars.trgs_pass = new std::vector<int>;
  treeVars.trgs_prescale = new std::vector<int>;
  
  treeVars.muons_pt = new std::vector<float>;
  treeVars.muons_eta = new std::vector<float>;
  treeVars.muons_phi = new std::vector<float>;
  treeVars.muons_energy = new std::vector<float>;
  treeVars.muons_charge = new std::vector<int>;
  treeVars.muons_dxy = new std::vector<float>;
  treeVars.muons_dxyErr = new std::vector<float>;
  treeVars.muons_dz = new std::vector<float>;
  treeVars.muons_dzErr = new std::vector<float>;
  treeVars.muons_isLoose = new std::vector<float>;
  treeVars.muons_isMedium = new std::vector<float>;
  treeVars.muons_isTight = new std::vector<float>;
  treeVars.muons_pfIsoChargedHadron = new std::vector<float>;
  treeVars.muons_pfIsoChargedParticle = new std::vector<float>;
  treeVars.muons_pfIsoNeutralHadron = new std::vector<float>;
  treeVars.muons_pfIsoPhoton = new std::vector<float>;
  treeVars.muons_pfIsoPU = new std::vector<float>;
  treeVars.muons_trackerLayersWithMeasurement = new std::vector<int>;
  
  treeVars.pfCandidates_pt = new std::vector<float>;
  treeVars.pfCandidates_eta = new std::vector<float>;
  treeVars.pfCandidates_phi = new std::vector<float>;
  treeVars.pfCandidates_energy = new std::vector<float>;
  treeVars.pfCandidates_dxy = new std::vector<float>;
  treeVars.pfCandidates_dz = new std::vector<float>;
  treeVars.pfCandidates_charge = new std::vector<int>;
  treeVars.pfCandidates_pdgId = new std::vector<int>;
  treeVars.pfCandidates_fromPV = new std::vector<int>;
  treeVars.pfCandidates_pvAssociationQuality = new std::vector<int>;
  treeVars.pfCandidates_isHighPurity = new std::vector<int>;
  
  //tree -> SetBranchStatus("*",0);
  
  if( chain_gen != NULL )
  {
    chain_gen -> SetBranchStatus("reso_pt",    1); chain_gen -> SetBranchAddress("reso_pt",    &treeVars.reso_pt);
    chain_gen -> SetBranchStatus("reso_eta",   1); chain_gen -> SetBranchAddress("reso_eta",   &treeVars.reso_eta);
    chain_gen -> SetBranchStatus("reso_phi",   1); chain_gen -> SetBranchAddress("reso_phi",   &treeVars.reso_phi);
    chain_gen -> SetBranchStatus("reso_energy",1); chain_gen -> SetBranchAddress("reso_energy",&treeVars.reso_energy);
    chain_gen -> SetBranchStatus("reso_charge",1); chain_gen -> SetBranchAddress("reso_charge",&treeVars.reso_charge);
    chain_gen -> SetBranchStatus("reso_pdgId", 1); chain_gen -> SetBranchAddress("reso_pdgId", &treeVars.reso_pdgId);
    
    chain_gen -> SetBranchStatus("resoDau1_n",     1); chain_gen -> SetBranchAddress("resoDau1_n",     &treeVars.resoDau1_n);
    chain_gen -> SetBranchStatus("resoDau1_pt",    1); chain_gen -> SetBranchAddress("resoDau1_pt",    &treeVars.resoDau1_pt);
    chain_gen -> SetBranchStatus("resoDau1_eta",   1); chain_gen -> SetBranchAddress("resoDau1_eta",   &treeVars.resoDau1_eta);
    chain_gen -> SetBranchStatus("resoDau1_phi",   1); chain_gen -> SetBranchAddress("resoDau1_phi",   &treeVars.resoDau1_phi);
    chain_gen -> SetBranchStatus("resoDau1_energy",1); chain_gen -> SetBranchAddress("resoDau1_energy",&treeVars.resoDau1_energy);
    chain_gen -> SetBranchStatus("resoDau1_charge",1); chain_gen -> SetBranchAddress("resoDau1_charge",&treeVars.resoDau1_charge);
    chain_gen -> SetBranchStatus("resoDau1_pdgId", 1); chain_gen -> SetBranchAddress("resoDau1_pdgId", &treeVars.resoDau1_pdgId);
    
    chain_gen -> SetBranchStatus("resoDau2_n",     1); chain_gen -> SetBranchAddress("resoDau2_n",     &treeVars.resoDau2_n);
    chain_gen -> SetBranchStatus("resoDau2_pt",    1); chain_gen -> SetBranchAddress("resoDau2_pt",    &treeVars.resoDau2_pt);
    chain_gen -> SetBranchStatus("resoDau2_eta",   1); chain_gen -> SetBranchAddress("resoDau2_eta",   &treeVars.resoDau2_eta);
    chain_gen -> SetBranchStatus("resoDau2_phi",   1); chain_gen -> SetBranchAddress("resoDau2_phi",   &treeVars.resoDau2_phi);
    chain_gen -> SetBranchStatus("resoDau2_energy",1); chain_gen -> SetBranchAddress("resoDau2_energy",&treeVars.resoDau2_energy);
    chain_gen -> SetBranchStatus("resoDau2_charge",1); chain_gen -> SetBranchAddress("resoDau2_charge",&treeVars.resoDau2_charge);
    chain_gen -> SetBranchStatus("resoDau2_pdgId", 1); chain_gen -> SetBranchAddress("resoDau2_pdgId", &treeVars.resoDau2_pdgId);
  }

  if( chain_reco != NULL )
  {
    chain_reco -> SetBranchStatus("trgs_n",       1); chain_reco -> SetBranchAddress("trgs_n",       &treeVars.trgs_n);
    chain_reco -> SetBranchStatus("trgs_name",    1); chain_reco -> SetBranchAddress("trgs_name",    &treeVars.trgs_name);
    chain_reco -> SetBranchStatus("trgs_pass",    1); chain_reco -> SetBranchAddress("trgs_pass",    &treeVars.trgs_pass);
    chain_reco -> SetBranchStatus("trgs_prescale",1); chain_reco -> SetBranchAddress("trgs_prescale",&treeVars.trgs_prescale);
    
    chain_reco -> SetBranchStatus("muons_n",       1); chain_reco -> SetBranchAddress("muons_n",       &treeVars.muons_n);
    chain_reco -> SetBranchStatus("muons_pt",      1); chain_reco -> SetBranchAddress("muons_pt",      &treeVars.muons_pt);
    chain_reco -> SetBranchStatus("muons_eta",     1); chain_reco -> SetBranchAddress("muons_eta",     &treeVars.muons_eta);
    chain_reco -> SetBranchStatus("muons_phi",     1); chain_reco -> SetBranchAddress("muons_phi",     &treeVars.muons_phi);
    chain_reco -> SetBranchStatus("muons_energy",  1); chain_reco -> SetBranchAddress("muons_energy",  &treeVars.muons_energy);
    chain_reco -> SetBranchStatus("muons_charge",  1); chain_reco -> SetBranchAddress("muons_charge",  &treeVars.muons_charge);
    chain_reco -> SetBranchStatus("muons_dxy",     1); chain_reco -> SetBranchAddress("muons_dxy",     &treeVars.muons_dxy);
    chain_reco -> SetBranchStatus("muons_dxyErr",  1); chain_reco -> SetBranchAddress("muons_dxyErr",  &treeVars.muons_dxyErr);
    chain_reco -> SetBranchStatus("muons_dz",      1); chain_reco -> SetBranchAddress("muons_dz",      &treeVars.muons_dz);
    chain_reco -> SetBranchStatus("muons_dzErr",   1); chain_reco -> SetBranchAddress("muons_dzErr",   &treeVars.muons_dzErr);
    chain_reco -> SetBranchStatus("muons_isLoose", 1); chain_reco -> SetBranchAddress("muons_isLoose", &treeVars.muons_isLoose);
    chain_reco -> SetBranchStatus("muons_isMedium",1); chain_reco -> SetBranchAddress("muons_isMedium",&treeVars.muons_isMedium);
    chain_reco -> SetBranchStatus("muons_isTight", 1); chain_reco -> SetBranchAddress("muons_isTight", &treeVars.muons_isTight);
    chain_reco -> SetBranchStatus("muons_pfIsoChargedHadron",  1); chain_reco -> SetBranchAddress("muons_pfIsoChargedHadron",  &treeVars.muons_pfIsoChargedHadron);
    chain_reco -> SetBranchStatus("muons_pfIsoChargedParticle",1); chain_reco -> SetBranchAddress("muons_pfIsoChargedParticle",&treeVars.muons_pfIsoChargedParticle);
    chain_reco -> SetBranchStatus("muons_pfIsoNeutralHadron",  1); chain_reco -> SetBranchAddress("muons_pfIsoNeutralHadron",  &treeVars.muons_pfIsoNeutralHadron);
    chain_reco -> SetBranchStatus("muons_pfIsoPhoton",         1); chain_reco -> SetBranchAddress("muons_pfIsoPhoton",         &treeVars.muons_pfIsoPhoton);
    chain_reco -> SetBranchStatus("muons_pfIsoPU",             1); chain_reco -> SetBranchAddress("muons_pfIsoPU",             &treeVars.muons_pfIsoPU);
    chain_reco -> SetBranchStatus("muons_trackerLayersWithMeasurement",1); chain_reco -> SetBranchAddress("muons_trackerLayersWithMeasurement",&treeVars.muons_trackerLayersWithMeasurement);
    
    chain_reco -> SetBranchStatus("pfCandidates_n",     1); chain_reco -> SetBranchAddress("pfCandidates_n",     &treeVars.pfCandidates_n);
    chain_reco -> SetBranchStatus("pfCandidates_pt",    1); chain_reco -> SetBranchAddress("pfCandidates_pt",    &treeVars.pfCandidates_pt);
    chain_reco -> SetBranchStatus("pfCandidates_eta",   1); chain_reco -> SetBranchAddress("pfCandidates_eta",   &treeVars.pfCandidates_eta);
    chain_reco -> SetBranchStatus("pfCandidates_phi",   1); chain_reco -> SetBranchAddress("pfCandidates_phi",   &treeVars.pfCandidates_phi);
    chain_reco -> SetBranchStatus("pfCandidates_energy",1); chain_reco -> SetBranchAddress("pfCandidates_energy",&treeVars.pfCandidates_energy);
    chain_reco -> SetBranchStatus("pfCandidates_charge",1); chain_reco -> SetBranchAddress("pfCandidates_charge",&treeVars.pfCandidates_charge);
    chain_reco -> SetBranchStatus("pfCandidates_dxy",   1); chain_reco -> SetBranchAddress("pfCandidates_dxy",   &treeVars.pfCandidates_dxy);
    chain_reco -> SetBranchStatus("pfCandidates_dz",    1); chain_reco -> SetBranchAddress("pfCandidates_dz",    &treeVars.pfCandidates_dz);
    chain_reco -> SetBranchStatus("pfCandidates_pdgId", 1); chain_reco -> SetBranchAddress("pfCandidates_pdgId", &treeVars.pfCandidates_pdgId);
    chain_reco -> SetBranchStatus("pfCandidates_fromPV",1); chain_reco -> SetBranchAddress("pfCandidates_fromPV", &treeVars.pfCandidates_fromPV);
    chain_reco -> SetBranchStatus("pfCandidates_pvAssociationQuality",1); chain_reco -> SetBranchAddress("pfCandidates_pvAssociationQuality", &treeVars.pfCandidates_pvAssociationQuality);
    chain_reco -> SetBranchStatus("pfCandidates_isHighPurity",1); chain_reco -> SetBranchAddress("pfCandidates_isHighPurity", &treeVars.pfCandidates_isHighPurity);
  }
}
