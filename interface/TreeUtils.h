#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"



/*** tree variables ***/
struct TreeVars
{
  std::vector<float>* reso_pt;
  std::vector<float>* reso_eta;
  std::vector<float>* reso_phi;
  std::vector<float>* reso_energy;
  std::vector<int>* reso_charge;
  std::vector<int>* reso_pdgId;
  
  std::vector<int>* resoDau1_n;
  std::vector<std::vector<float> >* resoDau1_pt;
  std::vector<std::vector<float> >* resoDau1_eta;
  std::vector<std::vector<float> >* resoDau1_phi;
  std::vector<std::vector<float> >* resoDau1_energy;
  std::vector<std::vector<int> >* resoDau1_charge;
  std::vector<std::vector<int> >* resoDau1_pdgId;
  
  std::vector<std::vector<int> >* resoDau2_n;
  std::vector<std::vector<std::vector<float> > >* resoDau2_pt;
  std::vector<std::vector<std::vector<float> > >* resoDau2_eta;
  std::vector<std::vector<std::vector<float> > >* resoDau2_phi;
  std::vector<std::vector<std::vector<float> > >* resoDau2_energy;
  std::vector<std::vector<std::vector<int> > >* resoDau2_charge;
  std::vector<std::vector<std::vector<int> > >* resoDau2_pdgId;

  int trgs_n;
  std::vector<std::string>* trgs_name;
  std::vector<int>* trgs_pass;
  std::vector<int>* trgs_prescale;
  
  int muons_n;
  std::vector<float>* muons_pt;
  std::vector<float>* muons_eta;
  std::vector<float>* muons_phi;
  std::vector<float>* muons_energy;
  std::vector<int>* muons_charge;
  std::vector<float>* muons_dxy;
  std::vector<float>* muons_dxyErr;
  std::vector<float>* muons_dz;
  std::vector<float>* muons_dzErr;
  std::vector<float>* muons_isLoose;
  std::vector<float>* muons_isMedium;
  std::vector<float>* muons_isTight;
  std::vector<float>* muons_pfIsoChargedHadron;
  std::vector<float>* muons_pfIsoChargedParticle;
  std::vector<float>* muons_pfIsoNeutralHadron;
  std::vector<float>* muons_pfIsoPhoton;
  std::vector<float>* muons_pfIsoPU;
  std::vector<int>* muons_trackerLayersWithMeasurement;
  
  int pfCandidates_n;
  std::vector<float>* pfCandidates_pt;
  std::vector<float>* pfCandidates_eta;
  std::vector<float>* pfCandidates_phi;
  std::vector<float>* pfCandidates_energy;
  std::vector<float>* pfCandidates_dxy;
  std::vector<float>* pfCandidates_dz;
  std::vector<int>* pfCandidates_charge;
  std::vector<int>* pfCandidates_pdgId;
  std::vector<int>* pfCandidates_fromPV;
  std::vector<int>* pfCandidates_pvAssociationQuality;
  std::vector<int>* pfCandidates_isHighPurity;

  float mu1_iso;
  float mu2_iso;

  float Ds_iso;
};

void InitTreeVars(TChain* chain_reco, TChain* chain_gen, TreeVars& treeVars);

#endif
