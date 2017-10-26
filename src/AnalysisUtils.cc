#include "interface/AnalysisUtils.h"



StdHistoSet::StdHistoSet(const std::string& label, TFile* outFile):
  label_(label),
  outFile_(outFile)
{
  outFile_ -> cd();
  outFile_ -> mkdir(label_.c_str());
  outFile_ -> cd(label_.c_str());
  
  h1s_[Form("%s_h1_W_pt",  label_.c_str())] = new TH1F(Form("%s_h1_W_pt",label_.c_str()),  "",1000,0.,500.);
  h1s_[Form("%s_h1_W_eta", label_.c_str())] = new TH1F(Form("%s_h1_W_eta",label_.c_str()), "",200,0.,10.);
  h1s_[Form("%s_h1_W_mass",label_.c_str())] = new TH1F(Form("%s_h1_W_mass",label_.c_str()),"",1000,0.,500.);
  
  h1s_[Form("%s_h1_JPsi_pt",  label_.c_str())] = new TH1F(Form("%s_h1_JPsi_pt",label_.c_str()),  "",10000,0.,500.);
  h1s_[Form("%s_h1_JPsi_eta", label_.c_str())] = new TH1F(Form("%s_h1_JPsi_eta",label_.c_str()), "",200,0.,10.);
  h1s_[Form("%s_h1_JPsi_mass",label_.c_str())] = new TH1F(Form("%s_h1_JPsi_mass",label_.c_str()),"",1000,0.,5.);

  h1s_[Form("%s_h1_Ds_pt",  label_.c_str())] = new TH1F(Form("%s_h1_Ds_pt",label_.c_str()),"",10000,0.,500.);
  h1s_[Form("%s_h1_Ds_eta", label_.c_str())] = new TH1F(Form("%s_h1_Ds_eta",label_.c_str()),"",200,0.,10.);
  h1s_[Form("%s_h1_Ds_mass",label_.c_str())] = new TH1F(Form("%s_h1_Ds_mass",label_.c_str()),"",1000,0.,5.);
  
  h1s_[Form("%s_h1_mu_pt1",      label_.c_str())] = new TH1F(Form("%s_h1_mu_pt1",label_.c_str()),      "",10000,0.,500.);
  h1s_[Form("%s_h1_mu_pt2",      label_.c_str())] = new TH1F(Form("%s_h1_mu_pt2",label_.c_str()),      "",10000,0.,500.);
  h1s_[Form("%s_h1_mu_eta",      label_.c_str())] = new TH1F(Form("%s_h1_mu_eta",label_.c_str()),      "",200,0.,10.);
  h1s_[Form("%s_h1_mu_DEta",     label_.c_str())] = new TH1F(Form("%s_h1_mu_DEta",label_.c_str()),     "",1000,0.,10.);
  h1s_[Form("%s_h1_mu_DPhi",     label_.c_str())] = new TH1F(Form("%s_h1_mu_DPhi",label_.c_str()),     "",1000,0.,3.15);
  h1s_[Form("%s_h1_mu_DR",       label_.c_str())] = new TH1F(Form("%s_h1_mu_DR",label_.c_str()),       "",1000,0.,10.);
  h1s_[Form("%s_h1_mu_dxy",      label_.c_str())] = new TH1F(Form("%s_h1_mu_dxy",label_.c_str()),      "",10000,0.,1.);
  h1s_[Form("%s_h1_mu_dxyPull",  label_.c_str())] = new TH1F(Form("%s_h1_mu_dxyPull",label_.c_str()),  "",10000,0.,10.);
  h1s_[Form("%s_h1_mu_dz",       label_.c_str())] = new TH1F(Form("%s_h1_mu_dz",label_.c_str()),       "",100000,0.,10.);
  h1s_[Form("%s_h1_mu_dzPull",   label_.c_str())] = new TH1F(Form("%s_h1_mu_dzPull",label_.c_str()),   "",10000,0.,10.);
  h1s_[Form("%s_h1_mu_isLoose1", label_.c_str())] = new TH1F(Form("%s_h1_mu_isLoose1",label_.c_str()), "",2,0.,2.);
  h1s_[Form("%s_h1_mu_isLoose2", label_.c_str())] = new TH1F(Form("%s_h1_mu_isLoose2",label_.c_str()), "",2,0.,2.);
  h1s_[Form("%s_h1_mu_isMedium1",label_.c_str())] = new TH1F(Form("%s_h1_mu_isMedium1",label_.c_str()),"",2,0.,2.);
  h1s_[Form("%s_h1_mu_isMedium2",label_.c_str())] = new TH1F(Form("%s_h1_mu_isMedium2",label_.c_str()),"",2,0.,2.);
  h1s_[Form("%s_h1_mu_isTight1", label_.c_str())] = new TH1F(Form("%s_h1_mu_isTight1",label_.c_str()), "",2,0.,2.);
  h1s_[Form("%s_h1_mu_isTight2", label_.c_str())] = new TH1F(Form("%s_h1_mu_isTight2",label_.c_str()), "",2,0.,2.);
  h1s_[Form("%s_h1_mu_relIso1",  label_.c_str())] = new TH1F(Form("%s_h1_mu_relIso1",label_.c_str()),  "",10000,0.,10.);
  h1s_[Form("%s_h1_mu_relIso2",  label_.c_str())] = new TH1F(Form("%s_h1_mu_relIso2",label_.c_str()),  "",10000,0.,10.);
  
  h1s_[Form("%s_h1_kepi_pt1", label_.c_str())] = new TH1F(Form("%s_h1_kepi_pt1",label_.c_str()), "",10000,0.,500.);
  h1s_[Form("%s_h1_kepi_pt2", label_.c_str())] = new TH1F(Form("%s_h1_kepi_pt2",label_.c_str()), "",10000,0.,500.);
  h1s_[Form("%s_h1_kepi_pt3", label_.c_str())] = new TH1F(Form("%s_h1_kepi_pt3",label_.c_str()), "",10000,0.,500.);
  h1s_[Form("%s_h1_kepi_eta", label_.c_str())] = new TH1F(Form("%s_h1_kepi_eta",label_.c_str()), "",200,0.,10.);
  h1s_[Form("%s_h1_kepi_DEta",label_.c_str())] = new TH1F(Form("%s_h1_kepi_DEta",label_.c_str()),"",1000,0.,10.);
  h1s_[Form("%s_h1_kepi_DPhi",label_.c_str())] = new TH1F(Form("%s_h1_kepi_DPhi",label_.c_str()),"",1000,0.,3.15);
  h1s_[Form("%s_h1_kepi_DR",  label_.c_str())] = new TH1F(Form("%s_h1_kepi_DR",label_.c_str()),  "",1000,0.,10.);
  h1s_[Form("%s_h1_kepi_dxy", label_.c_str())] = new TH1F(Form("%s_h1_kepi_dxy",label_.c_str()), "",10000,0.,1.);
  h1s_[Form("%s_h1_kepi_dz",  label_.c_str())] = new TH1F(Form("%s_h1_kepi_dz",label_.c_str()),  "",10000,0.,1.);
  h1s_[Form("%s_h1_kepi_fromPV",label_.c_str())] = new TH1F(Form("%s_h1_kepi_fromPV",label_.c_str()),"",5,0.,5.);
  h1s_[Form("%s_h1_kepi_PVAssociationQuality",label_.c_str())] = new TH1F(Form("%s_h1_kepi_PVAssociationQuality",label_.c_str()),"",8,0.,8.);
  h1s_[Form("%s_h1_kepi_isHighPurity",label_.c_str())] = new TH1F(Form("%s_h1_kepi_isHighPurity",label_.c_str()),"",2,0.,2.);

  h1s_[Form("%s_h1_Ds_relIso",label_.c_str())] = new TH1F(Form("%s_h1_Ds_relIso",label_.c_str()),"",10000,0.,10.);
    
  h1s_[Form("%s_h1_JPsi_Ds_DEta",label_.c_str())] = new TH1F(Form("%s_h1_JPsi_Ds_DEta",label_.c_str()),"",1000,0.,10.);
  h1s_[Form("%s_h1_JPsi_Ds_DPhi",label_.c_str())] = new TH1F(Form("%s_h1_JPsi_Ds_DPhi",label_.c_str()),"",1000,0.,3.15);
  h1s_[Form("%s_h1_JPsi_Ds_DR",  label_.c_str())] = new TH1F(Form("%s_h1_JPsi_Ds_DR",label_.c_str()),  "",1000,0.,10.);
  
  h1s_[Form("%s_h1_non_mu_dxy",     label_.c_str())] = new TH1F(Form("%s_h1_non_mu_dxy",label_.c_str()),     "",10000,0.,1.);
  h1s_[Form("%s_h1_non_mu_dz",      label_.c_str())] = new TH1F(Form("%s_h1_non_mu_dz",label_.c_str()),      "",10000,0.,1.);
  h1s_[Form("%s_h1_non_mu_isLoose", label_.c_str())] = new TH1F(Form("%s_h1_non_mu_isLoose",label_.c_str()), "",2,0.,2.);
  h1s_[Form("%s_h1_non_mu_isMedium",label_.c_str())] = new TH1F(Form("%s_h1_non_mu_isMedium",label_.c_str()),"",2,0.,2.);
  h1s_[Form("%s_h1_non_mu_isTight", label_.c_str())] = new TH1F(Form("%s_h1_non_mu_isTight",label_.c_str()), "",2,0.,2.);
  
  h1s_[Form("%s_h1_non_kepi_dxy", label_.c_str())] = new TH1F(Form("%s_h1_non_kepi_dxy",label_.c_str()), "",10000,0.,1.);
  h1s_[Form("%s_h1_non_kepi_dz",  label_.c_str())] = new TH1F(Form("%s_h1_non_kepi_dz",label_.c_str()),  "",10000,0.,1.);  
  h1s_[Form("%s_h1_non_kepi_fromPV",label_.c_str())] = new TH1F(Form("%s_h1_non_kepi_fromPV",label_.c_str()),"",5,0.,5.);
  h1s_[Form("%s_h1_non_kepi_PVAssociationQuality",label_.c_str())] = new TH1F(Form("%s_h1_non_kepi_PVAssociationQuality",label_.c_str()),"",5,0.,5.);
  h1s_[Form("%s_h1_non_kepi_isHighPurity",label_.c_str())] = new TH1F(Form("%s_h1_non_kepi_isHighPurity",label_.c_str()),"",2,0.,2.);  

  for(std::map<std::string,TH1F*>::const_iterator mapIt = h1s_.begin(); mapIt != h1s_.end(); ++mapIt)
    (mapIt->second) -> Sumw2();
  
  outFile_ -> cd();
}

StdHistoSet::~StdHistoSet()
{
  for(std::map<std::string,TH1F*>::const_iterator mapIt = h1s_.begin(); mapIt != h1s_.end(); ++mapIt)
    delete mapIt -> second;
}

void StdHistoSet::FillHistos(const float& weight,
                             std::vector<particle>& mu, std::vector<particle>& kepi, TreeVars* tv)
{
  particle JPsi;
  JPsi.charge = 0;
  for(unsigned int ii = 0; ii < mu.size(); ++ii)
  {
    JPsi.v += mu.at(ii).v;
    JPsi.charge += mu.at(ii).charge;
  }
  JPsi.pdgId = 443;
  
  particle Ds;
  Ds.charge = 0;
  for(unsigned int ii = 0; ii < kepi.size(); ++ii)
  {
    Ds.v += kepi.at(ii).v;
    Ds.charge += kepi.at(ii).charge;
  }
  Ds.pdgId = Ds.charge*431;
  
  particle W;
  W.v = JPsi.v + Ds.v;
  W.charge = JPsi.charge + Ds.charge;
  W.pdgId = W.charge*21;
  
  h1s_[Form("%s_h1_W_pt",label_.c_str())] -> Fill( W.v.Pt(),weight );
  h1s_[Form("%s_h1_W_eta",label_.c_str())] -> Fill( W.v.Eta(),weight );
  h1s_[Form("%s_h1_W_mass",label_.c_str())] -> Fill( W.v.M(),weight );
  
  h1s_[Form("%s_h1_JPsi_pt",label_.c_str())] -> Fill( JPsi.v.Pt(),weight );
  h1s_[Form("%s_h1_JPsi_eta",label_.c_str())] -> Fill( JPsi.v.Eta(),weight );
  h1s_[Form("%s_h1_JPsi_mass",label_.c_str())] -> Fill( JPsi.v.M(),weight );
  
  h1s_[Form("%s_h1_Ds_pt",label_.c_str())] -> Fill( Ds.v.Pt(),weight );
  h1s_[Form("%s_h1_Ds_eta",label_.c_str())] -> Fill( Ds.v.Eta(),weight );
  h1s_[Form("%s_h1_Ds_mass",label_.c_str())] -> Fill( Ds.v.M(),weight );
  
  h1s_[Form("%s_h1_mu_pt1", label_.c_str())] -> Fill( mu.at(0).v.Pt(),weight );
  h1s_[Form("%s_h1_mu_pt2", label_.c_str())] -> Fill( mu.at(1).v.Pt(),weight );
  h1s_[Form("%s_h1_mu_eta", label_.c_str())] -> Fill( mu.at(0).v.Eta(),weight );
  h1s_[Form("%s_h1_mu_eta", label_.c_str())] -> Fill( mu.at(1).v.Eta(),weight );
  h1s_[Form("%s_h1_mu_DEta",label_.c_str())] -> Fill( DeltaEta(mu.at(0).v.Eta(),mu.at(1).v.Eta()),weight );
  h1s_[Form("%s_h1_mu_DPhi",label_.c_str())] -> Fill( DeltaPhi(mu.at(0).v.Phi(),mu.at(1).v.Phi()),weight );
  h1s_[Form("%s_h1_mu_DR",  label_.c_str())] -> Fill( DeltaR(mu.at(0).v.Eta(),mu.at(0).v.Phi(),mu.at(1).v.Eta(),mu.at(1).v.Phi()),weight );
  
  h1s_[Form("%s_h1_kepi_pt1", label_.c_str())] -> Fill( kepi.at(0).v.Pt(),weight );
  h1s_[Form("%s_h1_kepi_eta", label_.c_str())] -> Fill( kepi.at(0).v.Eta(),weight );
  if( kepi.size() > 1 )
  {
    h1s_[Form("%s_h1_kepi_pt2", label_.c_str())] -> Fill( kepi.at(1).v.Pt(),weight );
    h1s_[Form("%s_h1_kepi_eta", label_.c_str())] -> Fill( kepi.at(1).v.Eta(),weight );
    h1s_[Form("%s_h1_kepi_DEta",label_.c_str())] -> Fill( DeltaEta(kepi.at(0).v.Eta(),kepi.at(1).v.Eta()),weight );
    h1s_[Form("%s_h1_kepi_DPhi",label_.c_str())] -> Fill( DeltaPhi(kepi.at(0).v.Phi(),kepi.at(1).v.Phi()),weight );
    h1s_[Form("%s_h1_kepi_DR",  label_.c_str())] -> Fill( DeltaR(kepi.at(0).v.Eta(),kepi.at(0).v.Phi(),kepi.at(1).v.Eta(),kepi.at(1).v.Phi()),weight );
  }
  if( kepi.size() > 2 )
  {
    h1s_[Form("%s_h1_kepi_pt3", label_.c_str())] -> Fill( kepi.at(2).v.Pt(),weight );
    h1s_[Form("%s_h1_kepi_eta", label_.c_str())] -> Fill( kepi.at(2).v.Eta(),weight );  
    h1s_[Form("%s_h1_kepi_DEta",label_.c_str())] -> Fill( DeltaEta(kepi.at(0).v.Eta(),kepi.at(2).v.Eta()),weight );
    h1s_[Form("%s_h1_kepi_DPhi",label_.c_str())] -> Fill( DeltaPhi(kepi.at(0).v.Phi(),kepi.at(2).v.Phi()),weight );
    h1s_[Form("%s_h1_kepi_DR",  label_.c_str())] -> Fill( DeltaR(kepi.at(0).v.Eta(),kepi.at(0).v.Phi(),kepi.at(2).v.Eta(),kepi.at(2).v.Phi()),weight );
    h1s_[Form("%s_h1_kepi_DEta",label_.c_str())] -> Fill( DeltaEta(kepi.at(1).v.Eta(),kepi.at(2).v.Eta()),weight );
    h1s_[Form("%s_h1_kepi_DPhi",label_.c_str())] -> Fill( DeltaPhi(kepi.at(1).v.Phi(),kepi.at(2).v.Phi()),weight );
    h1s_[Form("%s_h1_kepi_DR",  label_.c_str())] -> Fill( DeltaR(kepi.at(1).v.Eta(),kepi.at(1).v.Phi(),kepi.at(2).v.Eta(),kepi.at(2).v.Phi()),weight );
  }
  h1s_[Form("%s_h1_JPsi_Ds_DEta",label_.c_str())] -> Fill( DeltaEta(JPsi.v.Eta(),Ds.v.Eta()),weight );
  h1s_[Form("%s_h1_JPsi_Ds_DPhi",label_.c_str())] -> Fill( DeltaPhi(JPsi.v.Phi(),Ds.v.Phi()),weight );
  h1s_[Form("%s_h1_JPsi_Ds_DR",  label_.c_str())] -> Fill( DeltaR(JPsi.v.Eta(),JPsi.v.Phi(),Ds.v.Eta(),Ds.v.Phi()),weight );
  
  if( tv )
  {
    h1s_[Form("%s_h1_mu_dxy",      label_.c_str())] -> Fill( fabs(tv->muons_dxy->at(mu.at(0).it)),weight );
    h1s_[Form("%s_h1_mu_dxy",      label_.c_str())] -> Fill( fabs(tv->muons_dxy->at(mu.at(1).it)),weight );
    h1s_[Form("%s_h1_mu_dxyPull",  label_.c_str())] -> Fill( fabs(tv->muons_dxy->at(mu.at(0).it)/tv->muons_dxyErr->at(mu.at(0).it)),weight );
    h1s_[Form("%s_h1_mu_dxyPull",  label_.c_str())] -> Fill( fabs(tv->muons_dxy->at(mu.at(1).it)/tv->muons_dxyErr->at(mu.at(1).it)),weight );
    h1s_[Form("%s_h1_mu_dz",       label_.c_str())] -> Fill( fabs(tv->muons_dz->at(mu.at(0).it)),weight );
    h1s_[Form("%s_h1_mu_dz",       label_.c_str())] -> Fill( fabs(tv->muons_dz->at(mu.at(1).it)),weight );
    h1s_[Form("%s_h1_mu_dzPull",   label_.c_str())] -> Fill( fabs(tv->muons_dz->at(mu.at(0).it)/tv->muons_dzErr->at(mu.at(0).it)),weight );
    h1s_[Form("%s_h1_mu_dzPull",   label_.c_str())] -> Fill( fabs(tv->muons_dz->at(mu.at(1).it)/tv->muons_dzErr->at(mu.at(1).it)),weight );
    h1s_[Form("%s_h1_mu_isLoose1", label_.c_str())] -> Fill( tv->muons_isLoose->at(mu.at(0).it),weight );
    h1s_[Form("%s_h1_mu_isLoose2", label_.c_str())] -> Fill( tv->muons_isLoose->at(mu.at(1).it),weight );
    h1s_[Form("%s_h1_mu_isMedium1",label_.c_str())] -> Fill( tv->muons_isMedium->at(mu.at(0).it),weight );
    h1s_[Form("%s_h1_mu_isMedium2",label_.c_str())] -> Fill( tv->muons_isMedium->at(mu.at(1).it),weight );
    h1s_[Form("%s_h1_mu_isTight1", label_.c_str())] -> Fill( tv->muons_isTight->at(mu.at(0).it),weight );
    h1s_[Form("%s_h1_mu_isTight2", label_.c_str())] -> Fill( tv->muons_isTight->at(mu.at(1).it),weight );
    
    h1s_[Form("%s_h1_mu_relIso1",  label_.c_str())] -> Fill( tv->mu1_iso/mu.at(0).v.Pt(),weight );
    h1s_[Form("%s_h1_mu_relIso2",  label_.c_str())] -> Fill( tv->mu2_iso/mu.at(1).v.Pt(),weight );
    
    h1s_[Form("%s_h1_kepi_dxy",label_.c_str())] -> Fill( tv->pfCandidates_dxy->at(kepi.at(0).it),weight );
    h1s_[Form("%s_h1_kepi_dxy",label_.c_str())] -> Fill( tv->pfCandidates_dxy->at(kepi.at(1).it),weight );
    h1s_[Form("%s_h1_kepi_dxy",label_.c_str())] -> Fill( tv->pfCandidates_dxy->at(kepi.at(2).it),weight );
    h1s_[Form("%s_h1_kepi_dz", label_.c_str())] -> Fill( tv->pfCandidates_dz->at(kepi.at(0).it),weight );
    h1s_[Form("%s_h1_kepi_dz", label_.c_str())] -> Fill( tv->pfCandidates_dz->at(kepi.at(1).it),weight );
    h1s_[Form("%s_h1_kepi_dz", label_.c_str())] -> Fill( tv->pfCandidates_dz->at(kepi.at(2).it),weight );
    h1s_[Form("%s_h1_kepi_fromPV", label_.c_str())] -> Fill( tv->pfCandidates_fromPV->at(kepi.at(0).it),weight );
    h1s_[Form("%s_h1_kepi_fromPV", label_.c_str())] -> Fill( tv->pfCandidates_fromPV->at(kepi.at(1).it),weight );
    h1s_[Form("%s_h1_kepi_fromPV", label_.c_str())] -> Fill( tv->pfCandidates_fromPV->at(kepi.at(2).it),weight );
    h1s_[Form("%s_h1_kepi_PVAssociationQuality", label_.c_str())] -> Fill( tv->pfCandidates_pvAssociationQuality->at(kepi.at(0).it),weight );
    h1s_[Form("%s_h1_kepi_PVAssociationQuality", label_.c_str())] -> Fill( tv->pfCandidates_pvAssociationQuality->at(kepi.at(1).it),weight );
    h1s_[Form("%s_h1_kepi_PVAssociationQuality", label_.c_str())] -> Fill( tv->pfCandidates_pvAssociationQuality->at(kepi.at(2).it),weight );
    h1s_[Form("%s_h1_kepi_isHighPurity", label_.c_str())] -> Fill( tv->pfCandidates_isHighPurity->at(kepi.at(0).it),weight );
    h1s_[Form("%s_h1_kepi_isHighPurity", label_.c_str())] -> Fill( tv->pfCandidates_isHighPurity->at(kepi.at(1).it),weight );
    h1s_[Form("%s_h1_kepi_isHighPurity", label_.c_str())] -> Fill( tv->pfCandidates_isHighPurity->at(kepi.at(2).it),weight );
    
    h1s_[Form("%s_h1_Ds_relIso",label_.c_str())] -> Fill( tv->Ds_iso/Ds.v.Pt(),weight );
    
    for(int ii = 0; ii < tv->muons_n; ++ii)
    {
      if( ii == mu.at(0).it ) continue;
      if( ii == mu.at(1).it ) continue;
      
      h1s_[Form("%s_h1_non_mu_dxy",     label_.c_str())] -> Fill( tv->muons_dxy->at(ii),weight );
      h1s_[Form("%s_h1_non_mu_dz",      label_.c_str())] -> Fill( tv->muons_dz->at(ii),weight );
      h1s_[Form("%s_h1_non_mu_isLoose", label_.c_str())] -> Fill( tv->muons_isLoose->at(ii),weight );
      h1s_[Form("%s_h1_non_mu_isMedium",label_.c_str())] -> Fill( tv->muons_isMedium->at(ii),weight );
      h1s_[Form("%s_h1_non_mu_isTight", label_.c_str())] -> Fill( tv->muons_isTight->at(ii),weight );
    }

    for(int ii = 0; ii < tv->pfCandidates_n; ++ii)
    {
      if( ii == kepi.at(0).it ) continue;
      if( ii == kepi.at(1).it ) continue;
      if( ii == kepi.at(2).it ) continue;
      
      h1s_[Form("%s_h1_non_kepi_dxy",         label_.c_str())] -> Fill( tv->pfCandidates_dxy->at(ii),weight );
      h1s_[Form("%s_h1_non_kepi_dz",          label_.c_str())] -> Fill( tv->pfCandidates_dz->at(ii),weight );
      h1s_[Form("%s_h1_non_kepi_fromPV",      label_.c_str())] -> Fill( tv->pfCandidates_pvAssociationQuality->at(ii),weight );
      h1s_[Form("%s_h1_non_kepi_PVAssociationQuality",      label_.c_str())] -> Fill( tv->pfCandidates_pvAssociationQuality->at(ii),weight );
      h1s_[Form("%s_h1_non_kepi_isHighPurity",label_.c_str())] -> Fill( tv->pfCandidates_isHighPurity->at(ii),weight );      
    }
  }
}

float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) + 
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}



bool IsMatching(const particle& p1, const particle p2,
                const float& DRMax, const float& ptRatioMax)
{
  float DR = DeltaR(p1.v.Eta(),p1.v.Phi(),p2.v.Eta(),p2.v.Phi());

  bool isMatching = true;
  if( DR > DRMax ) isMatching = false;
  if( fabs(p1.v.Pt()/p2.v.Pt()-1.) > ptRatioMax ) isMatching = false;  

  return isMatching;
}

bool IsMatching(const std::vector<particle>& vp1, const std::vector<particle>& vp2,
                const float& DRMax, const float& ptRatioMax,
                const bool& verbosity)
{
  if( vp1.size() != vp2.size() )
  {
    std::cerr << "!!! IsMatching()::ERROR: the two vectors should have the same size !!!" << std::endl;
    return false;
  }
  
  std::vector<int> trial;
  std::vector<int> bestMatch;
  for(unsigned int ii = 0; ii < vp1.size(); ++ii)
  {
    trial.push_back( ii );
    bestMatch.push_back(-1);
  }
  
  float DRSumMin = 999999.;
  
  do {
    float DRSum = 0.;
    
    if( verbosity ) std::cout << "trying ";
    for(unsigned int ii = 0; ii < vp2.size(); ++ii)
    {
      float DR = DeltaR(vp1.at(trial.at(ii)).v.Eta(),vp1.at(trial.at(ii)).v.Phi(),vp2.at(ii).v.Eta(),vp2.at(ii).v.Phi());
      DRSum += DR;
      if( verbosity ) std::cout << trial.at(ii) << "-" << ii << " (" << DR << ")   ";
    }
    if( verbosity ) std::cout << std::endl;

    if( DRSum < DRSumMin )
    {
      DRSumMin = DRSum;
      bestMatch = trial;
    }
  } while ( std::next_permutation(trial.begin(),trial.end()) );


  bool isMatching = true;
  for(unsigned int ii = 0; ii < bestMatch.size(); ++ii)
  {
    float DR = DeltaR(vp1.at(bestMatch.at(ii)).v.Eta(),vp1.at(bestMatch.at(ii)).v.Phi(),vp2.at(ii).v.Eta(),vp2.at(ii).v.Phi());
    if( DR > DRMax ) isMatching = false;
    if( fabs(vp1.at(bestMatch.at(ii)).v.Pt()/vp2.at(ii).v.Pt()-1.) > ptRatioMax ) isMatching = false;
  }

  if( verbosity ) std::cout << (isMatching ? "MATCHING" : "NOT MATCHING") << std::endl;
  return isMatching;
}



int GetBestMatch(const particle& p, std::vector<particle>& vec, std::vector<int>* vetoVec)
{
  int bestMatch = -1;
  float DRMin = 999999.;
  
  for(unsigned int ii = 0; ii < vec.size(); ++ii)
  {
    bool skip = false;
    if( vetoVec)
      for(unsigned int jj = 0; jj < vetoVec->size(); ++jj)
        if( int(ii) == vetoVec->at(jj) ) skip = true;
    if( skip ) continue;
    
    float DR = DeltaR(vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),p.v.Eta(),p.v.Phi());
    if( DR < DRMin )
    {
      DRMin = DR;
      bestMatch = int(ii);
    }
  }

  return bestMatch;
}



std::vector<particle> AssignMassDs(std::vector<particle>& vec)
{
  int charge_Ds = vec.at(0).charge + vec.at(1).charge + vec.at(2).charge;
  
  std::vector<particle> temp1;
  std::vector<particle> temp2;

  int it = 0;
  for(int ii = 0; ii < 3; ++ii)
  {
    temp1.push_back(vec.at(ii));
    temp2.push_back(vec.at(ii));
    
    if( vec.at(ii).charge != charge_Ds )
    {
      temp1.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
      temp2.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
    }
    else
    {
      if( it == 0 )
      {
        temp1.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
        temp2.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_PI);
        ++it;
      }
      else
      {
        temp1.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_PI);
        temp2.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
      }
    }
  }

  float mass1 = (temp1.at(0).v + temp1.at(1).v + temp1.at(2).v).M();
  float mass2 = (temp2.at(0).v + temp2.at(1).v + temp2.at(2).v).M();

  std::vector<particle> ret;
  ret = fabs(mass1-MASS_DS) < fabs(mass2-MASS_DS) ? temp1 : temp2;
  return ret;
}



std::vector<std::vector<particle> > AssignMass(std::vector<particle>& vec)
{
  int charge_Ds = vec.at(0).charge + vec.at(1).charge + vec.at(2).charge;
  
  std::vector<particle> temp1;
  std::vector<particle> temp2;
  
  int it = 0;
  for(int ii = 0; ii < 3; ++ii)
  {
    temp1.push_back(vec.at(ii));
    temp2.push_back(vec.at(ii));
    
    if( vec.at(ii).charge != charge_Ds )
    {
      temp1.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
      temp2.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
    }
    else
    {
      if( it == 0 )
      {
        temp1.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
        temp2.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_PI);
        ++it;
      }
      else
      {
        temp1.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_PI);
        temp2.at(ii).v.SetPtEtaPhiM(vec.at(ii).v.Pt(),vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),MASS_K);
      }
    }
  }
  
  std::vector<std::vector<particle> > ret;
  ret.push_back(temp1);
  ret.push_back(temp2);
  return ret;
}



void PrintEvent(std::vector<particle>& mu, std::vector<particle>& kepi)
{
  particle JPsi;
  JPsi.charge = 0;
  for(unsigned int ii = 0; ii < mu.size(); ++ii)
  {
    JPsi.v += mu.at(ii).v;
    JPsi.charge += mu.at(ii).charge;
  }
  JPsi.pdgId = 443;
  
  particle Ds;
  Ds.charge = 0;
  for(unsigned int ii = 0; ii < kepi.size(); ++ii)
  {
    Ds.v += kepi.at(ii).v;
    Ds.charge += kepi.at(ii).charge;
  }
  Ds.pdgId = Ds.charge*431;
  
  particle W;
  W.v = JPsi.v + Ds.v;
  W.charge = JPsi.charge + Ds.charge;
  W.pdgId = W.charge*21;
  
  std::cout << "W: " << W << std::endl;
  std::cout << ">>> J/Psi: " << JPsi << std::endl;
  for(unsigned int ii = 0; ii < mu.size(); ++ii)
    std::cout << ">>>>>>   mu: " << mu.at(ii) << std::endl;
  std::cout << ">>>   D_s: " << Ds << std::endl;
  for(unsigned int ii = 0; ii < kepi.size(); ++ii)
    std::cout << ">>>>>> kepi: " << kepi.at(ii) << std::endl;
}

std::ostream& operator<<(std::ostream& os, const TLorentzVector& v)
{
  os << "(pt="    << std::fixed << std::setw(6) << std::setprecision(2) << v.Pt()
     << ", eta="  << std::fixed << std::setw(5) << std::setprecision(2) << v.Eta()
     << ", phi="  << std::fixed << std::setw(5) << std::setprecision(2) << v.Phi()
     << ", mass=" << std::fixed << std::setw(7) << std::setprecision(3) << v.M()
     << ")";

  return os;
}

std::ostream& operator<<(std::ostream& os, const particle& p)
{
  os << p.v
     << "   charge: " << std::fixed << std::setw(2) << p.charge
     << "   pdgId: "  << std::fixed << std::setw(6) << p.pdgId
     << " (" << GetParticleName(p.pdgId) << ")";
  
  return os;
}
