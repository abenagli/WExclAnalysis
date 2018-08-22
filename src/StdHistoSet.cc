#include "interface/StdHistoSet.h"



StdHistoSet::StdHistoSet(const std::string& label, TFile* outFile):
  label_(label),
  outFile_(outFile)
{
  outFile_ -> cd();
  outFile_ -> mkdir(label_.c_str());
  outFile_ -> cd(label_.c_str());
  
  h1s_[Form("%s_h1_H_pt",  label_.c_str())] = new TH1F(Form("%s_h1_H_pt",label_.c_str()),  "",1000,0.,500.);
  h1s_[Form("%s_h1_H_eta", label_.c_str())] = new TH1F(Form("%s_h1_H_eta",label_.c_str()), "", 200,0., 10.);
  h1s_[Form("%s_h1_H_mass",label_.c_str())] = new TH1F(Form("%s_h1_H_mass",label_.c_str()),"",1000,0.,500.);
  
  h1s_[Form("%s_h1_mu_pt1",      label_.c_str())] = new TH1F(Form("%s_h1_mu_pt1",label_.c_str()),      "",10000,0.,500.00);
  h1s_[Form("%s_h1_mu_pt2",      label_.c_str())] = new TH1F(Form("%s_h1_mu_pt2",label_.c_str()),      "",10000,0.,500.00);
  h1s_[Form("%s_h1_mu_eta",      label_.c_str())] = new TH1F(Form("%s_h1_mu_eta",label_.c_str()),      "",  200,0., 10.00);
  h1s_[Form("%s_h1_mu_DEta",     label_.c_str())] = new TH1F(Form("%s_h1_mu_DEta",label_.c_str()),     "", 1000,0., 10.00);
  h1s_[Form("%s_h1_mu_DPhi",     label_.c_str())] = new TH1F(Form("%s_h1_mu_DPhi",label_.c_str()),     "", 1000,0.,  3.15);
  h1s_[Form("%s_h1_mu_DR",       label_.c_str())] = new TH1F(Form("%s_h1_mu_DR",label_.c_str()),       "", 1000,0., 10.00);
  h1s_[Form("%s_h1_mu_dxy",      label_.c_str())] = new TH1F(Form("%s_h1_mu_dxy",label_.c_str()),      "",10000,0.,  1.00);
  h1s_[Form("%s_h1_mu_dxyPull",  label_.c_str())] = new TH1F(Form("%s_h1_mu_dxyPull",label_.c_str()),  "",10000,0., 10.00);
  h1s_[Form("%s_h1_mu_dz",       label_.c_str())] = new TH1F(Form("%s_h1_mu_dz",label_.c_str()),       "",10000,0.,  1.00);
  h1s_[Form("%s_h1_mu_dzPull",   label_.c_str())] = new TH1F(Form("%s_h1_mu_dzPull",label_.c_str()),   "",10000,0.,  1.00);
  h1s_[Form("%s_h1_mu_isLoose1", label_.c_str())] = new TH1F(Form("%s_h1_mu_isLoose1",label_.c_str()), "",    2,0.,  2.00);
  h1s_[Form("%s_h1_mu_isLoose2", label_.c_str())] = new TH1F(Form("%s_h1_mu_isLoose2",label_.c_str()), "",    2,0.,  2.00);
  h1s_[Form("%s_h1_mu_isMedium1",label_.c_str())] = new TH1F(Form("%s_h1_mu_isMedium1",label_.c_str()),"",    2,0.,  2.00);
  h1s_[Form("%s_h1_mu_isMedium2",label_.c_str())] = new TH1F(Form("%s_h1_mu_isMedium2",label_.c_str()),"",    2,0.,  2.00);
  h1s_[Form("%s_h1_mu_isTight1", label_.c_str())] = new TH1F(Form("%s_h1_mu_isTight1",label_.c_str()), "",    2,0.,  2.00);
  h1s_[Form("%s_h1_mu_isTight2", label_.c_str())] = new TH1F(Form("%s_h1_mu_isTight2",label_.c_str()), "",    2,0.,  2.00);
  h1s_[Form("%s_h1_mu_relIso1",  label_.c_str())] = new TH1F(Form("%s_h1_mu_relIso1",label_.c_str()),  "",10000,0., 10.00);
  h1s_[Form("%s_h1_mu_relIso2",  label_.c_str())] = new TH1F(Form("%s_h1_mu_relIso2",label_.c_str()),  "",10000,0., 10.00);
  
  h1s_[Form("%s_h1_jets_20GeV_n",      label_.c_str())] = new TH1F(Form("%s_h1_jets_20GeV_n",      label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_25GeV_n",      label_.c_str())] = new TH1F(Form("%s_h1_jets_25GeV_n",      label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_30GeV_n",      label_.c_str())] = new TH1F(Form("%s_h1_jets_30GeV_n",      label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_20GeV_bTagL_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_20GeV_bTagL_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_25GeV_bTagL_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_25GeV_bTagL_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_30GeV_bTagL_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_30GeV_bTagL_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_20GeV_bTagM_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_20GeV_bTagM_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_25GeV_bTagM_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_25GeV_bTagM_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_30GeV_bTagM_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_30GeV_bTagM_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_20GeV_bTagT_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_20GeV_bTagT_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_25GeV_bTagT_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_25GeV_bTagT_n",label_.c_str()),"",10,-0.5,9.5);
  h1s_[Form("%s_h1_jets_30GeV_bTagT_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_30GeV_bTagT_n",label_.c_str()),"",10,-0.5,9.5);
  
  // h1s_[Form("%s_h1_jets_puppi_20GeV_n",      label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_20GeV_n",      label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_25GeV_n",      label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_25GeV_n",      label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_30GeV_n",      label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_30GeV_n",      label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_20GeV_bTagL_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_20GeV_bTagL_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_25GeV_bTagL_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_25GeV_bTagL_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_30GeV_bTagL_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_30GeV_bTagL_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_20GeV_bTagM_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_20GeV_bTagM_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_25GeV_bTagM_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_25GeV_bTagM_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_30GeV_bTagM_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_30GeV_bTagM_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_20GeV_bTagT_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_20GeV_bTagT_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_25GeV_bTagT_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_25GeV_bTagT_n",label_.c_str()),"",20,-0.5,9.5);
  // h1s_[Form("%s_h1_jets_puppi_30GeV_bTagT_n",label_.c_str())] = new TH1F(Form("%s_h1_jets_puppi_30GeV_bTagT_n",label_.c_str()),"",20,-0.5,9.5);
  
  h1s_[Form("%s_h1_met_pt", label_.c_str())] = new TH1F(Form("%s_h1_met_pt", label_.c_str()), "",10000,0.,500.);
  h1s_[Form("%s_h1_met_sig",label_.c_str())] = new TH1F(Form("%s_h1_met_sig",label_.c_str()), "",10000,0.,100.);
  
  // h1s_[Form("%s_h1_met_puppi_pt",  label_.c_str())] = new TH1F(Form("%s_h1_met_puppi_pt",label_.c_str()), "",10000,0.,500.);
  // h1s_[Form("%s_h1_met_puppi_sig", label_.c_str())] = new TH1F(Form("%s_h1_met_puppi_sig",label_.c_str()),"",10000,0.,100.);
  
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
                             std::vector<particle>& mu, TreeVars* tv)
{
  particle H;
  H.charge = 0;
  for(unsigned int ii = 0; ii < mu.size(); ++ii)
  {
    H.v += mu.at(ii).v;
    H.charge += mu.at(ii).charge;
  }
  H.pdgId = 25;
  
  h1s_[Form("%s_h1_H_pt",  label_.c_str())] -> Fill( H.v.Pt(), weight );
  h1s_[Form("%s_h1_H_eta", label_.c_str())] -> Fill( H.v.Eta(),weight );
  h1s_[Form("%s_h1_H_mass",label_.c_str())] -> Fill( H.v.M(),  weight );
  
  h1s_[Form("%s_h1_mu_pt1", label_.c_str())] -> Fill( mu.at(0).v.Pt(), weight );
  h1s_[Form("%s_h1_mu_pt2", label_.c_str())] -> Fill( mu.at(1).v.Pt(), weight );
  h1s_[Form("%s_h1_mu_eta", label_.c_str())] -> Fill( mu.at(0).v.Eta(),weight );
  h1s_[Form("%s_h1_mu_eta", label_.c_str())] -> Fill( mu.at(1).v.Eta(),weight );
  h1s_[Form("%s_h1_mu_DEta",label_.c_str())] -> Fill( DeltaEta(mu.at(0).v.Eta(),mu.at(1).v.Eta()),weight );
  h1s_[Form("%s_h1_mu_DPhi",label_.c_str())] -> Fill( DeltaPhi(mu.at(0).v.Phi(),mu.at(1).v.Phi()),weight );
  h1s_[Form("%s_h1_mu_DR",  label_.c_str())] -> Fill( DeltaR(mu.at(0).v.Eta(),mu.at(0).v.Phi(),mu.at(1).v.Eta(),mu.at(1).v.Phi()),weight );
  
  if( tv )
  {
    float mu1_iso = tv->muons_pfIsoChargedHadron->at(mu.at(0).it) +
      std::max(0.,tv->muons_pfIsoNeutralHadron->at(mu.at(0).it)+tv->muons_pfIsoPhoton->at(mu.at(0).it)-0.5*tv->muons_pfIsoPU->at(mu.at(0).it));
    float mu2_iso = tv->muons_pfIsoChargedHadron->at(mu.at(1).it) +
      std::max(0.,tv->muons_pfIsoNeutralHadron->at(mu.at(1).it)+tv->muons_pfIsoPhoton->at(mu.at(1).it)-0.5*tv->muons_pfIsoPU->at(mu.at(1).it));
    
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
    h1s_[Form("%s_h1_mu_relIso1",  label_.c_str())] -> Fill( mu1_iso/mu.at(0).v.Pt(),weight );
    h1s_[Form("%s_h1_mu_relIso2",  label_.c_str())] -> Fill( mu2_iso/mu.at(1).v.Pt(),weight );
    
    h1s_[Form("%s_h1_met_pt", label_.c_str())] -> Fill( tv->met_pt, weight );
    h1s_[Form("%s_h1_met_sig",label_.c_str())] -> Fill( tv->met_sig,weight );
    // h1s_[Form("%s_h1_met_puppi_pt", label_.c_str())] -> Fill( tv->met_puppi_pt, weight );
    // h1s_[Form("%s_h1_met_puppi_sig",label_.c_str())] -> Fill( tv->met_puppi_sig,weight );
    
    int jets_20GeV_n = 0;
    int jets_25GeV_n = 0;
    int jets_30GeV_n = 0;
    int jets_20GeV_bTagL_n = 0;
    int jets_25GeV_bTagL_n = 0;
    int jets_30GeV_bTagL_n = 0;
    int jets_20GeV_bTagM_n = 0;
    int jets_25GeV_bTagM_n = 0;
    int jets_30GeV_bTagM_n = 0;
    int jets_20GeV_bTagT_n = 0;
    int jets_25GeV_bTagT_n = 0;
    int jets_30GeV_bTagT_n = 0;
    for(unsigned int ii = 0; ii < tv->jets_pt->size(); ++ii)
    {
      if( tv->jets_pt->at(ii) < 20. ) continue;
      if( fabs(tv->jets_eta->at(ii)) > 3.0 ) continue;
      
      bool skipJet = false;
      for(unsigned int jj = 0; jj < mu.size(); ++jj)
        if( DeltaR(tv->jets_eta->at(ii),tv->jets_phi->at(ii),mu.at(jj).v.Eta(),mu.at(jj).v.Phi()) < 0.4 ) { skipJet = true; break; }
      if( skipJet ) continue;
      
      if( tv->jets_pt->at(ii) >= 20. )
      {
        ++jets_20GeV_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.460 ) ++jets_20GeV_bTagL_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.800 ) ++jets_20GeV_bTagM_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.935 ) ++jets_20GeV_bTagT_n;
      }
      if( tv->jets_pt->at(ii) >= 25. )
      {
        ++jets_25GeV_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.460 ) ++jets_25GeV_bTagL_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.800 ) ++jets_25GeV_bTagM_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.935 ) ++jets_25GeV_bTagT_n;
      }
      if( tv->jets_pt->at(ii) >= 30. )
      {
        ++jets_30GeV_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.460 ) ++jets_30GeV_bTagL_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.800 ) ++jets_30GeV_bTagM_n;
        if( tv->jets_bTag->at(ii).at(0) > 0.935 ) ++jets_30GeV_bTagT_n;
      }
    }
    h1s_[Form("%s_h1_jets_20GeV_n",label_.c_str())] -> Fill( jets_20GeV_n);
    h1s_[Form("%s_h1_jets_25GeV_n",label_.c_str())] -> Fill( jets_25GeV_n);
    h1s_[Form("%s_h1_jets_30GeV_n",label_.c_str())] -> Fill( jets_30GeV_n);
    h1s_[Form("%s_h1_jets_20GeV_bTagL_n",label_.c_str())] -> Fill( jets_20GeV_bTagL_n);
    h1s_[Form("%s_h1_jets_25GeV_bTagL_n",label_.c_str())] -> Fill( jets_25GeV_bTagL_n);
    h1s_[Form("%s_h1_jets_30GeV_bTagL_n",label_.c_str())] -> Fill( jets_30GeV_bTagL_n);
    h1s_[Form("%s_h1_jets_20GeV_bTagM_n",label_.c_str())] -> Fill( jets_20GeV_bTagM_n);
    h1s_[Form("%s_h1_jets_25GeV_bTagM_n",label_.c_str())] -> Fill( jets_25GeV_bTagM_n);
    h1s_[Form("%s_h1_jets_30GeV_bTagM_n",label_.c_str())] -> Fill( jets_30GeV_bTagM_n);
    h1s_[Form("%s_h1_jets_20GeV_bTagT_n",label_.c_str())] -> Fill( jets_20GeV_bTagT_n);
    h1s_[Form("%s_h1_jets_25GeV_bTagT_n",label_.c_str())] -> Fill( jets_25GeV_bTagT_n);
    h1s_[Form("%s_h1_jets_30GeV_bTagT_n",label_.c_str())] -> Fill( jets_30GeV_bTagT_n);
    
    // int jets_puppi_20GeV_n = 0;
    // int jets_puppi_25GeV_n = 0;
    // int jets_puppi_30GeV_n = 0;
    // int jets_puppi_20GeV_bTagL_n = 0;
    // int jets_puppi_25GeV_bTagL_n = 0;
    // int jets_puppi_30GeV_bTagL_n = 0;
    // int jets_puppi_20GeV_bTagM_n = 0;
    // int jets_puppi_25GeV_bTagM_n = 0;
    // int jets_puppi_30GeV_bTagM_n = 0;
    // int jets_puppi_20GeV_bTagT_n = 0;
    // int jets_puppi_25GeV_bTagT_n = 0;
    // int jets_puppi_30GeV_bTagT_n = 0;
    // for(unsigned int ii = 0; ii < tv->jets_puppi_pt->size(); ++ii)
    // {
    //   if( tv->jets_puppi_pt->at(ii) < 20. ) continue;
    //   if( fabs(tv->jets_puppi_eta->at(ii)) > 3.0 ) continue;
      
    //   bool skipJet = false;
    //   for(unsigned int jj = 0; jj < mu.size(); ++jj)
    //     if( DeltaR(tv->jets_puppi_eta->at(ii),tv->jets_puppi_phi->at(ii),mu.at(jj).v.Eta(),mu.at(jj).v.Phi()) < 0.4 ) { skipJet = true; break; }
    //   if( skipJet ) continue;
      
    //   if( tv->jets_puppi_pt->at(ii) >= 20. )
    //   {
    //     ++jets_puppi_20GeV_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.460 ) ++jets_puppi_20GeV_bTagL_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.800 ) ++jets_puppi_20GeV_bTagM_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.935 ) ++jets_puppi_20GeV_bTagT_n;
    //   }
    //   if( tv->jets_puppi_pt->at(ii) >= 25. )
    //   {
    //     ++jets_puppi_25GeV_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.460 ) ++jets_puppi_25GeV_bTagL_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.800 ) ++jets_puppi_25GeV_bTagM_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.935 ) ++jets_puppi_25GeV_bTagT_n;
    //   }
    //   if( tv->jets_puppi_pt->at(ii) >= 30. )
    //   {
    //     ++jets_puppi_30GeV_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.460 ) ++jets_puppi_30GeV_bTagL_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.800 ) ++jets_puppi_30GeV_bTagM_n;
    //     if( tv->jets_puppi_bTag->at(ii).at(0) > 0.935 ) ++jets_puppi_30GeV_bTagT_n;
    //   }
    // }
    // h1s_[Form("%s_h1_jets_puppi_20GeV_n",label_.c_str())] -> Fill( jets_puppi_20GeV_n);
    // h1s_[Form("%s_h1_jets_puppi_25GeV_n",label_.c_str())] -> Fill( jets_puppi_25GeV_n);
    // h1s_[Form("%s_h1_jets_puppi_30GeV_n",label_.c_str())] -> Fill( jets_puppi_30GeV_n);
    // h1s_[Form("%s_h1_jets_puppi_20GeV_bTagL_n",label_.c_str())] -> Fill( jets_puppi_20GeV_bTagL_n);
    // h1s_[Form("%s_h1_jets_puppi_25GeV_bTagL_n",label_.c_str())] -> Fill( jets_puppi_25GeV_bTagL_n);
    // h1s_[Form("%s_h1_jets_puppi_30GeV_bTagL_n",label_.c_str())] -> Fill( jets_puppi_30GeV_bTagL_n);
    // h1s_[Form("%s_h1_jets_puppi_20GeV_bTagM_n",label_.c_str())] -> Fill( jets_puppi_20GeV_bTagM_n);
    // h1s_[Form("%s_h1_jets_puppi_25GeV_bTagM_n",label_.c_str())] -> Fill( jets_puppi_25GeV_bTagM_n);
    // h1s_[Form("%s_h1_jets_puppi_30GeV_bTagM_n",label_.c_str())] -> Fill( jets_puppi_30GeV_bTagM_n);
    // h1s_[Form("%s_h1_jets_puppi_20GeV_bTagT_n",label_.c_str())] -> Fill( jets_puppi_20GeV_bTagT_n);
    // h1s_[Form("%s_h1_jets_puppi_25GeV_bTagT_n",label_.c_str())] -> Fill( jets_puppi_25GeV_bTagT_n);
    // h1s_[Form("%s_h1_jets_puppi_30GeV_bTagT_n",label_.c_str())] -> Fill( jets_puppi_30GeV_bTagT_n);
  }
}
