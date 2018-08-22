#include "interface/PUReweighting.h"



PUReweighting::PUReweighting(TH1F* h_pileup_data, TH1F* h_pileup_mc)
{
  if( h_pileup_data->GetNbinsX() != h_pileup_mc->GetNbinsX() )
  {
    std::cerr << ">>>PUReweighting::PUReweighting::Warning: different number of bins in data and mc histograms detected!" << std::endl;
  }
  
  h_pileup_data -> Scale(1./h_pileup_data->Integral());
  h_pileup_mc -> Scale(1./h_pileup_mc->Integral());
  
  for(int bin = 1; bin <= h_pileup_data->GetNbinsX(); ++bin)
    pileup_data_[int(h_pileup_data->GetBinCenter(bin))] = h_pileup_data -> GetBinContent(bin);
  
  for(int bin = 1; bin <= h_pileup_mc->GetNbinsX(); ++bin)
    pileup_mc_[int(h_pileup_mc->GetBinCenter(bin))] = h_pileup_mc -> GetBinContent(bin);
}



PUReweighting::~PUReweighting()
{}



float PUReweighting::GetPUWeight(const float& trueNumInteractions)
{
  float val = 1. / pileup_mc_[trueNumInteractions] * pileup_data_[trueNumInteractions];
  if( isinf(val) )
  {
    std::cerr << "!!! PUReweighting::GetPUWeight error:" << std::endl;
    std::cerr << "   trueNumInteractions = " << trueNumInteractions << std::endl;
    std::cerr << "   pileup_mc_[trueNumInteractions] = " << pileup_mc_[trueNumInteractions] << std::endl;
    std::cerr << "   pileup_data_[trueNumInteractions] = " << pileup_data_[trueNumInteractions] << std::endl;
    std::cerr << "PUReweighting::GetPUWeight error !!!" << std::endl;
  }
  
  return val;
}



TH1F* PUReweighting::GetWeightsHistogram()
{
  TH1F* h_PUWeights = new TH1F("pileup_weights","",100,-0.5,99.5);
  for(int bin = 1; bin <= 100; ++bin)
    if( pileup_mc_[int(h_PUWeights->GetBinCenter(bin))] > 0. )
      h_PUWeights -> SetBinContent(bin,GetPUWeight(h_PUWeights->GetBinCenter(bin)));
  
  return h_PUWeights;
}
