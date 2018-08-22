#include <iostream>
#include <map>

#include "TFile.h"
#include "TH1F.h"



class PUReweighting
{
public:
  PUReweighting(TH1F* h_pileup_data, TH1F* h_pileup_mc);
  ~PUReweighting();
  
  float GetPUWeight(const float& trueNumInteractions);
  
  TH1F* GetWeightsHistogram();
  
private:
  std::map<int,float> pileup_data_;
  std::map<int,float> pileup_mc_;
};
