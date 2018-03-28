#include "interface/AnalysisUtils.h"

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"



class StdHistoSet
{
public:
  StdHistoSet(const std::string& label, TFile* outFile);
  ~StdHistoSet();
  
  void FillHistos(const float& weight,
                  std::vector<particle>& mu, TreeVars* tv = NULL);

private:
  std::string label_;
  TFile* outFile_;

  std::map<std::string,TH1F*> h1s_;
};
