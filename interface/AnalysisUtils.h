#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "interface/TreeUtils.h"
#include "interface/ParticleNames.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#define PI 3.14159265359
#define MASS_W 80.385
#define MASS_JPSI 3.096900
#define MASS_DS 1.96828
#define MASS_K 0.493677
#define MASS_PI 0.13957018
#define MASS_PI0 0.1349766



struct particle
{
  TLorentzVector v;
  int charge;
  int pdgId;
  int it;
};

struct PtSort
{
  inline bool operator()(const particle& p1, const particle& p2)
  {
    return (p1.v.Pt() > p2.v.Pt());
  }
};

struct EtaSort
{
  inline bool operator()(const particle& p1, const particle& p2)
  {
    return (fabs(p1.v.Eta()) > fabs(p2.v.Eta()));
  }
};

struct NameSort
{
  inline bool operator()(const particle& p1, const particle& p2)
  {
    return (GetParticleName(p1.pdgId) < GetParticleName(p2.pdgId));
  }
};


struct FindPair
{
FindPair(const std::string& key)
: m_key(key) {}
  std::string m_key;
  bool operator()
    ( const std::pair<std::string,int>& p )
    {
      return( p.first == m_key);
    }
};

struct PairSort
{
  inline bool operator()(const std::pair<std::string,int>& p1, const std::pair<std::string,int>& p2)
  {
    return (p1.second > p2.second);
  }
};

class StdHistoSet
{
public:
  StdHistoSet(const std::string& label, TFile* outFile);
  ~StdHistoSet();

  void FillHistos(const float& weight,
                  std::vector<particle>& mu, std::vector<particle>& kepi, TreeVars* tv = NULL);
  
private:
  std::string label_;
  TFile* outFile_;
  
  std::map<std::string,TH1F*> h1s_;
};


float DeltaEta(const float& eta1, const float& eta2);
float DeltaPhi(const float& phi1, const float& phi2);
float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);

bool IsMatching(const particle& p1, const particle p2,
                const float& DRMax = 0.05, const float& ptRatioMax = 0.1);
bool IsMatching(const std::vector<particle>& vp1, const std::vector<particle>& vp2,
                const float& DRMax = 0.05, const float& ptRatioMax = 0.1,
                const bool& verbosity = false);

int GetBestMatch(const particle& p, std::vector<particle>& vec, std::vector<int>* vetoVec = NULL);


std::vector<std::vector<particle> > AssignMass(std::vector<particle>& v);
std::vector<particle> AssignMassDs(std::vector<particle>& v);

void PrintEvent(std::vector<particle>& mu, std::vector<particle>& kepi);
std::ostream& operator<<(std::ostream& os, const TLorentzVector& v);
std::ostream& operator<<(std::ostream& os, const particle& p);

#endif
