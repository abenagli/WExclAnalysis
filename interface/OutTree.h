#ifndef _OUT__TREE
#define _OUT__TREE

#include "DynamicTTree/interface/DynamicTTreeBase.h"



//---Define the TTree branches
#define DYNAMIC_TREE_NAME OutTree

#define DATA_TABLE                    \
  DATA(float, weight)                 \
  DATA(float, weight_PU)              \
  DATA(float, weight_MC)              \
  DATA(float, H_pt)                   \
  DATA(float, H_eta)                  \
  DATA(float, H_phi)                  \
  DATA(float, H_mass)                 \
  DATA(float, mu1_pt)                 \
  DATA(float, mu2_pt)                 \
  DATA(float, mu1_eta)                \
  DATA(float, mu2_eta)                \
  DATA(float, mu1_phi)                \
  DATA(float, mu2_phi)                \
  DATA(float, mu1_dxy)                \
  DATA(float, mu2_dxy)                \
  DATA(float, mu1_dxyPull)            \
  DATA(float, mu2_dxyPull)            \
  DATA(float, mu1_dz)                 \
  DATA(float, mu2_dz)                 \
  DATA(float, mu1_dzPull)             \
  DATA(float, mu2_dzPull)             \
  DATA(float, mu1_relIso)             \
  DATA(float, mu2_relIso)             \
  DATA(int, mu1_isL)                  \
  DATA(int, mu2_isL)                  \
  DATA(int, mu1_isM)                  \
  DATA(int, mu2_isM)                  \
  DATA(int, mu1_isT)                  \
  DATA(int, mu2_isT)                  \
  DATA(float, mu_Deta)                \
  DATA(float, mu_Dphi)                \
  DATA(float, mu_DR)                  \
  DATA(int, jets_20GeV_n)             \
  DATA(int, jets_25GeV_n)             \
  DATA(int, jets_30GeV_n)             \
  DATA(int, jets_20GeV_bTagL_n)       \
  DATA(int, jets_25GeV_bTagL_n)       \
  DATA(int, jets_30GeV_bTagL_n)       \
  DATA(int, jets_20GeV_bTagM_n)       \
  DATA(int, jets_25GeV_bTagM_n)       \
  DATA(int, jets_30GeV_bTagM_n)       \
  DATA(int, jets_20GeV_bTagT_n)       \
  DATA(int, jets_25GeV_bTagT_n)       \
  DATA(int, jets_30GeV_bTagT_n)       \
  DATA(int, jets_puppi_20GeV_n)       \
  DATA(int, jets_puppi_25GeV_n)       \
  DATA(int, jets_puppi_30GeV_n)       \
  DATA(int, jets_puppi_20GeV_bTagL_n) \
  DATA(int, jets_puppi_25GeV_bTagL_n) \
  DATA(int, jets_puppi_30GeV_bTagL_n) \
  DATA(int, jets_puppi_20GeV_bTagM_n) \
  DATA(int, jets_puppi_25GeV_bTagM_n) \
  DATA(int, jets_puppi_30GeV_bTagM_n) \
  DATA(int, jets_puppi_20GeV_bTagT_n) \
  DATA(int, jets_puppi_25GeV_bTagT_n) \
  DATA(int, jets_puppi_30GeV_bTagT_n) \
  DATA(float, met_pt)                 \
  DATA(float, met_phi)                \
  DATA(float, met_sig)                \
  DATA(float, met_puppi_pt)           \
  DATA(float, met_puppi_phi)          \
  DATA(float, met_puppi_sig)

#include "DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
