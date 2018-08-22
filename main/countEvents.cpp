#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/SetTDRStyle.h"
#include "interface/Plotter.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> countEvents.cpp::usage:   " << argv[0] << " inputFileList" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::string inputFileList(argv[1]);
  
  TH1F* h_pileup = NULL;
  std::ifstream list(inputFileList.c_str(),std::ios::in);
  std::string fileName;
  while(1)
  {
    getline(list,fileName,'\n');
    if( !list.good() ) break;
    
    TFile* inFile = TFile::Open(fileName.c_str(),"READ");
    
    if( h_pileup == NULL )
      h_pileup = (TH1F*)( inFile->Get("DumpPU/pileup"));
    else
    {
      TH1F* h_temp = (TH1F*)( inFile->Get("DumpPU/pileup"));
      h_pileup -> Add( h_temp );
    } 
  }
  
  std::cout << ">>> nTotEvents: " << std::fixed << std::setprecision(0) << h_pileup->GetEntries() << std::endl;
  
  return 0;
}
