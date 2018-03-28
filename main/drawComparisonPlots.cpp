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
#include "TCanvas.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> drawComparisonPlots.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //--- open files and get trees
  std::vector<std::string> vars = opts.GetOpt<std::vector<std::string> >("Input.vars");
  std::string outputFolder = opts.GetOpt<std::string>("Output.outputFolder");
  
  setTDRStyle();
  
  Plotter myPlotter(opts,vars,outputFolder);
  myPlotter.DrawPlots();
  
  /*

  
  TH1F* h_data;
  TH1F* h_mc;
  TH1F* h_data_extra;
  TH1F* h_mc_extra;
  
  TCanvas* c;
  
  float yMax;
  float yMin;
  
  
  for(unsigned int ii = 0; ii < vars.size(); ++ii)
  {
    std::string varName = vars.at(ii);
    
    h_data = new TH1F(Form("h_data_%s",varName.c_str()),"",200,0.,200.);
    h_mc   = new TH1F(Form("h_mc_%s",varName.c_str()),  "",200,0.,200.);
    
    h_data_extra = new TH1F(Form("h_data_extra_%s",varName.c_str()),"",202,-1.,201.);
    h_mc_extra   = new TH1F(Form("h_mc_extra_%s",  varName.c_str()),"",202,-1.,201.);
    
    t_data -> Draw( Form("%s >>h_data_%s",varName.c_str(),varName.c_str()), "", "goff" );
    t_mc   -> Draw( Form("%s >>h_mc_%s",varName.c_str(),varName.c_str()), "", "goff" );
    
    h_data_extra -> Sumw2();
    h_mc_extra   -> Sumw2();
    for(int bin = 0; bin <= h_data->GetNbinsX()+1; ++bin)
    {
      h_data_extra -> SetBinContent(bin+1,h_data->GetBinContent(bin));
      h_data_extra -> SetBinError(bin+1,h_data->GetBinError(bin));
    }
    for(int bin = 0; bin <= h_mc->GetNbinsX()+1; ++bin)
    {
      h_mc_extra -> SetBinContent(bin+1,h_mc->GetBinContent(bin));
      h_mc_extra -> SetBinError(bin+1,h_mc->GetBinError(bin));
    }
    
    h_data_extra -> Scale(1./h_data_extra->Integral());
    h_mc_extra -> Scale(1./h_mc_extra->Integral());
    
    yMax = -999999999.;
    yMin = +999999999.;
    for(int bin = 1; bin <= h_data_extra->GetNbinsX(); ++bin)
    {
      if( h_data_extra->GetBinContent(bin) > yMax ) yMax = h_data_extra->GetBinContent(bin);
      if( h_data_extra->GetBinContent(bin) < yMin && h_data_extra->GetBinContent(bin) != 0. ) yMin = h_data_extra->GetBinContent(bin);
    }
    for(int bin = 1; bin <= h_mc_extra->GetNbinsX(); ++bin)
    {
      if( h_mc_extra->GetBinContent(bin) > yMax ) yMax = h_mc_extra->GetBinContent(bin);
      if( h_mc_extra->GetBinContent(bin) < yMin && h_mc_extra->GetBinContent(bin) != 0. ) yMin = h_mc_extra->GetBinContent(bin);
    }
    
    
    c = new TCanvas("c","c");
    
    h_data_extra -> SetMarkerColor(kBlack);
    h_data_extra -> SetTitle(";mass_{#mu^{+}#mu^{-} KK#pi} (GeV);event fraction");
    h_data_extra -> SetMaximum(1.25*yMax);
    
    h_mc_extra -> SetLineColor(kRed);
    h_mc_extra -> SetLineWidth(2);
    
    h_data_extra -> Draw("P");
    h_mc_extra -> Draw("hist,same");
    
    c -> Print(Form("%s/c_%s.pdf",outputFolder.c_str(),varName.c_str()));
    c -> Print(Form("%s/c_%s.png",outputFolder.c_str(),varName.c_str()));
    
    delete c;
    
    
    c = new TCanvas("c","c");
    
    h_data_extra -> SetMarkerColor(kBlack);
    h_data_extra -> SetTitle(";mass_{#mu^{+}#mu^{-} KK#pi} (GeV);event fraction");
    h_data_extra -> SetMaximum(5.*yMax);
    h_data_extra -> SetMinimum(yMin/5.);
    
    h_mc_extra -> SetLineColor(kRed);
    h_mc_extra -> SetLineWidth(2);
    
    h_data_extra -> Draw("P");
    h_mc_extra -> Draw("hist,same");
    
    gPad -> SetLogy();
    
    c -> Print(Form("%s/c_%s_log.pdf",outputFolder.c_str(),varName.c_str()));
    c -> Print(Form("%s/c_%s_log.png",outputFolder.c_str(),varName.c_str()));
    
    delete c;
  }
  
  */
  
  
  return 0;
}
