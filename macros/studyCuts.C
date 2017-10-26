#define PI 3.14159265359
#define CONF68 0.682689492137

Double_t (*pBound)(Double_t,Double_t,Double_t,Bool_t) = &TEfficiency::ClopperPearson; // default method



void studyCuts_DPhi()
{
  std::vector<float> DPhiCutVals;
  DPhiCutVals.push_back(0.);
  DPhiCutVals.push_back(0.25);
  DPhiCutVals.push_back(0.5);
  DPhiCutVals.push_back(0.75);
  DPhiCutVals.push_back(1.);
  DPhiCutVals.push_back(1.25);
  DPhiCutVals.push_back(1.5);
  DPhiCutVals.push_back(1.75);
  DPhiCutVals.push_back(2.);
  DPhiCutVals.push_back(2.25);
  DPhiCutVals.push_back(2.5);
  DPhiCutVals.push_back(2.75);
  DPhiCutVals.push_back(3.);
  
  TGraphAsymmErrors* g1_eff_DPhiCut = new TGraphAsymmErrors();
  TGraphAsymmErrors* g1_pur_DPhiCut = new TGraphAsymmErrors();
  TGraphAsymmErrors* g1_effpur_DPhiCut = new TGraphAsymmErrors();
  
  int point = 0;
  for(int ii = 0; ii < DPhiCutVals.size(); ++ii)
  {
    float DPhiCut = DPhiCutVals.at(ii);
    TFile* inFile = TFile::Open(Form("plots/plots_KKPi_PU__DRMin_Ds_JPsi_%.2f__DRMax_Ds_%.2f.root",DPhiCut,0.4),"READ");
    TH1F* histo = (TH1F*)( inFile->Get("h1_nEvents") );

    
    float t = histo->GetBinContent(7);
    float p = histo->GetBinContent(10);
    float eff = 0.;
    if(t) eff = p/t;
    
    float eLow = pBound(t,p,CONF68,false);
    float eHig = pBound(t,p,CONF68,true);
    
    g1_eff_DPhiCut -> SetPoint(point,DPhiCut,eff);
    g1_eff_DPhiCut -> SetPointError(point,0.,0.,eff-eLow,eHig-eff);

    
    t = histo->GetBinContent(10);
    p = histo->GetBinContent(11);
    eff = 0.;
    if(t) eff = p/t;
    
    eLow = pBound(t,p,CONF68,false);
    eHig = pBound(t,p,CONF68,true);
    
    g1_pur_DPhiCut -> SetPoint(point,DPhiCut,eff);
    g1_pur_DPhiCut -> SetPointError(point,0.,0.,eff-eLow,eHig-eff);
    
    
    t = histo->GetBinContent(7);
    p = histo->GetBinContent(11);
    eff = 0.;
    if(t) eff = p/t;
    
    eLow = pBound(t,p,CONF68,false);
    eHig = pBound(t,p,CONF68,true);
    
    g1_effpur_DPhiCut -> SetPoint(point,DPhiCut,eff);
    g1_effpur_DPhiCut -> SetPointError(point,0.,0.,eff-eLow,eHig-eff);

    
    ++point;
  }
  

  TCanvas* c1 = new TCanvas("c1","c1",1321,398);
  c1 -> Divide(3,1);
  c1 -> cd(1);
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-0.1,0.,PI+0.1,1.) );
  hPad -> SetTitle(Form(";#Delta#phi_{lead. track -- J/Psi} (rad.);selection efficiency"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  gPad -> SetGridy();
  g1_eff_DPhiCut -> SetMarkerSize(0.7);
  g1_eff_DPhiCut -> Draw("PL,same");
  c1 -> cd(2);
  hPad = (TH1F*)( gPad->DrawFrame(-0.1,0.,PI+0.1,1.) );
  hPad -> SetTitle(Form(";#Delta#phi_{lead. track -- J/Psi} (rad.);gen. match purity"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  gPad -> SetGridy();
  g1_pur_DPhiCut -> SetMarkerSize(0.7);
  g1_pur_DPhiCut -> Draw("PL,same");  
  c1 -> cd(3);
  hPad = (TH1F*)( gPad->DrawFrame(-0.1,0.,PI+0.1,1.) );
  hPad -> SetTitle(Form(";#Delta#phi_{lead. track -- J/Psi} (rad.);eff. #times purity"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  gPad -> SetGridy();
  g1_effpur_DPhiCut -> SetMarkerSize(0.7);
  g1_effpur_DPhiCut -> Draw("PL,same");  
}



void studyCuts_DR()
{
  std::vector<float> DRCutVals;
  DRCutVals.push_back(0.05);
  DRCutVals.push_back(0.06);
  DRCutVals.push_back(0.07);
  DRCutVals.push_back(0.08);
  DRCutVals.push_back(0.09);
  DRCutVals.push_back(0.10);
  DRCutVals.push_back(0.12);
  DRCutVals.push_back(0.14);
  DRCutVals.push_back(0.16);
  DRCutVals.push_back(0.18);
  DRCutVals.push_back(0.20);
  DRCutVals.push_back(0.25);
  DRCutVals.push_back(0.30);
  DRCutVals.push_back(0.35);
  DRCutVals.push_back(0.40);
  DRCutVals.push_back(0.50);
  DRCutVals.push_back(0.60);
  DRCutVals.push_back(0.70);
  DRCutVals.push_back(0.80);
  
  TGraphAsymmErrors* g1_eff_DRCut = new TGraphAsymmErrors();
  TGraphAsymmErrors* g1_pur_DRCut = new TGraphAsymmErrors();
  TGraphAsymmErrors* g1_effpur_DRCut = new TGraphAsymmErrors();
  
  int point = 0;
  for(int ii = 0; ii < DRCutVals.size(); ++ii)
  {
    float DRCut = DRCutVals.at(ii);
    TFile* inFile = TFile::Open(Form("plots/plots_KKPi__DRMin_Ds_JPsi_%.2f__DRMax_Ds_%.2f.root",0.,DRCut),"READ");
    TH1F* histo = (TH1F*)( inFile->Get("h1_nEvents") );

    
    float t = histo->GetBinContent(7);
    float p = histo->GetBinContent(10);
    float eff = 0.;
    if(t) eff = p/t;
    
    float eLow = pBound(t,p,CONF68,false);
    float eHig = pBound(t,p,CONF68,true);
    
    g1_eff_DRCut -> SetPoint(point,DRCut,eff);
    g1_eff_DRCut -> SetPointError(point,0.,0.,eff-eLow,eHig-eff);

    
    t = histo->GetBinContent(10);
    p = histo->GetBinContent(11);
    eff = 0.;
    if(t) eff = p/t;
    
    eLow = pBound(t,p,CONF68,false);
    eHig = pBound(t,p,CONF68,true);
    
    g1_pur_DRCut -> SetPoint(point,DRCut,eff);
    g1_pur_DRCut -> SetPointError(point,0.,0.,eff-eLow,eHig-eff);
    
    
    t = histo->GetBinContent(7);
    p = histo->GetBinContent(11);
    eff = 0.;
    if(t) eff = p/t;
    
    eLow = pBound(t,p,CONF68,false);
    eHig = pBound(t,p,CONF68,true);
    
    g1_effpur_DRCut -> SetPoint(point,DRCut,eff);
    g1_effpur_DRCut -> SetPointError(point,0.,0.,eff-eLow,eHig-eff);

    
    ++point;
  }
  

  TCanvas* c1 = new TCanvas("c1","c1",1321,398);
  c1 -> Divide(3,1);
  c1 -> cd(1);
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.1,1.) );
  hPad -> SetTitle(Form(";#DeltaR_{lead. track -- other tracks} (rad.);selection efficiency"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  gPad -> SetGridy();
  g1_eff_DRCut -> SetMarkerSize(0.7);
  g1_eff_DRCut -> Draw("PL,same");
  c1 -> cd(2);
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.,1.) );
  hPad -> SetTitle(Form(";#DeltaR_{lead. track -- other tracks} (rad.);gen. match purity"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  gPad -> SetGridy();
  g1_pur_DRCut -> SetMarkerSize(0.7);
  g1_pur_DRCut -> Draw("PL,same");  
  c1 -> cd(3);
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.,1.) );
  hPad -> SetTitle(Form(";#DeltaR_{lead. track -- other tracks} (rad.);eff. #times purity"));
  hPad -> SetLineColor(kWhite);
  hPad -> Draw();
  gPad -> SetGridy();
  g1_effpur_DRCut -> SetMarkerSize(0.7);
  g1_effpur_DRCut -> Draw("PL,same");  
}
