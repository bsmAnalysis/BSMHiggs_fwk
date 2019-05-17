// -*- C++ -*-
// Package:     BTaggingEffAnalyzer
//
// Based on: https://github.com/rappoccio/usercode/blob/Dev_53x/EDSHyFT/plugins/BTaggingEffAnalyzer.cc
// Created:  10 May 2019


#include <iostream>
#include <map>
//#include <filesystem>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TMath.h"

#define YEAR_2017

//namespace fs = std::filesystem;
using namespace std;

int main(int argc, char* argv[])
{
  if(argc<2) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    exit(0);
  }
  
  float CSVLooseWP = 0.5426;  float CSVMediumWP = 0.800;  float CSVTightWP = 0.935;
  float DeepCSVLooseWP = 0.2219;  float DeepCSVMediumWP = 0.6324; float DeepCSVTightWP = 0.8958;
#ifdef YEAR_2017
//2017 Btag Recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
  CSVLooseWP = 0.5803; CSVMediumWP = 0.8838; CSVTightWP = 0.9693;
  DeepCSVLooseWP = 0.1522; DeepCSVMediumWP = 0.4941; DeepCSVTightWP = 0.8001;
#endif
  float LooseWP = DeepCSVLooseWP; float MediumWP = DeepCSVMediumWP; float TightWP = DeepCSVTightWP;

  const int     ptNBins = 400;
  const double  ptMin = 0.;
  const double  ptMax = 4000.;
  const int     etaNBins = 60;
  const double  etaMin = -3.;
  const double  etaMax = 3.;

  gSystem->Load( "libFWCoreFWLite" );

  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString url = runProcess.getParameter<std::string>("input");
  TString outdir = runProcess.getParameter<std::string>("outdir");
  TString dtag= runProcess.getParameter<std::string>("tag");
  int evStart     = runProcess.getParameter<int>("evStart");
  int evEnd       = runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");
  bool useDeepCSV = runProcess.getParameter<bool>("useDeepCSV");
  bool isMC = runProcess.getParameter<bool>("isMC");

  TString csvTag = "DeepCSV";
  if( !useDeepCSV ){
    LooseWP = CSVLooseWP; TightWP = CSVLooseWP; LooseWP = CSVLooseWP;
    csvTag = "CSV";
  }

  TString outFileUrl( dtag );
  TString outUrl( outdir );
//  outUrl += "/Btag";
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += outFileUrl + "_Btag.root";

// Initialization of Histograms
  SmartSelectionMonitor mon;
  TH2F *h2_BTaggingEff_Denom_b    = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_BTaggingEff_Denom_b"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_BTaggingEff_Denom_c    = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_BTaggingEff_Denom_c"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_BTaggingEff_Denom_udsg = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_BTaggingEff_Denom_udsg"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  
  TH2F *h2_LooseBTaggingEff_Num_b      = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_LooseBTaggingEff_Num_b"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_LooseBTaggingEff_Num_c      = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_LooseBTaggingEff_Num_c"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_LooseBTaggingEff_Num_udsg   = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_LooseBTaggingEff_Num_udsg"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );

  TH2F *h2_MediumBTaggingEff_Num_b      = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_MediumBTaggingEff_Num_b"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_MediumBTaggingEff_Num_c      = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_MediumBTaggingEff_Num_c"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_MediumBTaggingEff_Num_udsg   = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_MediumBTaggingEff_Num_udsg"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );

  TH2F *h2_TightBTaggingEff_Num_b      = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_TightBTaggingEff_Num_b"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_TightBTaggingEff_Num_c      = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_TightBTaggingEff_Num_c"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
  TH2F *h2_TightBTaggingEff_Num_udsg   = (TH2F*) mon.addHistogram( new TH2F (csvTag+TString("_TightBTaggingEff_Num_udsg"), ";p_{T} [GeV];#eta", ptNBins, ptMin, ptMax, etaNBins, etaMin, etaMax) );
// Get ready for the event loop
  DataEvtSummaryHandler summaryHandler_;
  
//  std::string path = "/eos/cms/store/user/yuanc/results_2019_03_08/MuonEG/crab_Data13TeV_MuEG2017E_31Mar2018_v1_0/190506_201903/0000/analysis_103.root";
//  for (const auto & entry : fs::directory_iterator(path)){
//  url = entry.path();
  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0) { return -1; printf("file is 0"); }

  if(file->IsZombie()) return -1;
  if( !summaryHandler_.attachToTree( (TTree *)file->Get(dirname) ) ) {
    file->Close();
    return -1;
  }
  const Int_t totalEntries = summaryHandler_.getEntries();
  if(evEnd<0 || evEnd>totalEntries ) evEnd=totalEntries;
  if(evStart > evEnd ) {
    file->Close();
    return -1;
  }

// Event Loop
  int treeStep = (evEnd-evStart)/50;
  if(treeStep==0)treeStep=1;
  DuplicatesChecker duplicatesChecker;
  printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning the ntuple :");

  for( int iev=evStart; iev<evEnd; iev++) {
    
    if((iev-evStart)%treeStep==0) { 
      printf("."); fflush(stdout);
    }
    summaryHandler_.getEntry(iev);
    DataEvtSummary_t &ev=summaryHandler_.getEvent();
    if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) continue;
    PhysicsEvent_t phys=getPhysicsEventFrom(ev); 

    PhysicsObjectJetCollection &corrJets = phys.jets;
    for(size_t ijet=0; ijet<corrJets.size(); ijet++) {
      
      int partonFlavor = corrJets[ijet].flavid;
      double btag_dsc = -1.;
      if( useDeepCSV ) btag_dsc = corrJets[ijet].btag1;
      else  btag_dsc = corrJets[ijet].btag0;
      //std::cout << "partonFlavor: " << partonFlavor << std::endl;
      if( fabs(partonFlavor)==5 ){
        h2_BTaggingEff_Denom_b->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>LooseWP ) h2_LooseBTaggingEff_Num_b->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>MediumWP ) h2_MediumBTaggingEff_Num_b->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>TightWP ) h2_TightBTaggingEff_Num_b->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
      }
      else if( fabs(partonFlavor)==4 ){
        h2_BTaggingEff_Denom_c->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>LooseWP ) h2_LooseBTaggingEff_Num_c->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>MediumWP ) h2_MediumBTaggingEff_Num_c->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>TightWP ) h2_TightBTaggingEff_Num_c->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
      }
      else{
        h2_BTaggingEff_Denom_udsg->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>LooseWP ) h2_LooseBTaggingEff_Num_udsg->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>MediumWP ) h2_MediumBTaggingEff_Num_udsg->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
        if( btag_dsc>TightWP ) h2_TightBTaggingEff_Num_udsg->Fill(corrJets[ijet].pt(), corrJets[ijet].eta());
      }
    }

  }
  file->Close();
// }

//save to the file
  int nTrial = 0;
  TFile *ofile=TFile::Open(outUrl, "recreate");
  while( !ofile->IsOpen() || ofile->IsZombie() ){
    if(nTrial > 3){
      printf("Output file open failed!");
      return -1;
    }
    nTrial++;
    usleep(1000000*nTrial);
    ofile=TFile::Open(outUrl, "update");
  }
  mon.Write();
  ofile->Close();


}
