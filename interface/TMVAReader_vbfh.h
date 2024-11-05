#ifndef tmvaReader_h
#define tmvaReader_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <vector>
#include <string>
#include <utility>      // std::pair, std::make_pair

#include "TMVA/Tools.h"
//#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#endif

class TMVAReader
{
 public:
  float higgs_m=0;
  float higgs_pt=0;
  float higgs_eta=0;
  float costheta0=0;
  float dR_higgs_q1=0;
  float dR_higgs_q2=0;
  float phi_qq_higgs=0;
  float avDeltaR=0;
  float minDeltaM=0;
  float MassOfbbj=0;
  float DeltaEta=0;
  float qq_m=0;
  float ProductEta=0;
  float qq_dphi=0;
  float qq_dR=0;
  float alpha_qq=0;
  float jet1_eta=0;
  float jet1_pt=0;
  float jet2_pt=0;
  float jet3_pt=0;
  float jet4_pt=0;
  float jet5_pt=0;
  float N_jet_eta_cut=0;
  float H_T=0;
  float H_z=0;
  float H_Tvec=0;
  float bjet1_pt=0;
  float bjet2_pt=0;
  float bjet3_pt=0;
  float bjet1_eta=0;
  float bjet2_eta=0;
  float bjet3_eta=0;
  float bjet1_btag=0;
  float bjet2_btag=0;
  float bjet3_btag=0;
  float N_bjet=0;
  float qjet1_pt=0;
  float qjet2_pt=0;
  float qjet1_eta=0;
  float qjet2_eta=0;
  float pt_rest=0;
  float E_rest=0;
  float N_untaggedjet=0;
  float MET_pt=0;
  float MET_phi=0;
  
  TMVA::Reader *myreader;

  void InitTMVAReader();
  void SetupMVAReader( std::string methodName, std::string modelPath );
  float GenReMVAReader(float thishiggs_m, float thishiggs_pt, float thishiggs_eta, float thiscostheta0, float thisdR_higgs_q1, float thisdR_higgs_q2, float thisphi_qq_higgs, float thisavDeltaR, float thisminDeltaM, float thisMassOfbbj, float thisDeltaEta, float thisqq_m, float thisProductEta, float thisqq_dphi, float thisqq_dR, float thisalpha_qq, float thisjet1_eta, float thisjet1_pt, float thisjet2_pt, float thisjet3_pt, float thisjet4_pt, float thisjet5_pt, float thisN_jet_eta_cut, float thisH_T, float thisH_z, float thisH_Tvec, float thisbjet1_pt, float thisbjet2_pt, float thisbjet3_pt, float thisbjet1_eta, float thisbjet2_eta, float thisbjet3_eta, float thisbjet1_btag, float thisbjet2_btag, float thisbjet3_btag, float thisN_bjet, float thisqjet1_pt, float thisqjet2_pt, float thisqjet1_eta, float thisqjet2_eta, float thispt_rest, float thisE_rest, float thisN_untaggedjet, float thisMET_pt, float thisMET_phi, std::string methodName);

  void CloseMVAReader();

 private:

};
#endif
