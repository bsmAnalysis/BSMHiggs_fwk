#ifndef mvahandler_h
#define mvahandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <vector>

#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"

#endif

struct MVAEvtContainer
{
  float higgs_m;
  float higgs_pt;
  float higgs_eta;
  float costheta0;
  float dR_higgs_q1;
  float dR_higgs_q2;
  float phi_qq_higgs;
  float avDeltaR;
  float minDeltaM;
  float MassOfbbj;
  float DeltaEta;
  float qq_m;
  float ProductEta;
  float qq_dR;
  float qq_dphi;
  float alpha_qq;
  float N_jet;
  float N_jet_eta_cut;
  float jet1_pt;
  float jet2_pt;
  float jet3_pt;
  float jet4_pt;
  float jet5_pt;
  float jet6_pt;
  float jet1_eta;
  float jet2_eta;
  float jet3_eta;
  float jet4_eta;
  float jet5_eta;
  float jet6_eta;
  float H_T;
  float H_z;
  float H_Tvec;
  float pt_rest;
  float E_rest;
  float N_bjet;
  float bjet1_pt;
  float bjet2_pt;
  float bjet3_pt;
  float bjet4_pt;
  float bjet1_eta;
  float bjet2_eta;
  float bjet3_eta;
  float bjet4_eta;
  float bjet1_btag;
  float bjet2_btag;
  float bjet3_btag;
  float bjet4_btag;
  float N_qjet;
  float qjet1_pt;
  float qjet2_pt;
  float qjet1_eta;
  float qjet2_eta;
  float N_untaggedjet;
  float MET_pt;
  float MET_phi;
  float weight;
};

class MVAHandler 
{
 public:
  
  MVAHandler();
  ~MVAHandler();

  // Current event
  MVAEvtContainer evSummary_;
  MVAEvtContainer &getEvent();

  // Read mode, from calculated var
  void resetStruct();
  void getEntry(float higgs_m,
		float higgs_pt,
		float higgs_eta,
		float costheta0,
		float dR_higgs_q1,
		float dR_higgs_q2,
		float phi_qq_higgs,
		float avDeltaR,
		float minDeltaM,
		float MassOfbbj,
		float DeltaEta,
		float qq_m,
		float ProductEta,
		float qq_dR,
		float qq_dphi,
		float alpha_qq,
		float N_jet,
		float N_jet_eta_cut,
		float jet1_pt,
		float jet2_pt,
		float jet3_pt,
		float jet4_pt,
		float jet5_pt,
		float jet6_pt,
		float jet1_eta,
		float jet2_eta,
		float jet3_eta,
		float jet4_eta,
		float jet5_eta,
		float jet6_eta,
		float H_T,
		float H_z,
		float H_Tvec,
		float pt_rest,
		float E_rest,
		float N_bjet,
		float bjet1_pt,
		float bjet2_pt,
		float bjet3_pt,
		float bjet4_pt,
		float bjet1_eta,
		float bjet2_eta,
		float bjet3_eta,
		float bjet4_eta,
		float bjet1_btag,
		float bjet2_btag,
		float bjet3_btag,
		float bjet4_btag,
		float N_qjet,
		float qjet1_pt,
		float qjet2_pt,
		float qjet1_eta,
		float qjet2_eta,
		float N_untaggedjet,
		float MET_pt,
		float MET_phi,
		float weight
		);
	       
  // Write mode, to mva tree
  
  TFile* MVAofile;
  TTree* t1;
  bool initTree();
  void fillTree();
  void writeTree(TString mvaout);
 private:
};
#endif
