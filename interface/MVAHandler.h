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
  //catagory type
  bool isSR1=false;
  bool isSR2=false;
  bool isSR3=false;
  float m4b=-1;
  float pt4b=-1;
   float m4b_SR1=-1;
  float pt4b_SR1=-1;
  float ptf1=-1;
  float sd_mass1=-1;
  float ht=-1;
  float met=-1;
  float ht_SR1=-1;
  float met_SR1=-1;
  float xbb1=-2;
  float xbbccqq1=-2;
  float ptb1=-1;
  float ptb1_SR1=-1;
  float ptb2=-1;
  float n_ad_j=-1;
  float n_ad_j_SR1=-1;
  float btag1=-1;
  float btag3=-1;
  float drjj=-1;
  float dphi_met_j=-1;
  float dilep_pt=-1;
  float drll=-1;
  float dphiHZ=-1;
  float dphi_met_l=-1;
  float btag3_SR1=-1;
  float drjj_SR1=-1;
  float dphi_met_j_SR1=-1;
  float dilep_pt_SR1=-1;
  float drll_SR1=-1;
  float dphiHZ_SR1=-1;
  float dphi_met_l_SR1=-1;
  float weight=0;
};

class MVAHandler 
{
 public:
  
  MVAHandler();
  ~MVAHandler();

  //current event
  MVAEvtContainer evSummary_;
  MVAEvtContainer &getEvent();

  //read mode, from calculated var
  void resetStruct();
  void getEntry(bool isSR1,
		bool isSR2,
		bool isSR3,
		float m4b,
		float pt4b,
		float m4b_SR1,
		float pt4b_SR1,
		float ptf1,
		float sd_mass1,
		float ht,
		float met,
		float ht_SR1,
		float met_SR1,
		float xbb1,
		float xbbccqq1,
		float ptb1,
		float ptb1_SR1,
		float ptb2,
		float n_ad_j,
		float n_ad_j_SR1,
		float btag1,
		float btag3,
		float drjj,
		float dphi_met_j,
		float dilep_pt,
		float drll,
		float dphiHZ,
		float dphi_met_l,
		float btag3_SR1,
		float drjj_SR1,
		float dphi_met_j_SR1,
		float dilep_pt_SR1,
		float drll_SR1,
		float dphiHZ_SR1,
		float dphi_met_l_SR1,
	        float weight);
	       
  //write mode, to mva tree
  
  TFile* MVAofile;
  TTree *t1;
  TTree *t2;
  TTree *t3;
  bool initTree();
  void fillTree();
  void writeTree(TString mvaout);
 private:
};
#endif
