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
  float Njets;
  float Nlepton;
  float HT;
  float H_pt;
  float W_pt;
  float lepton_pt;
  float bjet1_pt;
  float MET;
  float H_M;
  float MTW;
  float dphiWH;
  float dphiMetLepton;
  float dphiMetJetMin;
  float DRbbav;
  float DeltaMassbbMin;
  float mbbuntag;
  float btag1;
  float btag2;
  float btag3;
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
  void getEntry(bool isSR1,
		bool isSR2,
		float Njets,
		float Nlepton,
		float HT,
		float H_pt,
		float W_pt,
		float lepton_pt,
		float bjet1_pt,
		float MET,
		float H_M,
		float MTW,
		float dphiWH,
		float dphiMetLepton,
		float dphiMetJetMin,
		float DRbbav,
		float DeltaMassbbMin,
		float mbbuntag,
		float btag1,
		float btag2,
		float btag3,
		float weight
		);
	       
  // Write mode, to mva tree
  
  TFile* MVAofile;
  TTree* tSR0;
  TTree* tSR1;
  TTree* tSR2;
  bool initTree();
  void fillTree();
  void writeTree(TString mvaout);
 private:
};
#endif
