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
  bool is3b = false, is4b = false; 
  //W boson related only related var
  float WpT = -1.0;
  //Higgs boson related only var
  float Hmass = -1.0, HpT = -1.0;
  float bbdRAve = -1.0, bbdMMin = -1.0;
  float HHt = -1.0;
  //dr W and Higgs 
  float WHdR = -1.0;
  float lepPt = -1.0;
  float pfMET = -1.0; float MTw = -1.0;
  float ljDR = -1.0;
  //weight
  float weight = -1.0;
  //AUX
  int lheNJets = -1;
};

class MVAHandler 
{
 public:
  //
  MVAHandler();
  ~MVAHandler();

  //current event
  MVAEvtContainer evSummary_;
  MVAEvtContainer &getEvent();

  //read mode, from calculated var
  void resetStruct();
  void getEntry(
                bool is3b, bool is4b,
                float Wpt, //W only
                float Hmass, float HpT, float bbdRAve, float bbdMMin, float HHt, //Higgs only
                float WHdR, //W and H
		float lepPt, float pfMET, float MTw,
		float ljDR,
                float weight,
                int lheNJets
               );

  //write mode, to mva tree
  TFile* MVAofile;
  //the tree, 2 for 3b 4b separately
  TTree *to3b, *to4b, *toSignal;
  bool initTree(TString mvaout);
  void fillTree();
  void writeTree();
 private:
};
#endif
