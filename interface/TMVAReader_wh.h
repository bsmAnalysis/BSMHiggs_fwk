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
  float Njets=0;
  float HTjets=0;
  float H_pt=0;
  float W_pt=0;
  float lepton_pt=0;
  float bjet1_pt=0;
  float MET=0;
  float H_M=0;
  float MTW=0;
  float dphiWH=0;
  float dphiMetJetMin=0;
  float btag1=0;
  float btag2=0;
  float btag3=0;
  TMVA::Reader *myreader;

  void InitTMVAReader();
  void SetupMVAReader( std::string methodName, std::string modelPath );
  float GenReMVAReader(float thisNjets, float thisHTjets, float thisH_pt, float thisW_pt, float thislepton_pt, float thisbjet1_pt, float thisMET, float thisH_M, float thisMTW, float thisdphiWH, float thisdphiMetJeMin, float thisbtag1, float thisbtag2, float thisbtag3, std::string methodName);

  void CloseMVAReader();

 private:

};
#endif
