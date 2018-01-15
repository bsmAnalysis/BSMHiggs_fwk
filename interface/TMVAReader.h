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
  Float_t WpT = 80.0, Hmass = 100.0, HpT = 80.0, bbdRAve = 1.5, bbdMMin = 10.0, HHt = 150.0, WHdR = 3.0;
  TMVA::Reader *myreader;

  void InitTMVAReader();
  void SetupMVAReader( std::string methodName, std::string modelPath );
  float GenReMVAReader(
                       float thisWpT,
                       float thisHmass, float thisHpT, float thisbbdRAve, float thisbbdMMin, float thisHHt, 
                       float thisWHdR,
                       std::string methodName
                      );
  void CloseMVAReader();

 private:

};
#endif
