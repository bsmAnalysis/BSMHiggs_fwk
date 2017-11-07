#ifndef dataevtsummaryhandler_h
#define dataevtsummaryhandler_h

#include <string.h>
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
};

class MVAHandler 
{
 public:
  //
  MVAHandler();
  ~MVAHandler();

  //current event
  MVAEvtContainer evSummary_;
  MVAEvtContainer &getEvent() 
  {
    return evSummary_;
  }

  //read mode, from calculated var
  void getEntry(int ientry);
  //{
    //resetStruct();
    
    //if(ti) ti->GetEntry(ientry);
    //return ;
  //}
  void resetStruct();

  //write mode, to mva tree
  bool initTree(TTree *t3b, TTree *t4b);
  void fillTree();

 private:
  //the tree, 2 for 3b 4b separately
  TTree *to3b, *to4b;
};
#endif
