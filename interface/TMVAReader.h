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
  float m4b=120;
  float pt4b=20;
  float m4b_SR1=120;
  float pt4b_SR1=20;
  float ptf1=200;
  float sd_mass1=20;
  float ht=300;
  float met=50;
  float ht_SR1=300;
  float met_SR1=50;
  float xbb1=0;
  float xbbccqq1=0.5;
  float ptb1=40;
  float ptb2=30;
  float n_ad_j=2;
  float ptb1_SR1=40;
 
  float n_ad_j_SR1=2;
  float btag1=0.5;
  float btag3=0.7;
  float drjj=0.8;
  float dphi_met_j=0;
  float dilep_pt=80;
  float drll=0.8;
  float dphiHZ=0;
  float dphi_met_l=0;
  float btag3_SR1=0.7;
  float drjj_SR1=0.8;
  float dphi_met_j_SR1=0;
  float dilep_pt_SR1=80;
  float drll_SR1=0.8;
  float dphiHZ_SR1=0;
  float dphi_met_l_SR1=0;
         
  TMVA::Reader *myreader;

  void InitTMVAReader();
  void SetupMVAReader( std::string methodName, std::string modelPath );
  float GenReMVAReader(
		       float thisdilep_pt, float thisdrll,float thisdhiHZ,float thisdphi_met_l,  float thisdphi_met_j,
		       float thism4b,float thispt4b,float thismet,  float thisht,
		       float thisptf1,float thissd_mass1,float thisxbb1, float thisxbbccqq1 ,
		       float thisdrjj,float thisn_ad_j,float thisptb1,float thisbtag3,
		       std::string methodName);
  float GenReMVAReader(
                       
                       float thism4b,float thispt4b,float thismet,  float thisht,
                       float thisptf1,float thissd_mass1,float thisxbb1, float thisxbbccqq1 ,
                       float thisdrjj,float thisn_ad_j,float thisptb1,float thisbtag3,
                       std::string methodName);
  float GenReMVAReader(
				   float thisdilep_pt, float thisdrll,float thisdhiHZ,float thisdphi_met_l,  float thisdphi_met_j,
				   float thism4b,float thispt4b,float thismet,  float thisht,float thisdrjj,
				   float thisn_ad_j,float thisptb1,float thisptb2,float thisbtag1,float thisbtag3,
				   std::string methodName
				   );
  float GenReMVAReader(
		       
		       float thism4b,float thispt4b,float thismet,  float thisht,float thisdrjj,
		       float thisn_ad_j,float thisptb1,float thisptb2,float thisbtag1,float thisbtag3,
		       std::string methodName
		       );

  float GenReMVAReader(
				   float thisdilep_pt, float thisdrll,float thisdhiHZ,float thisdphi_met_l,  float thisdphi_met_j,
				   float thism4b,float thispt4b,float thismet,  float thisht,
				   float thisn_ad_j,float thisptb1,float thisptb2,float thisbtag1,float thisbtag3,
				   std::string methodName
				   );
  float GenReMVAReader(
		       float thism4b,float thispt4b,float thismet,  float thisht,
		       float thisn_ad_j,float thisptb1,float thisptb2,float thisbtag1,float thisbtag3,
		       std::string methodName
		       );
  void CloseMVAReader();

 private:

};
#endif
