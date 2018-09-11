/*************************************************************

Class Usage:

This class should only be used for upgrading and downgrading 
if a single operating point is used in an analysis. 

bool isBTagged = b-tag flag for jet
int pdgIdPart = parton id
float Btag_SF = MC/data scale factor for b/c-tagging efficiency
float Btag_eff = b/c-tagging efficiency in data
float Bmistag_SF = MC/data scale factor for mistag efficiency
float Bmistag_eff = mistag efficiency in data

Author: Michael Segala
Contact: michael.segala@gmail.com
Updated: Ulrich Heintz 12/23/2011

v 1.1

*************************************************************/


#ifndef BTagSFUtil_lite_h
#define BTagSFUtil_lite_h

#include <Riostream.h>
#include "TRandom3.h"
#include "TMath.h"


class BTagSFUtil{

 public:
    
  BTagSFUtil( int seed=0 );
  ~BTagSFUtil();
    
  void SetSeed(int seed=0);
  float getBTagEff(double pt, TString key);
  void modifyBTagsWithSF( bool& isBTagged, float Btag_SF = 0.98, float Btag_eff = 1.0);


 private:
  
  bool applySF(bool& isBTagged, float Btag_SF = 0.98, float Btag_eff = 1.0);
  
  TRandom3* rand_;

};


BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom3(seed);

}

BTagSFUtil::~BTagSFUtil() {

  delete rand_;

}

void BTagSFUtil::SetSeed( int seed ) {

  rand_ = new TRandom3(seed);

}

float BTagSFUtil::getBTagEff(double pt, TString key) {
  //  AN2017_018_v9: Appendix C -> Parameterization functions for the efficiency of the DeepCSV as a functions of the jet pt
  double btag_eff = 1.;

  if (key=="bLOOSE") {
    if (pt>20 && pt<100) btag_eff=0.491+0.0191*pt-0.0004172*pt*pt+4.893*0.000001*pt*pt*pt-3.266*0.00000001*pt*pt*pt*pt+
		      1.238*0.0000000001*pt*pt*pt*pt*pt-2.474*0.0000000000001*pt*pt*pt*pt*pt*pt+2.021*0.0000000000000001*pt*pt*pt*pt*pt*pt*pt;
    else if (pt>=100 && pt<300) btag_eff=0.912-0.0001846*pt+2.479*0.00001*pt*pt-1.417*0.0000001*pt*pt*pt+3.617*0.0000000001*pt*pt*pt*pt-
				  3.433*0.0000000000001*pt*pt*pt*pt*pt;
    else if (pt>=300 && pt<1000) btag_eff=0.892-0.00014*pt+1.011*0.0000001*pt*pt;
  } else if (key=="cLOOSE") {
    if (pt>20 && pt<220) btag_eff=0.421-0.000563*pt+5.096*0.000001*pt*pt-1.4*0.00000001*pt*pt*pt+1.709*0.00000000001*pt*pt*pt*pt-7.627*0.000000000000001*pt*pt*pt*pt*pt;
    else if (pt>=220 && pt<1000) btag_eff=0.383+0.000217*pt-1.997*0.00000001*pt*pt;
  } else if (key=="lLOOSE") {
    if (pt>20 && pt<250) btag_eff=0.2407-0.00593*pt+8.5*0.00001*pt*pt-5.658*0.0000001*pt*pt*pt+1.828*0.000000001*pt*pt*pt*pt-2.287*0.000000000001*pt*pt*pt*pt*pt;
    else if (pt>=250 && pt<1000) btag_eff=0.0541+0.00036*pt-7.392*0.00000001*pt*pt;
  }

  return btag_eff;

}

void BTagSFUtil::modifyBTagsWithSF(bool& isBTagged, float tag_SF, float tag_Eff){
  bool newBTag = isBTagged;
  newBTag = applySF(isBTagged, tag_SF, tag_Eff);
  isBTagged = newBTag;
}


bool BTagSFUtil::applySF(bool& isBTagged, float Btag_SF, float Btag_eff){
  
  bool newBTag = isBTagged;

  if (Btag_SF == 1) return newBTag; //no correction needed 

  //throw die
  float coin = rand_->Uniform(1.);    
  
  if(Btag_SF > 1){  // use this if SF>1
    if( !isBTagged ) {
      //fraction of jets that need to be upgraded
      float mistagPercent = (1.0 - Btag_SF) / (1.0 - (Btag_SF/Btag_eff) );
      //upgrade to tagged
      if( coin < mistagPercent ) {newBTag = true;}
    }
  }else{  // use this if SF<1      
    //downgrade tagged to untagged
    if( isBTagged && coin > Btag_SF ) {newBTag = false;}
  }
  return newBTag;
}


#endif
