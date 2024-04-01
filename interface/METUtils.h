//
//  METUtils.h
//

#ifndef ____METUtils__
#define ____METUtils__

#include <utility>
#include <vector>

#include "Math/LorentzVector.h"
#include "TVector2.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "JetMETCorrections/Modules/src/JetResolution.cc"

#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;


namespace METUtils {

  double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass);
  double response(LorentzVector &Z, LorentzVector &MET);
  PhysicsObject_Jet smearedJet(const LorentzVector &origJet, double genJetPt, Int_t yearBits, int mode=0);
  PhysicsObject_Jet smearedJet(JME::JetResolutionScaleFactor &jer_sf, const LorentzVector &origJet, double genJetPt, Int_t yearBits, int mode=0);

  enum UncertaintyVariations { JER, JER_UP, JER_DOWN,UMET_UP,UMET_DOWN,LES_UP,LES_DOWN}; 
  //  enum UncertaintyVariations { JER, JER_UP, JER_DOWN, JES_UP, JES_DOWN,UMET_UP,UMET_DOWN,LES_UP,LES_DOWN};
  void computeJetVariation(PhysicsObjectJetCollection& jets,
			   PhysicsObjectLeptonCollection &leptons,       
			   std::vector<PhysicsObjectJetCollection>& jetsVar, 
			   std::vector<JetCorrectionUncertainty*>& jecUnc,
			   Int_t yearBits);
  void computeJetVariation(JME::JetResolutionScaleFactor &jer_sf,
			   PhysicsObjectJetCollection& jets,
			   PhysicsObjectLeptonCollection &leptons,       
			   std::vector<PhysicsObjectJetCollection>& jetsVar, 
			   std::vector< std::vector<double> >& variedJECSFs, 
			   std::vector<JetCorrectionUncertainty*>& jecUnc,
			   Int_t yearBits);

  void computeVariation(PhysicsObjectJetCollection& jets,
                        PhysicsObjectLeptonCollection &leptons,
                        LorentzVector& met,
                        std::vector<PhysicsObjectJetCollection>& jetsVar,
                        LorentzVectorCollection& metsVar,
			std::vector<JetCorrectionUncertainty*>& jecUnc,
			Int_t yearBits);
			//                        JetCorrectionUncertainty *jecUnc);
  LorentzVector applyMETXYCorr(LorentzVector met, bool isMC, int nvtx);
}
#endif /* defined(____METUtils__) */
