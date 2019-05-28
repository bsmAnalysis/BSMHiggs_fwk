//
//  BSMPhysicsEvent.h
//

#ifndef ____BSMPhysicsEvent__
#define ____BSMPhysicsEvent__

//#define YEAR_2017

#include <stdio.h>
#include <vector>
#include "TH2F.h"

#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "DataFormats/Math/interface/deltaR.h"

enum PhysicsObjects   { MET=0,JET=1,TOP=6,ELECTRON=11, MUON=13, TAU=15, GLUON=21, PHOTON=22, Z=23, W=24};
enum LeptonChannels { UNKNOWN=0,MUMU=1,MU=2,EE=3,E=4,EMU=5,ETAU=6,MUTAU=7, GAMMA=22};

class PhysicsObject : public LorentzVector {
public :
 PhysicsObject(LorentzVector vec, Int_t id_, Int_t momid_, Int_t momidx_, Int_t status_):
  LorentzVector(vec), id(id_), momid(momid_), momidx(momidx_), status(status_) { }
    Int_t id;
    Int_t momid;
    Int_t momidx;
    Int_t status;
};


//
class PhysicsObject_Lepton : public LorentzVector {
public :
    PhysicsObject_Lepton(LorentzVector vec, Int_t id_):
        LorentzVector(vec), id(id_) { }
    void setLeptonIDInfo(bool isLooseMu_, bool isTightMu_, bool isMediumMu_, bool isSoftMu_, bool isHighPtMu_,
			 bool passIdMu_, bool passIdLooseMu_, bool passSoftMuon_, 
                         bool isElpassVeto_, bool isElpassLoose_, bool isElpassMedium_, bool isElpassTight_,
			 bool passIdEl_, bool passIdLooseEl_, 
			 bool isTauDM_) {
        isLooseMu = isLooseMu_;
        isTightMu = isTightMu_;
	isMediumMu = isMediumMu_;
        isSoftMu = isSoftMu_;
        isHighPtMu = isHighPtMu_;

	passIdMu = passIdMu_;
	passIdLooseMu = passIdLooseMu_;
	passSoftMuon = passSoftMuon_;
	//passIsoMu = passIsoMu_;

        isElpassVeto = isElpassVeto_;
        isElpassLoose = isElpassLoose_;
        isElpassMedium = isElpassMedium_;
        isElpassTight = isElpassTight_;

	passIdEl = passIdEl_;
	passIdLooseEl = passIdLooseEl_;
	//	passIsoEl = passIsoEl_;

	isTauDM = isTauDM_;
    }

    void setLeptonVar(float en_cor_en_, float en_EtaSC_, float en_R9_, int en_gainSeed_,
		      float mn_validMuonHits_, float mn_trkLayersWithMeasurement_, float mn_pixelLayersWithMeasurement_) {

      en_cor_en = en_cor_en_;
      en_EtaSC = en_EtaSC_;
      en_R9 = en_R9_;
      en_gainSeed = en_gainSeed_;

      mn_validMuonHits = mn_validMuonHits_;
      mn_trkLayersWithMeasurement = mn_trkLayersWithMeasurement_;
      mn_pixelLayersWithMeasurement = mn_pixelLayersWithMeasurement_;

    }

    void setLeptonIsoInfo(float mn_pileupIso_, float mn_chargedIso_, float mn_photonIso_, float mn_neutralHadIso_, bool passIsoMu_, float mn_relIso_, float mn_trkrelIso_,
                          float en_pileupIso_, float en_chargedIso_, float en_photonIso_, float en_neutralHadIso_, float en_relIsoWithEA_, bool passIsoEl_,
			  bool ta_IsLooseIso_, bool ta_IsMediumIso_, bool ta_IsTightIso_, float en_relIso_ ) {

        mn_pileupIso = mn_pileupIso_;
        mn_chargedIso = mn_chargedIso_;
        mn_photonIso = mn_photonIso_;
        mn_neutralHadIso = mn_neutralHadIso_;

	passIsoMu = passIsoMu_;          

        en_pileupIso = en_pileupIso_;
        en_chargedIso = en_chargedIso_;
        en_photonIso = en_photonIso_;
        en_neutralHadIso = en_neutralHadIso_;
        en_relIsoWithEA = en_relIsoWithEA_;

	passIsoEl = passIsoEl_; 

        mn_relIso = mn_relIso_;
        en_relIso = en_relIso_;
        mn_trkrelIso = mn_trkrelIso_;

	ta_IsLooseIso = ta_IsLooseIso_;
	ta_IsMediumIso = ta_IsMediumIso_;
	ta_IsTightIso = ta_IsTightIso_;
    }

#ifdef YEAR_2017
    void setLeptonScaleFac(float en_enSmearNrSigma_, float en_enScaleValue_,
                           float en_enScaleStatUp_, float en_enScaleStatDown_, float en_enScaleSystUp_, float en_enScaleSystDown_, float en_enScaleGainUp_, float en_enScaleGainDown_, float en_enSigmaRhoUp_,  float en_enSigmaRhoDown_, float en_enSigmaPhiDown_) {
    en_enSmearNrSigma = en_enSmearNrSigma_;
    en_enScaleValue   = en_enScaleValue_;
    en_enScaleStatUp  = en_enScaleStatUp_;
    en_enScaleStatDown= en_enScaleStatDown_;
    en_enScaleSystUp  = en_enScaleSystUp_;
    en_enScaleSystDown= en_enScaleSystDown_;
    en_enScaleGainUp  = en_enScaleGainUp_;
    en_enScaleGainDown= en_enScaleGainDown_;
    en_enSigmaRhoUp   = en_enSigmaRhoUp_;
    en_enSigmaRhoDown = en_enSigmaRhoDown_;
    en_enSigmaPhiDown = en_enSigmaPhiDown_;
    }
#endif

    float e_pfRelIsoDbeta() {
        return (en_chargedIso + TMath::Max(0., en_neutralHadIso + en_photonIso - 0.5 * en_pileupIso) )/pt();
    }

    float e_relIsoWithEA() {
	return en_relIsoWithEA;
    }

    float m_pfRelIsoDbeta() {
        return (mn_chargedIso + TMath::Max(0., mn_neutralHadIso + mn_photonIso - 0.5 * mn_pileupIso) )/pt();
    }

    Int_t id;

    float en_relIso, mn_relIso, mn_trkrelIso;
    bool isLooseMu, isTightMu, isMediumMu, isSoftMu, isHighPtMu;
    bool passIdMu, passIdLooseMu, passSoftMuon, passIsoMu;

    int en_gainSeed;
    float en_cor_en, en_EtaSC, en_R9;
    float mn_validMuonHits, mn_trkLayersWithMeasurement, mn_pixelLayersWithMeasurement;

    bool isElpassVeto, isElpassLoose, isElpassMedium, isElpassTight;
    bool passIdEl, passIdLooseEl, passIsoEl;

    bool isTauDM;
    float mn_pileupIso, mn_chargedIso, mn_photonIso, mn_neutralHadIso;
    float en_pileupIso, en_chargedIso, en_photonIso, en_neutralHadIso;
    float en_relIsoWithEA;
    bool ta_IsLooseIso, ta_IsMediumIso, ta_IsTightIso;

#ifdef YEAR_2017
    float en_enSmearNrSigma, en_enScaleValue;
    float en_enScaleStatUp, en_enScaleStatDown, en_enScaleSystUp, en_enScaleSystDown, en_enScaleGainUp, en_enScaleGainDown, en_enSigmaRhoUp,  en_enSigmaRhoDown, en_enSigmaPhiDown;
#endif
};



//
class PhysicsObject_Jet : public LorentzVector {
public :
 PhysicsObject_Jet(LorentzVector vec, Float_t pumva_, Bool_t isPFLoose_, Bool_t isPFTight_):
  LorentzVector(vec), pumva(pumva_), isPFLoose(isPFLoose_), isPFTight(isPFTight_) { }

  void setBtagInfo(Float_t btag0_, Float_t btag1_) { //, Float_t btag1_, Float_t btag2_, Float_t btag3_, Float_t btag4_, Float_t btag5_, Float_t btag6_, Float_t btag7_) {
        btag0=btag0_;
        btag1=btag1_;
	/*
        btag1=btag1_;
        btag2=btag2_;
        btag3=btag3_;
        btag4=btag4_;
        btag5=btag5_;
        btag6=btag6_;
        btag7=btag7_;
	*/
    }

    void setGenInfo(Int_t flavid_, Int_t partonid_, Int_t motherid_, Float_t parton_px_, Float_t parton_py_, Float_t parton_pz_, Float_t parton_en_, Float_t genPt_)
    {
	flavid=flavid_;
	partonid=partonid_;
	motherid=motherid_;
	parton_px=parton_px_;
	parton_py=parton_py_;
	parton_pz=parton_pz_;
	parton_en=parton_en_;
	genPt=genPt_;
    }

    Float_t btag0, btag1; //, btag1, btag2, btag3, btag4, btag5, btag6, btag7;
    Float_t pumva;
    Bool_t isPFLoose,isPFTight;
    Int_t flavid, partonid, motherid;
    Float_t parton_px, parton_py, parton_pz, parton_en;
    Float_t genPt;

};
/*
class PhysicsObject_FatJet : public LorentzVector {
 public:
 PhysicsObject_FatJet(LorentzVector vec) :
  LorentzVector(vec) { }

  void setBtagInfo(Float_t btag0_) {
    btag0=btag0_;
  }

  void setSubjetInfo(Float_t prunedM_, Float_t softdropM_, Float_t tau1_, Float_t tau2_, Float_t tau3_)
  {
    prunedM=prunedM_;
    softdropM=softdropM_;
    tau1=tau1_;
    tau2=tau2_;
    tau3=tau3_;
  }

  void setSubjets(Int_t nSubj_, Float_t subjet_px_[2], Float_t subjet_py_[2], Float_t subjet_pz_[2], Float_t subjet_en_[2])
  {
    nSubj=nSubj_;

    for (int i=0; i<2; i++) {
      LorentzVector subjet;
      subjet.SetPxPyPzE(subjet_px_[i], subjet_py_[i], subjet_pz_[i], subjet_en_[i]);
      if (subjet.Pt()>20) { subjets.push_back(subjet); }
    }

  }

  void setGenInfo(Int_t flavid_, Int_t partonid_, Int_t motherid_, Float_t parton_px_, Float_t parton_py_, Float_t parton_pz_, Float_t parton_en_)
  {
    flavid=flavid_; 
    partonid=partonid_;
    motherid=motherid_;
    parton_px=parton_px_; 
    parton_py=parton_py_;
    parton_pz=parton_pz_;
    parton_en=parton_en_;
  }

  Float_t btag0;
  Float_t prunedM, softdropM;
  Float_t tau1, tau2, tau3;

  Int_t nSubj;
  //  Float_t subjet_px[2], subjet_py[2], subjet_pz[2], subjet_en[2];
  std::vector<LorentzVector> subjets;

  Int_t flavid, partonid, motherid;
  Float_t parton_px, parton_py, parton_pz, parton_en;

};
*/

class PhysicsObject_SV : public LorentzVector {
 public:
 PhysicsObject_SV(LorentzVector vec, Float_t chi2_, Float_t ndof_) :
  LorentzVector(vec), chi2(chi2_), ndof(ndof_) { }

  void setSVinfo(Int_t ntrk_, Float_t dxy_, Float_t dxyz_, Float_t dxyz_signif_, Float_t cos_dxyz_p_)
  {
    ntrk=ntrk_;
    dxy=dxy_;
    dxyz=dxyz_;
    dxyz_signif=dxyz_signif_;
    cos_dxyz_p=cos_dxyz_p_;
  }

  void setGenInfo(Int_t sv_mc_nbh_moms_, Int_t sv_mc_nbh_daus_, Int_t sv_mc_mcbh_ind_)
  {
    sv_mc_nbh_moms=sv_mc_nbh_moms_;
    sv_mc_nbh_daus=sv_mc_nbh_daus_;
    sv_mc_mcbh_ind=sv_mc_mcbh_ind_;
  }
  
  Float_t chi2, ndof;

  Int_t ntrk;
  Float_t dxy, dxyz, dxyz_signif;
  Float_t cos_dxyz_p;

  Int_t sv_mc_nbh_moms, sv_mc_nbh_daus, sv_mc_mcbh_ind;

};


typedef std::vector<PhysicsObject>        PhysicsObjectCollection;
typedef std::vector<PhysicsObject_Lepton> PhysicsObjectLeptonCollection;
typedef std::vector<PhysicsObject_Jet>    PhysicsObjectJetCollection;
//typedef std::vector<PhysicsObject_FatJet> PhysicsObjectFatJetCollection;
typedef std::vector<PhysicsObject_SV>     PhysicsObjectSVCollection;


//
struct PhysicsEvent_t {
  int run,event,lumi;
  int nvtx;
  
  PhysicsObjectLeptonCollection leptons;
  PhysicsObjectJetCollection jets;
  //  PhysicsObjectFatJetCollection fatjets;
  PhysicsObjectSVCollection svs;
  LorentzVector met, metNoHF;
  LorentzVector imet[11];
  LorentzVectorCollection variedMet;

  PhysicsObjectCollection genparticles;
  PhysicsObjectCollection genneutrinos,genleptons,genHiggs,genpartons;
  PhysicsObjectCollection genjets;
};



//
PhysicsEvent_t getPhysicsEventFrom(DataEvtSummary_t &ev);
int getLeptonId(int id);
int getDileptonId(int id1, int id2);

bool isDYToLL(int id1, int id2);
bool isDYToTauTau(int id1, int id2);

float getSFfrom1DHist(double xval, TH1F* h_);
float getSFfrom2DHist(double xval, double yval, TH2F* h_);

float getNLOEWKZZWeight(double trailing_pt);
float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ);



#endif /* defined(____BSMPhysicsEvent__) */
