#ifndef dataevtsummaryhandler_h
#define dataevtsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

//#define YEAR_2017

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;

const Int_t MAXSB = 2;

//#define MAXSB 2
#define MAXPARTICLES 50
#define MAXMCPARTICLES 250
#define MAXLHEWEIGHTS 500

struct DataEvtSummary_t {

    Int_t run,lumi;
    Long64_t event;
    //Float_t puWeight;
    Float_t curAvgInstLumi,curIntegLumi;
    Bool_t hasTrigger;
    Int_t triggerType;

    //primary vertex
    Int_t nvtx;
    Float_t vtx_x, vtx_y, vtx_z;
    //Float_t fixedGridRhoAll,fixedGridRhoFastjetAll,fixedGridRhoFastjetAllCalo,fixedGridRhoFastjetCentralCalo,fixedGridRhoFastjetCentralChargedPileUp,fixedGridRhoFastjetCentralNeutral;

    //generator info
    Int_t ngenITpu,ngenOOTpu,ngenOOTpum1, ngenTruepu;
    Float_t pthat,genWeight, qscale, x1,x2;
    Int_t id1,id2;
    /*
    Float_t weight_QCDscale_muR1_muF1;
    Float_t weight_QCDscale_muR1_muF2;
    Float_t weight_QCDscale_muR1_muF0p5;
    Float_t weight_QCDscale_muR2_muF1;
    Float_t weight_QCDscale_muR2_muF2;
    Float_t weight_QCDscale_muR2_muF0p5;
    Float_t weight_QCDscale_muR0p5_muF1;
    Float_t weight_QCDscale_muR0p5_muF2;
    Float_t weight_QCDscale_muR0p5_muF0p5;
    */

    Int_t   npdfs;
    Float_t pdfWeights[MAXLHEWEIGHTS];
    Int_t   nalphaS;
    Float_t alphaSWeights[MAXLHEWEIGHTS];
    //lhe njet
    Int_t lheNJets;
  Float_t lheHt;

    //gen level event
    Int_t nmcparticles, nmcjparticles;
    Float_t mc_px[MAXMCPARTICLES],mc_py[MAXMCPARTICLES],mc_pz[MAXMCPARTICLES],mc_en[MAXMCPARTICLES];
    Int_t mc_id[MAXMCPARTICLES], mc_status[MAXMCPARTICLES], mc_mom[MAXMCPARTICLES],  mc_momidx[MAXMCPARTICLES];
    Float_t mcj_px[MAXMCPARTICLES],mcj_py[MAXMCPARTICLES],mcj_pz[MAXMCPARTICLES],mcj_en[MAXMCPARTICLES];  
    Int_t mcj_id[MAXMCPARTICLES], mcj_status[MAXMCPARTICLES], mcj_mom[MAXMCPARTICLES];

    //gen ground state B hadrons
    Int_t   mcbh ;
    Float_t mcbh_px[MAXMCPARTICLES], mcbh_py[MAXMCPARTICLES], mcbh_pz[MAXMCPARTICLES], mcbh_en[MAXMCPARTICLES] ;
    Int_t   mcbh_id[MAXMCPARTICLES] ;

    //muon
    Int_t mn;
    Float_t mn_px[MAXPARTICLES],mn_py[MAXPARTICLES],mn_pz[MAXPARTICLES],mn_en[MAXPARTICLES];
    Int_t mn_id[MAXPARTICLES], mn_type[MAXPARTICLES];

    Float_t mn_d0[MAXPARTICLES],mn_dZ[MAXPARTICLES],mn_ip3d[MAXPARTICLES],mn_ip3dsig[MAXPARTICLES];
    Bool_t mn_IsLoose[MAXPARTICLES],mn_IsMedium[MAXPARTICLES],mn_IsTight[MAXPARTICLES],mn_IsSoft[MAXPARTICLES],mn_IsHighPt[MAXPARTICLES];
    Float_t mn_pileupIsoR03[MAXPARTICLES],mn_chargedIsoR03[MAXPARTICLES],mn_photonIsoR03[MAXPARTICLES],mn_neutralHadIsoR03[MAXPARTICLES];
    Float_t mn_pileupIsoR04[MAXPARTICLES],mn_chargedIsoR04[MAXPARTICLES],mn_photonIsoR04[MAXPARTICLES],mn_neutralHadIsoR04[MAXPARTICLES];

    Bool_t mn_passId[MAXPARTICLES],mn_passIdLoose[MAXPARTICLES],mn_passSoftMuon[MAXPARTICLES],mn_passIso[MAXPARTICLES];
    Float_t mn_nMatches[MAXPARTICLES],mn_nMatchedStations[MAXPARTICLES],mn_validMuonHits[MAXPARTICLES],mn_innerTrackChi2[MAXPARTICLES],mn_trkLayersWithMeasurement[MAXPARTICLES],mn_pixelLayersWithMeasurement[MAXPARTICLES];
    Float_t mn_relIso[MAXPARTICLES], mn_trkrelIso[MAXPARTICLES];

    //electron
    Int_t en;
  Float_t en_px[MAXPARTICLES],en_py[MAXPARTICLES],en_pz[MAXPARTICLES],en_en[MAXPARTICLES],en_cor_en[MAXPARTICLES];
    Int_t en_id[MAXPARTICLES];
  Int_t en_gainSeed[MAXPARTICLES]; 
  //  Float_t en_scale_corr[MAXPARTICLES];
    //Float_t en_d0[MAXPARTICLES],en_dZ[MAXPARTICLES];
  Float_t en_EtaSC[MAXPARTICLES];

  Float_t en_dEtaIn[MAXPARTICLES],en_dPhiIn[MAXPARTICLES],en_hOverE[MAXPARTICLES],en_R9[MAXPARTICLES],en_sigmaIetaIeta[MAXPARTICLES],en_sigmaIetaIeta5x5[MAXPARTICLES],en_ooEmooP[MAXPARTICLES];
    Float_t en_pileupIso[MAXPARTICLES],en_chargedIso[MAXPARTICLES],en_photonIso[MAXPARTICLES],en_neutralHadIso[MAXPARTICLES];
    Float_t en_relIsoWithEA[MAXPARTICLES],en_relIsoWithDBeta[MAXPARTICLES],en_MissingHits[MAXPARTICLES],en_passConversionVeto[MAXPARTICLES];
    Bool_t en_passVeto[MAXPARTICLES],en_passLoose[MAXPARTICLES],en_passMedium[MAXPARTICLES],en_passTight[MAXPARTICLES],en_passHEEP[MAXPARTICLES];
  
    Bool_t en_passMVATrigMedium[MAXPARTICLES], en_passMVATrigTight[MAXPARTICLES];
    Float_t en_IDMVATrigValue[MAXPARTICLES];
    Int_t   en_IDMVATrigCategory[MAXPARTICLES];

    Int_t en_istrue[MAXPARTICLES];
    Bool_t en_passId[MAXPARTICLES],en_passIdLoose[MAXPARTICLES],en_passIso[MAXPARTICLES];
    Float_t en_relIso[MAXPARTICLES];

#ifdef YEAR_2017
    Float_t en_enSmearNrSigma[MAXPARTICLES],en_enScaleValue[MAXPARTICLES];
    Float_t en_enScaleStatUp[MAXPARTICLES],en_enScaleStatDown[MAXPARTICLES],en_enScaleSystUp[MAXPARTICLES],en_enScaleSystDown[MAXPARTICLES],en_enScaleGainUp[MAXPARTICLES],en_enScaleGainDown[MAXPARTICLES],en_enSigmaRhoUp[MAXPARTICLES],en_enSigmaRhoDown[MAXPARTICLES],en_enSigmaPhiDown[MAXPARTICLES];
#endif
    //tau
    Int_t ta;
    Float_t ta_px[MAXPARTICLES],ta_py[MAXPARTICLES],ta_pz[MAXPARTICLES],ta_en[MAXPARTICLES];
    Int_t ta_id[MAXPARTICLES];
    Bool_t ta_dm[MAXPARTICLES],ta_newdm[MAXPARTICLES];
    Bool_t ta_IsLooseIso[MAXPARTICLES],ta_IsMediumIso[MAXPARTICLES],ta_IsTightIso[MAXPARTICLES];
    Float_t ta_combIsoDBeta3Hits[MAXPARTICLES],ta_chargedIso[MAXPARTICLES],ta_neutralIso[MAXPARTICLES],ta_pileupIso[MAXPARTICLES];
    Bool_t ta_passEleVetoLoose[MAXPARTICLES],ta_passEleVetoMedium[MAXPARTICLES],ta_passEleVetoTight[MAXPARTICLES],ta_passMuVetoLoose3[MAXPARTICLES],ta_passMuVetoTight3[MAXPARTICLES];

    //jet (ak4PFJetsCHS)
    Int_t jet;
    Float_t jet_px[MAXPARTICLES],jet_py[MAXPARTICLES],jet_pz[MAXPARTICLES],jet_en[MAXPARTICLES];
    Float_t jet_btag0[MAXPARTICLES], jet_btag1[MAXPARTICLES];//,jet_btag1[MAXPARTICLES],jet_btag2[MAXPARTICLES],jet_btag3[MAXPARTICLES];
    //Float_t jet_btag4[MAXPARTICLES],jet_btag5[MAXPARTICLES],jet_btag6[MAXPARTICLES],jet_btag7[MAXPARTICLES];
    //Float_t jet_btag8[MAXPARTICLES],jet_btag9[MAXPARTICLES],jet_btag10[MAXPARTICLES];
    Float_t jet_mass[MAXPARTICLES],jet_area[MAXPARTICLES],jet_pu[MAXPARTICLES],jet_puId[MAXPARTICLES],jet_genpt[MAXPARTICLES];
    Bool_t jet_PFLoose[MAXPARTICLES], jet_PFTight[MAXPARTICLES];
    Int_t jet_partonFlavour[MAXPARTICLES], jet_hadronFlavour[MAXPARTICLES], jet_mother_id[MAXPARTICLES];
    Float_t jet_parton_px[MAXPARTICLES], jet_parton_py[MAXPARTICLES], jet_parton_pz[MAXPARTICLES], jet_parton_en[MAXPARTICLES];


    //sv : Inclusive Secondary Vertices from slimmedSecondaryVertices
    Int_t sv ;
    Float_t sv_px[MAXPARTICLES], sv_py[MAXPARTICLES], sv_pz[MAXPARTICLES], sv_en[MAXPARTICLES] ;
    Int_t   sv_ntrk[MAXPARTICLES] ;
    Float_t sv_dxy[MAXPARTICLES], sv_dxyz[MAXPARTICLES], sv_dxyz_signif[MAXPARTICLES] ;
    Float_t sv_cos_dxyz_p[MAXPARTICLES] ;
    Float_t sv_chi2[MAXPARTICLES], sv_ndof[MAXPARTICLES] ;
    Int_t   sv_mc_nbh_moms[MAXPARTICLES] ;
    Int_t   sv_mc_nbh_daus[MAXPARTICLES] ;
    Int_t   sv_mc_mcbh_ind[MAXPARTICLES] ;

    /*
    //jet (slimmedJetsPuppi)
    Int_t pjet;
    Float_t pjet_px[MAXPARTICLES],pjet_py[MAXPARTICLES],pjet_pz[MAXPARTICLES],pjet_en[MAXPARTICLES];
    Float_t pjet_btag0[MAXPARTICLES],pjet_btag1[MAXPARTICLES],pjet_btag2[MAXPARTICLES],pjet_btag3[MAXPARTICLES];
    Float_t pjet_btag4[MAXPARTICLES],pjet_btag5[MAXPARTICLES],pjet_btag6[MAXPARTICLES],pjet_btag7[MAXPARTICLES];
    Float_t pjet_btag8[MAXPARTICLES],pjet_btag9[MAXPARTICLES],pjet_btag10[MAXPARTICLES];
    Float_t pjet_genpt[MAXPARTICLES];
    */

    //fjet (ak8PFJetsCHS)
  /*
    Int_t fjet;
    Float_t fjet_px[MAXPARTICLES],fjet_py[MAXPARTICLES],fjet_pz[MAXPARTICLES],fjet_en[MAXPARTICLES];
    Float_t fjet_btag0[MAXPARTICLES], fjet_btag1[MAXPARTICLES];
    Float_t fjet_prunedM[MAXPARTICLES], fjet_softdropM[MAXPARTICLES]; //fjet_trimmedM[MAXPARTICLES],fjet_filteredM[MAXPARTICLES];
    Float_t fjet_tau1[MAXPARTICLES],fjet_tau2[MAXPARTICLES],fjet_tau3[MAXPARTICLES];
    Float_t fjet_genpt[MAXPARTICLES];
    Int_t fjet_partonFlavour[MAXPARTICLES], fjet_hadronFlavour[MAXPARTICLES], fjet_mother_id[MAXPARTICLES];
    Float_t fjet_parton_px[MAXPARTICLES], fjet_parton_py[MAXPARTICLES], fjet_parton_pz[MAXPARTICLES], fjet_parton_en[MAXPARTICLES];
    Int_t fjet_subjet_count[MAXPARTICLES];
    Float_t fjet_subjets_px[MAXPARTICLES][MAXSB], fjet_subjets_py[MAXPARTICLES][MAXSB], fjet_subjets_pz[MAXPARTICLES][MAXSB], fjet_subjets_en[MAXPARTICLES][MAXSB];
  */
    //met
  Float_t imet_pt[11], imet_phi[11];
    Float_t met_pt,met_phi,met_sumMET;
    Float_t metNoHF_pt,metNoHF_phi,metNoHF_sumMET;
    Float_t metPuppi_pt,metPuppi_phi,metPuppi_sumMET;
    Float_t rawpfmet_pt,rawpfmet_phi,rawpfmet_sumMET;
    Float_t rawcalomet_pt,rawcalomet_phi,rawcalomet_sumMET;
    /*
    Bool_t flag_HBHENoiseFilter,flag_CSCTightHaloFilter,flag_hcalLaserEventFilter,flag_EcalDeadCellTriggerPrimitiveFilter,flag_goodVertices;
    Bool_t flag_HBHENoiseIsoFilter;
    Bool_t flag_EcalDeadCellBoundaryEnergyFilter;
    Bool_t flag_trackingFailureFilter,flag_eeBadScFilter,flag_ecalLaserCorrFilter,flag_trkPOGFilters,flag_trkPOG_manystripclus53X,flag_trkPOG_toomanystripclus53X;
    Bool_t flag_trkPOG_logErrorTooManyClusters,flag_METFilters;
    */
};

class DataEvtSummaryHandler {
  public:
    //
    DataEvtSummaryHandler();
    ~DataEvtSummaryHandler();

    //current event
    DataEvtSummary_t evSummary_;
    DataEvtSummary_t &getEvent() {
        return evSummary_;
    }

    //write mode
    bool initTree(TTree *t);
    void fillTree();

    //read mode
    bool attachToTree(TTree *t);
    int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
    void getEntry(int ientry) {
        resetStruct();
        if(t_) t_->GetEntry(ientry);
    }

    void resetStruct();

  private:
    //the tree
    TTree *t_;
};

#endif
