#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#define YEAR_2017
using namespace std;

//
DataEvtSummaryHandler::DataEvtSummaryHandler()
{
}

//
bool DataEvtSummaryHandler::initTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;

    //event info
    //    t_->Branch("run",           &evSummary_.run,            "run/I");
    //  t_->Branch("lumi",          &evSummary_.lumi,           "lumi/I");
    //  t_->Branch("event",         &evSummary_.event,          "event/L");
    //    t_->Branch("puWeight",      &evSummary_.puWeight,       "puWeight/F"); 
    /* tmp
    t_->Branch("curAvgInstLumi",&evSummary_.curAvgInstLumi, "curAvgInstLumi/F");
    t_->Branch("curIntegLumi",  &evSummary_.curIntegLumi,   "curIntegLumi/F");
    */
    t_->Branch("hasTrigger",    &evSummary_.hasTrigger,     "hasTrigger/O");
    t_->Branch("triggerType",   &evSummary_.triggerType,    "triggerType/I");

    //primary vertex
    //   t_->Branch("nvtx",          &evSummary_.nvtx,           "nvtx/I");
    //  t_->Branch("vtx_x",         &evSummary_.vtx_x,          "vtx_x/F");
    // t_->Branch("vtx_y",         &evSummary_.vtx_y,          "vtx_y/F");
    // t_->Branch("vtx_z",         &evSummary_.vtx_z,          "vtx_z/F");

    /*
    t_->Branch("fixedGridRhoAll",                           &evSummary_.fixedGridRhoAll,                        "fixedGridRhoAll/F");
    t_->Branch("fixedGridRhoFastjetAll",                    &evSummary_.fixedGridRhoFastjetAll,                 "fixedGridRhoFastjetAll/F");
    t_->Branch("fixedGridRhoFastjetAllCalo",                &evSummary_.fixedGridRhoFastjetAllCalo,             "fixedGridRhoFastjetAllCalo/F");
    t_->Branch("fixedGridRhoFastjetCentralCalo",            &evSummary_.fixedGridRhoFastjetCentralCalo,         "fixedGridRhoFastjetCentralCalo/F");
    t_->Branch("fixedGridRhoFastjetCentralChargedPileUp",   &evSummary_.fixedGridRhoFastjetCentralChargedPileUp,"fixedGridRhoFastjetCentralChargedPileUp/F");
    t_->Branch("fixedGridRhoFastjetCentralNeutral",         &evSummary_.fixedGridRhoFastjetCentralNeutral,      "fixedGridRhoFastjetCentralNeutral/F");
    */

    //generator level info
     
    //    t_->Branch("ngenITpu",      &evSummary_.ngenITpu,       "ngenITpu/I");
    //    t_->Branch("ngenOOTpu",     &evSummary_.ngenOOTpu,      "ngenOOTpu/I");
    //    t_->Branch("ngenOOTpum1",   &evSummary_.ngenOOTpum1,    "ngenOOTpum1/I");
    //    t_->Branch("ngenTruepu",    &evSummary_.ngenTruepu,     "ngenTruepu/I");
    /*r
    t_->Branch("pthat",         &evSummary_.pthat,          "pthat/F");
    t_->Branch("genWeight",     &evSummary_.genWeight,      "genWeight/F");
    t_->Branch("qscale",        &evSummary_.qscale,         "qscale/F");
    t_->Branch("x1",            &evSummary_.x1,             "x1/F");
    t_->Branch("x2",            &evSummary_.x2,             "x2/F");
    t_->Branch("id1",           &evSummary_.id1,            "id1/I");
    t_->Branch("id2",           &evSummary_.id2,            "id2/I");
	*/
    t_->Branch("lheHt",      &evSummary_.lheHt,       "lheHt/F");     
    t_->Branch("lheNJets",      &evSummary_.lheNJets,       "lheNJets/I");
    /*
    t_->Branch("weight_QCDscale_muR1_muF1",     &evSummary_.weight_QCDscale_muR1_muF1,      "weight_QCDscale_muR1_muF1/F");
    t_->Branch("weight_QCDscale_muR1_muF2",     &evSummary_.weight_QCDscale_muR1_muF2,      "weight_QCDscale_muR1_muF2/F");
    t_->Branch("weight_QCDscale_muR1_muF0p5",     &evSummary_.weight_QCDscale_muR1_muF0p5,      "weight_QCDscale_muR1_muF0p5/F");
    t_->Branch("weight_QCDscale_muR2_muF1",     &evSummary_.weight_QCDscale_muR2_muF1,      "weight_QCDscale_muR2_muF1/F");
    t_->Branch("weight_QCDscale_muR2_muF2",     &evSummary_.weight_QCDscale_muR2_muF2,      "weight_QCDscale_muR2_muF2/F");
    t_->Branch("weight_QCDscale_muR2_muF0p5",     &evSummary_.weight_QCDscale_muR2_muF0p5,      "weight_QCDscale_muR2_muF0p5/F");
    t_->Branch("weight_QCDscale_muR0p5_muF1",     &evSummary_.weight_QCDscale_muR0p5_muF1,      "weight_QCDscale_muR0p5_muF1/F");
    t_->Branch("weight_QCDscale_muR0p5_muF2",     &evSummary_.weight_QCDscale_muR0p5_muF2,      "weight_QCDscale_muR0p5_muF2/F");
    t_->Branch("weight_QCDscale_muR0p5_muF0p5",     &evSummary_.weight_QCDscale_muR0p5_muF0p5,      "weight_QCDscale_muR0p5_muF0p5/F");
    */

    /* t_->Branch("npdfs", &evSummary_.npdfs, "npdfs/I");
    t_->Branch("pdfWeights", evSummary_.pdfWeights, "pdfWeights[npdfs]/F");
    t_->Branch("nalphaS", &evSummary_.nalphaS, "nalphaS/I");
    t_->Branch("alphaSWeights", evSummary_.alphaSWeights, "alphaSWeights[nalphaS]/F"); */
    
    //mc truth
    /*t_->Branch("nmcparticles",  &evSummary_.nmcparticles,   "nmcparticles/I");
    t_->Branch("mc_px",         evSummary_.mc_px,           "mc_px[nmcparticles]/F");
    t_->Branch("mc_py",         evSummary_.mc_py,           "mc_py[nmcparticles]/F");
    t_->Branch("mc_pz",         evSummary_.mc_pz,           "mc_pz[nmcparticles]/F");
    t_->Branch("mc_en",         evSummary_.mc_en,           "mc_en[nmcparticles]/F");
    t_->Branch("mc_id",         evSummary_.mc_id,           "mc_id[nmcparticles]/I");
    t_->Branch("mc_status",     evSummary_.mc_status,       "mc_status[nmcparticles]/I");
    t_->Branch("mc_mom",        evSummary_.mc_mom,          "mc_mom[nmcparticles]/I");
    t_->Branch("mc_momidx",     evSummary_.mc_momidx,       "mc_momidx[nmcparticles]/I");*/
	
    //gen ground state B hadrons
//    t_->Branch("mcbh",          &evSummary_.mcbh,           "mcbh/I");
//    t_->Branch("mcbh_px",       evSummary_.mcbh_px,         "mcbh_px[mcbh]/F");
//    t_->Branch("mcbh_py",       evSummary_.mcbh_py,         "mcbh_py[mcbh]/F");
//    t_->Branch("mcbh_pz",       evSummary_.mcbh_pz,         "mcbh_pz[mcbh]/F");
//    t_->Branch("mcbh_en",       evSummary_.mcbh_en,         "mcbh_en[mcbh]/F");
//    t_->Branch("mcbh_id",       evSummary_.mcbh_id,         "mcbh_en[mcbh]/I");

    /* tmp
    // gen jets
    t_->Branch("nmcjparticles",  &evSummary_.nmcjparticles,   "nmcjparticles/I"); 
    t_->Branch("mcj_px",         evSummary_.mcj_px,           "mcj_px[nmcjparticles]/F"); 
    t_->Branch("mcj_py",         evSummary_.mcj_py,           "mcj_py[nmcjparticles]/F"); 
    t_->Branch("mcj_pz",         evSummary_.mcj_pz,           "mcj_pz[nmcjparticles]/F"); 
    t_->Branch("mcj_en",         evSummary_.mcj_en,           "mcj_en[nmcjparticles]/F"); 
    t_->Branch("mcj_id",         evSummary_.mcj_id,           "mcj_id[nmcjparticles]/I"); 
    t_->Branch("mcj_status",     evSummary_.mcj_status,       "mcj_status[nmcjparticles]/I"); 
    t_->Branch("mcj_mom",        evSummary_.mcj_mom,          "mcj_mom[nmcjparticles]/I"); 
    */
    //muon
    t_->Branch("mn",                    &evSummary_.mn,                     "mn/I");
    t_->Branch("mn_px",                 evSummary_.mn_px,                   "mn_px[mn]/F");
    t_->Branch("mn_py",                 evSummary_.mn_py,                   "mn_py[mn]/F");
    t_->Branch("mn_pz",                 evSummary_.mn_pz,                   "mn_pz[mn]/F");
    t_->Branch("mn_en",                 evSummary_.mn_en,                   "mn_en[mn]/F");
    t_->Branch("mn_id",                 evSummary_.mn_id,                   "mn_id[mn]/I");
    t_->Branch("mn_type",               evSummary_.mn_type,                 "mn_type[mn]/I");
    /* tmp
    t_->Branch("mn_d0",                 evSummary_.mn_d0,                   "mn_d0[mn]/F");
    t_->Branch("mn_dZ",                 evSummary_.mn_dZ,                   "mn_dZ[mn]/F");
    t_->Branch("mn_ip3d",               evSummary_.mn_ip3d,                 "mn_ip3d[mn]/F");
    t_->Branch("mn_ip3dsig",            evSummary_.mn_ip3dsig,              "mn_ip3dsig[mn]/F");
    t_->Branch("mn_IsLoose",            evSummary_.mn_IsLoose,              "mn_IsLoose[mn]/O");
    t_->Branch("mn_IsMedium",           evSummary_.mn_IsMedium,              "mn_IsMedium[mn]/O");
    t_->Branch("mn_IsTight",            evSummary_.mn_IsTight,              "mn_IsTight[mn]/O");
    t_->Branch("mn_IsSoft",             evSummary_.mn_IsSoft,               "mn_IsSoft[mn]/O");
    t_->Branch("mn_IsHighPt",           evSummary_.mn_IsHighPt,             "mn_IsHighPt[mn]/O");

    t_->Branch("mn_pileupIsoR03",       evSummary_.mn_pileupIsoR03,         "mn_pileupIsoR03[mn]/F");
    t_->Branch("mn_chargedIsoR03",      evSummary_.mn_chargedIsoR03,        "mn_chargedIsoR03[mn]/F");
    t_->Branch("mn_photonIsoR03",       evSummary_.mn_photonIsoR03,         "mn_photonIsoR03[mn]/F");
    t_->Branch("mn_neutralHadIsoR03",   evSummary_.mn_neutralHadIsoR03,     "mn_neutralHadIsoR03[mn]/F");
    t_->Branch("mn_pileupIsoR04",       evSummary_.mn_pileupIsoR04,         "mn_pileupIsoR04[mn]/F");
    t_->Branch("mn_chargedIsoR04",      evSummary_.mn_chargedIsoR04,        "mn_chargedIsoR04[mn]/F");
    t_->Branch("mn_photonIsoR04",       evSummary_.mn_photonIsoR04,         "mn_photonIsoR04[mn]/F");
    t_->Branch("mn_neutralHadIsoR04",   evSummary_.mn_neutralHadIsoR04,     "mn_neutralHadIsoR04[mn]/F");
    */
    t_->Branch("mn_passId",             evSummary_.mn_passId,               "mn_passId[mn]/O");
    t_->Branch("mn_passIdLoose",        evSummary_.mn_passIdLoose,          "mn_passIdLoose[mn]/O");
    ///t_->Branch("mn_passSoftMuon",       evSummary_.mn_passSoftMuon,         "mn_passSoftMuon[mn]/O");
    t_->Branch("mn_passIso",            evSummary_.mn_passIso,              "mn_passIso[mn]/O");
    // t_->Branch("mn_relIso",             evSummary_.mn_relIso,               "mn_relIso[mn]/F");
    //t_->Branch("mn_trkrelIso",          evSummary_.mn_trkrelIso,            "mn_trkrelIso[mn]/F");

    //    t_->Branch("mn_nMatches",                   evSummary_.mn_nMatches,                     "mn_nMatches[mn]/F");
    //t_->Branch("mn_nMatchedStations",           evSummary_.mn_nMatchedStations,             "mn_nMatchedStations[mn]/F");
    ///t_->Branch("mn_validMuonHits",              evSummary_.mn_validMuonHits,                "mn_validMuonHits[mn]/F");
    //    t_->Branch("mn_innerTrackChi2",             evSummary_.mn_innerTrackChi2,               "mn_innerTrackChi2[mn]/F");
    ///t_->Branch("mn_trkLayersWithMeasurement",   evSummary_.mn_trkLayersWithMeasurement,     "mn_trkLayersWithMeasurement[mn]/F");
    ///t_->Branch("mn_pixelLayersWithMeasurement", evSummary_.mn_pixelLayersWithMeasurement,   "mn_pixelLayersWithMeasurement[mn]/F");


    //electron
    t_->Branch("en",                    &evSummary_.en,                     "en/I");
    t_->Branch("en_px",                 evSummary_.en_px,                   "en_px[en]/F");
    t_->Branch("en_py",                 evSummary_.en_py,                   "en_py[en]/F");
    t_->Branch("en_pz",                 evSummary_.en_pz,                   "en_pz[en]/F");
    t_->Branch("en_en",                 evSummary_.en_en,                   "en_en[en]/F");
    t_->Branch("en_id",                 evSummary_.en_id,                   "en_id[en]/I");
    ///t_->Branch("en_cor_en",                 evSummary_.en_cor_en,                   "en_cor_en[en]/F"); 
    //    t_->Branch("en_scale_corr",                 evSummary_.en_scale_corr,                   "en_scale_corr[en]/F");  
    ///t_->Branch("en_EtaSC",              evSummary_.en_EtaSC,                "en_EtaSC[en]/F");  
    ///t_->Branch("en_R9",                 evSummary_.en_R9,                   "en_R9[en]/F"); 
    ///t_->Branch("en_gainSeed",                 evSummary_.en_gainSeed,                   "en_gainSeed[en]/I");     
    /*
    t_->Branch("en_d0",                 evSummary_.en_d0,                   "en_d0[en]/F");
    t_->Branch("en_dZ",                 evSummary_.en_dZ,                   "en_dZ[en]/F");
    t_->Branch("en_EtaSC",              evSummary_.en_EtaSC,                "en_EtaSC[en]/F");
    t_->Branch("en_PhiSC",              evSummary_.en_PhiSC,                "en_PhiSC[en]/F");
    t_->Branch("en_EnSC",               evSummary_.en_EnSC,                 "en_EnSC[en]/F");
    t_->Branch("en_dEtaIn",             evSummary_.en_dEtaIn,               "en_dEtaIn[en]/F");
    t_->Branch("en_dPhiIn",             evSummary_.en_dPhiIn,               "en_dPhiIn[en]/F");
    t_->Branch("en_hOverE",             evSummary_.en_hOverE,               "en_hOverE[en]/F");
    t_->Branch("en_R9",                 evSummary_.en_R9,                   "en_R9[en]/F");
    t_->Branch("en_sigmaIetaIeta",      evSummary_.en_sigmaIetaIeta,        "en_sigmaIetaIeta[en]/F");
    t_->Branch("en_sigmaIetaIeta5x5",   evSummary_.en_sigmaIetaIeta5x5,     "en_sigmaIetaIeta5x5[en]/F");
    t_->Branch("en_ooEmooP",            evSummary_.en_ooEmooP,              "en_ooEmooP[en]/F");
    */
    /* tmp
    t_->Branch("en_pileupIso",          evSummary_.en_pileupIso,            "en_pileupIso[en]/F");
    t_->Branch("en_chargedIso",         evSummary_.en_chargedIso,           "en_chargedIso[en]/F");
    t_->Branch("en_photonIso",          evSummary_.en_photonIso,            "en_photonIso[en]/F");
    t_->Branch("en_neutralHadIso",      evSummary_.en_neutralHadIso,        "en_neutralHadIso[en]/F");
    t_->Branch("en_relIsoWithEA",       evSummary_.en_relIsoWithEA,         "en_relIsoWithEA[en]/F");
    t_->Branch("en_relIsoWithDBeta",    evSummary_.en_relIsoWithDBeta,      "en_relIsoWithDBeta[en]/F");
    t_->Branch("en_MissingHits",        evSummary_.en_MissingHits,          "en_MissingHits[en]/F");
    t_->Branch("en_passConversionVeto", evSummary_.en_passConversionVeto,   "en_passConversionVeto[en]/F");
    t_->Branch("en_passVeto",           evSummary_.en_passVeto,             "en_passVeto[en]/O");
    t_->Branch("en_passLoose",          evSummary_.en_passLoose,            "en_passLoose[en]/O");
    t_->Branch("en_passMedium",         evSummary_.en_passMedium,           "en_passMedium[en]/O");
    t_->Branch("en_passTight",          evSummary_.en_passTight,            "en_passTight[en]/O");
    t_->Branch("en_passHEEP",           evSummary_.en_passHEEP,             "en_passHEEP[en]/O");

    t_->Branch("en_passMVATrigMedium",  evSummary_.en_passMVATrigMedium,    "en_passMVATrigMedium[en]/O");
    t_->Branch("en_passMVATrigTight",   evSummary_.en_passMVATrigTight,     "en_passMVATrigTight[en]/O");
    t_->Branch("en_IDMVATrigValue",     evSummary_.en_IDMVATrigValue,       "en_IDMVATrigValue[en]/F");
    t_->Branch("en_IDMVATrigCategory",  evSummary_.en_IDMVATrigCategory,    "en_IDMVATrigCategory[en]/I");
    t_->Branch("en_istrue",             evSummary_.en_istrue,               "en_istrue[en]/I");
    */
    t_->Branch("en_passId",             evSummary_.en_passId,               "en_passId[en]/O");
    t_->Branch("en_passIdLoose",        evSummary_.en_passIdLoose,          "en_passIdLoose[en]/O"); 
    t_->Branch("en_passIso",            evSummary_.en_passIso,              "en_passIso[en]/O"); 
    t_->Branch("en_relIso",             evSummary_.en_relIso,               "en_relIso[en]/F");
/*
#ifdef YEAR_2017
    t_->Branch("en_enSigmaValue",       evSummary_.en_enSigmaValue,         "en_enSigmaValue[en]/F");
    t_->Branch("en_enScaleValue",       evSummary_.en_enScaleValue,         "en_enScaleValue[en]/F");
    t_->Branch("en_enScaleStatUp",      evSummary_.en_enScaleStatUp,        "en_enScaleStatUp[en]/F");
    t_->Branch("en_enScaleStatDown",    evSummary_.en_enScaleStatDown,      "en_enScaleStatDown[en]/F");
    t_->Branch("en_enScaleSystUp",      evSummary_.en_enScaleSystUp,        "en_enScaleSystUp[en]/F");
    t_->Branch("en_enScaleSystDown",    evSummary_.en_enScaleSystDown,      "en_enScaleSystDown[en]/F");
    t_->Branch("en_enScaleGainUp",      evSummary_.en_enScaleGainUp,        "en_enScaleGainUp[en]/F");
    t_->Branch("en_enScaleGainDown",    evSummary_.en_enScaleGainDown,      "en_enScaleGainDown[en]/F");
    t_->Branch("en_enSigmaRhoUp",       evSummary_.en_enSigmaRhoUp,         "en_enSigmaRhoUp[en]/F");
    t_->Branch("en_enSigmaRhoDown",     evSummary_.en_enSigmaRhoDown,       "en_enSigmaRhoDown[en]/F");
    t_->Branch("en_enSigmaPhiDown",     evSummary_.en_enSigmaPhiDown,       "en_enSigmaPhiDown[en]/F");
#endif
*/
    /* tmp
    //tau
    t_->Branch("ta",                    &evSummary_.ta,                     "ta/I");
    t_->Branch("ta_px",                 evSummary_.ta_px,                   "ta_px[ta]/F");
    t_->Branch("ta_py",                 evSummary_.ta_py,                   "ta_py[ta]/F");
    t_->Branch("ta_pz",                 evSummary_.ta_pz,                   "ta_pz[ta]/F");
    t_->Branch("ta_en",                 evSummary_.ta_en,                   "ta_en[ta]/F");
    t_->Branch("ta_id",                 evSummary_.ta_id,                   "ta_id[ta]/I");
    t_->Branch("ta_dm",                 evSummary_.ta_dm,                   "ta_dm[ta]/O");
    t_->Branch("ta_newdm",              evSummary_.ta_newdm,                "ta_newdm[ta]/O");
    t_->Branch("ta_IsLooseIso",         evSummary_.ta_IsLooseIso,           "ta_IsLooseIso[ta]/O");
    t_->Branch("ta_IsMediumIso",        evSummary_.ta_IsMediumIso,          "ta_IsMediumIso[ta]/O");
    t_->Branch("ta_IsTightIso",         evSummary_.ta_IsTightIso,           "ta_IsTightIso[ta]/O");
    t_->Branch("ta_combIsoDBeta3Hits",  evSummary_.ta_combIsoDBeta3Hits,    "ta_combIsoDBeta3Hits[ta]/F");
    t_->Branch("ta_chargedIso",         evSummary_.ta_chargedIso,           "ta_chargedIso[ta]/F");
    t_->Branch("ta_neutralIso",         evSummary_.ta_neutralIso,           "ta_neutralIso[ta]/F");
    t_->Branch("ta_pileupIso",          evSummary_.ta_pileupIso,            "ta_pileupIso[ta]/F");
    t_->Branch("ta_passEleVetoLoose",   evSummary_.ta_passEleVetoLoose,     "ta_passEleVetoLoose[ta]/O");
    t_->Branch("ta_passEleVetoMedium",  evSummary_.ta_passEleVetoMedium,    "ta_passEleVetoMedium[ta]/O");
    t_->Branch("ta_passEleVetoTight",   evSummary_.ta_passEleVetoTight,     "ta_passEleVetoTight[ta]/O");
    t_->Branch("ta_passMuVetoLoose3",   evSummary_.ta_passMuVetoLoose3,     "ta_passMuVetoLoose3[ta]/O");
    t_->Branch("ta_passMuVetoTight3",   evSummary_.ta_passMuVetoTight3,     "ta_passMuVetoTight3[ta]/O");

    */

    //jet (ak4PFJetsCHS)
    t_->Branch("jet",                   &evSummary_.jet,                    "jet/I");
    t_->Branch("jet_px",                evSummary_.jet_px,                  "jet_px[jet]/F");
    t_->Branch("jet_py",                evSummary_.jet_py,                  "jet_py[jet]/F");
    t_->Branch("jet_pz",                evSummary_.jet_pz,                  "jet_pz[jet]/F");
    t_->Branch("jet_en",                evSummary_.jet_en,                  "jet_en[jet]/F");
    //t_->Branch("jet_btag0",             evSummary_.jet_btag0,               "jet_btag0[jet]/F");
    t_->Branch("jet_btag1",             evSummary_.jet_btag1,               "jet_btag1[jet]/F");
    /*
    t_->Branch("jet_btag1",             evSummary_.jet_btag1,               "jet_btag1[jet]/F");
    t_->Branch("jet_btag2",             evSummary_.jet_btag2,               "jet_btag2[jet]/F");
    t_->Branch("jet_btag3",             evSummary_.jet_btag3,               "jet_btag3[jet]/F");
    t_->Branch("jet_btag4",             evSummary_.jet_btag4,               "jet_btag4[jet]/F");
    t_->Branch("jet_btag5",             evSummary_.jet_btag5,               "jet_btag5[jet]/F");
    t_->Branch("jet_btag6",             evSummary_.jet_btag6,               "jet_btag6[jet]/F");
    t_->Branch("jet_btag7",             evSummary_.jet_btag7,               "jet_btag7[jet]/F");
    t_->Branch("jet_btag8",             evSummary_.jet_btag8,               "jet_btag8[jet]/F");
    t_->Branch("jet_btag9",             evSummary_.jet_btag9,               "jet_btag9[jet]/F");
    t_->Branch("jet_btag10",            evSummary_.jet_btag10,              "jet_btag10[jet]/F");
    */
    //t_->Branch("jet_mass",              evSummary_.jet_mass,                "jet_mass[jet]/F");
    /* ///
    t_->Branch("jet_area",              evSummary_.jet_area,                "jet_area[jet]/F");
    t_->Branch("jet_pu",                evSummary_.jet_pu,                  "jet_pu[jet]/F");
    t_->Branch("jet_puId",              evSummary_.jet_puId,                "jet_puId[jet]/F");
	*/
    //t_->Branch("jet_PFLoose",           evSummary_.jet_PFLoose,             "jet_PFLoose[jet]/O");
    t_->Branch("jet_PFTight",           evSummary_.jet_PFTight,             "jet_PFTight[jet]/O");
    
    t_->Branch("jet_partonFlavour",     evSummary_.jet_partonFlavour,       "jet_partonFlavour[jet]/I");
    /* t_->Branch("jet_hadronFlavour",     evSummary_.jet_hadronFlavour,       "jet_hadronFlavour[jet]/I");
    t_->Branch("jet_mother_id",         evSummary_.jet_mother_id,           "jet_mother_id[jet]/I");
    t_->Branch("jet_parton_px",         evSummary_.jet_parton_px,        "jet_parton_px[jet]/F");
    t_->Branch("jet_parton_py",         evSummary_.jet_parton_py,        "jet_parton_py[jet]/F");
    t_->Branch("jet_parton_pz",         evSummary_.jet_parton_pz,        "jet_parton_pz[jet]/F");
    t_->Branch("jet_parton_en",         evSummary_.jet_parton_en,        "jet_parton_en[jet]/F");
    t_->Branch("jet_genpt",             evSummary_.jet_genpt,               "jet_genpt[jet]/F");
    */

    /*
    //sv : Inclusive Secondary Vertices from slimmedSecondaryVertices
    t_->Branch("sv"              , &evSummary_.sv              , "sv/I" ) ;
    t_->Branch("sv_px"           , evSummary_.sv_px            , "sv_px[sv]/F" ) ;
    t_->Branch("sv_py"           , evSummary_.sv_py            , "sv_py[sv]/F" ) ;
    t_->Branch("sv_pz"           , evSummary_.sv_pz            , "sv_pz[sv]/F" ) ;
    t_->Branch("sv_en"           , evSummary_.sv_en            , "sv_en[sv]/F" ) ;
    t_->Branch("sv_ntrk"         , evSummary_.sv_ntrk          , "sv_ntrk[sv]/I" ) ;
    t_->Branch("sv_dxy"          , evSummary_.sv_dxy           , "sv_dxy[sv]/F" ) ;
    t_->Branch("sv_dxyz"         , evSummary_.sv_dxyz          , "sv_dxyz[sv]/F" ) ;
    t_->Branch("sv_dxyz_signif"  , evSummary_.sv_dxyz_signif   , "sv_dxyz_signif[sv]/F" ) ;
    t_->Branch("sv_cos_dxyz_p"   , evSummary_.sv_cos_dxyz_p    , "sv_cos_dxyz_p[sv]/F" ) ;
    t_->Branch("sv_chi2"         , evSummary_.sv_chi2          , "sv_chi2[sv]/F" ) ;
    t_->Branch("sv_ndof"         , evSummary_.sv_ndof          , "sv_ndof[sv]/F" ) ;
    t_->Branch("sv_mc_nbh_moms"  , evSummary_.sv_mc_nbh_moms   , "sv_mc_nbh_moms[sv]/I" ) ;
    t_->Branch("sv_mc_nbh_daus"  , evSummary_.sv_mc_nbh_daus   , "sv_mc_nbh_daus[sv]/I" ) ;
//    t_->Branch("sv_mc_mcbh_ind"  , evSummary_.sv_mc_mcbh_ind   , "sv_mc_mcbh_ind[sv]/I" ) ;
    */

    //fjet (ak8PFJetsCHS)
    //
    t_->Branch("fjet",                    &evSummary_.fjet,                   "fjet/I");
    //t_->Branch("fjet_pt",                 evSummary_.fjet_pt,                 "fjet_pt[fjet]/F");
    t_->Branch("fjet_px",                 evSummary_.fjet_px,                 "fjet_px[fjet]/F");
    t_->Branch("fjet_py",                 evSummary_.fjet_py,                 "fjet_py[fjet]/F");
    t_->Branch("fjet_pz",                 evSummary_.fjet_pz,                 "fjet_pz[fjet]/F");
    t_->Branch("fjet_en",                 evSummary_.fjet_en,                 "fjet_en[fjet]/F");
    /*
    t_->Branch("fjet_btag0",              evSummary_.fjet_btag0,              "fjet_btag0[fjet]/F");  
    t_->Branch("fjet_btag1",              evSummary_.fjet_btag1,              "fjet_btag1[fjet]/F");  
    t_->Branch("fjet_btag2",              evSummary_.fjet_btag2,              "fjet_btag2[fjet]/F");
    t_->Branch("fjet_btag3",              evSummary_.fjet_btag3,              "fjet_btag3[fjet]/F");
    t_->Branch("fjet_btag4",              evSummary_.fjet_btag4,              "fjet_btag4[fjet]/F");
    t_->Branch("fjet_btag5",              evSummary_.fjet_btag5,              "fjet_btag5[fjet]/F");
    t_->Branch("fjet_btag6",              evSummary_.fjet_btag6,              "fjet_btag6[fjet]/F");
    t_->Branch("fjet_btag7",              evSummary_.fjet_btag7,              "fjet_btag7[fjet]/F");
    t_->Branch("fjet_btag8",              evSummary_.fjet_btag8,              "fjet_btag8[fjet]/F");
    t_->Branch("fjet_btag9",              evSummary_.fjet_btag9,              "fjet_btag9[fjet]/F");
    */
    t_->Branch("fjet_btag10",             evSummary_.fjet_btag10,             "fjet_btag10[fjet]/F");
    t_->Branch("fjet_btag11",             evSummary_.fjet_btag11,             "fjet_btag11[fjet]/F");
    t_->Branch("fjet_btag12",             evSummary_.fjet_btag12,             "fjet_btag12[fjet]/F");
    t_->Branch("fjet_btag13",             evSummary_.fjet_btag13,             "fjet_btag13[fjet]/F");
    t_->Branch("fjet_btag14",             evSummary_.fjet_btag14,             "fjet_btag14[fjet]/F");
    t_->Branch("fjet_btag15",             evSummary_.fjet_btag15,             "fjet_btag15[fjet]/F");
    t_->Branch("fjet_btag16",             evSummary_.fjet_btag16,             "fjet_btag16[fjet]/F");
    t_->Branch("fjet_btag17",             evSummary_.fjet_btag17,             "fjet_btag17[fjet]/F");




    // t_->Branch("fjet_genpt",              evSummary_.fjet_genpt,              "fjet_genpt[fjet]/F");
    ///t_->Branch("fjet_prunedM",            evSummary_.fjet_prunedM,            "fjet_prunedM[fjet]/F");
    // t_->Branch("fjet_softdropM",          evSummary_.fjet_softdropM,          "fjet_softdropM[fjet]/F"); 
    // t_->Branch("fjet_trimmeen_dM",         evSummary_.fjet_trimmedM,           "fjet_trimmedM[fjet]/F");
    // t_->Branch("fjet_filteredM",        evSummary_.fjet_filteredM,          "fjet_filteredM[fjet]/F");
    /*t_->Branch("fjet_tau1",               evSummary_.fjet_tau1,               "fjet_tau1[fjet]/F");
    t_->Branch("fjet_tau2",               evSummary_.fjet_tau2,               "fjet_tau2[fjet]/F");
    t_->Branch("fjet_tau3",               evSummary_.fjet_tau3,               "fjet_tau3[fjet]/F");
    t_->Branch("fjet_tau4",               evSummary_.fjet_tau4,               "fjet_tau4[fjet]/F"); 
    
    t_->Branch("fjet_mother_id",          evSummary_.fjet_mother_id,          "fjet_mother_id[fjet]/I"); 
    t_->Branch("fjet_parton_px",          evSummary_.fjet_parton_px,          "fjet_parton_px[fjet]/F");
    t_->Branch("fjet_parton_py",          evSummary_.fjet_parton_py,          "fjet_parton_py[fjet]/F");
    t_->Branch("fjet_parton_pz",          evSummary_.fjet_parton_pz,          "fjet_parton_pz[fjet]/F");
    t_->Branch("fjet_parton_en",          evSummary_.fjet_parton_en,          "fjet_parton_en[fjet]/F");
    t_->Branch("fjet_partonFlavour",      evSummary_.fjet_partonFlavour,      "fjet_partonFlavour[fjet]/I");
    t_->Branch("fjet_hadronFlavour",      evSummary_.fjet_hadronFlavour,      "fjet_hadronFlavour[fjet]/I");*/
    
    // t_->Branch("fjet_chf",  		  evSummary_.fjet_chf,		      "fjet_chf[fjet]/F");
    // t_->Branch("fjet_nhf",                evSummary_.fjet_nhf,                "fjet_nhf[fjet]/F");
    // t_->Branch("fjet_phf",                evSummary_.fjet_phf,                "fjet_phf[fjet]/F");
    // t_->Branch("fjet_muf",                evSummary_.fjet_muf,                "fjet_muf[fjet]/F");
    // t_->Branch("fjet_elf",                evSummary_.fjet_elf,                "fjet_elf[fjet]/F");

    // t_->Branch("fjet_ecfB1N2",	    evSummary_.fjet_ecfB1N2,	       "fjet_ecfB1N2[fjet]/F");
    // t_->Branch("fjet_ecfB1N3",            evSummary_.fjet_ecfB1N3,           "fjet_ecfB1N3[fjet]/F");
    // t_->Branch("fjet_ecfB2N2",            evSummary_.fjet_ecfB2N2,           "fjet_ecfB2N2[fjet]/F");
    // t_->Branch("fjet_ecfB2N3",            evSummary_.fjet_ecfB2N3,           "fjet_ecfB2N3[fjet]/F");

    t_->Branch("fjet_subjet_count",       evSummary_.fjet_subjet_count,       "fjet_subjet_count[fjet]/I");
    /* t_->Branch("fjet_subjets_px",         &evSummary_.fjet_subjets_px,        "fjet_subjets_px[fjet][2]/F");
    t_->Branch("fjet_subjets_py",         &evSummary_.fjet_subjets_py,        "fjet_subjets_py[fjet][2]/F");
    t_->Branch("fjet_subjets_pz",         &evSummary_.fjet_subjets_pz,        "fjet_subjets_pz[fjet][2]/F"); 
    t_->Branch("fjet_subjets_en",         &evSummary_.fjet_subjets_en,        "fjet_subjets_en[fjet][2]/F");
    t_->Branch("fjet_subjets_btag",       &evSummary_.fjet_subjets_btag,      "fjet_subjets_btag[fjet][2]/F");*/ 
    ///t_->Branch("fjet_subjets_partonFlavour",       &evSummary_.fjet_subjets_partonFlavour,      "fjet_subjets_partonFlavour[fjet][2]/I");
    ///t_->Branch("fjet_subjets_hadronFlavour",       &evSummary_.fjet_subjets_hadronFlavour,      "fjet_subjets_hadronFlavour[fjet][2]/I");



    //
    //met
    ///t_->Branch("imet_pt",                &evSummary_.imet_pt,                 "imet_pt[11]/F");    
    ///t_->Branch("imet_phi",               &evSummary_.imet_phi,                "imet_phi[11]/F");      

    t_->Branch("met_pt",                &evSummary_.met_pt,                 "met_pt/F");
    t_->Branch("met_phi",               &evSummary_.met_phi,                "met_phi/F");
    t_->Branch("met_sumMET",            &evSummary_.met_sumMET,             "met_sumMET/F");
    /*
    t_->Branch("rawpfmet_pt",                &evSummary_.rawpfmet_pt,                 "rawpfmet_pt/F");
    t_->Branch("rawpfmet_phi",               &evSummary_.rawpfmet_phi,                "rawpfmet_phi/F");
    t_->Branch("rawpfmet_sumMET",            &evSummary_.rawpfmet_sumMET,             "rawpfmet_sumMET/F");
    t_->Branch("rawcalomet_pt",              &evSummary_.rawcalomet_pt,               "rawcalomet_pt/F");
    t_->Branch("rawcalomet_phi",             &evSummary_.rawcalomet_phi,              "rawcalomet_phi/F");
    t_->Branch("rawcalomet_sumMET",          &evSummary_.rawcalomet_sumMET,           "rawcalomet_sumMET/F");
    t_->Branch("metNoHF_pt",                 &evSummary_.metNoHF_pt,                  "metNoHF_pt/F");
    t_->Branch("metNoHF_phi",                &evSummary_.metNoHF_phi,                 "metNoHF_phi/F");
    t_->Branch("metNoHF_sumMET",             &evSummary_.metNoHF_sumMET,              "metNoHF_sumMET/F");
    t_->Branch("metPuppi_pt",                &evSummary_.metPuppi_pt,                 "metPuppi_pt/F");
    t_->Branch("metPuppi_phi",               &evSummary_.metPuppi_phi,                "metPuppi_phi/F");
    t_->Branch("metPuppi_sumMET",            &evSummary_.metPuppi_sumMET,             "metPuppi_sumMET/F");
    */

    return true;
}

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    //event info
    /*
    t_->SetBranchAddress("run",             &evSummary_.run);
    t_->SetBranchAddress("lumi",            &evSummary_.lumi);
    t_->SetBranchAddress("event",           &evSummary_.event);
    */
    //t_->SetBranchAddress("puWeight",        &evSummary_.puWeight);  
    /* tmp
    t_->SetBranchAddress("curAvgInstLumi",  &evSummary_.curAvgInstLumi);
    t_->SetBranchAddress("curIntegLumi",    &evSummary_.curIntegLumi);
    */
    t_->SetBranchAddress("hasTrigger",      &evSummary_.hasTrigger);
    t_->SetBranchAddress("triggerType",     &evSummary_.triggerType);

    //primary vertex
    /*
    t_->SetBranchAddress("nvtx",            &evSummary_.nvtx);
    t_->SetBranchAddress("vtx_x",           &evSummary_.vtx_x);
    t_->SetBranchAddress("vtx_y",           &evSummary_.vtx_y);
    t_->SetBranchAddress("vtx_z",           &evSummary_.vtx_z);
    */
    /*
    t_->SetBranchAddress("fixedGridRhoAll",                           &evSummary_.fixedGridRhoAll);
    t_->SetBranchAddress("fixedGridRhoFastjetAll",                    &evSummary_.fixedGridRhoFastjetAll);
    t_->SetBranchAddress("fixedGridRhoFastjetAllCalo",                &evSummary_.fixedGridRhoFastjetAllCalo);
    t_->SetBranchAddress("fixedGridRhoFastjetCentralCalo",            &evSummary_.fixedGridRhoFastjetCentralCalo);
    t_->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp",   &evSummary_.fixedGridRhoFastjetCentralChargedPileUp);
    t_->SetBranchAddress("fixedGridRhoFastjetCentralNeutral",         &evSummary_.fixedGridRhoFastjetCentralNeutral);
    */
    //generator level info
    /*
    t_->SetBranchAddress("ngenITpu",        &evSummary_.ngenITpu);
    t_->SetBranchAddress("ngenOOTpu",       &evSummary_.ngenOOTpu);
    t_->SetBranchAddress("ngenOOTpum1",     &evSummary_.ngenOOTpum1);
    t_->SetBranchAddress("ngenTruepu",      &evSummary_.ngenTruepu);
    */
    /*
    t_->SetBranchAddress("pthat",           &evSummary_.pthat);
    t_->SetBranchAddress("genWeight",       &evSummary_.genWeight);
    t_->SetBranchAddress("qscale",          &evSummary_.qscale);
    t_->SetBranchAddress("x1",              &evSummary_.x1);
    t_->SetBranchAddress("x2",              &evSummary_.x2);
    t_->SetBranchAddress("id1",             &evSummary_.id1);
    t_->SetBranchAddress("id2",             &evSummary_.id2);
	*/
    t_->SetBranchAddress("lheHt",           &evSummary_.lheHt);  
    t_->SetBranchAddress("lheNJets",        &evSummary_.lheNJets);
    /*
    t_->SetBranchAddress("weight_QCDscale_muR1_muF1",       &evSummary_.weight_QCDscale_muR1_muF1);
    t_->SetBranchAddress("weight_QCDscale_muR1_muF2",       &evSummary_.weight_QCDscale_muR1_muF2);
    t_->SetBranchAddress("weight_QCDscale_muR1_muF0p5",       &evSummary_.weight_QCDscale_muR1_muF0p5);
    t_->SetBranchAddress("weight_QCDscale_muR2_muF1",       &evSummary_.weight_QCDscale_muR2_muF1);
    t_->SetBranchAddress("weight_QCDscale_muR2_muF2",       &evSummary_.weight_QCDscale_muR2_muF2);
    t_->SetBranchAddress("weight_QCDscale_muR2_muF0p5",       &evSummary_.weight_QCDscale_muR2_muF0p5);
    t_->SetBranchAddress("weight_QCDscale_muR0p5_muF1",       &evSummary_.weight_QCDscale_muR0p5_muF1);
    t_->SetBranchAddress("weight_QCDscale_muR0p5_muF2",       &evSummary_.weight_QCDscale_muR0p5_muF2);
    t_->SetBranchAddress("weight_QCDscale_muR0p5_muF0p5",       &evSummary_.weight_QCDscale_muR0p5_muF0p5);
    */

    /* t_->SetBranchAddress("npdfs", &evSummary_.npdfs);
    t_->SetBranchAddress("pdfWeights", evSummary_.pdfWeights);
    t_->SetBranchAddress("nalphaS", &evSummary_.nalphaS);
    t_->SetBranchAddress("alphaSWeights", evSummary_.alphaSWeights); */

	
    //mc truth
    /*t_->SetBranchAddress("nmcparticles",    &evSummary_.nmcparticles);
    t_->SetBranchAddress("mc_px",           evSummary_.mc_px);
    t_->SetBranchAddress("mc_py",           evSummary_.mc_py);
    t_->SetBranchAddress("mc_pz",           evSummary_.mc_pz);
    t_->SetBranchAddress("mc_en",           evSummary_.mc_en);
    t_->SetBranchAddress("mc_id",           evSummary_.mc_id);
    t_->SetBranchAddress("mc_status",       evSummary_.mc_status);
    t_->SetBranchAddress("mc_mom",          evSummary_.mc_mom);
    t_->SetBranchAddress("mc_momidx",       evSummary_.mc_momidx);*/
	

    //gen ground state B hadrons
//    t_->SetBranchAddress("mcbh",          &evSummary_.mcbh   );
//    t_->SetBranchAddress("mcbh_px",       evSummary_.mcbh_px );
//    t_->SetBranchAddress("mcbh_py",       evSummary_.mcbh_py );
//    t_->SetBranchAddress("mcbh_pz",       evSummary_.mcbh_pz );
//    t_->SetBranchAddress("mcbh_en",       evSummary_.mcbh_en );
//    t_->SetBranchAddress("mcbh_id",       evSummary_.mcbh_id );




    /* tmp
    // gen jets
    t_->SetBranchAddress("nmcjparticles",    &evSummary_.nmcjparticles); 
    t_->SetBranchAddress("mcj_px",           evSummary_.mcj_px); 
    t_->SetBranchAddress("mcj_py",           evSummary_.mcj_py); 
    t_->SetBranchAddress("mcj_pz",           evSummary_.mcj_pz); 
    t_->SetBranchAddress("mcj_en",           evSummary_.mcj_en); 
    t_->SetBranchAddress("mcj_id",           evSummary_.mcj_id); 
    t_->SetBranchAddress("mcj_status",       evSummary_.mcj_status); 
    t_->SetBranchAddress("mcj_mom",          evSummary_.mcj_mom); 
    */
    //muon
    t_->SetBranchAddress("mn",                      &evSummary_.mn);
    t_->SetBranchAddress("mn_px",                   evSummary_.mn_px);
    t_->SetBranchAddress("mn_py",                   evSummary_.mn_py);
    t_->SetBranchAddress("mn_pz",                   evSummary_.mn_pz);
    t_->SetBranchAddress("mn_en",                   evSummary_.mn_en);
    t_->SetBranchAddress("mn_id",                   evSummary_.mn_id);
    t_->SetBranchAddress("mn_type",                 evSummary_.mn_type);
    /* tmp
    t_->SetBranchAddress("mn_d0",                   evSummary_.mn_d0);
    t_->SetBranchAddress("mn_dZ",                   evSummary_.mn_dZ);
    t_->SetBranchAddress("mn_ip3d",                 evSummary_.mn_ip3d);
    t_->SetBranchAddress("mn_ip3dsig",              evSummary_.mn_ip3dsig);
    t_->SetBranchAddress("mn_IsLoose",              evSummary_.mn_IsLoose);
    t_->SetBranchAddress("mn_IsMedium",              evSummary_.mn_IsMedium);
    t_->SetBranchAddress("mn_IsTight",              evSummary_.mn_IsTight);
    t_->SetBranchAddress("mn_IsSoft",               evSummary_.mn_IsSoft);
    t_->SetBranchAddress("mn_IsHighPt",             evSummary_.mn_IsHighPt);

    t_->SetBranchAddress("mn_pileupIsoR03",         evSummary_.mn_pileupIsoR03);
    t_->SetBranchAddress("mn_chargedIsoR03",        evSummary_.mn_chargedIsoR03);
    t_->SetBranchAddress("mn_photonIsoR03",         evSummary_.mn_photonIsoR03);
    t_->SetBranchAddress("mn_neutralHadIsoR03",     evSummary_.mn_neutralHadIsoR03);
    t_->SetBranchAddress("mn_pileupIsoR04",         evSummary_.mn_pileupIsoR04);
    t_->SetBranchAddress("mn_chargedIsoR04",        evSummary_.mn_chargedIsoR04);
    t_->SetBranchAddress("mn_photonIsoR04",         evSummary_.mn_photonIsoR04);
    t_->SetBranchAddress("mn_neutralHadIsoR04",     evSummary_.mn_neutralHadIsoR04);
    */
    t_->SetBranchAddress("mn_passId",  evSummary_.mn_passId);
    t_->SetBranchAddress("mn_passIdLoose",  evSummary_.mn_passIdLoose);
    ///t_->SetBranchAddress("mn_passSoftMuon",  evSummary_.mn_passSoftMuon);
    t_->SetBranchAddress("mn_passIso",  evSummary_.mn_passIso);
    //    t_->SetBranchAddress("mn_relIso",  evSummary_.mn_relIso);
    //    t_->SetBranchAddress("mn_trkrelIso",  evSummary_.mn_trkrelIso);
    //    t_->SetBranchAddress("mn_nMatches",                   evSummary_.mn_nMatches);
    //t_->SetBranchAddress("mn_nMatchedStations",           evSummary_.mn_nMatchedStations);
    ///t_->SetBranchAddress("mn_validMuonHits",              evSummary_.mn_validMuonHits);
    //    t_->SetBranchAddress("mn_innerTrackChi2",             evSummary_.mn_innerTrackChi2);
    ///t_->SetBranchAddress("mn_trkLayersWithMeasurement",   evSummary_.mn_trkLayersWithMeasurement);
    ///t_->SetBranchAddress("mn_pixelLayersWithMeasurement", evSummary_.mn_pixelLayersWithMeasurement);


    //electron
    t_->SetBranchAddress("en",                      &evSummary_.en);
    t_->SetBranchAddress("en_px",                   evSummary_.en_px);
    t_->SetBranchAddress("en_py",                   evSummary_.en_py);
    t_->SetBranchAddress("en_pz",                   evSummary_.en_pz);
    t_->SetBranchAddress("en_en",                   evSummary_.en_en);
    t_->SetBranchAddress("en_id",                   evSummary_.en_id);
    ///t_->SetBranchAddress("en_cor_en",                   evSummary_.en_cor_en); 
    //    t_->SetBranchAddress("en_scale_corr",           evSummary_.en_scale_corr);  
    ///t_->SetBranchAddress("en_EtaSC",                evSummary_.en_EtaSC);           
    ///t_->SetBranchAddress("en_R9",                   evSummary_.en_R9);    
    ///t_->SetBranchAddress("en_gainSeed",                   evSummary_.en_gainSeed); 

    /*
    t_->SetBranchAddress("en_d0",                   evSummary_.en_d0);
    t_->SetBranchAddress("en_dZ",                   evSummary_.en_dZ);
    t_->SetBranchAddress("en_EtaSC",                evSummary_.en_EtaSC);
    t_->SetBranchAddress("en_PhiSC",                evSummary_.en_PhiSC);
    t_->SetBranchAddress("en_EnSC",                 evSummary_.en_EnSC);
    t_->SetBranchAddress("en_dEtaIn",               evSummary_.en_dEtaIn);
    t_->SetBranchAddress("en_dPhiIn",               evSummary_.en_dPhiIn);
    t_->SetBranchAddress("en_hOverE",               evSummary_.en_hOverE);
    t_->SetBranchAddress("en_R9",                   evSummary_.en_R9);
    t_->SetBranchAddress("en_sigmaIetaIeta",        evSummary_.en_sigmaIetaIeta);
    t_->SetBranchAddress("en_sigmaIetaIeta5x5",     evSummary_.en_sigmaIetaIeta5x5);
    t_->SetBranchAddress("en_ooEmooP",              evSummary_.en_ooEmooP);
    */
    /* tmp
    t_->SetBranchAddress("en_pileupIso",            evSummary_.en_pileupIso);
    t_->SetBranchAddress("en_chargedIso",           evSummary_.en_chargedIso);
    t_->SetBranchAddress("en_photonIso",            evSummary_.en_photonIso);
    t_->SetBranchAddress("en_neutralHadIso",        evSummary_.en_neutralHadIso);
    t_->SetBranchAddress("en_relIsoWithEA",         evSummary_.en_relIsoWithEA);
    t_->SetBranchAddress("en_relIsoWithDBeta",      evSummary_.en_relIsoWithDBeta);
    t_->SetBranchAddress("en_MissingHits",          evSummary_.en_MissingHits);
    t_->SetBranchAddress("en_passConversionVeto",   evSummary_.en_passConversionVeto);
    t_->SetBranchAddress("en_passVeto",             evSummary_.en_passVeto);
    t_->SetBranchAddress("en_passLoose",            evSummary_.en_passLoose);
    t_->SetBranchAddress("en_passMedium",           evSummary_.en_passMedium);
    t_->SetBranchAddress("en_passTight",            evSummary_.en_passTight);
    t_->SetBranchAddress("en_passHEEP",             evSummary_.en_passHEEP);
    
    t_->SetBranchAddress("en_passMVATrigMedium",    evSummary_.en_passMVATrigMedium);
    t_->SetBranchAddress("en_passMVATrigTight",     evSummary_.en_passMVATrigTight);
    t_->SetBranchAddress("en_IDMVATrigValue",       evSummary_.en_IDMVATrigValue);
    t_->SetBranchAddress("en_IDMVATrigCategory",    evSummary_.en_IDMVATrigCategory);
    t_->SetBranchAddress("en_istrue",               evSummary_.en_istrue);
    */
    t_->SetBranchAddress("en_passId", evSummary_.en_passId);
    t_->SetBranchAddress("en_passIdLoose", evSummary_.en_passIdLoose);
    t_->SetBranchAddress("en_passIso", evSummary_.en_passIso);
    t_->SetBranchAddress("en_relIso", evSummary_.en_relIso);
	/*
#ifdef YEAR_2017
    t_->SetBranchAddress("en_enSigmaValue",       evSummary_.en_enSigmaValue);
    t_->SetBranchAddress("en_enScaleValue",       evSummary_.en_enScaleValue);
    t_->SetBranchAddress("en_enScaleStatUp",      evSummary_.en_enScaleStatUp);
    t_->SetBranchAddress("en_enScaleStatDown",    evSummary_.en_enScaleStatDown);
    t_->SetBranchAddress("en_enScaleSystUp",      evSummary_.en_enScaleSystUp);
    t_->SetBranchAddress("en_enScaleSystDown",    evSummary_.en_enScaleSystDown);
    t_->SetBranchAddress("en_enScaleGainUp",      evSummary_.en_enScaleGainUp);
    t_->SetBranchAddress("en_enScaleGainDown",    evSummary_.en_enScaleGainDown);
    t_->SetBranchAddress("en_enSigmaRhoUp",       evSummary_.en_enSigmaRhoUp);
    t_->SetBranchAddress("en_enSigmaRhoDown",     evSummary_.en_enSigmaRhoDown);
    t_->SetBranchAddress("en_enSigmaPhiDown",     evSummary_.en_enSigmaPhiDown);
#endif
	*/
    /* tmp
    //tau
    t_->SetBranchAddress("ta",                      &evSummary_.ta);
    t_->SetBranchAddress("ta_px",                   evSummary_.ta_px);
    t_->SetBranchAddress("ta_py",                   evSummary_.ta_py);
    t_->SetBranchAddress("ta_pz",                   evSummary_.ta_pz);
    t_->SetBranchAddress("ta_en",                   evSummary_.ta_en);
    t_->SetBranchAddress("ta_id",                   evSummary_.ta_id);
    t_->SetBranchAddress("ta_dm",                   evSummary_.ta_dm);
    t_->SetBranchAddress("ta_newdm",                evSummary_.ta_newdm);
    t_->SetBranchAddress("ta_IsLooseIso",           evSummary_.ta_IsLooseIso);
    t_->SetBranchAddress("ta_IsMediumIso",          evSummary_.ta_IsMediumIso);
    t_->SetBranchAddress("ta_IsTightIso",           evSummary_.ta_IsTightIso);
    t_->SetBranchAddress("ta_combIsoDBeta3Hits",    evSummary_.ta_combIsoDBeta3Hits);
    t_->SetBranchAddress("ta_chargedIso",           evSummary_.ta_chargedIso);
    t_->SetBranchAddress("ta_neutralIso",           evSummary_.ta_neutralIso);
    t_->SetBranchAddress("ta_pileupIso",            evSummary_.ta_pileupIso);
    t_->SetBranchAddress("ta_passEleVetoLoose",     evSummary_.ta_passEleVetoLoose);
    t_->SetBranchAddress("ta_passEleVetoMedium",    evSummary_.ta_passEleVetoMedium);
    t_->SetBranchAddress("ta_passEleVetoTight",     evSummary_.ta_passEleVetoTight);
    t_->SetBranchAddress("ta_passMuVetoLoose3",     evSummary_.ta_passMuVetoLoose3);
    t_->SetBranchAddress("ta_passMuVetoTight3",     evSummary_.ta_passMuVetoTight3);

    */

    //jet (ak4PFJetsCHS)
    t_->SetBranchAddress("jet",                     &evSummary_.jet);
    t_->SetBranchAddress("jet_px",                  evSummary_.jet_px);
    t_->SetBranchAddress("jet_py",                  evSummary_.jet_py);
    t_->SetBranchAddress("jet_pz",                  evSummary_.jet_pz);
    t_->SetBranchAddress("jet_en",                  evSummary_.jet_en);
    // t_->SetBranchAddress("jet_btag0",               evSummary_.jet_btag0);
    t_->SetBranchAddress("jet_btag1",               evSummary_.jet_btag1);
    /*
    t_->SetBranchAddress("jet_btag1",               evSummary_.jet_btag1);
    t_->SetBranchAddress("jet_btag2",               evSummary_.jet_btag2);
    t_->SetBranchAddress("jet_btag3",               evSummary_.jet_btag3);
    t_->SetBranchAddress("jet_btag4",               evSummary_.jet_btag4);
    t_->SetBranchAddress("jet_btag5",               evSummary_.jet_btag5);
    t_->SetBranchAddress("jet_btag6",               evSummary_.jet_btag6);
    t_->SetBranchAddress("jet_btag7",               evSummary_.jet_btag7);
    t_->SetBranchAddress("jet_btag8",               evSummary_.jet_btag8);
    t_->SetBranchAddress("jet_btag9",               evSummary_.jet_btag9);
    t_->SetBranchAddress("jet_btag10",               evSummary_.jet_btag10);
    */
    // t_->SetBranchAddress("jet_mass",                evSummary_.jet_mass);
    /*///
    t_->SetBranchAddress("jet_area",                evSummary_.jet_area);
    t_->SetBranchAddress("jet_pu",                  evSummary_.jet_pu);
    t_->SetBranchAddress("jet_puId",                evSummary_.jet_puId);
	*/
    // t_->SetBranchAddress("jet_PFLoose",             evSummary_.jet_PFLoose);
    t_->SetBranchAddress("jet_PFTight",             evSummary_.jet_PFTight);
    /*
    t_->SetBranchAddress("jet_mother_id",           evSummary_.jet_mother_id);
    t_->SetBranchAddress("jet_partonFlavour",       evSummary_.jet_partonFlavour);
    t_->SetBranchAddress("jet_hadronFlavour",       evSummary_.jet_hadronFlavour);
    t_->SetBranchAddress("jet_parton_px",           evSummary_.jet_parton_px);
    t_->SetBranchAddress("jet_parton_py",           evSummary_.jet_parton_py);
    t_->SetBranchAddress("jet_parton_pz",           evSummary_.jet_parton_pz);
    t_->SetBranchAddress("jet_parton_en",           evSummary_.jet_parton_en);
    t_->SetBranchAddress("jet_genpt",               evSummary_.jet_genpt);
	*/

    /*   
    //sv : Inclusive Secondary Vertices from slimmedSecondaryVertices
    t_->SetBranchAddress("sv"              , &evSummary_.sv               ) ;
    t_->SetBranchAddress("sv_px"           , evSummary_.sv_px             ) ;
    t_->SetBranchAddress("sv_py"           , evSummary_.sv_py             ) ;
    t_->SetBranchAddress("sv_pz"           , evSummary_.sv_pz             ) ;
    t_->SetBranchAddress("sv_en"           , evSummary_.sv_en             ) ;
    t_->SetBranchAddress("sv_ntrk"         , evSummary_.sv_ntrk           ) ;
    t_->SetBranchAddress("sv_dxy"          , evSummary_.sv_dxy            ) ;
    t_->SetBranchAddress("sv_dxyz"         , evSummary_.sv_dxyz           ) ;
    t_->SetBranchAddress("sv_dxyz_signif"  , evSummary_.sv_dxyz_signif    ) ;
    t_->SetBranchAddress("sv_cos_dxyz_p"   , evSummary_.sv_cos_dxyz_p     ) ;
    t_->SetBranchAddress("sv_chi2"         , evSummary_.sv_chi2           ) ;
    t_->SetBranchAddress("sv_ndof"         , evSummary_.sv_ndof           ) ;
    t_->SetBranchAddress("sv_mc_nbh_moms"  , evSummary_.sv_mc_nbh_moms    ) ;
    t_->SetBranchAddress("sv_mc_nbh_daus"  , evSummary_.sv_mc_nbh_daus    ) ;
//    t_->SetBranchAddress("sv_mc_mcbh_ind"  , evSummary_.sv_mc_mcbh_ind    ) ;
    */

    /*
    //pjet: slimmedJetsPuppi
    t_->SetBranchAddress("pjet",                     &evSummary_.pjet);
    t_->SetBranchAddress("pjet_px",                  evSummary_.pjet_px);
    t_->SetBranchAddress("pjet_py",                  evSummary_.pjet_py);
    t_->SetBranchAddress("pjet_pz",                  evSummary_.pjet_pz);
    t_->SetBranchAddress("pjet_en",                  evSummary_.pjet_en);
    t_->SetBranchAddress("pjet_genpt",               evSummary_.pjet_genpt);
    t_->SetBranchAddress("pjet_btag0",               evSummary_.pjet_btag0);
    t_->SetBranchAddress("pjet_btag1",               evSummary_.pjet_btag1);
    t_->SetBranchAddress("pjet_btag2",               evSummary_.pjet_btag2);
    t_->SetBranchAddress("pjet_btag3",               evSummary_.pjet_btag3);
    t_->SetBranchAddress("pjet_btag4",               evSummary_.pjet_btag4);
    t_->SetBranchAddress("pjet_btag5",               evSummary_.pjet_btag5);
    t_->SetBranchAddress("pjet_btag6",               evSummary_.pjet_btag6);
    t_->SetBranchAddress("pjet_btag7",               evSummary_.pjet_btag7);
    t_->SetBranchAddress("pjet_btag8",               evSummary_.pjet_btag8);
    t_->SetBranchAddress("pjet_btag9",               evSummary_.pjet_btag9);
    t_->SetBranchAddress("pjet_btag10",               evSummary_.pjet_btag10);
    */

    //fjet (ak8PFJetsCHS)
    
    t_->SetBranchAddress("fjet",                    &evSummary_.fjet);
    t_->SetBranchAddress("fjet_px",                 evSummary_.fjet_px);
    t_->SetBranchAddress("fjet_py",                 evSummary_.fjet_py);
    t_->SetBranchAddress("fjet_pz",                 evSummary_.fjet_pz);
    t_->SetBranchAddress("fjet_en",                 evSummary_.fjet_en);
    /*
    t_->SetBranchAddress("fjet_btag0",              evSummary_.fjet_btag0);     
    t_->SetBranchAddress("fjet_btag1",              evSummary_.fjet_btag1);
    t_->SetBranchAddress("fjet_btag2",              evSummary_.fjet_btag2);
    t_->SetBranchAddress("fjet_btag3",              evSummary_.fjet_btag3);
    t_->SetBranchAddress("fjet_btag4",              evSummary_.fjet_btag4);
    t_->SetBranchAddress("fjet_btag5",              evSummary_.fjet_btag5);
    t_->SetBranchAddress("fjet_btag6",              evSummary_.fjet_btag6);
    t_->SetBranchAddress("fjet_btag7",              evSummary_.fjet_btag7);
    t_->SetBranchAddress("fjet_btag8",              evSummary_.fjet_btag8);
    t_->SetBranchAddress("fjet_btag9",              evSummary_.fjet_btag9);
    */
    t_->SetBranchAddress("fjet_btag10",             evSummary_.fjet_btag10);
    t_->SetBranchAddress("fjet_btag11",             evSummary_.fjet_btag11);
    t_->SetBranchAddress("fjet_btag12",             evSummary_.fjet_btag12);
    t_->SetBranchAddress("fjet_btag13",             evSummary_.fjet_btag13);
    t_->SetBranchAddress("fjet_btag14",             evSummary_.fjet_btag14);
    t_->SetBranchAddress("fjet_btag15",             evSummary_.fjet_btag15);
    t_->SetBranchAddress("fjet_btag16",             evSummary_.fjet_btag16);
    t_->SetBranchAddress("fjet_btag17",             evSummary_.fjet_btag17);


    //    t_->SetBranchAddress("fjet_genpt",              evSummary_.fjet_genpt);
    // t_->SetBranchAddress("fjet_prunedM",            evSummary_.fjet_prunedM);
    // t_->SetBranchAddress("fjet_softdropM",          evSummary_.fjet_softdropM);
    // t_->SetBranchAddress("fjet_trimmedM",           evSummary_.fjet_trimmedM);
    // t_->SetBranchAddress("fjet_filteredM",          evSummary_.fjet_filteredM);
    // t_->SetBranchAddress("fjet_tau1",               evSummary_.fjet_tau1);
    // t_->SetBranchAddress("fjet_tau2",               evSummary_.fjet_tau2);
    // t_->SetBranchAddress("fjet_tau3",               evSummary_.fjet_tau3);
    // t_->SetBranchAddress("fjet_tau4",               evSummary_.fjet_tau4);
    /*
    t_->SetBranchAddress("fjet_mother_id",          evSummary_.fjet_mother_id);
    t_->SetBranchAddress("fjet_partonFlavour",      evSummary_.fjet_partonFlavour);
    t_->SetBranchAddress("fjet_hadronFlavour",      evSummary_.fjet_hadronFlavour); 
    t_->SetBranchAddress("fjet_parton_px",          evSummary_.fjet_parton_px);
    t_->SetBranchAddress("fjet_parton_py",          evSummary_.fjet_parton_py);
    t_->SetBranchAddress("fjet_parton_pz",          evSummary_.fjet_parton_pz);
    t_->SetBranchAddress("fjet_parton_en",          evSummary_.fjet_parton_en);
    
    t_->SetBranchAddress("fjet_chf",                evSummary_.fjet_chf);
    t_->SetBranchAddress("fjet_nhf",                evSummary_.fjet_nhf);
    t_->SetBranchAddress("fjet_phf",                evSummary_.fjet_phf);
    t_->SetBranchAddress("fjet_muf",                evSummary_.fjet_muf);
    t_->SetBranchAddress("fjet_elf",                evSummary_.fjet_elf);*/

    //t_->SetBranchAddress("fjet_ecfB1N2",            evSummary_.fjet_ecfB1N2);
    //t_->SetBranchAddress("fjet_ecfB1N3",            evSummary_.fjet_ecfB1N3);
    //t_->SetBranchAddress("fjet_ecfB2N2",            evSummary_.fjet_ecfB2N2);
    //t_->SetBranchAddress("fjet_ecfB2N3",            evSummary_.fjet_ecfB2N3);
    
    t_->SetBranchAddress("fjet_subjet_count",       evSummary_.fjet_subjet_count);
    // t_->SetBranchAddress("fjet_subjets_px",         &evSummary_.fjet_subjets_px);
    // t_->SetBranchAddress("fjet_subjets_py",         &evSummary_.fjet_subjets_py);
    // t_->SetBranchAddress("fjet_subjets_pz",         &evSummary_.fjet_subjets_pz);
    // t_->SetBranchAddress("fjet_subjets_en",         &evSummary_.fjet_subjets_en);
    // t_->SetBranchAddress("fjet_subjets_btag",       &evSummary_.fjet_subjets_btag);
    // t_->SetBranchAddress("fjet_subjets_partonFlavour",       &evSummary_.fjet_subjets_partonFlavour);
    // t_->SetBranchAddress("fjet_subjets_hadronFlavour",       &evSummary_.fjet_subjets_hadronFlavour);


    //met
    ///t_->SetBranchAddress("imet_pt",                  &evSummary_.imet_pt);
    ///t_->SetBranchAddress("imet_phi",                 &evSummary_.imet_phi); 

    t_->SetBranchAddress("met_pt",                  &evSummary_.met_pt);
    t_->SetBranchAddress("met_phi",                 &evSummary_.met_phi);
    t_->SetBranchAddress("met_sumMET",              &evSummary_.met_sumMET);
    /* tmp
    t_->SetBranchAddress("rawpfmet_pt",                  &evSummary_.rawpfmet_pt);
    t_->SetBranchAddress("rawpfmet_phi",                 &evSummary_.rawpfmet_phi);
    t_->SetBranchAddress("rawpfmet_sumMET",              &evSummary_.rawpfmet_sumMET);
    t_->SetBranchAddress("rawcalomet_pt",                  &evSummary_.rawcalomet_pt);
    t_->SetBranchAddress("rawcalomet_phi",                 &evSummary_.rawcalomet_phi);
    t_->SetBranchAddress("rawcalomet_sumMET",              &evSummary_.rawcalomet_sumMET);
    t_->SetBranchAddress("metNoHF_pt",                  &evSummary_.metNoHF_pt);
    t_->SetBranchAddress("metNoHF_phi",                 &evSummary_.metNoHF_phi);
    t_->SetBranchAddress("metNoHF_sumMET",              &evSummary_.metNoHF_sumMET);
    t_->SetBranchAddress("metPuppi_pt",                  &evSummary_.metPuppi_pt);
    t_->SetBranchAddress("metPuppi_phi",                 &evSummary_.metPuppi_phi);
    t_->SetBranchAddress("metPuppi_sumMET",              &evSummary_.metPuppi_sumMET);
    */

    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nmcparticles=0;
    evSummary_.nmcjparticles=0;
    evSummary_.run=0;
    evSummary_.lumi=0;
    evSummary_.event=0;
    evSummary_.mn=0;
    evSummary_.en=0;
    //evSummary_.ta=0;
    evSummary_.jet=0;
    // evSummary_.sv=0;
    ///// evSummary_.pjet=0;
    evSummary_.fjet=0;
    // evSummary_.npdfs=0;
    // evSummary_.nalphaS=0;

    evSummary_.ngenITpu=0;
    evSummary_.ngenOOTpu=0;
    evSummary_.ngenOOTpum1=0;
    evSummary_.ngenTruepu=0;
    evSummary_.genWeight=0;
    evSummary_.qscale=0;
    evSummary_.x1=0;
    evSummary_.x2=0;
    evSummary_.id1=0;
    evSummary_.id2=0;
    
    evSummary_.lheNJets=0;
    evSummary_.pthat=0;
    /*
    evSummary_.weight_QCDscale_muR1_muF1=0;
    evSummary_.weight_QCDscale_muR1_muF2=0;
    evSummary_.weight_QCDscale_muR1_muF0p5=0;
    evSummary_.weight_QCDscale_muR2_muF1=0;
    evSummary_.weight_QCDscale_muR2_muF2=0;
    evSummary_.weight_QCDscale_muR2_muF0p5=0;
    evSummary_.weight_QCDscale_muR0p5_muF1=0;
    evSummary_.weight_QCDscale_muR0p5_muF2=0;
    evSummary_.weight_QCDscale_muR0p5_muF0p5=0;
    */
    return ;
}

//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
    return ;
}

//
DataEvtSummaryHandler::~DataEvtSummaryHandler()
{
}
