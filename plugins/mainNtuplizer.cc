// -*- C++ -*-
//
// Package:    bsmhiggs_fwk
// Class:      mainNtuplizer
// 
/**\class mainNtuplizer mainNtuplizer.cc bwmhiggs_fwk/plugins/mainNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Georgia Karapostoli
//         Created:  Sun, 17 Sep 2017 08:27:19 GMT
//
//
#define YEAR_2017
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#ifdef YEAR_2017
#include "FWCore/Framework/interface/Run.h"
#endif

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Top pt reweighting:
//#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

// BSM Higgs code:
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
//#include "UserCode/bsmhiggs_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/bsmhiggs_fwk/interface/PDFInfo.h"
#include "UserCode/bsmhiggs_fwk/interface/rochcor2016.h"
#include "UserCode/bsmhiggs_fwk/interface/muresolution_run2.h"
#include "UserCode/bsmhiggs_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/bsmhiggs_fwk/interface/BtagUncertaintyComputer.h"

#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/EwkCorrections.h"
#include "UserCode/bsmhiggs_fwk/interface/ZZatNNLO.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h" 
#include <Math/VectorUtil.h>
#include "TRandom3.h"
//
// class declaration
//
using namespace edm;
using namespace reco;
using namespace pat;
using namespace std;
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
//public edm::one::EDAnalyzer<edm::one::SharedResources>  {

class mainNtuplizer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
  public:
      explicit mainNtuplizer(const edm::ParameterSet&);
      ~mainNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
      bool isAncestor( const reco::Candidate& ancestor, const reco::Candidate& daughter );

   private:

    edm::EDGetTokenT<reco::VertexCollection> vtxTag_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotTag_;
    edm::EDGetTokenT<pat::MuonCollection> muonTag_;
    edm::EDGetTokenT<pat::ElectronCollection> electronTag_;
  //  edm::EDGetTokenT<edm::View<pat::Electron> > electronTag_;

    // edm::EDGetTokenT<pat::TauCollection> tauTag_;
    // edm::EDGetTokenT<pat::PhotonCollection> photonTag_;
    edm::EDGetTokenT<pat::JetCollection> jetTag_;
//    edm::EDGetTokenT<pat::JetCollection> jetPuppiTag_;
    edm::EDGetTokenT<pat::JetCollection> fatjetTag_;
    edm::EDGetTokenT<pat::METCollection> metTag_;
    edm::EDGetTokenT<pat::METCollection> metTagData_;
    edm::EDGetTokenT<pat::METCollection> metNoHFTag_;
    edm::EDGetTokenT<pat::METCollection> metPuppiTag_;

    edm::EDGetTokenT<edm::TriggerResults> metFilterBitsTag_;

    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svTag_;
  
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenTag_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoTag_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoTag_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenTag_;
  // edm::EDGetTokenT<edm::View<reco::GenJet> > genjetTag_;
  edm::InputTag lheRunInfoTag_;
  edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
    edm::EDGetTokenT<double> rhoFastjetAllTag_;

  //  EnergyScaleCorrection_class eScaler_;     
  edm::EDGetTokenT<EcalRecHitCollection> reducedEBRecHitCollectionToken_;
  edm::EDGetTokenT<EcalRecHitCollection> reducedEERecHitCollectionToken_;

    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

    std::vector<std::string> DoubleMuTrigs_, DoubleEleTrigs_, SingleMuTrigs_, SingleEleTrigs_, MuEGTrigs_;// DoubleTauTrigs_;

    string proc_;
    bool isMC_;
    double xsec_;
    int mctruthmode_;
    bool verbose_;
  
    DataEvtSummaryHandler summaryHandler_;

  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iRun, edm::EventSetup const&) override {};
  void endJob() override;


  // ----------member data ---------------------------
    float curAvgInstLumi_;
    float curIntegLumi_;
  
    int firstPdfWeight;
    int lastPdfWeight;
    int firstAlphasWeight;
    int lastAlphasWeight;

  // Some histograms
  TH1F * h_nevents, *h_negevents, *h_posevents;
  TH1F * h_pileup, * h_pileuptrue;
  TH1F * h_sumWeights, * h_sumScaleWeights , * h_sumPdfWeights ,* h_sumAlphasWeights; 
  TH1F * h_metFilter;
  
  bool isMC_ttbar;
  /*
  std::string bit_string_stat = "001";
  std::string bit_string_syst = "010";
  std::string bit_string_gain = "100";
  std::bitset<6> bit_stat(bit_string_stat);
  std::bitset<6> bit_syst(bit_string_syst);
  std::bitset<6> bit_gain(bit_string_gain);
  */
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
mainNtuplizer::mainNtuplizer(const edm::ParameterSet& iConfig):
  vtxTag_(		consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("verticesTag"))		),
    beamSpotTag_(	consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))			),
    muonTag_(		consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonsTag"))			),
    electronTag_(	consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronsTag"))		),
  //    electronTag_(	consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronsTag"))	),
//  tauTag_(		consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tausTag"))			),
//   photonTag_(		consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonsTag"))		),
    jetTag_(		consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsTag"))			),
//    jetPuppiTag_(       consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppiTag"))               ),
    fatjetTag_(		consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjetsTag"))			),
    metTag_(		consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsTag"))			),
    metTagData_(		consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsTagData"))			),
    metNoHFTag_(        consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsNoHFTag"))                ),
    metPuppiTag_(       consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsPuppiTag"))               ),
    metFilterBitsTag_(	consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBitsTag"))		),
    svTag_(		consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("svTag"))			),
    prunedGenTag_(	consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedTag"))	),
    puInfoTag_(         consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoTag"))     ),
    genInfoTag_(        consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoTag"))                ),
    packedGenTag_(	consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedTag"))	),
  //  genjetTag_(		consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsTag"))		),
  lheRunInfoTag_(     iConfig.getParameter<edm::InputTag>("lheInfo")                                                  ),
  lheRunInfoToken_(   consumes<LHERunInfoProduct,edm::InRun>(lheRunInfoTag_)						),
    rhoFastjetAllTag_(  	consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll")) 			),
  //  eScaler_(iConfig.getParameter<std::string>("correctionFile") ),
	   reducedEBRecHitCollectionToken_(     consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEcalRecHitsEB")) ),
	   reducedEERecHitCollectionToken_(     consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEcalRecHitsEE"))  ),
    triggerBits_(	consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))			),
    triggerObjects_(	consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(	consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))		),
    DoubleMuTrigs_(	iConfig.getParameter<std::vector<std::string> >("DoubleMuTrigs")				),
    DoubleEleTrigs_(    iConfig.getParameter<std::vector<std::string> >("DoubleEleTrigs")                               ),
    SingleMuTrigs_(	iConfig.getParameter<std::vector<std::string> >("SingleMuTrigs")				),
    SingleEleTrigs_(    iConfig.getParameter<std::vector<std::string> >("SingleEleTrigs")				),
    MuEGTrigs_(		iConfig.getParameter<std::vector<std::string> >("MuEGTrigs")					),
    proc_(	iConfig.getParameter<std::string>("dtag")							),
    isMC_(		iConfig.getParameter<bool>("isMC")								),
    xsec_(  iConfig.getParameter<double>("xsec")                                                                          ),
    mctruthmode_( iConfig.getParameter<int>("mctruthmode")                                                                ),
    verbose_(		iConfig.getParameter<bool>("verbose")								),
  // std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");
    //rhoAllTag_(			consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))				),
    //rhoFastjetAllCaloTag_( 	consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo")) 		),
    //rhoFastjetCentralCaloTag_(  consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo")) 		),
    //rhoFastjetCentralChargedPileUpTag_(  consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp")) ),
    //rhoFastjetCentralNeutralTag_(  	consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))	 ),
  //    eleMediumIdMapTokenTrig_(	consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMapTrig"))	),
  // eleTightIdMapTokenTrig_(	consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMapTrig"))	),
    curAvgInstLumi_(0),
    curIntegLumi_(0)
{
   //now do what ever initialization is needed

  consumesMany<LHEEventProduct>();

  edm::Service<TFileService> fs;
  
  summaryHandler_.initTree(  fs->make<TTree>("data","Event Summary") );
  TFileDirectory baseDir=fs->mkdir(iConfig.getParameter<std::string>("dtag"));

  if ( verbose_ ) {
     printf("  Verbose set to true.  Will print info for each event.\n") ;
  } else {
     printf("  Verbose set to false.  Will be quiet.\n") ;
  }

  
  printf("Definition of plots\n");
  h_nevents = fs->make< TH1F>("nevents",";nevents; nevents",1,-0.5,0.5);
  h_negevents = fs->make< TH1F>("n_negevents",";n_negevents; n_negevents",1,-0.5,0.5);
  h_posevents = fs->make< TH1F>("n_posevents",";n_posevents; n_posevents",1,-0.5,0.5);
  h_pileup = fs->make< TH1F>("pileup", ";Pileup; Events",100,-0.5,99.5);
  //mon_.addHistogram(new TH1F("integlumi", ";Integrated luminosity ; Events",100,0,1e5);
  //mon_.addHistogram(new TH1F("instlumi", ";Max average inst. luminosity; Events",100,0,1e5);
  h_pileuptrue = fs->make< TH1F>("pileuptrue", ";True pileup; Events",100,-0.5,99.5);
  
  h_sumWeights = fs->make< TH1F>("sumWeights",";;sumWeights;",1,-0.5,0.5);
  h_sumScaleWeights = fs->make< TH1F>("sumScaleWeights",";;sumScaleWeights;",9,-0.5,8.5);
  h_sumPdfWeights = fs->make< TH1F>("sumPdfWeights",";;sumPdfWeights;",100,-0.5,99.5);
  h_sumAlphasWeights = fs->make< TH1F>("sumAlphasWeights",";;sumAlphasWeights;",2,-0.5,1.5);
  
  h_metFilter = fs->make<TH1F>( "metFilter",";metEventflow",20,0,20);
  h_metFilter->GetXaxis()->SetBinLabel(1,"raw");
  h_metFilter->GetXaxis()->SetBinLabel(2,"globalTightHalo2016Filter");
  h_metFilter->GetXaxis()->SetBinLabel(3,"goodVertices");
  h_metFilter->GetXaxis()->SetBinLabel(4,"eeBadScFilter");
  h_metFilter->GetXaxis()->SetBinLabel(5,"EcalDeadCellTriggerPrimitiveFilter");
  h_metFilter->GetXaxis()->SetBinLabel(6,"HBHENoiseFilter");
  h_metFilter->GetXaxis()->SetBinLabel(7,"HBHENoiseIsoFilter");
  h_metFilter->GetXaxis()->SetBinLabel(8,"BadPFMuonFilter");
  h_metFilter->GetXaxis()->SetBinLabel(9,"BadChargedCandidateFilte");
  h_metFilter->GetXaxis()->SetBinLabel(10,"badMuonHIPFilter");
  h_metFilter->GetXaxis()->SetBinLabel(11,"duplicateMuonHIPFilter");

  // patUtils::MetFilter metFilter;
  //MET CORRection level
  //  pat::MET::METCorrectionLevel metcor = pat::MET::METCorrectionLevel::Type1XY;
  
  // Use for Top pt re-weighting
  isMC_ttbar = isMC_ && (string(proc_.c_str()).find("TeV_TTJets") != string::npos);

  // Energy scale and resolution residuals
  // Remember we are using the MC-based electron calibration in miniAOD
  /*
  if (isMC_) {
    eScaler_.doScale=false;
    eScaler_.doSmearings=true;
  } else {
    eScaler_.doScale=true;
    eScaler_.doSmearings=false;
  }
  */
  //  usesResource("TFileService");

}


mainNtuplizer::~mainNtuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
mainNtuplizer::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //mon_.fillHisto("nevents","all",1.,1.); //increment event count
   h_nevents->Fill(0.);
   
   summaryHandler_.resetStruct();
   //event summary to be filled
   DataEvtSummary_t &ev=summaryHandler_.getEvent();

   //float weight = xsecWeight;

   // PU weights
   //float puWeight_(1.0);
   //   ev.puWeight = puWeight_;

   ev.run = event.eventAuxiliary().run() ;
   ev.lumi = event.eventAuxiliary().luminosityBlock() ;
   ev.event = event.eventAuxiliary().event() ;

   if ( verbose_ ) { printf("\n\n ================= Run %u , lumi %u , event %lld\n\n", ev.run, ev.lumi, ev.event ) ; }

   ev.mcbh = 0 ;
   std::vector<reco::GenParticle> b_hadrons ;
   
   if (isMC_) {
     
     edm::Handle< std::vector<PileupSummaryInfo> > puInfoH;
     event.getByToken(puInfoTag_,puInfoH);
     //   puInfoH.getByLabel(event, "slimmedAddPileupInfo");
     
     int npuOOT(0),npuIT(0),npuOOTm1(0);
     float truePU(0);
     if(puInfoH.isValid()) {
       for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++) {
	 if(it->getBunchCrossing()==0) {
	   npuIT += it->getPU_NumInteractions();
	   truePU = it->getTrueNumInteractions();
	 } else                          npuOOT += it->getPU_NumInteractions();
	 if(it->getBunchCrossing()<0)  npuOOTm1 += it->getPU_NumInteractions();
	 
       }
     }
     ev.ngenITpu=npuIT;
     ev.ngenOOTpu=npuOOT;
     ev.ngenOOTpum1=npuOOTm1;
     ev.ngenTruepu=truePU;
     h_pileup->Fill(ev.ngenITpu,1);
     h_pileuptrue->Fill(truePU);
     //     mon_.fillHisto("pileup","all",ev.ngenITpu,1);
     //     mon_.fillHisto("pileuptrue","all",truePU,1);
     
     if ( verbose_ ) { printf("  MC : Npu= %3d, truePU = %5.1f\n", npuIT, truePU ) ; }
     
     //retrieve pdf info
     GenEventInfoProduct eventInfo;
     edm::Handle< GenEventInfoProduct > genEventInfoProd;
     event.getByToken(genInfoTag_, genEventInfoProd);
     if(genEventInfoProd.isValid()){ eventInfo = *genEventInfoProd;}                
     
     ev.genWeight = eventInfo.weight();
     float SignGenWeight=1;
     if(ev.genWeight<0) SignGenWeight=-1;

     h_sumWeights->Fill(0.,SignGenWeight);
     //  mon_.fillHisto("sumWeights","all",0.,SignGenWeight);
     ev.qscale = eventInfo.qScale();
     if(eventInfo.pdf()) {
       ev.x1  = eventInfo.pdf()->x.first;
       ev.x2  = eventInfo.pdf()->x.second;
       ev.id1 = eventInfo.pdf()->id.first;
       ev.id2 = eventInfo.pdf()->id.second;
     }
     if(eventInfo.binningValues().size()>0) ev.pthat = eventInfo.binningValues()[0];
     
     if(ev.genWeight<0) h_negevents->Fill(0.); //mon_.fillHisto("n_negevents","all",1.0,1); //increment negative event count
     if(ev.genWeight>0) h_posevents->Fill(0.); // mon_.fillHisto("n_posevents","all",1.0,1); //increment positive event count
     
     //scale variations
    
     //https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
     vector< edm::Handle<LHEEventProduct> > EvtHandles; 
     event.getManyByType(EvtHandles);
     ev.npdfs=0;
     ev.nalphaS=0;
     ev.lheNJets=0;
     ev.lheHt=0.;           
     if(EvtHandles.size()>0)
     {
       edm::Handle<LHEEventProduct> EvtHandle = EvtHandles.front();
       if(EvtHandle.isValid() && EvtHandle->weights().size()>0)
       {       
         //fill pdf+alpha_s variation weights
         if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(EvtHandle->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100)
         { 
           //fill pdf variation weights after converting with mc2hessian transformation
           //std::array<double, 100> inpdfweights;
           for (int iwgt=firstPdfWeight; iwgt<=lastPdfWeight; ++iwgt)
           {
             ev.pdfWeights[ev.npdfs] = SignGenWeight * EvtHandle->weights()[iwgt].wgt/EvtHandle->originalXWGTUP();
             h_sumPdfWeights->Fill(double(ev.npdfs), ev.pdfWeights[ev.npdfs]);
             //mon_.fillHisto("sumPdfWeights","all",double(ev.npdfs), ev.pdfWeights[ev.npdfs]);
             ev.npdfs++;
           }
   
           //fill alpha_s variation weights
           if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(EvtHandle->weights().size()))
           {
             for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt)
             {
               ev.alphaSWeights[ev.nalphaS] = SignGenWeight * EvtHandle->weights()[iwgt].wgt/EvtHandle->originalXWGTUP();
               h_sumAlphasWeights->Fill(double(ev.nalphaS), ev.alphaSWeights[ev.nalphaS]);
               //mon_.fillHisto("sumAlphasWeights","all",double(ev.nalphaS), ev.alphaSWeights[ev.nalphaS]);
               ev.nalphaS++;
             }
           }   
         } // pdf variation weights END

         //Add lhe njets into DataEvtSummaryHandler
         const lhef::HEPEUP& lheEvent = EvtHandle->hepeup();
         std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

	 float genparticles_lheHt(0.);

         size_t numParticles = lheParticles.size();
         for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle )
         {
           int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
           int status = lheEvent.ISTUP[idxParticle];
           if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) )
           { // quarks and gluons
	     genparticles_lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
             ev.lheNJets++;
           }
         }
	 ev.lheHt = genparticles_lheHt;

       }// EvtHandle.isValid
     }
     
     //
     // gen particles
     //
     // reco::GenParticleCollection gen;
     Handle<edm::View<reco::GenParticle> > pruned;
     event.getByToken(prunedGenTag_,pruned);

     std::vector<TLorentzVector> chLeptons;       
     
     if ( verbose_ ) { printf("\n\n Gen particles:\n" ) ; }
     
     //Look for mother particle and Fill gen variables
     ev.nmcparticles = 0;  
   
     //for(unsigned int igen=0; igen<gen.size(); igen++){
     int igen(0);
     for (auto & it : *pruned) {
       
       {
	 int pdgId = it.pdgId();
	 if (  abs( pdgId ) == 511
	       || abs( pdgId ) == 521
	       || abs( pdgId ) == 531
	       || abs( pdgId ) == 541
	       || abs( pdgId ) == 5122
	       || abs( pdgId ) == 5112
	       || abs( pdgId ) == 5222
	       || abs( pdgId ) == 5132
	       || abs( pdgId ) == 5232
	       || abs( pdgId ) == 5332 ) {
	   b_hadrons.emplace_back( it ) ;
	   ev.mcbh_id[ev.mcbh] = pdgId ;
	   ev.mcbh_px[ev.mcbh] = it.px() ;
	   ev.mcbh_py[ev.mcbh] = it.py() ;
	   ev.mcbh_pz[ev.mcbh] = it.pz() ;
	   ev.mcbh_en[ev.mcbh] = it.energy() ;
	   ev.mcbh ++ ;
	 }
       }
       
       // Specific info for Top pt re-weighting
       if (isMC_ttbar) {
	 //'isLastCopy' definition of the parton-level top quark (after radiation and before decay) in Pythia8
	 // It is crucial that the pT value used when you calculate the weights is derived from the 'isLastCopy' definition
	 if (it.isLastCopy()) {
	   if (abs(it.pdgId()) == 6) {

	     if (verbose_) { printf("  Top : ID=%6d, m=%5.1f, momID=%6d : pt=%6.1f, status=%d\n",  
				   it.pdgId(),
				   it.mass(),
				   findFirstMotherWithDifferentID(&it)->pdgId(),
				   it.pt(),
				   it.status()
				   );
	     }
	     
	     ev.mc_px[ev.nmcparticles] = it.px();
	     ev.mc_py[ev.nmcparticles] = it.py(); 
	     ev.mc_pz[ev.nmcparticles] = it.pz(); 
	     ev.mc_en[ev.nmcparticles] = it.energy(); 
	     ev.mc_id[ev.nmcparticles] = it.pdgId(); 
	     ev.mc_mom[ev.nmcparticles] = findFirstMotherWithDifferentID(&it)->pdgId(); 

	     ev.mc_momidx[ev.nmcparticles] = -999;
	     ev.mc_status[ev.nmcparticles] = it.status(); 

	     ev.nmcparticles++;    
	   }
	   //	   printf("GenParticle isLastCopy iwth pdgId=%d\n\n",it.pdgId());
	 }
       }

       //       if(!it.isHardProcess()) continue; 
       if(it.isHardProcess()){

       //find the ID of the first mother that has a different ID than the particle itself
	 const reco::Candidate* mom = findFirstMotherWithDifferentID(&it);
	 
	 if (mom) {
	   int pid = it.pdgId();
	   
	   ev.mc_px[ev.nmcparticles] = it.px();
	   ev.mc_py[ev.nmcparticles] = it.py();
	   ev.mc_pz[ev.nmcparticles] = it.pz();
	   ev.mc_en[ev.nmcparticles] = it.energy();
	   ev.mc_id[ev.nmcparticles] = it.pdgId();
	   ev.mc_mom[ev.nmcparticles] = mom->pdgId();
	   
	   // loop over genParticles to find the mom index
	   int idx=0; int idxx=0;
	   for (auto & ig : *pruned) {
	     //	 for(unsigned int ig=0; ig<gen.size(); ig++){ 
	     if(!ig.isHardProcess()) continue; 
	     
	     const reco::Candidate* imom = findFirstMotherWithDifferentID(&ig);
	     if (imom) {
	       if ( mom->p4() == ig.p4() && idxx==0) {
		 idxx=idx; 
	       }
	       idx++;      
	     }
	   }
	   
	   ev.mc_momidx[ev.nmcparticles] = idxx; 
	   ev.mc_status[ev.nmcparticles] = it.status();
	   
	   TLorentzVector p4( it.px(), it.py(), it.pz(), it.energy() );
	   if(abs(pid)==11 || abs(pid)==13 || abs(pid)==15) {
	     chLeptons.push_back(p4);
	   }
	   ev.nmcparticles++;
	   
	   if ( verbose_ ) {
	     printf("  %3d : ID=%6d, m=%5.1f, momID=%6d : pt=%6.1f, eta=%7.3f, phi=%7.3f, status=%d\n",
		    igen,
		    it.pdgId(),
		    it.mass(),
		    mom->pdgId(),
		    it.pt(),
		    it.eta(),
		    it.phi(),
		    it.status()
		    ) ;
	   }
	   
	 } // has mom?
	 igen++;
       } // if is hardProcess
     } // igen
     
     if ( verbose_ ) {
       printf("\n\n ---- ground state B hadrons:\n") ;
       int bi(0);
       for (auto & it : b_hadrons) {
       //       for ( int bi=0; bi<b_hadrons.size(); bi++ ) {
	 printf( " %2d %p : ID=%6d : m=%6.2f : pt=%6.1f, eta=%7.3f, phi=%7.3f\n",
                 bi, it, it.pdgId(), it.mass(), it.pt(), it.eta(), it.phi()) ;
	 bi++;
       } // bi
       printf("\n\n") ;
     }
     
     // //
     // // gen jets
     // //
     // ev.nmcjparticles = 0;  
     
     // reco::GenJetCollection genJets;
     // edm::Handle< reco::GenJetCollection > genJetsHandle;
     // event.getByToken(genjetTag_, genJetsHandle);
     // if(genJetsHandle.isValid()){ genJets = *genJetsHandle;}
     
     // std::vector<TLorentzVector> jets;
     // //   for(size_t j=0; j<genJets.size(); j++) {
     // for (auto & it : genJets ) {
     //   const reco::GenJet genJet = it;
       
     //   TLorentzVector p4( genJet.px(), genJet.py(), genJet.pz(), genJet.energy() );
     //   if(p4.Pt()<10 || fabs(p4.Eta())>2.5) continue;
       
     //   bool matchesLepton(false);
     //   for(size_t i=0; i<chLeptons.size(); i++) {
     // 	 float dR=p4.DeltaR(chLeptons[i]);
     // 	 if(dR>0.4) continue;
     // 	 matchesLepton=true;
     // 	 break;
     //   }
     //   if(matchesLepton) continue;
       
     //   jets.push_back(p4);
     //   ev.mcj_px[ev.nmcjparticles]=genJet.px();
     //   ev.mcj_py[ev.nmcjparticles]=genJet.py();
     //   ev.mcj_pz[ev.nmcjparticles]=genJet.pz();
     //   ev.mcj_en[ev.nmcjparticles]=genJet.energy();
     //   ev.mcj_status[ev.nmcjparticles]=0; // special for genjet
     //   ev.mcj_id[ev.nmcjparticles]=1; // special for genjet
     //   ev.mcj_mom[ev.nmcjparticles] = 0;
     //   ev.nmcjparticles++;
     // }
     
   } // end MC
   
 //
       // Trigger
       //
   edm::Handle<edm::TriggerResults> triggerBits; 
   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects; 
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   event.getByToken(triggerBits_, triggerBits);
   event.getByToken(triggerObjects_, triggerObjects);
   event.getByToken(triggerPrescales_, triggerPrescales);

   bool debug_=false;
   if (debug_) {
     
     const edm::TriggerNames &names = event.triggerNames(*triggerBits); 
     
     std::cout << "\n === TRIGGER PATHS === " << std::endl; 
     for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) { 
       std::cout << "Trigger " << names.triggerName(i) <<
	 ", prescale " << triggerPrescales->getPrescaleForIndex(i) << 
	 ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
		 << std::endl; 
     }        
   }
       
   //apply trigger and require compatibilitiy of the event with the PD
   edm::TriggerResultsByName tr(nullptr,nullptr);
#ifdef YEAR_2017
   tr = event.triggerResultsByName(*triggerBits);
#else
   tr = event.triggerResultsByName("HLT");
#endif
   if(!tr.isValid()  )return;
   
   // float triggerPrescale(1.0),triggerThreshold(0), triggerThresholdHigh(99999);
   bool mumuTrigger(true); bool muTrigger(true); bool muTrigger2(true);
   bool eeTrigger(true); bool eTrigger(true); bool emuTrigger(true);
   bool highPTeTrigger(true);

#ifdef YEAR_2017
//https://indico.cern.ch/event/682891/contributions/2810364/attachments/1570825/2820752/20171206_CMSWeek_MuonHLTReport_KPLee_v3_4.pdf
   mumuTrigger        = utils::passTriggerPatterns(tr,"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
   muTrigger          = utils::passTriggerPatterns(tr,"HLT_IsoMu24_v*","HLT_IsoMu24_eta2p1_v*","HLT_IsoMu27_v*");
// For 2017, muTrigger2 is always set to True
//   muTrigger2         = utils::passTriggerPatterns(tr, "HLT_Mu50_v*", "HLT_TkMu50_v*");
//   https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#E/gamma%20Trigger%20Recomendations
   eeTrigger          = utils::passTriggerPatterns(tr,"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_DoubleEle33_CaloIdL_MW_v*","HLT_DoubleEle25_CaloIdL_MW_v*");
   eTrigger           = utils::passTriggerPatterns(tr,"HLT_Ele32_WPTight_Gsf_v*","HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*");
   emuTrigger         = utils::passTriggerPatterns(tr,"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*") || utils::passTriggerPatterns(tr,"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
#else
   mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*" , "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
   //   muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu22_v*","HLT_IsoTkMu22_v*", "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*");
   muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*");
   muTrigger2         = utils::passTriggerPatterns(tr, "HLT_Mu50_v*", "HLT_TkMu50_v*");
   eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DoubleEle33_CaloIdL_v*");
   //   highPTeeTrigger    = utils::passTriggerPatterns(tr, "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"); //"HLT_ECALHT800_v*");
   highPTeTrigger    = utils::passTriggerPatterns(tr, "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*");
   eTrigger           = utils::passTriggerPatterns(tr, "HLT_Ele27_WPTight_Gsf_v*") ;
   emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*" , "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*") || utils::passTriggerPatterns(tr,"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
#endif 
   
   ev.hasTrigger  = ( mumuTrigger||muTrigger||muTrigger2||eeTrigger||highPTeTrigger||eTrigger||emuTrigger );
   //ev.hasTrigger  = ( muTrigger||eTrigger||emuTrigger ); 
   
   ev.triggerType = ( mumuTrigger  << 0 )
     | ( muTrigger  << 1 )
     | ( eeTrigger << 2 )
     | ( highPTeTrigger << 3 )
     | ( eTrigger << 4 )
     | ( emuTrigger << 5 ) 
     | ( muTrigger2 << 6 );
   
   //if(!isMC__ && !ev.hasTrigger) return; // skip the event if no trigger, only for Data
   if(!ev.hasTrigger) return; // skip the event if no trigger, for both Data and MC
   
   //##############################################   EVENT PASSED THE TRIGGER   ######################################
   //met filters
    edm::Handle<edm::TriggerResults> metFilterBits;
    event.getByToken(metFilterBitsTag_, metFilterBits);
    const edm::TriggerNames &metNames = event.triggerNames(*metFilterBits);
    bool passMETFilters(true);
    for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
        if(strcmp(metNames.triggerName(i).c_str(), 	 "Flag_goodVertices") == 0){
            passMETFilters &= metFilterBits->accept(i);
	    //	}else if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0){
	    //           passMETFilters &= metFilterBits->accept(i);
	}else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0){
            passMETFilters &= metFilterBits->accept(i);
        }else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0){
            passMETFilters &= metFilterBits->accept(i);
        }else if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalTightHalo2016Filter") == 0){
            passMETFilters &= metFilterBits->accept(i);
        }else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0){
            passMETFilters &= metFilterBits->accept(i);
	}else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadChargedCandidateFilter") == 0){
            passMETFilters &= metFilterBits->accept(i);
	}else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonFilter") == 0){
            passMETFilters &= metFilterBits->accept(i);
	}
	
    }
  
    if(!passMETFilters) return;
   //   if( metFilterValue!=0 ) continue;	 //Note this must also be applied on MC
   
   // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
    //   if (!filterbadPFMuon || !filterbadChCandidate) return;
   //##############################################   EVENT PASSED MET FILTER   #######################################
   
  
    //load all the objects we will need to access
       reco::VertexCollection vtx;
       edm::Handle< reco::VertexCollection > vtxHandle;
       event.getByToken(vtxTag_, vtxHandle);
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}
       
       double rho = 0;
       edm::Handle< double > rhoHandle;
       event.getByToken(rhoFastjetAllTag_,rhoHandle);
       if(rhoHandle.isValid()){ rho = *rhoHandle;}

       if (vtx.empty()) return; // skip the event if no PV found
       const reco::Vertex &PV = vtx.front();

       ev.vtx_x = PV.x();
       ev.vtx_y = PV.y();
       ev.vtx_z = PV.z();

       ev.nvtx = 0;
       //select good vertices
       for(unsigned int i = 0; i < vtx.size(); i++) {
	 if(vtx.at(i).isValid() && !vtx.at(i).isFake()) ev.nvtx++;
       }
       if(ev.nvtx == 0) return;


       pat::MuonCollection muons;
       edm::Handle< pat::MuonCollection > muonsHandle;
       event.getByToken(muonTag_, muonsHandle);
       //       muonsHandle.getByLabel(event, "slimmedMuons");
       if(muonsHandle.isValid()){ muons = *muonsHandle;}
       
       ev.mn=0;
       //       for (std::vector<pat::Muon >::const_iterator mu = muons.begin(); mu!=muons.end(); mu++) 
       for(pat::Muon &mu : muons) {
	 if(mu.pt() < 5.) continue;
	 ev.mn_px[ev.mn] = mu.px();
	 ev.mn_py[ev.mn] = mu.py();
	 ev.mn_pz[ev.mn] = mu.pz();
	 ev.mn_en[ev.mn] = mu.energy();
	 ev.mn_id[ev.mn] = 13*mu.charge();
	 /*
	 ev.mn_d0[ev.mn] = -mu.muonBestTrack()->dxy(PV.position());
	 ev.mn_dZ[ev.mn] = mu.muonBestTrack()->dz(PV.position());
	 ev.mn_ip3d[ev.mn] = mu.dB(pat::Muon::PV3D);
	 ev.mn_ip3dsig[ev.mn] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);

	 ev.mn_IsLoose[ev.mn] = mu.isLooseMuon();
	 ev.mn_IsMedium[ev.mn] = mu.isMediumMuon();
	 ev.mn_IsTight[ev.mn] = mu.isTightMuon(PV);
	 ev.mn_IsSoft[ev.mn] = mu.isSoftMuon(PV);
	 ev.mn_IsHighPt[ev.mn] = mu.isHighPtMuon(PV);

	 ev.mn_pileupIsoR03[ev.mn] = mu.pfIsolationR03().sumPUPt;
	 ev.mn_chargedIsoR03[ev.mn] = mu.pfIsolationR03().sumChargedHadronPt;
	 ev.mn_photonIsoR03[ev.mn] = mu.pfIsolationR03().sumPhotonEt;
	 ev.mn_neutralHadIsoR03[ev.mn] = mu.pfIsolationR03().sumNeutralHadronEt;

	 ev.mn_pileupIsoR04[ev.mn] = mu.pfIsolationR04().sumPUPt;
	 ev.mn_chargedIsoR04[ev.mn] = mu.pfIsolationR04().sumChargedHadronPt;
	 ev.mn_photonIsoR04[ev.mn] = mu.pfIsolationR04().sumPhotonEt;
	 ev.mn_neutralHadIsoR04[ev.mn] = mu.pfIsolationR04().sumNeutralHadronEt;
	 
	 ev.mn_nMatches[ev.mn]                   = mu.numberOfMatches();
	 ev.mn_nMatchedStations[ev.mn]           = mu.numberOfMatchedStations();
	 ev.mn_validMuonHits[ev.mn]              = mu.isGlobalMuon() ? mu.globalTrack().hitPattern().numberOfValidMuonHits() : 0.;
	 ev.mn_innerTrackChi2[ev.mn]             = mu.isTrackerMuon() ? mu.innerTrack().normalizedChi2() : 0.;
	 */
	 ev.mn_validMuonHits[ev.mn]              = mu.isGlobalMuon() ? mu.globalTrack()->hitPattern().numberOfValidMuonHits() : 0.;    
	 ev.mn_trkLayersWithMeasurement[ev.mn]   = mu.isTrackerMuon() ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : 0.;
	 ev.mn_pixelLayersWithMeasurement[ev.mn] = mu.isTrackerMuon() ? mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() : 0.;

         float relIso_mu = -1, trkrelIso = -1;
	 ev.mn_passId[ev.mn]  = patUtils::passId(mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::ICHEP16Cut);
	 ev.mn_passIdLoose[ev.mn] = patUtils::passId(mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::ICHEP16Cut);
	 ev.mn_passSoftMuon[ev.mn] = patUtils::passId(mu, vtx[0], patUtils::llvvMuonId::Soft, patUtils::CutVersion::ICHEP16Cut);
	 ev.mn_passIso[ev.mn] = patUtils::passIso(mu, patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut, &relIso_mu, &trkrelIso);
         ev.mn_relIso[ev.mn] = relIso_mu;
         ev.mn_trkrelIso[ev.mn] = trkrelIso;

	 ev.mn_type[ev.mn]   = (mu.isMuon() << 0)
	   | (mu.isGlobalMuon() << 1)
	   | (mu.isTrackerMuon() << 2)
	   | (mu.isStandAloneMuon()<< 3)
	   | (mu.isCaloMuon() << 4)
	   | (mu.isPFMuon() << 5)
	   | (mu.isRPCMuon()<< 6);
	 ev.mn++;
       } // mu


       pat::ElectronCollection electrons;
       edm::Handle< pat::ElectronCollection > electronsHandle;
       event.getByToken(electronTag_, electronsHandle);
       if(electronsHandle.isValid()){ electrons = *electronsHandle;}

       // to find Gain Seed 
       edm::Handle<EcalRecHitCollection> recHitCollectionEBHandle;
       edm::Handle<EcalRecHitCollection> recHitCollectionEEHandle;
       event.getByToken(reducedEBRecHitCollectionToken_, recHitCollectionEBHandle); // "reducedEgamma","reducedEBRecHits");
       event.getByToken(reducedEERecHitCollectionToken_, recHitCollectionEEHandle); //"reducedEgamma","reducedEERecHits");

       ev.en=0;
 
       for (pat::Electron &el : electrons) {
       //       for( View<pat::ElectronCollection>::const_iterator el = electrons.begin(); el != electrons.end(); el++ ) 
	 float pt_ = el.pt();
	 if (pt_ < 10.) continue;

	 // Kinematics
	 ev.en_px[ev.en] = el.px();
	 ev.en_py[ev.en] = el.py();
	 ev.en_pz[ev.en] = el.pz();
	 ev.en_en[ev.en] = el.energy();
	 ev.en_id[ev.en] = 11*el.charge();

	 //	 float aeta = std::abs(el.superCluster()->eta());
	 float cor_en = el.correctedEcalEnergy() ;           

	 ev.en_cor_en[ev.en] = cor_en;
	 ev.en_EtaSC[ev.en] = el.superCluster()->eta(); 
	 ev.en_R9[ev.en] = el.full5x5_r9();

	 /*
	 //Isolation
	 GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
	 ev.en_pileupIso[ev.en] = pfIso.sumPUPt;
	 ev.en_chargedIso[ev.en] = pfIso.sumChargedHadronPt;
	 ev.en_photonIso[ev.en] = pfIso.sumPhotonEt;
	 ev.en_neutralHadIso[ev.en] = pfIso.sumNeutralHadronEt;
	 */
         float relIso_el = -1;
	 ev.en_passId[ev.en] = patUtils::passId(el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut);
	 ev.en_passIdLoose[ev.en] = patUtils::passId(el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::ICHEP16Cut);
	 ev.en_passIso[ev.en] = patUtils::passIso(el, patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, &relIso_el, rho) ;
         ev.en_relIso[ev.en] = relIso_el;

	 const EcalRecHitCollection* recHits = (el.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product();
	 unsigned int gainSeed = patUtils::GainSeed(el,recHits);

	 ev.en_gainSeed[ev.en] = gainSeed;

	 /*
	 if(isMC_){
	   double sigma= eScaler_.getSmearingSigma(event.eventAuxiliary().run(),el.isEB(),el.full5x5_r9(), el.superCluster()->eta(), el.et(),gainSeed,0,0);
	   //Put the last two inputs at 0,0 for the nominal value of sigma
	   //Now smear the MC energy
	   TRandom3 *rgen_ = new TRandom3(0);
	   ev.en_scale_corr[ev.en] = (float)rgen_->Gaus(1, sigma) ;
	   //E_new=E_old*(rgen_->Gaus(1, sigma)) ;
	 } else {
	   ev.en_scale_corr[ev.en] = (float)eScaler_.ScaleCorrection(event.eventAuxiliary().run(),el.isEB(),el.full5x5_r9(), el.superCluster()->eta(), el.et(),gainSeed);
	   //At this point, the new data energy will be:
	   //	   E_new=E_old*(scale_corr); 
	 }
	 */
	 ev.en++;
       } // el


          //
       // jet selection (ak4PFJetsCHS)
       //
       pat::JetCollection jets;
       edm::Handle< pat::JetCollection > jetsHandle;
       event.getByToken(jetTag_, jetsHandle);
       if(jetsHandle.isValid()){ jets = *jetsHandle;}
       
       ev.jet=0;
       //       PFJetIDSelectionFunctor looseJetIdSelector(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::LOOSE);
       //       PFJetIDSelectionFunctor tightJetIdSelector(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::TIGHT);
       //       pat::strbitset hasLooseId =looseJetIdSelector.getBitTemplate();   
       //       bool passPFloose = patUtils::passPFJetID("Loose", jet); //looseJetIdSelector.getBitTemplate();
       //       pat::strbitset hasTightId = tightJetIdSelector.getBitTemplate();

       if ( verbose_ ) printf("\n\n ----- Reconstructed jets : %lu\n", jets.size() ) ;

       //for (std::vector<pat::Jet >::const_iterator j = jets.begin(); j!=jets.end(); j++) 
       int ijet(0) ;
       for (pat::Jet &j : jets) {
	 if(j.pt() < 20) continue;

	 //jet id
	 //	 hasLooseId.set(false);
	 //hasTightId.set(false);
	 //bool passLooseId(looseJetIdSelector( j, hasLooseId ));
	 //bool passTightId(tightJetIdSelector( j, hasTightId ));
	 ev.jet_PFLoose[ev.jet] = patUtils::passPFJetID("Loose", j);
	 ev.jet_PFTight[ev.jet] = patUtils::passPFJetID("Tight", j);

	 ev.jet_px[ev.jet] = j.px(); //correctedP4(0).px();
	 ev.jet_py[ev.jet] = j.py(); //correctedP4(0).py();
	 ev.jet_pz[ev.jet] = j.pz(); //correctedP4(0).pz();
	 ev.jet_en[ev.jet] = j.energy(); //correctedP4(0).energy();

	 ev.jet_btag0[ev.jet] = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
#ifdef YEAR_2017
   ev.jet_btag1[ev.jet] = j.bDiscriminator("pfDeepCSVJetTags:probb") + j.bDiscriminator("pfDeepCSVJetTags:probbb");
#else
   ev.jet_btag1[ev.jet] = j.bDiscriminator("deepFlavourJetTags:probb") + j.bDiscriminator("deepFlavourJetTags:probbb");
#endif
   // ev.jet_btag1[ev.jet] = j.bDiscriminator("pfJetBProbabilityBJetTags");
	 // ev.jet_btag2[ev.jet] = j.bDiscriminator("pfJetProbabilityBJetTags");
	 // ev.jet_btag3[ev.jet] = j.bDiscriminator("pfTrackCountingHighPurBJetTags");
	 // ev.jet_btag4[ev.jet] = j.bDiscriminator("pfTrackCountingHighEffBJetTags");
	 // ev.jet_btag5[ev.jet] = j.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
	 // ev.jet_btag6[ev.jet] = j.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags");
	 // ev.jet_btag7[ev.jet] = j.bDiscriminator("combinedSecondaryVertexBJetTags");
	 // ev.jet_btag8[ev.jet] = j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
	 // ev.jet_btag9[ev.jet] = j.bDiscriminator("pfCombinedSecondaryVertexSoftLeptonBJetTags");
	 // ev.jet_btag10[ev.jet] = j.bDiscriminator("pfCombinedMVABJetTags");

	 ev.jet_mass[ev.jet] = j.mass(); //correctedP4(0).M();
	 ev.jet_area[ev.jet] = j.jetArea();
	 ev.jet_pu[ev.jet] = j.pileup();
	 ev.jet_puId[ev.jet] = j.userFloat("pileupJetId:fullDiscriminant");

         if ( verbose_ ) {
            printf("    %2d : pt=%6.1f, eta=%7.3f, phi=%7.3f : ID=%s%s, bCSV=%7.3f, PUID=%7.3f\n",
                ijet, j.pt(), j.eta(), j.phi(),
                (patUtils::passPFJetID("Loose", j)?"L":" "),
                (patUtils::passPFJetID("Tight", j)?"T":" "),
                j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                j.userFloat("pileupJetId:fullDiscriminant")
             ) ;
         }
	 

	 ev.jet_mother_id[ev.jet] = 0;
	 
	 ev.jet_parton_px[ev.jet] = 0.; 
	 ev.jet_parton_py[ev.jet] = 0.; 
	 ev.jet_parton_pz[ev.jet] = 0.; 
	 ev.jet_parton_en[ev.jet] = 0.; 

	 ev.jet_partonFlavour[ev.jet] = 0;
	 ev.jet_hadronFlavour[ev.jet] = 0;
	 ev.jet_genpt[ev.jet] = 0.;
	 
	 if (isMC_) {

	   const reco::GenParticle *pJet = j.genParton();
	   if (pJet) {
	     const reco::Candidate* mom = findFirstMotherWithDifferentID(&(*pJet));
	     if (mom) {
	       ev.jet_mother_id[ev.jet] = mom->pdgId();

	       ev.jet_parton_px[ev.jet] = pJet->px();
	       ev.jet_parton_py[ev.jet] = pJet->py(); 
	       ev.jet_parton_pz[ev.jet] = pJet->pz(); 
	       ev.jet_parton_en[ev.jet] = pJet->energy();
	       
	     }
	   }

	   ev.jet_partonFlavour[ev.jet] = j.partonFlavour();
	   ev.jet_hadronFlavour[ev.jet] = j.hadronFlavour();
	   const reco::GenJet *gJet=j.genJet(); 
	   if(gJet) ev.jet_genpt[ev.jet] = gJet->pt();
	   else     ev.jet_genpt[ev.jet] = 0.;
	 } // isMC_
	 
	 ev.jet++;
         ijet++ ;
       }

       //
       // jet selection (AK8Jets)
       //

       pat::JetCollection fatjets; 
       edm::Handle< pat::JetCollection > jetsAK8Handle; 
       event.getByToken(fatjetTag_, jetsAK8Handle);
       if(jetsAK8Handle.isValid()){ fatjets = *jetsAK8Handle;}   

       ev.fjet=0;
       
       if ( verbose_ ) printf("\n\n ----- Reconstructed AK8 jets : %lu\n", fatjets.size() ) ;  

       int ifjet(0) ; 
       for (const pat::Jet &j : fatjets) {
	 ev.fjet_px[ev.fjet] = j.correctedP4(0).px();
	 ev.fjet_py[ev.fjet] = j.correctedP4(0).py();
	 ev.fjet_pz[ev.fjet] = j.correctedP4(0).pz();
	 ev.fjet_en[ev.fjet] = j.correctedP4(0).energy();

	 ev.fjet_btag0[ev.fjet] = j.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
	 ev.fjet_btag1[ev.fjet] = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

	 if ( verbose_ ) {  
	   printf("\n    %3d : pt=%6.1f, eta=%7.3f, phi=%7.3f : boostedSV=%7.3f, combinedSV=%7.3f\n", 
		  ifjet, j.pt(), j.eta(), j.phi(),
		  j.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"),
      j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")
		  ) ;
         }  
	 
#ifdef YEAR_2017
   ev.fjet_prunedM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass");
   ev.fjet_softdropM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass");
   ev.fjet_tau1[ev.fjet] =  (float) j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1");
   ev.fjet_tau2[ev.fjet] =  (float) j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2");
   ev.fjet_tau3[ev.fjet] =  (float) j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3");
#else
	 ev.fjet_prunedM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSPrunedMass");
	 ev.fjet_softdropM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSSoftDropMass");
	 //	 ev.fjet_filteredM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
	 ev.fjet_tau1[ev.fjet] =  (float) j.userFloat("NjettinessAK8CHS:tau1");
	 ev.fjet_tau2[ev.fjet] =  (float) j.userFloat("NjettinessAK8CHS:tau2");
	 ev.fjet_tau3[ev.fjet] =  (float) j.userFloat("NjettinessAK8CHS:tau3");
#endif

	 // Add soft drop subjets

	 if ( verbose_ ) {
	   printf("\n   This AK8 jet has N = %3d subjet collections\n",j.nSubjetCollections());

	   std::vector<std::string> const & sjNames = j.subjetCollectionNames();
	   for (auto const & it : sjNames) {
	     printf("   Subjet collection name %s , ", it.c_str());
	   } 
	   printf("\n") ; 
	 }

#ifdef YEAR_2017
   auto const & sdSubjets = j.subjets("SoftDropPuppi");  
#else
   auto const & sdSubjets = j.subjets("SoftDrop"); 
#endif
	 //The Soft Drop Subjets are stored in positions 0  in the subjet collection list.
	 int count_subj(0);

	 std::vector<TLorentzVector> softdrop_subjets; 
	 for ( auto const & it : sdSubjets ) {

	   TLorentzVector softdrop_subjet;
	   softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(), it->correctedP4(0).eta(), it->correctedP4(0).phi(), it->correctedP4(0).mass()); 
	   
	   softdrop_subjets.push_back(softdrop_subjet);

	   count_subj++;
	   if ( verbose_ ) {
	     
	     printf("\n    Soft drop subjet  %3d : pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f",
		    count_subj,
		    softdrop_subjet.Pt(),
		    softdrop_subjet.Eta(),
		    softdrop_subjet.Phi(),
		    softdrop_subjet.M()
		    );

	   }
	 } // subjets

	 for (int i=0; i<2; i++) { // store up to 4 subjets for each AK8 jet ?
	   ev.fjet_subjets_px[ev.fjet][i] = 0.;
	   ev.fjet_subjets_py[ev.fjet][i] = 0.;
	   ev.fjet_subjets_pz[ev.fjet][i] = 0.;
	   ev.fjet_subjets_en[ev.fjet][i] = 0.;
	 }
       
	 int csb(0);
	 for ( auto & it : softdrop_subjets ) {
	   if (csb>1) break; // up to 2 subjets  
	   //	   if (it.Pt()>20.) { // only store subjets above 20 GeV ?
	   ev.fjet_subjets_px[ev.fjet][csb] = it.Px();
	   ev.fjet_subjets_py[ev.fjet][csb] = it.Py();
	   ev.fjet_subjets_pz[ev.fjet][csb] = it.Pz();
	   ev.fjet_subjets_en[ev.fjet][csb] = it.E();
	   
	     csb++;
	 }
	 ev.fjet_subjet_count[ev.fjet] = csb;


	 ev.fjet_partonFlavour[ev.fjet] = 0;
	 ev.fjet_hadronFlavour[ev.fjet] = 0;

	 ev.fjet_mother_id[ev.fjet] = 0; 

	 ev.fjet_parton_px[ev.fjet] = 0.; 
	 ev.fjet_parton_py[ev.fjet] = 0.; 
	 ev.fjet_parton_pz[ev.fjet] = 0.; 
	 ev.fjet_parton_en[ev.fjet] = 0.; 
	 
	 if (isMC_) {
#ifndef YEAR_2017
	   const reco::GenParticle *pJet = j.genParton();
	   if (pJet) {
	     const reco::Candidate* mom = findFirstMotherWithDifferentID(&(*pJet));
	     if (mom) {
	       ev.fjet_mother_id[ev.fjet] = mom->pdgId(); 

	       ev.fjet_parton_px[ev.fjet] = pJet->px();
	       ev.fjet_parton_py[ev.fjet] = pJet->py();
	       ev.fjet_parton_pz[ev.fjet] = pJet->pz();
	       ev.fjet_parton_en[ev.fjet] = pJet->energy();
	     }
	   }

#endif
	   ev.fjet_partonFlavour[ev.fjet] = j.partonFlavour();
	   ev.fjet_hadronFlavour[ev.fjet] = j.hadronFlavour();

	   const reco::GenJet *gJet=j.genJet();
	   if(gJet) ev.fjet_genpt[ev.fjet] = gJet->pt();
	   else     ev.fjet_genpt[ev.fjet] = 0.;
	   
	 } // isMC_

	 ev.fjet++;
	 ifjet++;
       }
       
     pat::METCollection mets;
     edm::Handle< pat::METCollection > metsHandle;
     // if(isMC_)
     event.getByToken(metTag_, metsHandle);
     //     if(!isMC_)event.getByToken(metTagData_, metsHandle);
     if(metsHandle.isValid()){ mets = *metsHandle;}
     pat::MET met = mets[0];

     //       const pat::MET &met = mets.front();

       //MET CORRection level
       pat::MET::METCorrectionLevel metcor = pat::MET::METCorrectionLevel::Type1XY;
       //recompute MET with variation
       //       LorentzVector imet[];
       ev.imet_pt[0] = met.corP4(metcor).pt();
       ev.imet_phi[0] = met.corP4(metcor).phi();
       ev.imet_pt[1] = met.shiftedP4(pat::MET::METUncertainty::JetEnUp           , metcor).pt();
       ev.imet_phi[1] = met.shiftedP4(pat::MET::METUncertainty::JetEnUp           , metcor).phi();   
       ev.imet_pt[2] = met.shiftedP4(pat::MET::METUncertainty::JetEnDown         , metcor).pt();
       ev.imet_phi[2] = met.shiftedP4(pat::MET::METUncertainty::JetEnDown         , metcor).phi();    
       ev.imet_pt[3] = met.shiftedP4(pat::MET::METUncertainty::JetResUp          , metcor).pt();
       ev.imet_phi[3] = met.shiftedP4(pat::MET::METUncertainty::JetResUp          , metcor).phi();  
       ev.imet_pt[4] = met.shiftedP4(pat::MET::METUncertainty::JetResDown        , metcor).pt();
       ev.imet_phi[4] = met.shiftedP4(pat::MET::METUncertainty::JetResDown        , metcor).phi();  
       ev.imet_pt[5] = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp   , metcor).pt();
       ev.imet_phi[5] = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp   , metcor).phi();   
       ev.imet_pt[6] = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown , metcor).pt();
       ev.imet_phi[6] = met.shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown , metcor).phi();  
       ev.imet_pt[7] = met.shiftedP4(pat::MET::METUncertainty::MuonEnUp          , metcor).pt();
       ev.imet_phi[7] = met.shiftedP4(pat::MET::METUncertainty::MuonEnUp          , metcor).phi(); 
       ev.imet_pt[8] = met.shiftedP4(pat::MET::METUncertainty::MuonEnDown        , metcor).pt();
       ev.imet_phi[8] = met.shiftedP4(pat::MET::METUncertainty::MuonEnDown        , metcor).phi();     
       ev.imet_pt[9] = met.shiftedP4(pat::MET::METUncertainty::ElectronEnUp      , metcor).pt();
       ev.imet_phi[9] = met.shiftedP4(pat::MET::METUncertainty::ElectronEnUp      , metcor).phi(); 
       ev.imet_pt[10] = met.shiftedP4(pat::MET::METUncertainty::ElectronEnDown , metcor).pt();
       ev.imet_phi[10] = met.shiftedP4(pat::MET::METUncertainty::ElectronEnDown , metcor).phi();

       //PF type-1 ETmiss
       ev.met_pt = met.pt();
       ev.met_phi = met.phi();
       ev.met_sumMET = met.sumEt();

       // raw PF ETmiss
       ev.rawpfmet_pt = met.uncorPt();
       ev.rawpfmet_phi = met.uncorPhi();
       ev.rawpfmet_sumMET = met.uncorSumEt();

       // raw calo ETmiss
       ev.rawcalomet_pt = met.caloMETPt();
       ev.rawcalomet_phi = met.caloMETPhi();
       ev.rawcalomet_sumMET = met.caloMETSumEt();
       /*
       // type1 PF MET but excluding HF
       edm::Handle<pat::METCollection> metsNoHF;
       event.getByToken(metNoHFTag_, metsNoHF);
       if(metsNoHF.isValid()) {
	 const pat::MET &metNoHF = metsNoHF->front();
	 ev.metNoHF_pt = metNoHF.pt();
	 ev.metNoHF_phi = metNoHF.phi();
	 ev.metNoHF_sumMET = metNoHF.sumEt();
       }
       */

       //-- Inclusive Secondary Vertices

       reco::VertexCompositePtrCandidateCollection sec_vert ;
       edm::Handle< reco::VertexCompositePtrCandidateCollection > svHandle ;
       event.getByToken(svTag_, svHandle);
       if ( svHandle.isValid() ) { sec_vert = *svHandle ;} 
       else { printf("\n\n *** bad handle for reco::VertexCompositePtrCandidateCollection\n\n") ; gSystem -> Exit(-1) ; }

       ev.sv = sec_vert.size() ;
       if ( verbose_ ) printf("\n\n\n ---- Inclusive Secondary Vertices:\n" ) ;
       for ( unsigned int isv=0; isv<sec_vert.size(); isv++ ) {

          if (verbose_ ) {
            printf(" %3d :   x,y,z = %9.5f, %9.5f, %9.5f :  Ntrk = %2lu : chi2 = %7.3f, Ndof = %5.2f\n",
             isv,
             sec_vert[isv].position().x(),
             sec_vert[isv].position().y(),
             sec_vert[isv].position().z(),
             sec_vert[isv].numberOfDaughters(),
             sec_vert[isv].vertexChi2(),
             sec_vert[isv].vertexNdof()
             ) ;
          }

          ev.sv_chi2[isv] = sec_vert[isv].vertexChi2() ;
          ev.sv_ndof[isv] = sec_vert[isv].vertexNdof() ;

          ev.sv_dxy[isv]  = sqrt( pow((sec_vert[isv].position().x() - PV.x()),2) + pow((sec_vert[isv].position().y() - PV.y()),2) ) ;
          ev.sv_dxyz[isv] = sqrt( pow((sec_vert[isv].position().x() - PV.x()),2) + pow((sec_vert[isv].position().y() - PV.y()),2) + pow((sec_vert[isv].position().z() - PV.z()),2) ) ;

          if ( verbose_ ) {
             printf("      dx,dy,dz = %9.5f, %9.5f, %9.5f :  dxy = %9.5f , dxyz = %9.5f\n",
               sec_vert[isv].position().x() - PV.x(),
               sec_vert[isv].position().y() - PV.y(),
               sec_vert[isv].position().z() - PV.z(),
               ev.sv_dxy[isv],
               ev.sv_dxyz[isv]
             ) ;
          }


          TLorentzVector sv_p4 ;

          for ( unsigned int id=0; id<sec_vert[isv].numberOfDaughters(); id++ ) {

             reco::CandidatePtr dau = sec_vert[isv].daughterPtr(id) ;
             TLorentzVector svd_p4( dau->px(), dau->py(), dau->pz(), dau->energy() ) ;
             sv_p4 += svd_p4 ;

             if ( verbose_ ) {
                printf("      trk %2d :  pt=%6.1f, eta=%7.3f, phi = %7.3f\n",
                   id,
                   dau->pt(),
                   dau->eta(),
                   dau->phi()
                ) ;
             }

          } // id

          ev.sv_ntrk[isv] = sec_vert[isv].numberOfDaughters() ;

          ev.sv_px[isv] = sv_p4.Px() ;
          ev.sv_py[isv] = sv_p4.Py() ;
          ev.sv_pz[isv] = sv_p4.Pz() ;
          ev.sv_en[isv] = sv_p4.E() ;

          GlobalVector dxyz( sec_vert[isv].position().x() - PV.x(), sec_vert[isv].position().y() - PV.y(), sec_vert[isv].position().z() - PV.z() ) ;
          GlobalVector sv_p3( sv_p4.Px(), sv_p4.Py(), sv_p4.Pz() ) ;
          ev.sv_cos_dxyz_p[isv] = -2. ;
          if ( sv_p3.mag() * dxyz.mag() > 0 ) {
             ev.sv_cos_dxyz_p[isv] = sv_p3.dot( dxyz ) / ( sv_p3.mag() * dxyz.mag() ) ;
          }
          double sv_pt = sqrt( pow( sv_p3.x(), 2. ) + pow( sv_p3.y(), 2. ) ) ;
          if ( verbose_ ) {
             printf("  secondary vertex pt = %6.1f, mass = %6.2f\n", sv_pt, sv_p4.M() ) ;
             printf("  cos(pv,sv) = %6.3f\n", ev.sv_cos_dxyz_p[isv] ) ;
          }

          const reco::Vertex sv( sec_vert[isv].position(), sec_vert[isv].error() ) ;
          Measurement1D projected_flight_length = reco::SecondaryVertex::computeDist3d( PV, sv, sv_p3, true ) ;
          if ( verbose_ ) {
             printf("  projected flight length: val = %9.5f , err = %9.5f , signif = %9.5f\n",
                 projected_flight_length.value(), projected_flight_length.error(), projected_flight_length.significance() ) ;
          }
          ev.sv_dxyz_signif[isv] = projected_flight_length.significance() ;


          ev.sv_mc_nbh_moms[isv] = 0 ;
          ev.sv_mc_nbh_daus[isv] = 0 ;
          ev.sv_mc_mcbh_ind[isv] = -1 ;

          if (isMC_) {

             //--- go through daughters and find which are from B hadrons

             std::vector<int> bhmomind ;
             std::vector<int> bhmom_ndau ;

             for ( unsigned int id=0; id<sec_vert[isv].numberOfDaughters(); id++ ) {

                reco::CandidatePtr dau = sec_vert[isv].daughterPtr(id) ;

                //--- find closest match to a charged particle in packedGenParticles
		
                edm::View<pat::PackedGenParticle> packed_gen;
                edm::Handle< edm::View<pat::PackedGenParticle> > packed_genHandle;
		event.getByToken(packedGenTag_,packed_genHandle);
                if(packed_genHandle.isValid()){ packed_gen = *packed_genHandle; } else { printf("\n\n *** bad handle for pat::PackedGenParticleCollection\n\n") ; gSystem->Exit(-1) ; }

                double minDr = 9999. ;
                int match_ipgp = -1 ;

                for ( unsigned int ipgp=0; ipgp<packed_gen.size(); ipgp++ ) {

                   if ( packed_gen[ipgp].charge() == 0 ) continue ;

                   double gp_eta = packed_gen[ipgp].eta() ;
                   double gp_phi = packed_gen[ipgp].phi() ;

                   double deta = fabs( dau->eta() - gp_eta ) ;
                   double dphi = fabs( dau->phi() - gp_phi ) ;
                   if ( dphi > 3.14159265 ) dphi -= 2*3.14159265 ;
                   if ( dphi <-3.14159265 ) dphi += 2*3.14159265 ;
                   double dr = sqrt( dphi*dphi + deta*deta ) ;
                   if ( dr < minDr ) {
                      minDr = dr ;
                      match_ipgp = ipgp ;
                   }

                } // ipgp

                if ( minDr < 0.02 && match_ipgp >= 0 ) {
                   if ( verbose_ ) printf("  trk %2d matches pgp %3d , dr = %.4f", id, match_ipgp, minDr ) ; 
		   int bi(0);
		   for (auto & it : b_hadrons) {
		   //   for ( int bi=0; bi<b_hadrons.size(); bi++ ) {
                      if ( isAncestor( it, packed_gen[match_ipgp] ) ) {
                         if ( verbose_ ) printf(" : is daughter of B had %d (pt=%5.1f, eta=%7.3f, phi=%7.3f)",
                          bi, it.pt(), it.eta(), it.phi()) ;
                         ev.sv_mc_nbh_daus[isv] ++ ;
                         if ( ev.sv_mc_mcbh_ind[isv] < 0 && ev.sv_mc_nbh_moms[isv] == 0 ) {
                            ev.sv_mc_nbh_moms[isv] = 1 ;
                            ev.sv_mc_mcbh_ind[isv] = bi ;
                            bhmomind.emplace_back( bi ) ;
                            bhmom_ndau.emplace_back( 1 ) ;
                         } else if ( ev.sv_mc_nbh_moms[isv] >= 1 ) {
                            bool found(false) ;
                            for ( unsigned int i=0; i<bhmomind.size(); i++ ) {
                               if ( bhmomind[i] == bi ) {
                                  found = true ;
                                  bhmom_ndau[i] ++ ;
                                  break ;
                               }
                            } // i
                            if ( !found ) {
                               ev.sv_mc_nbh_moms[isv] ++ ;
                               bhmomind.emplace_back( bi ) ;
                               bhmom_ndau.emplace_back( 1 ) ;
                            }
                         }
                      }
		      bi++;
                   } // bi
                   if ( verbose_ ) printf("\n") ;
                } // found a reasonable match?
                

             } // id

             if ( bhmom_ndau.size() > 1 ) {
                int max(0) ;
                float maxpt(0.) ;
                for ( unsigned int i = 0 ; i < bhmom_ndau.size() ; i++ ) {
                   if ( bhmom_ndau[i] > max ) {
                      max = bhmom_ndau[i] ;
                      maxpt = sqrt( pow(ev.mcbh_px[bhmomind[i]], 2. ) + pow(ev.mcbh_py[bhmomind[i]], 2. ) ) ;
                      ev.sv_mc_mcbh_ind[isv] = bhmomind[i] ;
                   } else if ( bhmom_ndau[i] == max ) {
                      float thispt = sqrt( pow(ev.mcbh_px[bhmomind[i]], 2. ) + pow(ev.mcbh_py[bhmomind[i]], 2. ) ) ;
                      if ( thispt > maxpt ) {
                         max = bhmom_ndau[i] ;
                         maxpt = thispt ;
                         ev.sv_mc_mcbh_ind[isv] = bhmomind[i] ;
                      }
                   }
                } // i
             }

             if ( verbose_ ) {
                if ( ev.sv_mc_nbh_moms[isv] > 0 ) {
                   printf( "         MC b hadron : Nmoms = %d , Ndaus = %d ,  Index of mom with most tracks = %d\n",
                    ev.sv_mc_nbh_moms[isv], ev.sv_mc_nbh_daus[isv], ev.sv_mc_mcbh_ind[isv] ) ;
                } else {
                   printf( "         MC b hadron : Nmoms = 0\n" ) ;
                }
             }

          } // MC?



          if ( verbose_ ) printf("\n") ;

       } // isv

        // Fill Tree
       summaryHandler_.fillTree();
 
}

//========================================================================

const reco::Candidate* mainNtuplizer::findFirstMotherWithDifferentID(const reco::Candidate *particle) 
{ 
  
  if( particle == 0 ) {
    printf("ERROR! null candidate pointer, this should never happen\n"); 
    return 0; 
  }

  // go deeper into recursion 
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) { 
    if (particle->pdgId() == particle->mother(0)->pdgId()) { 
      return findFirstMotherWithDifferentID(particle->mother(0)); 
    } else { 
      return particle->mother(0); 
    } 
  }
  return 0; 
}         

//========================================================================

bool mainNtuplizer::isAncestor( const reco::Candidate& ancestor, const reco::Candidate& daughter ) {

   ////printf(" in isAncestor:  ancestor ID = %d     ,     daughter ID = %d\n", ancestor.pdgId(),  daughter.pdgId() ) ;

  //-- totally stupid way that works...
   if ( ancestor.pdgId() == daughter.pdgId() ) {
      if ( fabs( ancestor.pt() - daughter.pt() ) < 0.1
         && fabs( ancestor.phi() - daughter.phi() ) < 0.01
         && fabs( ancestor.eta() - daughter.eta() ) < 0.01 )
         //printf(" these two are the same.\n" ) ;
         return true ;
   }

   for (size_t i=0; i< daughter.numberOfMothers(); i++) {
      if ( isAncestor( ancestor, *(daughter.mother(i)) ) ) {
         return true ;
      }
   } // i
   return false ;

} // isAncestor


// // ------------ method called once each job just before starting event loop  ------------
void 
mainNtuplizer::beginJob() 
{
}


void 
mainNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{

  if (isMC_) {

    // edm::Handle<LHERunInfoProduct> lheRunInfo; 
    // typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    
    // iRun.getByLabel( "externalLHEProducer", lheRunInfo );
    // LHERunInfoProduct myLHERunInfoProduct = *(lheRunInfo.product());
    
    // for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    //   std::cout << iter->tag() << std::endl;
    //   std::vector<std::string> lines = iter->lines();
    //   for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
    // 	std::cout << lines.at(iLine);
    //   }
    // }
    edm::Handle<LHERunInfoProduct> lheRunInfo;
    iRun.getByLabel(lheRunInfoTag_,lheRunInfo);
    
    if (lheRunInfo.isValid()) {
      
      int pdfidx = lheRunInfo->heprup().PDFSUP.first;
      
      if (pdfidx == 263000) {
	if(verbose_) cout << "NNPDF30_lo_as_0130 (nf5) for LO madgraph samples" << endl;
	//      pdfweightshelper.Init(100,60,edm::FileInPath("llvvAnalysis/DMAnalysis/data/NNPDF30_lo_as_0130_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;
      } else if (pdfidx == 263400) {
	if(verbose_) cout << "NNPdf30_lo_as_0130_nf4 for LO madgraph samples" << endl;
	//      pdfweightshelper.Init(100,60,edm::FileInPath("llvvAnalysis/DMAnalysis/data/NNPDF30_lo_as_0130_nf_4_hessian_60.csv"));
	firstPdfWeight = 111;
	lastPdfWeight = 210;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;
      } else if (pdfidx == 260000 || pdfidx == -1) {
	if(verbose_) cout << "NNPdf30_nlo_as_0118 (nf5) for NLO powheg samples" << endl;
	//(work around apparent bug in current powheg samples by catching "-1" as well)
	//      pdfweightshelper.Init(100,60,edm::FileInPath("llvvAnalysis/DMAnalysis/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	firstPdfWeight = 9;
	lastPdfWeight = 108;
	firstAlphasWeight = 109;
	lastAlphasWeight = 110;
      } else if (pdfidx == 292200) {
	if(verbose_) cout << "NNPdf30_nlo_as_0118 (nf5) with built-in alphas variations for NLO aMC@NLO samples" << endl;
	//      pdfweightshelper.Init(100,60,edm::FileInPath("llvvAnalysis/DMAnalysis/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	firstPdfWeight = 9;
	lastPdfWeight = 108;
	firstAlphasWeight = 109;
	lastAlphasWeight = 110;
      } else if (pdfidx == 292000) {
	if(verbose_) cout << "NNPdf30_nlo_as_0118_nf4 with built-in alphas variations for NLO aMC@NLO samples" << endl;
	//      pdfweightshelper.Init(100,60,edm::FileInPath("llvvAnalysis/DMAnalysis/data/NNPDF30_nlo_as_0118_nf_4_hessian_60.csv"));
	firstPdfWeight = 9;
	lastPdfWeight = 108;
	firstAlphasWeight = 109;
	lastAlphasWeight = 110;
      } else {
	firstPdfWeight = -1;
	lastPdfWeight = -1;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;
      }
      
    } else {
      firstPdfWeight = -1;
      lastPdfWeight = -1;
      firstAlphasWeight = -1;
      lastAlphasWeight = -1;
    }
    
  } else {
    firstPdfWeight = -1;
    lastPdfWeight = -1;
    firstAlphasWeight = -1;
    lastAlphasWeight = -1;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
mainNtuplizer::endJob() 
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
mainNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(mainNtuplizer);
