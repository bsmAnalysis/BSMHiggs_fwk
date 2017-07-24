//#####################################################
// ########## Author: Georgia Karapostoli #############
// ####################################################

#include <iostream>
#include <boost/shared_ptr.hpp>

#include <memory>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

//#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
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

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/PhotonEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
//#include "UserCode/bsmhiggs_fwk/interface/HiggsUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"

#include "UserCode/bsmhiggs_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/bsmhiggs_fwk/interface/TMVAUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/bsmhiggs_fwk/interface/PhotonEfficiencySF.h"
#include "UserCode/bsmhiggs_fwk/interface/PDFInfo.h"
#include "UserCode/bsmhiggs_fwk/interface/rochcor2015.h"
#include "UserCode/bsmhiggs_fwk/interface/rochcor2016.h"
#include "UserCode/bsmhiggs_fwk/interface/muresolution_run2.h"
#include "UserCode/bsmhiggs_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/bsmhiggs_fwk/interface/BtagUncertaintyComputer.h"

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
#include <Math/VectorUtil.h>
#include "TRandom3.h"

#include <time.h>

using namespace std;


const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle) 
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


int main(int argc, char* argv[])
{

  //##############################################
  //########    GLOBAL INITIALIZATION     ########
  //##############################################

  // check arguments
  if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  FWLiteEnabler::enable();

  DataEvtSummaryHandler summaryHandler_;

  // configure the process
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  int mctruthmode=runProcess.getParameter<int>("mctruthmode");
  TString dtag=runProcess.getParameter<std::string>("dtag");

  TString suffix=runProcess.getParameter<std::string>("suffix");
  std::vector<std::string> urls=runProcess.getUntrackedParameter<std::vector<std::string> >("input");
  TString outUrl = runProcess.getParameter<std::string>("outfile");

  //##############################################
  //########    INITIATING TREE      #############
  //##############################################

  fwlite::TFileService fs = fwlite::TFileService("test.root");//outUrl.Data());

  TFileDirectory baseDir=fs.mkdir(runProcess.getParameter<std::string>("dtag"));          
  summaryHandler_.initTree(  fs.make<TTree>("data","Event Summary") );   

  //##############################################
  //########    INITIATING HISTOGRAMS     ########
  //##############################################
  SmartSelectionMonitor mon_;
  printf("Definition of plots");

  mon_.addHistogram(new TH1F("nevents",";nevents; nevents",1,-0.5,0.5));
  mon_.addHistogram(new TH1F("n_negevents",";n_negevents; n_negevents",1,-0.5,0.5));
  mon_.addHistogram(new TH1F("n_posevents",";n_posevents; n_posevents",1,-0.5,0.5));
  mon_.addHistogram(new TH1F("pileup", ";Pileup; Events",100,-0.5,99.5));
  //mon_.addHistogram(new TH1F("integlumi", ";Integrated luminosity ; Events",100,0,1e5);
  //mon_.addHistogram(new TH1F("instlumi", ";Max average inst. luminosity; Events",100,0,1e5);
  mon_.addHistogram(new TH1F("pileuptrue", ";True pileup; Events",100,-0.5,99.5));

  mon_.addHistogram(new TH1F("sumWeights",";;sumWeights;",1,-0.5,0.5));
  mon_.addHistogram(new TH1F("sumScaleWeights",";;sumScaleWeights;",9,-0.5,8.5));
  mon_.addHistogram(new TH1F("sumPdfWeights",";;sumPdfWeights;",100,-0.5,99.5));
  mon_.addHistogram(new TH1F("sumAlphasWeights",";;sumAlphasWeights;",2,-0.5,1.5));

  TH1F *hm=(TH1F*) mon_.addHistogram( new TH1F( "metFilter",";metEventflow",20,0,20) );
  hm->GetXaxis()->SetBinLabel(1,"raw");
  hm->GetXaxis()->SetBinLabel(2,"globalTightHalo2016Filter");
  hm->GetXaxis()->SetBinLabel(3,"goodVertices");
  hm->GetXaxis()->SetBinLabel(4,"eeBadScFilter");
  hm->GetXaxis()->SetBinLabel(5,"EcalDeadCellTriggerPrimitiveFilter");
  hm->GetXaxis()->SetBinLabel(6,"HBHENoiseFilter");
  hm->GetXaxis()->SetBinLabel(7,"HBHENoiseIsoFilter");
  hm->GetXaxis()->SetBinLabel(8,"BadPFMuonFilter");
  hm->GetXaxis()->SetBinLabel(9,"BadChargedCandidateFilte");
  hm->GetXaxis()->SetBinLabel(10,"badMuonHIPFilter");
  hm->GetXaxis()->SetBinLabel(11,"duplicateMuonHIPFilter");
  

  //MC normalization (to 1/pb)
  double xsecWeight = 1.0;
  string debugText = "";

  float curAvgInstLumi_;
  float curIntegLumi_;

  int firstPdfWeight=-1;
  int lastPdfWeight=-1;
  int firstAlphasWeight=-1;
  int lastAlphasWeight=-1;

  //good lumi MASK                                                                                                              
  lumiUtils::GoodLumiFilter goodLumiFilter(runProcess.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess", std::vector<edm::LuminosityBlockRange>())); 

  patUtils::MetFilter metFilter;

  //##############################################
  //########           EVENT LOOP         ########
  //##############################################
  //loop on all the events

  printf("Progressing Bar           :0%%       20%%       40%%       60%%       80%%       100%%\n");

  for(unsigned int f=0;f<urls.size();f++){
     TFile* file = TFile::Open(urls[f].c_str() );
     fwlite::Event event(file);
     printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)urls.size());
     int iev=0;
     int treeStep(event.size()/50);
     if(treeStep==0){ treeStep = 1;}
     for(event.toBegin(); !event.atEnd(); ++event){ 
       iev++;
       if(iev%treeStep==0){printf(".");fflush(stdout);}
       
       mon_.fillHisto("nevents","all",1.0,0); //increment event count

       summaryHandler_.resetStruct();
       //event summary to be filled
       DataEvtSummary_t &ev=summaryHandler_.getEvent();

       float weight = xsecWeight;

       //##############################################   EVENT LOOP STARTS   ##############################################
       //if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }
       
       //Skip bad lumi
       if(!isMC && !goodLumiFilter.isGoodLumi(event.eventAuxiliary().run(),event.eventAuxiliary().luminosityBlock()))continue;

       if (isMC) {

	 ev.nmcparticles = 0;
	 
	 fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
	 puInfoH.getByLabel(event, "slimmedAddPileupInfo");

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
	 mon_.fillHisto("pileup","all",ev.ngenITpu,0);
	 mon_.fillHisto("pileuptrue","all",truePU,0);

	 //retrieve pdf info
	 GenEventInfoProduct eventInfo;
	 fwlite::Handle< GenEventInfoProduct > genEventInfoProd;
	 genEventInfoProd.getByLabel(event, "generator");
	 if(genEventInfoProd.isValid()){ eventInfo = *genEventInfoProd;}                

	 ev.genWeight = eventInfo.weight();
	 float SignGenWeight=1;
	 if(ev.genWeight<0) SignGenWeight=-1;
	 
	 mon_.fillHisto("sumWeights","all",0.,SignGenWeight);
	 ev.qscale = eventInfo.qScale();
	 if(eventInfo.pdf()) {
	   ev.x1  = eventInfo.pdf()->x.first;
	   ev.x2  = eventInfo.pdf()->x.second;
	   ev.id1 = eventInfo.pdf()->id.first;
	   ev.id2 = eventInfo.pdf()->id.second;
	 }
	 if(eventInfo.binningValues().size()>0) ev.pthat = eventInfo.binningValues()[0];

	 if(ev.genWeight<0) mon_.fillHisto("n_negevents","all",1.0,0); //increment negative event count
	 if(ev.genWeight>0) mon_.fillHisto("n_posevents","all",1.0,0); //increment positive event count

	 //scale variations
	 //https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
	 fwlite::Handle<LHEEventProduct> EvtHandle; 
	 EvtHandle.getByLabel( event, "externalLHEProducer");
	 ev.npdfs=0;
	 ev.nalphaS=0;
	 //if(EvtHandles.size()>0) {
	 if(EvtHandle.isValid() && EvtHandle->weights().size()>0) {

	     //fill pdf+alpha_s variation weights
	   if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(EvtHandle->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {

	     //fill pdf variation weights after converting with mc2hessian transformation
	     //std::array<double, 100> inpdfweights;
	     for (int iwgt=firstPdfWeight; iwgt<=lastPdfWeight; ++iwgt) {
	       ev.pdfWeights[ev.npdfs] = SignGenWeight * EvtHandle->weights()[iwgt].wgt/EvtHandle->originalXWGTUP();
	       mon_.fillHisto("sumPdfWeights","all",double(ev.npdfs), ev.pdfWeights[ev.npdfs]);
	       ev.npdfs++;
	     }
	     
	     //fill alpha_s variation weights
	     if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(EvtHandle->weights().size())) {
	       for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
		 ev.alphaSWeights[ev.nalphaS] = SignGenWeight * EvtHandle->weights()[iwgt].wgt/EvtHandle->originalXWGTUP();
		 mon_.fillHisto("sumAlphasWeights","all",double(ev.nalphaS), ev.alphaSWeights[ev.nalphaS]);
		 ev.nalphaS++;
	       }
	     }
	     
	   } // pdf variation weights END
	   
	 }// EvtHandle.isValid

	 //
	 // gen particles
	 //
	 reco::GenParticleCollection gen;
	 fwlite::Handle< reco::GenParticleCollection > genHandle;
	 genHandle.getByLabel(event, "prunedGenParticles");
	 if(genHandle.isValid()){ gen = *genHandle;}
	 
	 std::vector<TLorentzVector> chLeptons;       
	 //Look for mother particle and Fill gen variables
	 for(unsigned int igen=0; igen<gen.size(); igen++){ 
	   if(!gen[igen].isHardProcess()) continue; 
	   
	     //find the ID of the first mother that has a different ID than the particle itself
	   const reco::Candidate* mom = findFirstMotherWithDifferentID(&gen[igen]);
	     
	   if (mom) {
	     int pid = gen[igen].pdgId();
	     
	     ev.mc_px[ev.nmcparticles] = gen[igen].px();
	     ev.mc_py[ev.nmcparticles] = gen[igen].py();
	     ev.mc_pz[ev.nmcparticles] = gen[igen].pz();
	     ev.mc_en[ev.nmcparticles] = gen[igen].energy();
	     ev.mc_id[ev.nmcparticles] = gen[igen].pdgId();
	     ev.mc_mom[ev.nmcparticles] = mom->pdgId();

	     // loop over genParticles to find the mom index
	     int idx=0; int idxx=0;
	     for(unsigned int ig=0; ig<gen.size(); ig++){ 
	       if(!gen[ig].isHardProcess()) continue; 

	       const reco::Candidate* imom = findFirstMotherWithDifferentID(&gen[ig]);
	       if (imom) {
		 if ( mom->p4() == gen[ig].p4() && idxx==0) {
		   idxx=idx; 
		 }
		 idx++;      
	       }
	     }

	     ev.mc_momidx[ev.nmcparticles] = idxx; 
	     ev.mc_status[ev.nmcparticles] = gen[igen].status();
	     
	     TLorentzVector p4( gen[igen].px(), gen[igen].py(), gen[igen].pz(), gen[igen].energy() );
	     if(abs(pid)==11 || abs(pid)==13 || abs(pid)==15) {
	       chLeptons.push_back(p4);
	     }
	     ev.nmcparticles++;
	   }
	   
	 }
	 
	 //
	 // gen jets
	 //
	 ev.nmcjparticles = 0;  
	 
	 reco::GenJetCollection genJets;
	 fwlite::Handle< reco::GenJetCollection > genJetsHandle;
	 genJetsHandle.getByLabel(event, "slimmedGenJets");
	 if(genJetsHandle.isValid()){ genJets = *genJetsHandle;}

	 std::vector<TLorentzVector> jets;
	 for(size_t j=0; j<genJets.size(); j++) {

	   const reco::GenJet genJet = genJets[j];

	   TLorentzVector p4( genJet.px(), genJet.py(), genJet.pz(), genJet.energy() );
	   if(p4.Pt()<10 || fabs(p4.Eta())>2.5) continue;
	   
	   bool matchesLepton(false);
	   for(size_t i=0; i<chLeptons.size(); i++) {
	     float dR=p4.DeltaR(chLeptons[i]);
	     if(dR>0.4) continue;
	     matchesLepton=true;
	     break;
	   }
	   if(matchesLepton) continue;
	   
	   jets.push_back(p4);
	   ev.mcj_px[ev.nmcjparticles]=genJet.px();
	   ev.mcj_py[ev.nmcjparticles]=genJet.py();
	   ev.mcj_pz[ev.nmcjparticles]=genJet.pz();
	   ev.mcj_en[ev.nmcjparticles]=genJet.energy();
	   ev.mcj_status[ev.nmcjparticles]=0; // special for genjet
	   ev.mcj_id[ev.nmcjparticles]=1; // special for genjet
	   ev.mcj_mom[ev.nmcjparticles] = 0;
	   ev.nmcjparticles++;
	 }

       } // end MC

       //
       // Trigger
       //
       fwlite::Handle<edm::TriggerResults> triggerBits; 
       fwlite::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects; 
       fwlite::Handle<pat::PackedTriggerPrescales> triggerPrescales; 
       triggerBits.getByLabel(event,"TriggerResults::HLT"); 
       triggerObjects.getByLabel(event,"selectedPatTrigger"); 
       triggerPrescales.getByLabel(event,"patTrigger"); 

       bool verbose_=false;
       if (verbose_) {

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
       tr = event.triggerResultsByName("HLT");
       if(!tr.isValid()  )return false;
       
       float triggerPrescale(1.0),triggerThreshold(0), triggerThresholdHigh(99999);
       bool mumuTrigger(true); bool muTrigger(true);	bool eeTrigger(true); bool eTrigger(true); bool emuTrigger(true);
       bool highPTeeTrigger(true);

       int metFilterValue = 0;
       
       bool filterbadPFMuon = true;
       bool filterbadChCandidate = true;
       bool filterbadMuonHIP = true;
       bool filterduplicateMuonHIP = true;
       std::unique_ptr<std::vector<reco::Muon*>> outbadMuon(new std::vector<reco::Muon*>());
       std::unique_ptr<std::vector<reco::Muon*>> outduplicateMuon(new std::vector<reco::Muon*>());
	 
       mumuTrigger        = utils::passTriggerPatterns(tr, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*" , "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
       muTrigger          = utils::passTriggerPatterns(tr, "HLT_IsoMu22_v*","HLT_IsoTkMu22_v*", "HLT_IsoMu24_v*", "HLT_IsoTkMu24_v*");
       eeTrigger          = utils::passTriggerPatterns(tr, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_DoubleEle33_CaloIdL_v*");
       highPTeeTrigger    = utils::passTriggerPatterns(tr, "HLT_ECALHT800_v*");
       eTrigger           = utils::passTriggerPatterns(tr, "HLT_Ele27_eta2p1_WPLoose_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*") ;
       emuTrigger         = utils::passTriggerPatterns(tr, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*" , "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*") || utils::passTriggerPatterns(tr,"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
	 
       metFilterValue = metFilter.passMetFilterInt( event );
	 
       // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
       filterbadChCandidate = metFilter.passBadChargedCandidateFilter(event); if (!filterbadChCandidate) {  metFilterValue=9; }
       filterbadPFMuon = metFilter.passBadPFMuonFilter(event); if (!filterbadPFMuon) { metFilterValue=8; }
       filterbadMuonHIP = metFilter.BadGlobalMuonTaggerFilter(event,outbadMuon,false); if (!filterbadMuonHIP) { metFilterValue=10; }
       filterduplicateMuonHIP = metFilter.BadGlobalMuonTaggerFilter(event,outduplicateMuon,true); if (!filterduplicateMuonHIP) { metFilterValue=11; }
       
       mon_.fillHisto("metFilter", "all", metFilterValue, 0);
       
       ev.hasTrigger  = ( mumuTrigger||muTrigger||eeTrigger||highPTeeTrigger||eTrigger||emuTrigger );
       
       ev.triggerType = ( mumuTrigger  << 0 )
	 | ( muTrigger  << 1 )
	 | ( eeTrigger << 2 )
	 | ( highPTeeTrigger << 3 )
	 | ( eTrigger << 4 )
         | ( emuTrigger << 5 ) ;

       //if(!isMC_ && !ev.hasTrigger) return; // skip the event if no trigger, only for Data
       //       if(!ev.hasTrigger) return false; // skip the event if no trigger, for both Data and MC

       //##############################################   EVENT PASSED THE TRIGGER   ######################################
       if (metFilterValue==10 || metFilterValue==11) { metFilterValue=0; }
       if( metFilterValue!=0 ) continue;	 //Note this must also be applied on MC

       // Apply Bad Charged Hadron and Bad Muon Filters from MiniAOD (for Run II 2016 only )
       //	  if (!filterbadPFMuon || !filterbadChCandidate) continue;
       //##############################################   EVENT PASSED MET FILTER   #######################################
       
       
       //load all the objects we will need to access
       reco::VertexCollection vtx;
       fwlite::Handle< reco::VertexCollection > vtxHandle;
       vtxHandle.getByLabel(event, "offlineSlimmedPrimaryVertices");
       if(vtxHandle.isValid()){ vtx = *vtxHandle;}
       
       double rho = 0;
       fwlite::Handle< double > rhoHandle;
       rhoHandle.getByLabel(event, "fixedGridRhoFastjetAll");
       if(rhoHandle.isValid()){ rho = *rhoHandle;}


       if (vtx.empty()) return false; // skip the event if no PV found
       const reco::Vertex &PV = vtx.front();

       ev.vtx_x = PV.x();
       ev.vtx_y = PV.y();
       ev.vtx_z = PV.z();

       ev.nvtx = 0;
       //select good vertices
       for(unsigned int i = 0; i < vtx.size(); i++) {
	 if(vtx.at(i).isValid() && !vtx.at(i).isFake()) ev.nvtx++;
       }
       if(ev.nvtx == 0) return false;

       pat::MuonCollection muons;
       fwlite::Handle< pat::MuonCollection > muonsHandle;
       muonsHandle.getByLabel(event, "slimmedMuons");
       if(muonsHandle.isValid()){ muons = *muonsHandle;}
       
       ev.mn=0;
       //       for (std::vector<pat::Muon >::const_iterator mu = muons.begin(); mu!=muons.end(); mu++) {
       for(pat::Muon &mu : muons) {
	 if(mu.pt() < 3) continue;
	 ev.mn_px[ev.mn] = mu.px();
	 ev.mn_py[ev.mn] = mu.py();
	 ev.mn_pz[ev.mn] = mu.pz();
	 ev.mn_en[ev.mn] = mu.energy();
	 ev.mn_id[ev.mn] = 13*mu.charge();

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
	 /*
	 ev.mn_nMatches[ev.mn]                   = mu.numberOfMatches();
	 ev.mn_nMatchedStations[ev.mn]           = mu.numberOfMatchedStations();
	 ev.mn_validMuonHits[ev.mn]              = mu.isGlobalMuon() ? mu.globalTrack().hitPattern().numberOfValidMuonHits() : 0.;
	 ev.mn_innerTrackChi2[ev.mn]             = mu.isTrackerMuon() ? mu.innerTrack().normalizedChi2() : 0.;
	 ev.mn_trkLayersWithMeasurement[ev.mn]   = mu.track().hitPattern().trackerLayersWithMeasurement();
	 ev.mn_pixelLayersWithMeasurement[ev.mn] = mu.isTrackerMuon() ? mu.innerTrack().hitPattern().pixelLayersWithMeasurement() : 0.;
	 */
	 ev.mn_passId[ev.mn]  = patUtils::passId(mu, vtx[0], patUtils::llvvMuonId::Tight, patUtils::CutVersion::ICHEP16Cut);
	 ev.mn_passIdLoose[ev.mn] = patUtils::passId(mu, vtx[0], patUtils::llvvMuonId::Loose, patUtils::CutVersion::ICHEP16Cut);
	 ev.mn_passSoftMuon[ev.mn] = patUtils::passId(mu, vtx[0], patUtils::llvvMuonId::Soft, patUtils::CutVersion::ICHEP16Cut);
	 ev.mn_passIso[ev.mn] = patUtils::passIso(mu, patUtils::llvvMuonIso::Tight, patUtils::CutVersion::ICHEP16Cut);
	 

	 ev.mn_type[ev.mn]   = (mu.isMuon() << 0)
	   | (mu.isGlobalMuon() << 1)
	   | (mu.isTrackerMuon() << 2)
	   | (mu.isStandAloneMuon()<< 3)
	   | (mu.isCaloMuon() << 4)
	   | (mu.isPFMuon() << 5)
	   | (mu.isRPCMuon()<< 6);
	 ev.mn++;
       }


       pat::ElectronCollection electrons;
       fwlite::Handle< pat::ElectronCollection > electronsHandle;
       electronsHandle.getByLabel(event, "slimmedElectrons");
       if(electronsHandle.isValid()){ electrons = *electronsHandle;}

       ev.en=0;
 
       for (pat::Electron &el : electrons) {
       //       for( View<pat::ElectronCollection>::const_iterator el = electrons.begin(); el != electrons.end(); el++ ) {
	 float pt_ = el.pt();
	 if (pt_ < 5) continue;

	 // Kinematics
	 ev.en_px[ev.en] = el.px();
	 ev.en_py[ev.en] = el.py();
	 ev.en_pz[ev.en] = el.pz();
	 ev.en_en[ev.en] = el.energy();
	 ev.en_id[ev.en] = 11*el.charge();

	 /*
	 //Isolation
	 GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
	 ev.en_pileupIso[ev.en] = pfIso.sumPUPt;
	 ev.en_chargedIso[ev.en] = pfIso.sumChargedHadronPt;
	 ev.en_photonIso[ev.en] = pfIso.sumPhotonEt;
	 ev.en_neutralHadIso[ev.en] = pfIso.sumNeutralHadronEt;
	 */

	 ev.en_passId[ev.en] = patUtils::passId(el, vtx[0], patUtils::llvvElecId::Tight, patUtils::CutVersion::ICHEP16Cut);
	 ev.en_passIdLoose[ev.en] = patUtils::passId(el, vtx[0], patUtils::llvvElecId::Loose, patUtils::CutVersion::ICHEP16Cut);
	 ev.en_passIso[ev.en] = patUtils::passIso(el, patUtils::llvvElecIso::Tight, patUtils::CutVersion::ICHEP16Cut, rho) ;


	 ev.en++;
       }

       
       //
       // jet selection (ak4PFJetsCHS)
       //
       pat::JetCollection jets;
       fwlite::Handle< pat::JetCollection > jetsHandle;
       jetsHandle.getByLabel(event, "slimmedJets");
       if(jetsHandle.isValid()){ jets = *jetsHandle;}
       
       ev.jet=0;
       //       PFJetIDSelectionFunctor looseJetIdSelector(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::LOOSE);
       //       PFJetIDSelectionFunctor tightJetIdSelector(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::TIGHT);
       //       pat::strbitset hasLooseId =looseJetIdSelector.getBitTemplate();   
       //       bool passPFloose = patUtils::passPFJetID("Loose", jet); //looseJetIdSelector.getBitTemplate();
       //       pat::strbitset hasTightId = tightJetIdSelector.getBitTemplate();

       //for (std::vector<pat::Jet >::const_iterator j = jets.begin(); j!=jets.end(); j++) {
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
	 ev.jet_btag1[ev.jet] = j.bDiscriminator("pfJetBProbabilityBJetTags");
	 ev.jet_btag2[ev.jet] = j.bDiscriminator("pfJetProbabilityBJetTags");
	 ev.jet_btag3[ev.jet] = j.bDiscriminator("pfTrackCountingHighPurBJetTags");
	 ev.jet_btag4[ev.jet] = j.bDiscriminator("pfTrackCountingHighEffBJetTags");
	 ev.jet_btag5[ev.jet] = j.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
	 ev.jet_btag6[ev.jet] = j.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags");
	 ev.jet_btag7[ev.jet] = j.bDiscriminator("combinedSecondaryVertexBJetTags");
	 ev.jet_btag8[ev.jet] = j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
	 ev.jet_btag9[ev.jet] = j.bDiscriminator("pfCombinedSecondaryVertexSoftLeptonBJetTags");
	 ev.jet_btag10[ev.jet] = j.bDiscriminator("pfCombinedMVABJetTags");

	 ev.jet_mass[ev.jet] = j.mass(); //correctedP4(0).M();
	 ev.jet_area[ev.jet] = j.jetArea();
	 ev.jet_pu[ev.jet] = j.pileup();
	 ev.jet_puId[ev.jet] = j.userFloat("pileupJetId:fullDiscriminant");
	 ev.jet_partonFlavour[ev.jet] = j.partonFlavour();
	 ev.jet_hadronFlavour[ev.jet] = j.hadronFlavour();
	 const reco::GenJet *gJet=j.genJet();
	 if(gJet) ev.jet_genpt[ev.jet] = gJet->pt();
	 else     ev.jet_genpt[ev.jet] = 0;

	 ev.jet++;
       }

       //
       // jet selection (AK8Jets)
       //
       pat::JetCollection fatjets; 
       fwlite::Handle< pat::JetCollection > jetsAK8Handle; 
       jetsAK8Handle.getByLabel(event, "slimmedJetsAK8"); 
       if(jetsAK8Handle.isValid()){ fatjets = *jetsAK8Handle;}   

       ev.fjet=0;
       
       for (const pat::Jet &j : fatjets) {
	 ev.fjet_px[ev.fjet] = j.correctedP4(0).px();
	 ev.fjet_py[ev.fjet] = j.correctedP4(0).py();
	 ev.fjet_pz[ev.fjet] = j.correctedP4(0).pz();
	 ev.fjet_en[ev.fjet] = j.correctedP4(0).energy();

	 ev.fjet_btag0[ev.jet] = j.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
	 
	 const reco::GenJet *gJet=j.genJet();
	 if(gJet) ev.fjet_genpt[ev.fjet] = gJet->pt();
	 else     ev.fjet_genpt[ev.fjet] = 0;

	 ev.fjet_prunedM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSPrunedMass");
	 ev.fjet_softdropM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSSoftDropMass");
	 //	 ev.fjet_filteredM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
	 ev.fjet_tau1[ev.fjet] =  (float) j.userFloat("NjettinessAK8:tau1");
	 ev.fjet_tau2[ev.fjet] =  (float) j.userFloat("NjettinessAK8:tau2");
	 ev.fjet_tau3[ev.fjet] =  (float) j.userFloat("NjettinessAK8:tau3");

	 ev.fjet++;
       }

       
       /*
       pat::PhotonCollection photons;
       fwlite::Handle< pat::PhotonCollection > photonsHandle;
       photonsHandle.getByLabel(event, "slimmedPhotons");
       if(photonsHandle.isValid()){ photons = *photonsHandle;}
       */
       
       pat::METCollection mets;
       fwlite::Handle< pat::METCollection > metsHandle;
       if(!isMC)metsHandle.getByLabel(event, "slimmedMETsMuEGClean");
       if(isMC)metsHandle.getByLabel(event, "slimmedMETs");
       if(metsHandle.isValid()){ mets = *metsHandle;}
       //       pat::MET met = mets[0];

       const pat::MET &met = mets.front();

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

       summaryHandler_.fillTree();

     }
     /*
       pat::METCollection puppimets;
       fwlite::Handle< pat::METCollection > puppimetsHandle;
       puppimetsHandle.getByLabel(event, "slimmedMETsPuppi");
       if(puppimetsHandle.isValid()){ puppimets = *puppimetsHandle;}
       LorentzVector puppimet = puppimets[0].p4();
       
       fwlite::Handle<pat::PackedCandidateCollection> PFparticles;
       PFparticles.getByLabel(event, "packedPFCandidates");
     */
       
     printf("\n");
     delete file;

  }
  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //

  //scale all events by 1/N to avoid the initial loop to stupidly count the events
  //mon.Scale(1.0/totalNumEvent);

  TString terminationCmd = "";
  //save control plots to file
  printf("Results save in local directory and moved to %s\n", outUrl.Data());

  //save all to the file
  terminationCmd += TString("mv test.root ") + outUrl + ";";
 
  if(!isMC && debugText!=""){
     TString outTxtUrl= outUrl + ".txt";
     terminationCmd += TString("mv out.txt ") + outTxtUrl + ";";
     FILE* outTxtFile = fopen("out.txt", "w");
     fprintf(outTxtFile, "%s", debugText.c_str());
     printf("TextFile URL = %s\n",outTxtUrl.Data());
     if(outTxtFile)fclose(outTxtFile);
  }

  //Now that everything is done, dump the list of lumiBlock that we processed in this job
  if(!isMC){
     terminationCmd += TString("mv out.json ") + ((outUrl.ReplaceAll(".root",""))+".json") + ";";
     goodLumiFilter.FindLumiInFiles(urls);
     goodLumiFilter.DumpToJson("out.json");
  }

  system(terminationCmd.Data());
       
}

