#define YEAR_2017
#include <iostream>
#include <map>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakePyBind11ParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/MVAHandler_vbfh.h"
#include "UserCode/bsmhiggs_fwk/interface/TMVAReader_vbfh.h"
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/bsmhiggs_fwk/interface/PDFInfo.h"
#include "UserCode/bsmhiggs_fwk/interface/RoccoR.h"
#include "UserCode/bsmhiggs_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/bsmhiggs_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/bsmhiggs_fwk/interface/METUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/XYMETCorrection.h"
#include "UserCode/bsmhiggs_fwk/interface/EventCategory.h"
#include "UserCode/bsmhiggs_fwk/interface/statWgt.h"

#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "correction.h"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TROOT.h>
#include <TChain.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TPaveStats.h>

#include <limits>
#include <utility>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <cassert>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <random>

// Initialize struct for findBoostandRotation function
struct triplet {
  std::vector<double> b;
  TVector3 k;
  double theta;
};

// Define structure for storage of jet along with its corresponding btag value and matching category: "1" for b quark, "2" for outgoing quark, "-1" for none, "0" if no matching info is available (background)
struct jets_struct {
  TLorentzVector momentum;
  double btag;
  int category;
};

// Initialize functions 
double                   getDeltaR                (TLorentzVector vec_1,                          TLorentzVector vec_2                                                                         );
double                   getDeltaRyy              (TLorentzVector vec_1,                          TLorentzVector vec_2                                                                         );
bool                     match_jets               (TLorentzVector det_vec,                        std::vector<TLorentzVector> gen_vec,               double threshold                          );
std::pair<double, int>   getDeltaRmin             (TLorentzVector det_vec,                        std::vector<TLorentzVector> gen_vec                                                          );
triplet                  findBoostandRotation     (TLorentzVector p_mom                                                                                                                        );
int                      rand_int                 (int min,                                       int max                                                                                      );
double                   getDeltaMmin             (const std::vector<TLorentzVector>& vec_bjets,  double& unpairedMass                                                                         );
double                   getMassOfbbj             (const std::vector<jets_struct>& vec_nonqjets,  int iEvt                                                                                     );
double                   getDeltaMminNew          (const std::vector<TLorentzVector>& vec_bjets                                                                                                );
int                      EvtCatIndex              (const std::vector<jets_struct>& vec_bjets,     const std::vector<jets_struct>& vec_nonqjets,      const std::vector<jets_struct>& vec_qjets );

using namespace std;
namespace fs = boost::filesystem;

// Main function
int main(int argc, char* argv[])
{

  //#######################################################################################################################################//
  //---------------------------------------------------------GLOBAL INITIALIZATION---------------------------------------------------------//
  //#######################################################################################################################################//

  gSystem->Load("libTree.so");
  
  // Check arguments
  if(argc<2)
    {
      std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
      exit(0);
    }

  // Load framework libraries
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();
  MVAHandler myMVAHandler_;
  myMVAHandler_.initTree();

  // Configure the process
  const edm::ParameterSet &runProcess = edm::cmspybind11::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  
  // Will produce the input root trees to BDT training (optional)
  bool runMVA = runProcess.getParameter<bool>("runMVA");
  
  // For the Trigger Efficiency plots
  bool runTriggerStudies = runProcess.getParameter<bool>("runTriggerStudies");

  // For the detector and generator object matching
  bool runMatchingStudies = runProcess.getParameter<bool>("runMatchingStudies");
  
  // For the QCD control region
  bool runQCD = runProcess.getParameter<bool>("runQCD");

  bool   isMC    = runProcess.getParameter<bool>("isMC");
  double xsec    = runProcess.getParameter<double>("xsec");
  double nevts   = runProcess.getParameter<double>("nevts");

  TString proc   = runProcess.getParameter<std::string>("proc");
  TString dtag   = runProcess.getParameter<std::string>("tag");
  TString suffix = runProcess.getParameter<std::string>("suffix");

  bool is2018data  = (!isMC && dtag.Contains("2018"));
  bool is2018MC    = (isMC && dtag.Contains("2018"));
  
  bool verbose = runProcess.getParameter<bool>("verbose");

  TString url = runProcess.getParameter<std::string>("input");
  TString outFileUrl(dtag); 
  gSystem->BaseName( url );

  bool isQCD     = isMC && dtag.Contains("QCD");
  bool isMC_Wh   = isMC && (string(url.Data()).find("Wh_amass")   != string::npos); 
  bool isMC_Zh   = isMC && (string(url.Data()).find("Zh_amass")   != string::npos);
  bool isMC_VBFh = isMC && (string(url.Data()).find("VBFh_amass") != string::npos);
  bool isSignal  = isMC_VBFh;

  TString outdir = runProcess.getParameter<std::string>("outdir");
  TString outUrl( outdir );
  TString outTxtUrl = outUrl + "/TXT";
  gSystem->Exec("mkdir -p " + outUrl);
  gSystem->Exec("mkdir -p " + outTxtUrl);

  TString outTxtUrl_final= outTxtUrl + "/" + outFileUrl + "_FinalList.txt";
  FILE* outTxtFile_final = NULL;
  outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
  printf("TextFile URL = %s\n",outTxtUrl_final.Data());
  fprintf(outTxtFile_final,"run lumi event\n");

  //#######################################################################################################################################//
  //--------------------------------------------------------INITIALIZE HISTOGRAMS----------------------------------------------------------//
  //#######################################################################################################################################//

  SmartSelectionMonitor mon;

  // Event flow
  TH1F *h_event_flow = (TH1F*) mon.addHistogram ( new TH1F ("event_flow", "Cut Flow;Step;N_{Events}", 9, 0, 9) );
  h_event_flow->GetXaxis()->SetBinLabel(1,"0: Raw Events");
  h_event_flow->GetXaxis()->SetBinLabel(2,"1: N_leptons=0");
  h_event_flow->GetXaxis()->SetBinLabel(3,"2: VBF Triggers");
  h_event_flow->GetXaxis()->SetBinLabel(4,"3: N_jets>=5");
  h_event_flow->GetXaxis()->SetBinLabel(5,"4: Jet pT cuts");
  h_event_flow->GetXaxis()->SetBinLabel(6,"5: N_bjets>=2 & N_qjets>=2 & N_nonqjets >=3");
  h_event_flow->GetXaxis()->SetBinLabel(7,"6: SR/CR");
  h_event_flow->GetXaxis()->SetBinLabel(8,"7: Dh>2.5 and Mqq>250");
  h_event_flow->GetXaxis()->SetBinLabel(9,"8: Dphi<1.6");

  // Event Categorization
  TH1F *h_event_cat = (TH1F*) mon.addHistogram ( new TH1F ("event_cat", "Event Categorization;Category;N_{Events}", 12, 0, 12) );
  h_event_cat->GetXaxis()->SetBinLabel(1,  "(2b, 3j)");
  h_event_cat->GetXaxis()->SetBinLabel(2,  "(2b, 4j)");
  h_event_cat->GetXaxis()->SetBinLabel(3,  "(2b, >=5j)");
  h_event_cat->GetXaxis()->SetBinLabel(4,  "(3b, 3j)");
  h_event_cat->GetXaxis()->SetBinLabel(5,  "(3b, 4j)");
  h_event_cat->GetXaxis()->SetBinLabel(6,  "(3b, >=5j)");
  h_event_cat->GetXaxis()->SetBinLabel(7,  "(4b, 3j)");
  h_event_cat->GetXaxis()->SetBinLabel(8,  "(4b, 4j)");
  h_event_cat->GetXaxis()->SetBinLabel(9,  "(4b, >=5j)");
  h_event_cat->GetXaxis()->SetBinLabel(10, "(5b, 3j)");
  h_event_cat->GetXaxis()->SetBinLabel(11, "(5b, 4j)");
  h_event_cat->GetXaxis()->SetBinLabel(12, "(5b, >=5j)");

  // Particles/Objects kinematics (single)
  mon.addHistogram ( new TH1F ("pt", "Transverse Momentum;p_{T} [GeV];N_{Events}", 500, 0, 2000) );
  mon.addHistogram ( new TH1F ("eta", "Pseudorapidity;#eta [-];N_{Events}", 150, -6, 6) );
  mon.addHistogram ( new TH1F ("psi", "Rapidity;y [-];N_{Events}", 150, -6, 6) );
  mon.addHistogram ( new TH1F ("phi", "Azimuthal Angle;#phi [-];N_{Events}", 150, -TMath::Pi(), TMath::Pi()) );
  mon.addHistogram ( new TH1F ("M", "Mass;M [GeV];N_{Events}", 500, 0, 2000) );
  mon.addHistogram ( new TH1F ("E", "Energy;E [GeV];N_{Events}", 500, 0, 2000) );

  mon.addHistogram ( new TH1F ("theta", "Polar Angle;#theta [-];N_{Events}", 150, 0, TMath::Pi()) );
  mon.addHistogram ( new TH1F ("costheta", "Polar Angle;cos #theta [-];N_{Events}", 150, -1, 1) );

  // -//- (multi)
  mon.addHistogram ( new TH1F ("deltaR", "Angular Distance;#Delta R [-];N_{Events}", 150, 0, 16) );
  mon.addHistogram ( new TH1F ("deltaRyy", "Angular Distance;#Delta R [-];N_{Events}", 150, 0, 16) );
  mon.addHistogram ( new TH1F ("deltaEta", "Pseudorapidity Difference;#left|#Delta#eta#right| [-];N_{Events}", 150, 0, 10) );
  mon.addHistogram ( new TH1F ("deltaPhi", "Azimuthal Angle Difference;#left|#Delta#phi#right| [-];N_{Events}", 150, 0, 10) );
  mon.addHistogram ( new TH1F ("prodEta", "Pseudorapidity Product;#eta_{1} #times #eta_{2} [-];N_{Events}", 300, -25, 25) );
  mon.addHistogram ( new TH1F ("Mqq", "Mass;M [GeV];N_{Events}", 500, 0, 2000) );

  // Multiplicity
  mon.addHistogram ( new TH1F ("multi", ";N [-] ;N_{Events}", 18, 0, 18) );

  // B-tag discriminator
  mon.addHistogram ( new TH1F ( "btag", "b-tag discriminator;score [-];N_{Events}", 150, 0, 1) );
 
  // HT, Hz
  mon.addHistogram ( new TH1F ("HT", ";H_{T} [GeV];N_{Events}", 200, 0, 2000) );
  mon.addHistogram ( new TH1F ("HT_studies", ";H_{T} [GeV];N_{Events}", 500, 0, 3000) );
  mon.addHistogram ( new TH1F ("Hz", "Longitudinal Momentum;#sum_{i}{p_{z,{b_i}}} + #sum_{j}{p_{z,{q_j}}} [GeV];N_{Events}", 600, -6000, 6000) );
  mon.addHistogram ( new TH1F ("HT_vec", ";#frac{#sum_{i}{#vec{p_i}}_T}{#sum_{i}{p_{i,T}}} [-];N_{Events}", 150, 0, 2) );

  // Flavour
  mon.addHistogram( new TH1F( "flavour", "flavour", 6, 0, 6) );

  // Matching plots
  TH1F *h_match_bjets = (TH1F*) mon.addHistogram ( new TH1F ("match_bjets", "b-tagged jets;true category;N_{Events}", 3, 0, 3) );
  h_match_bjets->GetXaxis()->SetBinLabel(1, "b quark");
  h_match_bjets->GetXaxis()->SetBinLabel(2, "outgoing quark");
  h_match_bjets->GetXaxis()->SetBinLabel(3, "none");

  TH1F *h_match_qjets = (TH1F*) mon.addHistogram ( new TH1F ("match_qjets", "q-tagged jets;true category;N_{Events}", 3, 0, 3) );
  h_match_qjets->GetXaxis()->SetBinLabel(1, "b quark");
  h_match_qjets->GetXaxis()->SetBinLabel(2, "outgoing quark");
  h_match_qjets->GetXaxis()->SetBinLabel(3, "none");

  TH1F *h_match_untaggedjets = (TH1F*) mon.addHistogram ( new TH1F ("match_untaggedjets", "untagged jets;true category;N_{Events}", 3, 0, 3) );
  h_match_untaggedjets->GetXaxis()->SetBinLabel(1, "b quark");
  h_match_untaggedjets->GetXaxis()->SetBinLabel(2, "outgoing quark");
  h_match_untaggedjets->GetXaxis()->SetBinLabel(3, "none");

  // 2-d histograms 
  TH2F *h_match_eta = (TH2F*) mon.addHistogram ( new TH2F (TString("match_eta"), "Gen-Reco Matching vs Eta;#eta [-];Generator Match;N_{Events}", 150, -6, 6, 3, 0, 3) );
  h_match_eta->GetYaxis()->SetBinLabel(1, "b quark");
  h_match_eta->GetYaxis()->SetBinLabel(2, "outgoing quark");
  h_match_eta->GetYaxis()->SetBinLabel(3, "none");

  TH2F *h_match_pt = (TH2F*) mon.addHistogram ( new TH2F (TString("match_pt"), "Gen-Reco Matching vs p_{T};p_{T} [GeV];Generator Match;N_{Events}", 150, 0, 150, 3, 0, 3) );
  h_match_pt->GetYaxis()->SetBinLabel(1, "b quark");
  h_match_pt->GetYaxis()->SetBinLabel(2, "outgoing quark");
  h_match_pt->GetYaxis()->SetBinLabel(3, "none");

  mon.addHistogram ( new TH2F (TString("mindeltaR_vs_pt"),"Angular Distance vs p_{T};p_{T} [GeV];min#Delta R [-];N_{Events}", 150, 0, 500, 150, 0, 7) );
  mon.addHistogram ( new TH2F (TString("mindeltaR_vs_eta"), "Angular Distance vs Eta;#eta [-];min#Delta R [-];N_{Events}", 150, -6, 6, 150, 0, 7) );

  // 2-d histograms for jets and b-jets pt_gen vs pt_reco
  TH2F *h_pt_reco_gen = (TH2F*) mon.addHistogram ( new TH2F (TString("pt_reco_gen"), "Gen-Reco p_{T};p_{T,reco} [GeV];3p_{T,gen} [GeV];N_{Events}", 500, 0, 2000, 500, 0, 2000) );
  
  // BDT score  
  mon.addHistogram ( new TH1F ("BDT", "BDT;BDT Score [-];N_{Events}", 200, -1, 1) );

  // PileUp distribution
  mon.addHistogram ( new TH1F ("pileup", "Vertex Multiplicity;Vertices [-];N_{Events}", 100, 0, 100) );

  //#######################################################################################################################################//
  //----------------------------------------------------INITIALIZE KINEMATIC VARIABLES-----------------------------------------------------//
  //#######################################################################################################################################//

  // Category tag (Signal or Control Region)
  TString tag;
  
  // Higgs kinematics
  float higgs_m;
  float higgs_pt;
  float higgs_eta;
  float costheta0;

  //Higgs - VBF jets kinematics
  float dR_higgs_q1;
  float dR_higgs_q2;
  float phi_qq_higgs;

  // b pair kinematics
  float avDeltaR;
  float minDeltaM;
  float MassOfbbj;
  float minDeltaPhi;

  // Outgoing quark pair kinematics
  float DeltaEta;
  float ProductEta;
  float qq_m;
  float qq_dR;
  float qq_dphi;
  float alpha_qq;

  // Jet kinematics
  float N_jet;
  float N_jet_eta_cut;
  float jet1_pt;
  float jet2_pt;
  float jet3_pt;
  float jet4_pt;
  float jet5_pt;
  float jet6_pt;
  float jet1_eta;
  float jet2_eta;
  float jet3_eta;
  float jet4_eta;
  float jet5_eta;
  float jet6_eta;
  float H_T;
  float H_z;
  float H_Tvec;

  // bjet kinematics
  float N_bjet;
  float bjet1_pt;
  float bjet2_pt;
  float bjet3_pt;
  float bjet4_pt;
  float bjet1_eta;
  float bjet2_eta;
  float bjet3_eta;
  float bjet4_eta;
  float bjet1_btag;
  float bjet2_btag;
  float bjet3_btag;
  float bjet4_btag;

  // qjet kinematics
  float N_qjet;
  float qjet1_pt;
  float qjet2_pt;
  float qjet1_eta;
  float qjet2_eta;

  // Untagged jets
  float N_untaggedjet;
  float pt_rest;
  float E_rest;

  // MET
  float MET_pt;
  float MET_phi;
  
  //#######################################################################################################################################//

  //--------------------------------------------------------------Get Tree info------------------------------------------------------------//
  int evStart     = runProcess.getParameter<int>("evStart");
  int evEnd       = runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  // Open the file and get events tree
  DataEvtSummaryHandler summaryHandler_;

  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0)
    {
      return -1;
      printf("file is 0");
    }

  if(file->IsZombie()) return -1;
  
  if(!summaryHandler_.attachToTree( (TTree*)file->Get(dirname) ) )
    {
      file->Close();
      return -1;
    }

  // Check run range to compute scale factor (if not all entries are used)
  const Int_t totalEntries= summaryHandler_.getEntries();
  float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
  if(evEnd<0 || evEnd>summaryHandler_.getEntries() ) evEnd=totalEntries;
  
  if(evStart > evEnd )
    {
      file->Close();
      return -1;
    }

  //#######################################################################################################################################//
  //--------------------------------------------------------------TMVA READER--------------------------------------------------------------//
  //#######################################################################################################################################//

  std::string chpath = "VBFHaa4bMVA/";
  TMVAReader VBFhAnalysisReader;
  VBFhAnalysisReader.InitTMVAReader();
  std::string VBFhAnalysis_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"MVAnalysis_BDT.weights.xml";
  VBFhAnalysisReader.SetupMVAReader("VBFhAnalysisClass", VBFhAnalysis_xml_path);

  //#######################################################################################################################################//
  //----------------------------------------------------------PREPARE FOR THE LOOP---------------------------------------------------------//
  //#######################################################################################################################################//

  // Define Event Category Boolean
  bool isSR(false);

  // Initialize weight
  double weight(1.0);
  
  // QCD k factor (we need to calculate this in our own)
  float k_factor = 1;
  
  //##############################################//
  // X sec weighting per luminosity unit [pb^-1]
  //##############################################//

  double weight_xsec = isMC ? xsec / nevts : 1.;
    
  //##############################################//
  // Pile up weighting
  //##############################################//

  // Read pile up distribution (data)
  TString PU = runProcess.getParameter<std::string>("PU_Central");
  gSystem->ExpandPathName(PU);
  cout << "Loading PU weights: " << PU << endl;
  TFile *PU_File = TFile::Open(PU);
  TH1F* PU_intended = (TH1F *) PU_File->Get("pileup"); // Data pileup distribution
  
  // Initialize histograms
  TH1F* PU_generated=NULL; // MC pileup distribution 
  TH1F* PU_weight = new TH1F("hPUweight", "", 100, 0, 100);

  if (isMC) {
    // Read and normalize the MC pileup distribution
    TH1F* PU_generated = (TH1F*)file->Get("mainNtuplizer/pileuptrue");
    PU_intended->Scale(1.0 / PU_intended->Integral()); // Normalize data distribution
    PU_generated->Scale(1.0 / PU_generated->Integral()); // Normalize MC distribution

    // Calculate weights as the ratio of data to MC distributions
    PU_intended->Divide(PU_generated);

    // Fill PU_weight histogram with calculated weights
    for (int ibin = 0; ibin < 102; ++ibin) {
      float weight = PU_intended->GetBinContent(ibin);
      PU_weight->SetBinContent(ibin, weight);
      if (verbose) printf("pu = %3d has weight = %7.3f \n", ibin, weight);
    }
  }

  //##############################################//
  // Prepare for b-tagging
  //##############################################//
  
  // DeepJet working points
                              //efficiency:
  double looseWP  = 0.0490;   //65.1%
  double mediumWP = 0.2783;   //80.7%
  double tightWP  = 0.7100;   //91.5%

  // Read b tagging efficiencies
  TH3F* btagEffLoose_num    = new TH3F(), *btagEffMedium_num    = new TH3F(), *btagEffTight_num    = new TH3F(), *btagEff_denom    = new TH3F();
  TH3F* btagEffLoose        = new TH3F(), *btagEffMedium        = new TH3F(), *btagEffTight        = new TH3F();

  // pt bin edges: {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000}
  // eta bin edges: -3 to 3 with 60 total bins
  // flavour: 0 for b, 1 for c and 2 for u,d,s,g

  if (isMC) {
    // Retrieve numinator and denominator histograms for all wp's
    btagEffLoose_num  = (TH3F *)file->Get("mainNtuplizer/LooseBTaggingEff_Num");
    btagEffMedium_num = (TH3F *)file->Get("mainNtuplizer/MediumBTaggingEff_Num");
    btagEffTight_num  = (TH3F *)file->Get("mainNtuplizer/TightBTaggingEff_Num");
    btagEff_denom     = (TH3F *)file->Get("mainNtuplizer/BTaggingEff_Denom");

    // Divide for efficiency
    btagEffLoose = (TH3F*)btagEffLoose_num->Clone("btagEffLoose");
    btagEffLoose->Divide(btagEff_denom);
    btagEffMedium = (TH3F*)btagEffMedium_num->Clone("btagEffMedium");
    btagEffMedium->Divide(btagEff_denom);
    btagEffTight = (TH3F*)btagEffTight_num->Clone("btagEffTight");
    btagEffTight->Divide(btagEff_denom);
  }

  if (verbose){
    cout << endl;
    cout << "The efficiency for (pT, eta, flavour) = (160, 0.8, b) and Loose WP is " << getSFfrom3DHist(160., 0.8, 0., btagEffLoose) << endl;
    cout << endl;
  }

  // Read b-tagging scale factors
  fs::path sf_file_name = std::string(std::getenv("CMSSW_BASE"))+
    "/src/UserCode/bsmhiggs_fwk/data/weights/btagging.json.gz";
  auto cset = correction::CorrectionSet::from_file(sf_file_name.string());

  // For b and c flavour the "comb" SF's are recommended.
  // For light flavours the "incl" SF's are recommended.

  auto hf_sf = cset->at("deepJet_comb");
  auto lf_sf = cset->at("deepJet_incl");

  // Get the scale factor for parameters: ('systematic', 'working_point', 'flavor', 'abseta', 'pt')
  if (verbose){
    cout << endl;
    cout << "The SF for ('systematic', 'working_point', 'flavor', 'abseta', 'pt') = (central, M, b, 0.8, 160) is " << hf_sf->evaluate({"central", "M", 5, 0.8, 160.}) << endl;
    cout << endl;
  }
  
  
  //##############################################//
  // PileUpID Loose WP vs (pt, eta)
  //##############################################//
  const int PUptNBins = 4;
  const int PUetaNBins = 4;

  double PUptBins[PUptNBins + 1] = {10, 20, 30, 40, 50};  // pT bins
  double PUetaBins[PUetaNBins + 1] = {0.0, 2.5, 2.75, 3.0, 5.0};  // abs(eta) bins

  // Create TH2F histogram
  TH2F *PULooseWP = new TH2F("PULooseWP", "PU Loose WPs",
			   PUptNBins, PUptBins,  // X-axis (pT)
			   PUetaNBins, PUetaBins);   // Y-axis (abs(eta))

  // Table values corresponding to [pT bin][eta bin]
  double values[PUptNBins][PUetaNBins] = {
    {-0.95, -0.72, -0.68, -0.47},
    {-0.88, -0.55, -0.60, -0.43},
    {-0.63, -0.18, -0.43, -0.24},
    {-0.19,  0.22, -0.13, -0.03}
  };

  // Fill the TH2F histogram
  for (int i = 0; i < PUptNBins; i++) {
    for (int j = 0; j < PUetaNBins; j++) {
      PULooseWP->SetBinContent(i + 1, j + 1, values[i][j]);  // ROOT bins start from 1
    }
  }

  if (verbose){
    cout << endl;
    cout << "The PileUP Loose WP for (pT, abs(eta)) = (30, 0.8 ) is " << getSFfrom2DHist(30., 0.8, PULooseWP) << endl;
    cout << endl;
  }


  //#######################################################################################################################################//
  //--------------------------------------------------------------EVENT LOOP---------------------------------------------------------------//
  //#######################################################################################################################################//

  // Loop on all the events
  int treeStep = (evEnd-evStart)/50; 
  if(treeStep==0)treeStep=1;
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);

  for(int iev=evStart; iev<evEnd; iev++)
    {
      // Load the event content from tree
      summaryHandler_.getEntry(iev);
      DataEvtSummary_t &ev=summaryHandler_.getEvent();

      // Check for duplicates
      if(!isMC && duplicatesChecker.isDuplicate(ev.run, ev.lumi, ev.event))
	{
	  nDuplicates++;
	  cout << "nDuplicates: " << nDuplicates << endl;
	  continue;
	}

      // Initialize weight
      weight = 1.0;

      // X sec weighting
      weight *= weight_xsec;

      // Control Plot: Vertex Multiplicity BEFORE the pileUP reweighting
      if(isMC) mon.fillHisto("pileup", "before", ev.ngenTruepu, weight);

      // PileUp weighting
      float weight_PU(1.0);
      if(isMC) {
	weight_PU = getSFfrom1DHist(ev.ngenTruepu, PU_weight); //removed the +1
	if ( verbose ) printf("pu = %3d has weight = %7.3f \n", ev.ngenTruepu + 1, weight_PU);
	weight *= weight_PU;
      }
      
      // Control Plot: Vertex Multiplicity AFTER the pileUP reweighting
      if(isMC) mon.fillHisto("pileup", "after", ev.ngenTruepu, weight);

      // Initialize boolean for HEM veto
      bool afterRun319077(false);
      if(is2018data) afterRun319077 = (ev.run >= 319077);
      
      //#######################################################################################################################################//
      //--------------------------------------------------------GENERATOR LEVEL ANALYSIS-------------------------------------------------------//
      //#######################################################################################################################################//

      //--------------------------------------------------------------Step 1: Truth------------------------------------------------------------//

      // Define vector to store Higgs momentum
      std::vector<TLorentzVector> Higgs;
      
      // Define vector to store A Bosons momentum
      std::vector<TLorentzVector> A_Boson;
      std::vector<TLorentzVector> A_Boson_1;
      std::vector<TLorentzVector> A_Boson_2;

      // Define vector to store the outgoing q's momentum (and id)
      std::vector<TLorentzVector> out_quarks;
      std::vector<int> out_quarks_id;

      // Define vector to store the b quarks momentum
      std::vector<TLorentzVector> b_quarks_1;
      std::vector<TLorentzVector> b_quarks_2;
      std::vector<TLorentzVector> b_quarks;

      //----------------------------------------------------Step 2: Truth (after the acceptance cuts)------------------------------------------//

      // Define vector to store the outgoing q's momentum
      std::vector<TLorentzVector> out_quarks_acc;

      // Define vector to store the b quarks momentum
      std::vector<TLorentzVector> b_quarks_1_acc;
      std::vector<TLorentzVector> b_quarks_2_acc;

      //--------------------------------------------Step 3: Generator level (without using the mom's id)---------------------------------------//

      // Define vector in order to store quarks & gluons and bottoms in general
      std::vector<TLorentzVector> bottoms;
      std::vector<TLorentzVector> quarks_gluons;

      //----------perform this analysis only for the signal----------//
      if (isSignal)
	{
	  //----------------------------------------------------------------------------//
	  //---------------------------loop over MC particles---------------------------// 
	  //----------------------------------------------------------------------------//

	  for (int imc=0; imc<ev.nmcparticles; imc++)
	    {
	      TLorentzVector p_particle;
	      p_particle.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Printouts

	      if(iev<5 && verbose) 
		{
		  cout << " Particle no.  " << imc << "  is a  " << ev.mc_id[imc] << "  and has a mother at  " << ev.mc_momidx[imc] <<  "  ( "<<ev.mc_mom[imc] << " )" <<"  has status  " << ev.mc_status[imc] <<  "  and has a 4-vector:  p = (" << p_particle.Pt() << ", " << p_particle.Eta() << ", " << p_particle.Phi() << ", " << p_particle.M() << ") " << endl;
		  cout << endl;
		}

	      // Extract the info using the mothers etc. (TRUTH ANALYSIS)
	    
	      // Find the Higgs boson
	      if(ev.mc_id[imc]==25) Higgs.push_back(p_particle);	      

	      // Find the A bosons
	      if(ev.mc_id[imc]==36)
		{ 
		  A_Boson.push_back(p_particle);
		  // Catch A1 and A2 boson
		  if (imc == 5) A_Boson_1.push_back(p_particle);
		  else A_Boson_2.push_back(p_particle);	
		}
	  
	      // Find the outgoing quarks
	      if((abs(ev.mc_id[imc])>=1 && abs(ev.mc_id[imc])<=5) && ev.mc_mom[imc]!=36 && ev.mc_status[imc]==23)
		{
		  out_quarks.push_back(p_particle);
		  out_quarks_id.push_back(abs(ev.mc_id[imc]));

		  // Apply acceptance cut
		  if (p_particle.Pt()>20. && abs(p_particle.Eta())<=4.7) out_quarks_acc.push_back(p_particle);		  
		}

	      // Find the b quarks
	      if(abs(ev.mc_id[imc])==5 && ev.mc_mom[imc]==36 && ev.mc_status[imc]==23)
		{
		  b_quarks.push_back(p_particle);
		  if(ev.mc_momidx[imc]==5)
		    {
		      b_quarks_1.push_back(p_particle);
		      if (p_particle.Pt()>20. && abs(p_particle.Eta())<=4.7)  b_quarks_1_acc.push_back(p_particle);
		    }	      
		  if(ev.mc_momidx[imc]==6)
		    {
		      b_quarks_2.push_back(p_particle);
		      if (p_particle.Pt()>20. && abs(p_particle.Eta())<=4.7) b_quarks_2_acc.push_back(p_particle);
		    }
		}
	  
	      // Extract the info only using the type of the particle (and also applying the acceptance cuts) (GENERATOR LEVEL FINAL STATE ANALYSIS)

	      // Quarks - gluons vector
	      if(((abs(ev.mc_id[imc])>=1 && abs(ev.mc_id[imc])<=4) || (ev.mc_id[imc]==21)) && (ev.mc_status[imc]==23 || ev.mc_status[imc]==24))
		{
		  if (p_particle.Pt()>20. && abs(p_particle.Eta())<=4.7) quarks_gluons.push_back(p_particle);
		}
	
	      // Bottoms vector
	      if(abs(ev.mc_id[imc])==5 && (ev.mc_status[imc]==23 || ev.mc_status[imc]==24))
		{
		  if (p_particle.Pt()>20. && abs(p_particle.Eta())<=4.7) bottoms.push_back(p_particle);
		}	  	  
	    }      // END loop of MC particles


	  // Printouts
	  if(iev < 5 && verbose)
	    {
	      cout << "mass of outqoing quarks: " << out_quarks[0].M() << endl;
	      cout << endl;
	    }

	  //----------------------------------------------------------------------------//
	  //------------------------------Fill histogramms------------------------------// 
	  //----------------------------------------------------------------------------//

	  //-------------------------------Step 1: Truth--------------------------------//

	  // b quark multiplicity
	  mon.fillHisto("multi", "b_truth", b_quarks.size(), 1);

	  // Higgs kinematics
	  TLorentzVector pH = Higgs.back();

	  mon.fillHisto("pt", "H", pH.Pt(), 1);
	  mon.fillHisto("eta", "H", pH.Eta(), 1);
	  mon.fillHisto("psi", "H", pH.Rapidity(), 1);
	  mon.fillHisto("phi", "H", pH.Phi(), 1);
	  mon.fillHisto("M", "H", pH.M(), 1);

	  // Find boosting direction and rotation direction and angle (Higgs)
	  std::vector<double> b_H = findBoostandRotation(pH).b;
	  TVector3 k_H = findBoostandRotation(pH).k;
	  double theta_H = findBoostandRotation(pH).theta;

	  // A boson kinematics

	  for (size_t i = 0; i < A_Boson.size(); ++i)
	    {
	      TLorentzVector pA = A_Boson[i];
	      
	      mon.fillHisto("pt", "A", pA.Pt(), 1);
	      mon.fillHisto("eta", "A", pA.Eta(), 1);
	      mon.fillHisto("psi", "A", pA.Rapidity(), 1);
	      mon.fillHisto("phi", "A", pA.Phi(), 1);
	      mon.fillHisto("M", "A", pA.M(), 1);
	      mon.fillHisto("theta", "A", pA.Theta(), 1);
	      mon.fillHisto("costheta", "A", pA.CosTheta(), 1);

	      pA.Boost(-b_H[0], -b_H[1], -b_H[2]);
	      pA.Rotate(theta_H, k_H);

	      mon.fillHisto("costheta", "0A", pA.CosTheta(), 1);

	    }

	  //
	  // Outgoing quarks
	  //
	  
	  // Sort the quarks on pT 
	  std::sort(out_quarks.begin(), out_quarks.end(), [](const TLorentzVector &a, const TLorentzVector &b){
	      return a.Pt() > b.Pt();
	    });

	  mon.fillHisto("pt", "q1", out_quarks[0].Pt(), 1);
	  mon.fillHisto("eta", "q1", out_quarks[0].Eta(), 1);
	  mon.fillHisto("phi", "q1", out_quarks[0].Phi(), 1);
	  mon.fillHisto("M", "q1", out_quarks[0].M(), 1);

	  mon.fillHisto("pt", "q2", out_quarks[1].Pt(), 1);
	  mon.fillHisto("eta", "q2", out_quarks[1].Eta(), 1);
	  mon.fillHisto("phi", "q2", out_quarks[1].Phi(), 1);
	  mon.fillHisto("M", "q2", out_quarks[1].M(), 1);
	  
	  mon.fillHisto("deltaEta", "q1q2", abs(out_quarks[0].Eta()-out_quarks[1].Eta()), 1);
	  mon.fillHisto("prodEta", "q1q2", out_quarks[0].Eta()*out_quarks[1].Eta(), 1);
	  mon.fillHisto("Mqq", "q1q2", (out_quarks[0]+out_quarks[1]).M(), 1);

	  mon.fillHisto("flavour", "q", out_quarks_id[0], 1);
	  mon.fillHisto("flavour", "q", out_quarks_id[1], 1);

	  //
	  // B quarks
	  //
	  
	  // Pick a random b quark
	  int random_number = rand_int(0, 3);
	  TLorentzVector b_random = b_quarks[random_number];

	  // Boost and rotate it to the Higgs f.o.m
	  b_random.Boost(-b_H[0], -b_H[1], -b_H[2]);
	  b_random.Rotate(theta_H, k_H);
	  mon.fillHisto("costheta", "0b_random", b_random.CosTheta(), 1);

	  // Sorting the b quarks on pT
	  std::sort(b_quarks.begin(), b_quarks.end(), [](const TLorentzVector &a, const TLorentzVector &b){
	      return a.Pt() > b.Pt();
	    });

	  mon.fillHisto("pt", "b1", b_quarks[0].Pt(), 1);
	  mon.fillHisto("eta", "b1", b_quarks[0].Eta(), 1);
	  mon.fillHisto("phi", "b1", b_quarks[0].Phi(), 1);
	  mon.fillHisto("M", "b1", b_quarks[0].M(), 1);

	  mon.fillHisto("pt", "b2", b_quarks[1].Pt(), 1);
	  mon.fillHisto("eta", "b2", b_quarks[1].Eta(), 1);
	  mon.fillHisto("phi", "b2", b_quarks[1].Phi(), 1);
	  mon.fillHisto("M", "b2", b_quarks[1].M(), 1);

	  mon.fillHisto("pt", "b3", b_quarks[2].Pt(), 1);
	  mon.fillHisto("eta", "b3", b_quarks[2].Eta(), 1);
	  mon.fillHisto("phi", "b3", b_quarks[2].Phi(), 1);
	  mon.fillHisto("M", "b3", b_quarks[2].M(), 1);

	  mon.fillHisto("pt", "b4", b_quarks[3].Pt(), 1);
	  mon.fillHisto("eta", "b4", b_quarks[3].Eta(), 1);
	  mon.fillHisto("phi", "b4", b_quarks[3].Phi(), 1);
	  mon.fillHisto("M", "b4", b_quarks[3].M(), 1);

	  //
	  // bb pairs, AA pair and bbbb
	  //

	  mon.fillHisto("M", "bb1", (b_quarks_1[0]+b_quarks_1[1]).M(), 1);
	  mon.fillHisto("deltaR", "bb1", getDeltaR(b_quarks_1[0], b_quarks_1[1]), 1);
	  
	  mon.fillHisto("M", "bb2", (b_quarks_2[0]+b_quarks_2[1]).M(), 1);
	  mon.fillHisto("deltaR", "bb2", getDeltaR(b_quarks_2[0], b_quarks_2[1]), 1);
	  
	  mon.fillHisto("M", "bbbb", (b_quarks_1[0]+b_quarks_1[1]+b_quarks_2[0]+b_quarks_2[1]).M(), 1);
	  mon.fillHisto("pt", "bbbb", (b_quarks_1[0]+b_quarks_1[1]+b_quarks_2[0]+b_quarks_2[1]).Pt(), 1);

	  mon.fillHisto("deltaR", "AA", getDeltaR(A_Boson[0], A_Boson[1]), 1);
	  mon.fillHisto("deltaRyy", "AA", getDeltaRyy(A_Boson[0], A_Boson[1]), 1);

	  // Boost to A boson rest frame

	  TLorentzVector pA1 = A_Boson_1.back();
	  TLorentzVector pA2 = A_Boson_2.back();

	  // Find boosting direction and rotation direction and angle (A1, A2 Boson)
	  std::vector<double> b_A1 = findBoostandRotation(pA1).b;
	  TVector3 k_A1 = findBoostandRotation(pA1).k;
	  double theta_A1 = findBoostandRotation(pA1).theta;

	  std::vector<double> b_A2 = findBoostandRotation(pA2).b;
	  TVector3 k_A2 = findBoostandRotation(pA2).k;
	  double theta_A2 = findBoostandRotation(pA2).theta;

	  // Boost and rotate
	  for (size_t i = 0; i < b_quarks_1.size(); ++i)
	    {
	      TLorentzVector p_b = b_quarks_1[i];
	      p_b.Boost(-b_A1[0], -b_A1[1], -b_A1[2]);
	      p_b.Rotate(theta_A1, k_A1);
	      
	      mon.fillHisto("costheta", "0b", p_b.CosTheta(), 1);
	      mon.fillHisto("costheta", "0b_all", p_b.CosTheta(), 1);
	    }
	  
	  for (size_t i = 0; i < b_quarks_2.size(); ++i)
	    {
	      TLorentzVector p_b = b_quarks_2[i];
	      p_b.Boost(-b_A2[0], -b_A2[1], -b_A2[2]);
	      p_b.Rotate(theta_A2, k_A2);
	      mon.fillHisto("costheta", "0b_all", p_b.CosTheta(), 1);
	    }

	  std::vector<TLorentzVector> b_quarks_boosted;
	  for (size_t i = 0; i < b_quarks.size(); ++i)
	    {
	      TLorentzVector p_b = b_quarks[i];
	      p_b.Boost(-b_H[0], -b_H[1], -b_H[2]);
	      p_b.Rotate(theta_H, k_H);
	      mon.fillHisto("costheta", "0b_H", p_b.CosTheta(), 1);
	      b_quarks_boosted.push_back(p_b);
	    }
	  
	  mon.fillHisto("costheta", "0b1_H", b_quarks_boosted[0].CosTheta(), 1);
	  mon.fillHisto("costheta", "0b2_H", b_quarks_boosted[1].CosTheta(), 1);
	  mon.fillHisto("costheta", "0b3_H", b_quarks_boosted[2].CosTheta(), 1);
	  mon.fillHisto("costheta", "0b4_H", b_quarks_boosted[3].CosTheta(), 1);
	  mon.fillHisto("costheta", "0b_H_avg", (b_quarks_boosted[3].CosTheta()+b_quarks_boosted[2].CosTheta()+b_quarks_boosted[1].CosTheta()+b_quarks_boosted[0].CosTheta())/4, 1);
	 		
	  //----------------------Step 2: Truth (Acceptance Cuts)-----------------------//

	  // b quark multiplicity

	  mon.fillHisto("multi", "b_truth_acc", b_quarks_1_acc.size()+b_quarks_2_acc.size(), 1);

	  if ((b_quarks_1_acc.size()+b_quarks_2_acc.size())>=3 && out_quarks_acc.size()>=2)
	    {

	      //
	      // "Ηiggs" - bbb(b)
	      //

	      // Create a new TLorentzVector to store the sum

	      TLorentzVector p_bbb_b;

	      // Sum the corresponding TLorentzVectors from both vectors
	      for (size_t i = 0; i < b_quarks_1_acc.size(); ++i) p_bbb_b += b_quarks_1_acc[i]; 
	      for (size_t i = 0; i < b_quarks_2_acc.size(); ++i) p_bbb_b += b_quarks_2_acc[i];

	      mon.fillHisto("M", "bbbb_acc", p_bbb_b.M(), 1);
	      mon.fillHisto("pt", "bbbb_acc", p_bbb_b.Pt(), 1);

	      // Find boosting direction and rotation direction and angle (Higgs)
	      std::vector<double> b_bbb_b = findBoostandRotation(p_bbb_b).b;
	      TVector3 k_bbb_b = findBoostandRotation(p_bbb_b).k;
	      double theta_bbb_b = findBoostandRotation(p_bbb_b).theta;

	      //
	      // Outgoing quarks
	      //
	      
	      // Sort the quarks on pT in ascending order
	      std::sort(out_quarks_acc.begin(), out_quarks_acc.end(), [](const TLorentzVector &a, const TLorentzVector &b){
		  return a.Pt() > b.Pt();
		});

	      mon.fillHisto("pt", "q1_acc", out_quarks_acc[0].Pt(), 1);
	      mon.fillHisto("eta", "q1_acc", out_quarks_acc[0].Eta(), 1);
	      mon.fillHisto("phi", "q1_acc", out_quarks_acc[0].Phi(), 1);
	      mon.fillHisto("M", "q1_acc", out_quarks_acc[0].M(), 1);

	      mon.fillHisto("pt", "q2_acc", out_quarks_acc[1].Pt(), 1);
	      mon.fillHisto("eta", "q2_acc", out_quarks_acc[1].Eta(), 1);
	      mon.fillHisto("phi", "q2_acc", out_quarks_acc[1].Phi(), 1);
	      mon.fillHisto("M", "q2_acc", out_quarks_acc[1].M(), 1);

	      mon.fillHisto("deltaEta", "q1q2_acc", abs(out_quarks_acc[0].Eta()-out_quarks_acc[1].Eta()), 1);
	      mon.fillHisto("prodEta", "q1q2_acc", out_quarks_acc[0].Eta()*out_quarks_acc[1].Eta(), 1);
	      mon.fillHisto("Mqq", "q1q2_acc", (out_quarks_acc[0]+out_quarks_acc[1]).M(), 1);

	      //
	      // b quarks
	      //
	      
	      std::vector<TLorentzVector> bbb_b;
	      for (size_t i = 0; i < b_quarks_1_acc.size(); ++i) bbb_b.push_back(b_quarks_1_acc[i]); 
	      for (size_t i = 0; i < b_quarks_2_acc.size(); ++i) bbb_b.push_back(b_quarks_2_acc[i]);

	      if (bbb_b.size()>=4)
		{
		  int random_number_acc = rand_int(0, 3);

		  TLorentzVector b_random_acc = bbb_b[random_number_acc];
		  b_random_acc.Boost(-b_bbb_b[0], -b_bbb_b[1], -b_bbb_b[2]);
		  b_random_acc.Rotate(theta_bbb_b, k_bbb_b);

		  mon.fillHisto("costheta", "0b_random_acc", b_random_acc.CosTheta(), 1);

		  for (size_t i = 0; i<4; i++)
		    {
		      TLorentzVector p_b = bbb_b[i];
		      p_b.Boost(-b_bbb_b[0], -b_bbb_b[1], -b_bbb_b[2]);
		      p_b.Rotate(theta_bbb_b, k_bbb_b);

		      mon.fillHisto("costheta", "0b_4_acc", p_b.CosTheta(), 1);
		      mon.fillHisto("costheta", "0b_acc", p_b.CosTheta(), 1);
		    }
		}
	      else
		{
		  int random_number_acc = rand_int(0, 2);

		  TLorentzVector b_random_acc = bbb_b[random_number_acc];
		  b_random_acc.Boost(-b_bbb_b[0], -b_bbb_b[1], -b_bbb_b[2]);
		  b_random_acc.Rotate(theta_bbb_b, k_bbb_b);

		  mon.fillHisto("costheta", "0b_random_acc", b_random_acc.CosTheta(), 1);
		
		  for (size_t i = 0; i<3; i++)
		    {
		      TLorentzVector p_b = bbb_b[i];
		      p_b.Boost(-b_bbb_b[0], -b_bbb_b[1], -b_bbb_b[2]);
		      p_b.Rotate(theta_bbb_b, k_bbb_b);

		      mon.fillHisto("costheta", "0b_3_acc", p_b.CosTheta(), 1);
		      mon.fillHisto("costheta", "0b_acc", p_b.CosTheta(), 1);
		    }
		}

	      // 4b case: take each pair and transform it back to its parent frame
	      // 3b case: do it for only one pair

	      if (b_quarks_1_acc.size()==2)
		{
		  TLorentzVector p_bb1 = b_quarks_1_acc[0] + b_quarks_1_acc[1];
		
		  // Find boosting direction and rotation direction and angle (bb1)
		  std::vector<double> b_bb1 = findBoostandRotation(p_bb1).b;
		  TVector3 k_bb1 = findBoostandRotation(p_bb1).k;
		  double theta_bb1 = findBoostandRotation(p_bb1).theta;

		  for (size_t i = 0; i<b_quarks_1_acc.size(); i++)
		    {
		      TLorentzVector p_b = b_quarks_1_acc[i];
		      p_b.Boost(-b_bb1[0], -b_bb1[1], -b_bb1[2]);
		      p_b.Rotate(theta_bb1, k_bb1);
		      mon.fillHisto("costheta", "0b_all_acc", p_b.CosTheta(), 1);
		    }
		}	    
	      if (b_quarks_2_acc.size()==2)
		{
		  TLorentzVector p_bb2 = b_quarks_2_acc[0] + b_quarks_2_acc[1];
		
		  // Find boosting direction and rotation direction and angle (bb2)
		  std::vector<double> b_bb2 = findBoostandRotation(p_bb2).b;
		  TVector3 k_bb2 = findBoostandRotation(p_bb2).k;
		  double theta_bb2 = findBoostandRotation(p_bb2).theta;

		  for (size_t i = 0; i<b_quarks_2_acc.size(); i++)
		    {
		      TLorentzVector p_b = b_quarks_2_acc[i];
		      p_b.Boost(-b_bb2[0], -b_bb2[1], -b_bb2[2]);
		      p_b.Rotate(theta_bb2, k_bb2);
		      mon.fillHisto("costheta", "0b_all_acc", p_b.CosTheta(), 1);
		    }
		}

	      // Sorting the b quarks on pT
	      std::sort(bbb_b.begin(), bbb_b.end(), [](const TLorentzVector &a, const TLorentzVector &b){
		  return a.Pt() > b.Pt();
		});

	      mon.fillHisto("pt", "b1_acc", bbb_b[0].Pt(), 1);
	      mon.fillHisto("eta", "b1_acc", bbb_b[0].Eta(), 1);
	      mon.fillHisto("phi", "b1_acc", bbb_b[0].Phi(), 1);
	      mon.fillHisto("M", "b1_acc", bbb_b[0].M(), 1);

	      mon.fillHisto("pt", "b2_acc", bbb_b[1].Pt(), 1);
	      mon.fillHisto("eta", "b2_acc", bbb_b[1].Eta(), 1);
	      mon.fillHisto("phi", "b2_acc", bbb_b[1].Phi(), 1);
	      mon.fillHisto("M", "b2_acc", bbb_b[1].M(), 1);

	      mon.fillHisto("pt", "b3_acc", bbb_b[2].Pt(), 1);
	      mon.fillHisto("eta", "b3_acc", bbb_b[2].Eta(), 1);
	      mon.fillHisto("phi", "b3_acc", bbb_b[2].Phi(), 1);
	      mon.fillHisto("M", "b3_acc", bbb_b[2].M(), 1);
	      
	      if (bbb_b.size()==4)
		{
		  mon.fillHisto("pt", "b4_acc", bbb_b[3].Pt(), 1);
		  mon.fillHisto("eta", "b4_acc", bbb_b[3].Eta(), 1);
		  mon.fillHisto("phi", "b4_acc", bbb_b[3].Phi(), 1);
		  mon.fillHisto("M", "b4_acc", bbb_b[3].M(), 1);
		}
	  
	      // Multiplicity
	      mon.fillHisto("multi", "truth", out_quarks_acc.size()+b_quarks_1_acc.size()+b_quarks_2_acc.size(), 1);
	      
	    } // END if b_quarks >=3

	  
	  //-----------------------Step 3: Final State Analysis-------------------------//

	  // Multiplicity
	  mon.fillHisto("multi", "gen", quarks_gluons.size()+bottoms.size(), 1);

	  // Demand at least 2 quarks-gluons and at least 3 bottoms
	  if (quarks_gluons.size()>=2 && bottoms.size()>=3)
	    {
	    
	      //
	      // Outgoing quarks
	      //
	    
	      // Sort the quarks and gluons on pT 
	      std::sort(quarks_gluons.begin(), quarks_gluons.end(), [](const TLorentzVector &a, const TLorentzVector &b){
		  return a.Pt() > b.Pt();
		});

	      // Plot the first two quarks
	      mon.fillHisto("pt", "q1_gen", quarks_gluons[0].Pt(), 1);
	      mon.fillHisto("eta", "q1_gen", quarks_gluons[0].Eta(), 1);

	      mon.fillHisto("pt", "q2_gen", quarks_gluons[1].Pt(), 1);
	      mon.fillHisto("eta", "q2_gen", quarks_gluons[1].Eta(), 1);

	      mon.fillHisto("deltaEta", "q1q2_gen", abs(quarks_gluons[0].Eta()-quarks_gluons[1].Eta()), 1);
	      mon.fillHisto("prodEta", "q1q2_gen", quarks_gluons[0].Eta()*quarks_gluons[1].Eta(), 1);
	      mon.fillHisto("Mqq", "q1q2_gen", (quarks_gluons[0]+quarks_gluons[1]).M(), 1);

	      //
	      // Bottoms
	      //
	      
	      TLorentzVector bottoms_total;
	      for (const auto& vec : bottoms) bottoms_total += vec;

	      if (bottoms.size()>=4)
		{
		  TLorentzVector bottoms_total_4;
		  for (const auto& vec : bottoms) bottoms_total_4 += vec;
		  mon.fillHisto("M", "bbbb_4_gen", bottoms_total_4.M(), 1);
		}

	      if (bottoms.size()==3)
		{
		  TLorentzVector bottoms_total_3;
		  for (const auto& vec : bottoms) bottoms_total_3 += vec;
		  mon.fillHisto("M", "bbbb_3_gen", bottoms_total_3.M(), 1);
		}

	      mon.fillHisto("M", "bbbb_gen", bottoms_total.M(), 1);
	      mon.fillHisto("pt", "bbbb_gen", bottoms_total.Pt(), 1);
	      mon.fillHisto("eta", "bbbb_gen", bottoms_total.Eta(), 1);

	      // Find boosting direction and rotation direction and angle (reconstructed Higgs)

	      std::vector<double> b = findBoostandRotation(bottoms_total).b;
	      TVector3 k = findBoostandRotation(bottoms_total).k;
	      double theta = findBoostandRotation(bottoms_total).theta;
	   
	      if (bottoms.size()>=4)
		{
		  int random_number_gen = rand_int(0, 3);
		  TLorentzVector b_random_gen = bottoms[random_number_gen];

		  b_random_gen.Boost(-b[0], -b[1], -b[2]);
		  b_random_gen.Rotate(theta, k);

		  mon.fillHisto("costheta", "0b_random_gen", b_random_gen.CosTheta(), 1);
		  mon.fillHisto("costheta", "0b_random_4_gen", b_random_gen.CosTheta(), 1);
		}
	      else
		{
		  int random_number_gen = rand_int(0, 2);
		  TLorentzVector b_random_gen = bottoms[random_number_gen];

		  b_random_gen.Boost(-b[0], -b[1], -b[2]);
		  b_random_gen.Rotate(theta, k);

		  mon.fillHisto("costheta", "0b_random_gen", b_random_gen.CosTheta(), 1);
		  mon.fillHisto("costheta", "0b_random_3_gen", b_random_gen.CosTheta(), 1);
		}
	      
	      // Sorting the b quarks on pT 
	      std::sort(bottoms.begin(), bottoms.end(), [](const TLorentzVector &a, const TLorentzVector &b){
		  return a.Pt() > b.Pt();
		});

	      mon.fillHisto("pt", "b1_gen", bottoms[0].Pt(), 1);
	      mon.fillHisto("eta", "b1_gen", bottoms[0].Eta(), 1);

	      mon.fillHisto("pt", "b2_gen", bottoms[1].Pt(), 1);
	      mon.fillHisto("eta", "b2_gen", bottoms[1].Eta(), 1);

	      mon.fillHisto("pt", "b3_gen", bottoms[2].Pt(), 1);
	      mon.fillHisto("eta", "b3_gen", bottoms[2].Eta(), 1);

	      if (bottoms.size()>=4)
		{
		  mon.fillHisto("pt", "b4_gen", bottoms[3].Pt(), 1);
		  mon.fillHisto("eta", "b4_gen", bottoms[3].Eta(), 1);
		}	   
	      
	    } // END if bottoms >=3
	  
	} // END isSignal for the generator level analysis
      

      //#######################################################################################################################################//
      //--------------------------------------------------------DETECTOR LEVEL ANALYSIS--------------------------------------------------------//
      //#######################################################################################################################################//

      // Define boolean for matching between detector and generator objects
      bool runMatching(runMatchingStudies && isSignal);
	
      // Leptons and jets vectors
      std::vector<TLorentzVector> leptons;
      std::vector<TLorentzVector> muons;
      std::vector<TLorentzVector> electrons;
      
      std::vector<TLorentzVector> jets;
      std::vector<jets_struct> b_jets;
      std::vector<jets_struct> non_b_jets;
      std::vector<jets_struct> q_jets;
      std::vector<jets_struct> untagged_jets;
      std::vector<jets_struct> non_q_jets;

      //----------------------------------------------------------------------------//
      //-----------------------------------Triggers---------------------------------// 
      //----------------------------------------------------------------------------//

      bool hasvbfTrigger1  = (ev.triggerType >> 14 ) & 0x1;
      bool hasvbfTrigger2 = (ev.triggerType >> 15 ) & 0x1;
      bool hasjetsTrigger = (ev.triggerType >> 16) & 0x1;

      //----------------------------------------------------------------------------//
      //------------------------loop over electrons and muons-----------------------// 
      //----------------------------------------------------------------------------//

      // Muons
      for (int imn = 0; imn < ev.mn; imn++)
	{
	  TLorentzVector p_muon;
	  p_muon.SetPxPyPzE(ev.mn_px[imn], ev.mn_py[imn], ev.mn_pz[imn], ev.mn_en[imn]);
	  
	  // Acceptance cuts
	  if (p_muon.Pt()<20. || abs(p_muon.Eta())>2.4) continue;

	  // Id + isolation
	  if (ev.mn_passId[imn] && ev.mn_passIso[imn])
	    {
	      muons.push_back(p_muon);
	      leptons.push_back(p_muon);
	    }
	} // END muon loop
     
      // Printout
      if (iev <5 && verbose)
	{
	  for (int i = 0; i < muons.size(); i++)
	    {
	      printf("Muons have: pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f\n",    
		     muons[i].Pt(),    
		     muons[i].Eta(),  
		     muons[i].Phi(),
		     muons[i].M()  
		     );
	    }
	}

      // Electrons

      bool eleinHEM(false); // to use for HEM veto in 2018

      for (int ien = 0; ien < ev.en; ien++)
	{
	  TLorentzVector p_electron;
	  p_electron.SetPxPyPzE(ev.en_px[ien], ev.en_py[ien], ev.en_pz[ien], ev.en_en[ien]);
	  
	  // Acceptance cuts
	  if (p_electron.Pt()<20. || abs(p_electron.Eta())>2.4) continue;

	  // Id + isolation
	  if (ev.en_passId[ien] && ev.en_passIso[ien])
	    {
	      electrons.push_back(p_electron);
	      leptons.push_back(p_electron);

	      // Check if (at least one) electron is in HEM region
	      if(!eleinHEM && p_electron.Pt()>25. &&
		 -1.57<p_electron.Phi() && p_electron.Phi()<-0.87 &&
		 -3.0<p_electron.Eta() && p_electron.Eta()<-1.4) eleinHEM = true;
	    }
	} // END electrons loop
      
      // Printout
      if (iev <5 && verbose)
	{
	  for (int i = 0; i < electrons.size(); i++)
	    {
	      printf("Electrons have: pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f\n",    
		     electrons[i].Pt(),    
		     electrons[i].Eta(),  
		     electrons[i].Phi(),
		     electrons[i].M()  
		     );
	    }
	}

      //----------------------------------------------------------------------------//
      //--------------------------------loop over jets------------------------------// 
      //----------------------------------------------------------------------------//

      bool jetinHEM(false); // to use for HEM veto in 2018
      
      for (int ijet = 0; ijet < ev.jet; ijet++)
	{
	  TLorentzVector p_jet;
	  p_jet.SetPxPyPzE(ev.jet_px[ijet], ev.jet_py[ijet], ev.jet_pz[ijet], ev.jet_en[ijet]);

	  // Check if (at least one) jet is in HEM region
	  if(!jetinHEM && p_jet.Pt()>25. &&
	     p_jet.Eta()>-3.2 && p_jet.Eta()<-1.2 &&
	     p_jet.Phi()>-1.77 && p_jet.Phi()<-0.67) jetinHEM = true;

	  // Acceptance cuts
	  if (p_jet.Pt()<20. || abs(p_jet.Eta())>4.7) continue;
	  
	  // Id
	  if (!ev.jet_PFTight[ijet]) continue;

	  // Loose PUID
	  float LoosePUID = getSFfrom2DHist(p_jet.Pt(), abs(p_jet.Eta()), PULooseWP);
	  bool passLoosePUID = (ev.jet_puId[ijet] > LoosePUID) || (p_jet.Pt()>50);

	  if (!passLoosePUID) continue;
	  
	  // Cross cleaning with respect to electrons and muons
	  bool overlap = false;

	  for (int imn = 0; imn < muons.size(); imn++)
	    {
	      float dR_jet_mn = getDeltaR(p_jet, muons[imn]);
	      mon.fillHisto("deltaR", "jet_mn_before", dR_jet_mn, 1.);
	      if (dR_jet_mn < 0.4) overlap = true;	     
	      else mon.fillHisto("deltaR", "jet_mn_after", dR_jet_mn, 1.);	   
	    }
	  
	  for (int ien = 0; ien < electrons.size(); ien++)
	    {
	      float dR_jet_en = getDeltaR(p_jet, electrons[ien]);
	      mon.fillHisto("deltaR", "jet_en_before", dR_jet_en, 1.);
	      if (dR_jet_en < 0.4) overlap = true;
	      else mon.fillHisto("deltaR", "jet_en_after", dR_jet_en, 1.);		   
	    }
	  
	  if (!overlap)
	    {
	      jets.push_back(p_jet);
	      
	      // b tagging algorithm
	      double btag_value = ev.jet_btag1[ijet];

	      bool match_bjets = match_jets(p_jet, b_quarks, 0.3);
	      bool match_qjets = match_jets(p_jet, out_quarks, 0.5);
	      
	      // Configure b jets and non b jets
	      if (ev.jet_btag1[ijet] > mediumWP) //medium working point
		{
		  // Demand b jets to have eta less than 2.4
		  if (abs(p_jet.Eta())>2.4) continue;
		  
		  if (runMatching) // This analysis corresponds only to signal (it can be removed when the fine tuning is done) and it matches the tagged jets with the "true" ones
		    {
		      if (match_bjets && !match_qjets)  b_jets.push_back({p_jet, btag_value, 1});
		      else if (match_qjets && !match_bjets)   b_jets.push_back({p_jet, btag_value, 2});
		      else if (match_bjets && match_qjets)
			{
			  if (getDeltaRmin(p_jet, b_quarks).first < getDeltaRmin(p_jet, out_quarks).first) b_jets.push_back({p_jet, btag_value, 1});
			  else b_jets.push_back({p_jet, btag_value, 2});
			}
		      else b_jets.push_back({p_jet, btag_value, -1});
		    } // END matching
		  else b_jets.push_back({p_jet, btag_value, 0});
		}
	      else
		{
		  if (runMatching)
		    {
		      if (match_bjets && !match_qjets)  non_b_jets.push_back({p_jet, btag_value, 1});
		      else if (match_qjets && !match_bjets)   non_b_jets.push_back({p_jet, btag_value, 2});
		      else if (match_bjets && match_qjets)
			{
			  if (getDeltaRmin(p_jet, b_quarks).first < getDeltaRmin(p_jet, out_quarks).first) non_b_jets.push_back({p_jet, btag_value, 1});
			  else non_b_jets.push_back({p_jet, btag_value, 2});
			}
		      else non_b_jets.push_back({p_jet, btag_value, -1});
		    } // END matching
		  else non_b_jets.push_back({p_jet, btag_value, 0});
		  
		} // END b tagging
	      
	    } // END cross cleaning

	} // END jets loop
      
      // Printout
      if (iev <5 && verbose)
	{
	  for (int i = 0; i < jets.size(); i++)
	    {
	      printf("Jets have: pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f\n",    
		     jets[i].Pt(),    
		     jets[i].Eta(),  
		     jets[i].Phi(),
		     jets[i].M()  
		     );
	      cout << endl;
	    }
	}

      // Sort jets on pt
      std::sort(jets.begin(), jets.end(), [](const TLorentzVector &a, const TLorentzVector &b){
	  return a.Pt() > b.Pt();
	});

      // Sort non b jets on pt
      std::sort(non_b_jets.begin(), non_b_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.momentum.Pt() > b.momentum.Pt();
	});

      //----------Configure qjets and untagged jets from non b jets----------//

      if (non_b_jets.size()==1) q_jets.push_back(non_b_jets[0]);
      if (non_b_jets.size()>1)
	{
	  q_jets.push_back(non_b_jets[0]);
	  q_jets.push_back(non_b_jets[1]);
	}
      if (non_b_jets.size()>2)
	{
	  for (size_t i = 2; i < non_b_jets.size(); i++)
	    {
	      untagged_jets.push_back(non_b_jets[i]);
	    }
	}

      // Sort q jets (unnecessary)
      std::sort(q_jets.begin(), q_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.momentum.Pt() > b.momentum.Pt();
	});

      // Remove forward untagged jets!!
      for (int i = untagged_jets.size() - 1; i >= 0; --i) {
	if (abs(untagged_jets[i].momentum.Eta()) > 2.4) {
	  untagged_jets.erase(untagged_jets.begin() + i);
	}
      }
      
      // Sort untagged jets (based on btag)
      std::sort(untagged_jets.begin(), untagged_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.btag > b.btag;
	});      
      
      //----------Configure non q jets from b jets and untagged jets----------//
      non_q_jets.reserve(b_jets.size() + untagged_jets.size());
      non_q_jets.insert(non_q_jets.end(), b_jets.begin(), b_jets.end());
      non_q_jets.insert(non_q_jets.end(), untagged_jets.begin(), untagged_jets.end());

       // Sort non q jets on btag
      std::sort(non_q_jets.begin(), non_q_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.btag > b.btag;
	});

      //
      // 2d Histogramms for jet and bjet pt_gen vs pt_reco
      //

      if (runMatching){
	for (const auto& jet : jets){
	  std::vector<TLorentzVector> jets_gen;
	  jets_gen.reserve(bottoms.size() + quarks_gluons.size());
	  jets_gen.insert(jets_gen.end(), bottoms.begin(), bottoms.end());
	  jets_gen.insert(jets_gen.end(), quarks_gluons.begin(), quarks_gluons.end());
	  
	  if (getDeltaRmin(jet, jets_gen).first < 0.5){
	    int index = getDeltaRmin(jet, jets_gen).second;
	    mon.fillHisto("pt_reco_gen", "jet", jet.Pt(), jets_gen[index].Pt());
	    }
	}//END jets loop

	for (const auto& b_jet : b_jets){	  
	  if (getDeltaRmin(b_jet.momentum, bottoms).first < 0.5){
	    int index = getDeltaRmin(b_jet.momentum, bottoms).second;
	    mon.fillHisto("pt_reco_gen", "bjet", b_jet.momentum.Pt(), bottoms[index].Pt());
	  }
	}// END b_jets loop
	
      }// END runMatching
      
      //
      // DEFINE EVENT CATEGORY
      //
      
      int iEvt = EvtCatIndex(b_jets, non_q_jets, q_jets);

      // here, we denote non outgoing jets with j
      // Event Category : (<2b or <2q or <3j) (2b, 2q, 3j) (2b, 2q, 4j) (2b, 2q, >=5j) (3b, 2q, 3j) (3b, 2q, 4j) (3b, 2q, >=5j) (4b, 2q, 3j) (4b, 2q, 4j) (4b, 2q, >=5j) (>=5b, 2q, 3j) (>=5b, 2q, 4j) (>=5b, 2q, >=5j)
      // iEvt           :         0                1            2            3              4            5            6              7            8            9                10             11             12
      
      if(iEvt>3){
	isSR = true;
	tag = "SR_";
      }
      else tag = "CR_";

      // Fill HT on generator level and on detector level before any cut (study the behaviour of QCD background)
      if (isQCD) {
	mon.fillHisto("HT_studies", "gen", ev.lheHt, weight);

	Float_t H_T_draft_step_0 = 0.;
	for (const auto& p_jet : jets) H_T_draft_step_0 += p_jet.Pt();
      
	mon.fillHisto("HT_studies", "reco", H_T_draft_step_0, weight);
      }
      
      // Compare with generator

      if (runMatching)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets","step0", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets","step0", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets","step0", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets","step0", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets","step0", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets","step0", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets","step0", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets","step0", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets","step0", 2.5, 1);
	    }
	} // END if runMatching

      //----------------------------------------------------------------------------//
      //---------------------------EVENT SELECTION CRITERIA-------------------------// 
      //----------------------------------------------------------------------------//

      // Apply the HEM veto
      if(is2018data && afterRun319077 && (jetinHEM || eleinHEM)) continue;
      if(is2018MC && (jetinHEM || eleinHEM)) weight *= 0.35;

      // Apply QCD k factor
      if(isQCD) weight *= k_factor;

      mon.fillHisto("event_flow", "SR_weighted", 0.5, weight);
      mon.fillHisto("event_flow", "CR_weighted", 0.5, weight);

      mon.fillHisto("event_flow", "SR_unweighted", 0.5, 1);
      mon.fillHisto("event_flow", "CR_unweighted", 0.5, 1);

      // CONTROL PLOT: Lepton multiplicity 
      mon.fillHisto("multi", "el", electrons.size(), weight);
      mon.fillHisto("multi", "mn", muons.size(), weight);
      mon.fillHisto("multi", "lep", leptons.size(), weight);	   

      // STEP1: Recquire zero leptons
      if (leptons.size()!=0) continue;
      
      mon.fillHisto("event_flow", "SR_weighted", 1.5, weight);
      mon.fillHisto("event_flow", "CR_weighted", 1.5, weight);

      mon.fillHisto("event_flow", "SR_unweighted", 1.5, 1);
      mon.fillHisto("event_flow", "CR_unweighted", 1.5, 1);
      
      // Save jets' Pt before and after the trigger for the trigger sensitivity plot
      if(runTriggerStudies){
	
      if(jets.size()>0) mon.fillHisto("pt", "jet1_before_trigger", jets[0].Pt(), 1.);
      if(jets.size()>1) mon.fillHisto("pt", "jet2_before_trigger", jets[1].Pt(), 1.);
      if(jets.size()>2) mon.fillHisto("pt", "jet3_before_trigger", jets[2].Pt(), 1.);
      if(jets.size()>3) mon.fillHisto("pt", "jet4_before_trigger", jets[3].Pt(), 1.);
      if(jets.size()>4) mon.fillHisto("pt", "jet5_before_trigger", jets[4].Pt(), 1.);

      if(hasvbfTrigger1){
	if(jets.size()>0) mon.fillHisto("pt", "jet1_after_vbfTrigger1", jets[0].Pt(), 1.);
	if(jets.size()>1) mon.fillHisto("pt", "jet2_after_vbfTrigger1", jets[1].Pt(), 1.);
	if(jets.size()>2) mon.fillHisto("pt", "jet3_after_vbfTrigger1", jets[2].Pt(), 1.);
	if(jets.size()>3) mon.fillHisto("pt", "jet4_after_vbfTrigger1", jets[3].Pt(), 1.);
	if(jets.size()>4) mon.fillHisto("pt", "jet5_after_vbfTrigger1", jets[4].Pt(), 1.);
      }

      if(hasvbfTrigger2){
	if(jets.size()>0) mon.fillHisto("pt", "jet1_after_vbfTrigger2", jets[0].Pt(), 1.);
	if(jets.size()>1) mon.fillHisto("pt", "jet2_after_vbfTrigger2", jets[1].Pt(), 1.);
	if(jets.size()>2) mon.fillHisto("pt", "jet3_after_vbfTrigger2", jets[2].Pt(), 1.);
	if(jets.size()>3) mon.fillHisto("pt", "jet4_after_vbfTrigger2", jets[3].Pt(), 1.);
	if(jets.size()>4) mon.fillHisto("pt", "jet5_after_vbfTrigger2", jets[4].Pt(), 1.);
      }

      if(hasjetsTrigger){
	if(jets.size()>0) mon.fillHisto("pt", "jet1_after_jetsTrigger", jets[0].Pt(), 1.);
	if(jets.size()>1) mon.fillHisto("pt", "jet2_after_jetsTrigger", jets[1].Pt(), 1.);
	if(jets.size()>2) mon.fillHisto("pt", "jet3_after_jetsTrigger", jets[2].Pt(), 1.);
	if(jets.size()>3) mon.fillHisto("pt", "jet4_after_jetsTrigger", jets[3].Pt(), 1.);
	if(jets.size()>4) mon.fillHisto("pt", "jet5_after_jetsTrigger", jets[4].Pt(), 1.);
      }

      // for the (vbfTrigger1 OR vbfTrigger2) Mqq efficiency curve
      if(q_jets.size()>1){
	float M_qq = (q_jets[0].momentum+q_jets[1].momentum).M();
	mon.fillHisto("Mqq", "before_trigger", M_qq, 1.);
      }
      
      }

      // STEP2:
      // Pass the VBF trigger2
      if (!hasvbfTrigger2) continue;
      
      mon.fillHisto("event_flow", "SR_weighted", 2.5, weight);
      mon.fillHisto("event_flow", "CR_weighted", 2.5, weight);

      mon.fillHisto("event_flow", "SR_unweighted", 2.5, 1);
      mon.fillHisto("event_flow", "CR_unweighted", 2.5, 1);

      // for the vbfTrigger1 OR vbfTrigger2 Mqq efficiency curve
      if(runTriggerStudies){
	if(q_jets.size()>1){
	  float M_qq = (q_jets[0].momentum+q_jets[1].momentum).M();
	  mon.fillHisto("Mqq", "after_trigger", M_qq, 1.);
	}
      }

      // CONTROL PLOT: jet multiplicity
      mon.fillHisto("multi", "jet", jets.size(), weight);
        
      // STEP3: Recquire at least 5 jets
      if(jets.size()<5) continue;
      
      mon.fillHisto("event_flow", "SR_weighted", 3.5, weight);
      mon.fillHisto("event_flow", "CR_weighted", 3.5, weight);

      mon.fillHisto("event_flow", "SR_unweighted", 3.5, 1);
      mon.fillHisto("event_flow", "CR_unweighted", 3.5, 1);

      // Compare with generator

      if (runMatching)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets","step3", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets","step3", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets","step3", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets","step3", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets","step3", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets","step3", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets","step3", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets","step3", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets","step3", 2.5, 1);
	    }
	} // END if runMatching

      //STEP4: Apply pt thresholds for jets
      if(jets[0].Pt()<110 || jets[1].Pt()<90 || jets[2].Pt()<80 || jets[3].Pt()<30) continue;

      mon.fillHisto("event_flow", "SR_weighted", 4.5, weight);
      mon.fillHisto("event_flow", "CR_weighted", 4.5, weight);

      mon.fillHisto("event_flow", "SR_unweighted", 4.5, 1);
      mon.fillHisto("event_flow", "CR_unweighted", 4.5, 1);
      
      // CONTROL PLOT: b and q jet multiplicity
      mon.fillHisto("multi", "bjet", b_jets.size(), weight);
      mon.fillHisto("multi", "qjet", q_jets.size(), weight);
      mon.fillHisto("multi", "untaggedjets", untagged_jets.size(), weight);
      
      //STEP5: Recquire at least 2 q tagged jets AND at least 2 b tagged jets AND at least 3 non q tagged jets
      if(iEvt==0) continue;

      mon.fillHisto("event_flow", "SR_weighted", 5.5, weight);
      mon.fillHisto("event_flow", "CR_weighted", 5.5, weight);

      mon.fillHisto("event_flow", "SR_unweighted", 5.5, 1);
      mon.fillHisto("event_flow", "CR_unweighted", 5.5, 1);

      //STEP6: SR or CR
      mon.fillHisto("event_flow", tag + "weighted", 6.5, weight);
      mon.fillHisto("event_flow", tag + "unweighted", 6.5, 1);

      // CONTROL PLOT: b jet, q jet and untagged jet multiplicity
      mon.fillHisto("multi", tag + "bjet_step5", b_jets.size(), weight);
      mon.fillHisto("multi", tag + "qjet_step5", q_jets.size(), weight);
      mon.fillHisto("multi", tag + "untaggedjet_step5", untagged_jets.size(), weight);         

      // Compare with generator

      if (runMatching)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets", tag + "step5", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets", tag + "step5", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets", tag + "step5", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets", tag + "step5", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets", tag + "step5", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets", tag + "step5", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets", tag + "step5", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets", tag + "step5", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets", tag + "step5", 2.5, 1);
	    }
	} // END if runMatching
      
      //sort b jets on btag
      std::sort(b_jets.begin(), b_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.btag > b.btag;
	});

      // // CONTROL PLOT: btag values
      // mon.fillHisto("btag","bjet1_before", b_jets[0].btag, weight);
      // mon.fillHisto("btag","bjet2_before", b_jets[1].btag, weight);
      // mon.fillHisto("btag","bjet3_before", b_jets[2].btag, weight);
      
      // // STEP6:
      // // TMML (inactive)
      // if (b_jets[0].btag < mediumWP || b_jets[2].btag < mediumWP) continue;
      // //if (b_jets[0].btag >= tightWP) continue;
      // }

      // mon.fillHisto("event_flow", tag + "weighted", 6.5, weight);
      // mon.fillHisto("event_flow", tag + "unweighted", 6.5, 1);

      // // Compare with generator

      // if (runMatching)
      // 	{
      // 	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

      // 	  //b tagged
      // 	  for (size_t i = 0; i < b_jets.size(); i++)
      // 	    {
      // 	      int category = b_jets[i].category;
      // 	      if (category==1) mon.fillHisto("match_bjets","step6", 0.5, 1);
      // 	      if (category==2) mon.fillHisto("match_bjets","step6", 1.5, 1);
      // 	      if (category==-1) mon.fillHisto("match_bjets","step6", 2.5, 1);
      // 	    }
	  
      // 	  //q tagged
      // 	  for (size_t i = 0; i < q_jets.size(); i++)
      // 	    {
      // 	      int category = q_jets[i].category;
      // 	      if (category==1) mon.fillHisto("match_qjets","step6", 0.5, 1);
      // 	      if (category==2) mon.fillHisto("match_qjets","step6", 1.5, 1);
      // 	      if (category==-1) mon.fillHisto("match_qjets","step6", 2.5, 1);
      // 	    }
	  
      // 	  //untagged
      // 	  for (size_t i = 0; i < untagged_jets.size(); i++)
      // 	    {
      // 	      int category = untagged_jets[i].category;
      // 	      if (category==1) mon.fillHisto("match_untaggedjets","step6", 0.5, 1);
      // 	      if (category==2) mon.fillHisto("match_untaggedjets","step6", 1.5, 1);
      // 	      if (category==-1) mon.fillHisto("match_untaggedjets","step6", 2.5, 1);
      // 	    }
      // 	} // END if runMatching

      // Define vbf jets quantities
      DeltaEta = abs(q_jets[0].momentum.Eta()-q_jets[1].momentum.Eta());
      ProductEta = q_jets[0].momentum.Eta()*q_jets[1].momentum.Eta();
      qq_m = (q_jets[0].momentum+q_jets[1].momentum).M();

      // CONTROL PLOT: qq mass, eta difference and eta product plots
      mon.fillHisto("deltaEta", tag + "qjet1_qjet2_before", DeltaEta, weight);
      mon.fillHisto("prodEta", tag + "qjet1_qjet2_before", ProductEta, weight);
      mon.fillHisto("Mqq", tag + "qjet1_qjet2_before", qq_m, weight);

      //STEP6: Recquire Mqq more than 250 GeV and Dhqq more than 2.5
      
      if (DeltaEta<2.5 || qq_m<250) continue;

      double deltaPhibjet1bjet2 = abs(non_q_jets[0].momentum.Phi() - non_q_jets[1].momentum.Phi());
      double deltaEtabjet1bjet2 = abs(non_q_jets[0].momentum.Eta() - non_q_jets[1].momentum.Eta());
      double deltaRbjet1bjet2 = getDeltaR(non_q_jets[0].momentum, non_q_jets[1].momentum);
		   
      mon.fillHisto("event_flow", tag + "weighted", 7.5, weight);
      mon.fillHisto("event_flow", tag + "unweighted", 7.5, 1);

      // Compare with generator

      if (runMatching)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets", tag + "step6", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets", tag + "step6", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets", tag + "step6", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets", tag + "step6", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets", tag + "step6", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets", tag + "step6", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets", tag + "step6", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets", tag + "step6", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets", tag + "step6", 2.5, 1);
	    }

	  //untagged n.0 matchings 

	  if (untagged_jets.size()>0 && (abs(untagged_jets[0].momentum.Eta()) < 2.4 || abs(untagged_jets[0].momentum.Eta()) > 3.3)) 
	    {

	      if (untagged_jets[0].category == 1)
		{
		  mon.fillHisto(TString("match_eta"), tag + "untagged0", untagged_jets[0].momentum.Eta(), 0.5, 1.);
		  mon.fillHisto(TString("match_pt"), tag + "untagged0", untagged_jets[0].momentum.Pt(), 0.5, 1.);
		  mon.fillHisto(TString("match_untaggedjets"), tag + "0_step6", 0.5, 1);
		}
	      if (untagged_jets[0].category == 2)
		{
		  mon.fillHisto(TString("match_eta"), tag + "untagged0", untagged_jets[0].momentum.Eta(), 1.5, 1.);
		  mon.fillHisto(TString("match_pt"), tag + "untagged0", untagged_jets[0].momentum.Pt(), 1.5, 1.);
		  mon.fillHisto(TString("match_untaggedjets"), tag + "0_step6", 1.5, 1);
		}
	      if (untagged_jets[0].category == -1)
		{
		  mon.fillHisto(TString("match_eta"), tag + "untagged0", untagged_jets[0].momentum.Eta(), 2.5, 1.);
		  mon.fillHisto(TString("match_pt"), tag + "untagged0", untagged_jets[0].momentum.Pt(), 2.5, 1.);
		  mon.fillHisto(TString("match_untaggedjets"), tag + "0_step6", 2.5, 1);
		}		  
	    }

	  // Calculate the minimum angular distance between each jet and the collection of generated b's and q's (only for signal)

	  for (size_t i=0; i<jets.size(); i++)
	    {
	      TLorentzVector p_jet = jets[i];
	      
	      double dRmin_jet_b = getDeltaRmin(p_jet,b_quarks).first;
	      double dRmin_jet_q = getDeltaRmin(p_jet,out_quarks).first;

	      mon.fillHisto("deltaR", tag + "min_reco_b", dRmin_jet_b, 1);
	      mon.fillHisto("deltaR", tag + "min_reco_q", dRmin_jet_q, 1);
	    }

	  for (size_t i=0; i<b_jets.size(); i++)
	    {
	      TLorentzVector p_jet = b_jets[i].momentum;

	      double dRmin_bjet_b = getDeltaRmin(p_jet,b_quarks).first;
	      double dRmin_bjet_q = getDeltaRmin(p_jet,out_quarks).first;

	      int iter_b = getDeltaRmin(p_jet,b_quarks).second;

	      mon.fillHisto("deltaR", tag + "min_bjet_b", dRmin_bjet_b, 1);
	      mon.fillHisto("deltaR", tag + "min_bjet_q", dRmin_bjet_q, 1);
	      mon.fillHisto(TString("mindeltaR_vs_pt"), tag + "bjet_b", b_quarks[iter_b].Pt(), dRmin_bjet_b, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), tag + "bjet_b", b_quarks[iter_b].Eta(), dRmin_bjet_b, 1.);
	    }

	  for (size_t i=0; i<q_jets.size(); i++)
	    {
	      TLorentzVector p_jet = q_jets[i].momentum;
       
	      double dRmin_qjet_b =  getDeltaRmin(p_jet,b_quarks).first;
	      double dRmin_qjet_q = getDeltaRmin(p_jet,out_quarks).first;

	      int iter_q = getDeltaRmin(p_jet,out_quarks).second;

	      mon.fillHisto("deltaR", tag + "min_qjet_b", dRmin_qjet_b, 1);
	      mon.fillHisto("deltaR", tag + "min_qjet_q", dRmin_qjet_q, 1);
	      mon.fillHisto(TString("mindeltaR_vs_pt"), tag + "qjet_q", out_quarks[iter_q].Pt(), dRmin_qjet_q, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), tag + "qjet_q", out_quarks[iter_q].Eta(), dRmin_qjet_q, 1.);

	    }

	  if (untagged_jets.size()>0)
	    {		
	      TLorentzVector p_jet_0 = untagged_jets[0].momentum;
	      
	      double dRmin_0_b = getDeltaRmin(p_jet_0,b_quarks).first;
	      double dRmin_0_q = getDeltaRmin(p_jet_0,out_quarks).first;

	      int iter_q = getDeltaRmin(p_jet_0,out_quarks).second;
	      int iter_b = getDeltaRmin(p_jet_0,b_quarks).second;

	      mon.fillHisto(TString("mindeltaR_vs_pt"), tag + "0_b", b_quarks[iter_b].Pt(), dRmin_0_b, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_pt"), tag + "0_q", out_quarks[iter_q].Pt(), dRmin_0_q, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), tag + "0_b", b_quarks[iter_b].Eta(), dRmin_0_b, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), tag + "0_q", out_quarks[iter_q].Eta(), dRmin_0_q, 1.);
	    
	      mon.fillHisto("deltaR", tag + "min_0_b", dRmin_0_b, 1);
	      mon.fillHisto("deltaR", tag + "min_0_q", dRmin_0_q, 1);
	    }
	  
	} // END if runMatching

      // "Higgs" four vector and vector of four b's
      TLorentzVector bjet_total;
      std::vector<TLorentzVector> bjet_total_vec;

      if((iEvt==1) || (iEvt>=4 && iEvt<=6)){
	for (int i = 0; i < 3; ++i){
	  bjet_total += non_q_jets[i].momentum;
	  bjet_total_vec.push_back(non_q_jets[i].momentum);
	}
      }
      else{
	for (int i = 0; i < 4; ++i){
	  bjet_total += non_q_jets[i].momentum;
	  bjet_total_vec.push_back(non_q_jets[i].momentum);
	}
      }

      //
      // Compute the minimum deltaPhi bewtween all the possible b combinations
      //

      // Variable to store the result
      double MinDeltaPhibb(1e6);
      
      for (size_t i = 0; i < bjet_total_vec.size(); ++i) {
        for (size_t j = i + 1; j < bjet_total_vec.size(); ++j) {
	  double deltaPhi = abs(bjet_total_vec[i].Phi()-bjet_total_vec[j].Phi());
	  MinDeltaPhibb = MinDeltaPhibb < deltaPhi ? MinDeltaPhibb : deltaPhi;
        }
      }

      minDeltaPhi = MinDeltaPhibb;

      // CONTROL PLOT: minDeltaPhi
      mon.fillHisto("deltaPhi", tag + "bjet1_bjet2_min_before", minDeltaPhi, weight);

      //STEP7: Recquire deltaPhibb to be less than 1.6
      if(minDeltaPhi>1.6) continue;

      mon.fillHisto("event_flow", tag + "weighted", 8.5, weight);
      mon.fillHisto("event_flow", tag + "unweighted", 8.5, 1);

      //
      // bbb(b) ("higgs") kinematic variables
      //
      
      higgs_m = bjet_total.M();
      higgs_pt = bjet_total.Pt();
      higgs_eta = bjet_total.Eta();
 
      //
      // qq quantities
      //

      qq_dR = getDeltaR(q_jets[0].momentum, q_jets[1].momentum);
      qq_dphi = abs(q_jets[0].momentum.Phi() - q_jets[1].momentum.Phi());

      //
      // jet kinematics
      //

      jet1_pt = jets[0].Pt();
      jet2_pt = jets[1].Pt();
      jet3_pt = jets[2].Pt();
      jet4_pt = jets[3].Pt();
      jet5_pt = jets[4].Pt();

      jet1_eta = jets[0].Eta();
      jet2_eta = jets[1].Eta();
      jet3_eta = jets[2].Eta();
      jet4_eta = jets[3].Eta();
      jet5_eta = jets[4].Eta();

      if (jets.size()>=6)
	{
	  jet6_pt = jets[5].Pt();
	  jet6_eta = jets[5].Eta();
	}
      else
	{
	  jet6_pt = std::nan("");
	  jet6_eta = std::nan("");
	}

      //
      // bjet, qjet kinematics (pt,eta) and btag values
      //

      bjet1_pt = non_q_jets[0].momentum.Pt();
      bjet2_pt = non_q_jets[1].momentum.Pt();
      bjet3_pt = non_q_jets[2].momentum.Pt();
      
      bjet1_eta = non_q_jets[0].momentum.Eta();
      bjet2_eta = non_q_jets[1].momentum.Eta();
      bjet3_eta = non_q_jets[2].momentum.Eta();

      bjet1_btag = non_q_jets[0].btag;
      bjet2_btag = non_q_jets[1].btag;
      bjet3_btag = non_q_jets[2].btag;

      qjet1_pt  = q_jets[0].momentum.Pt();
      qjet1_eta = q_jets[0].momentum.Eta();
      qjet2_pt  = q_jets[1].momentum.Pt();
      qjet2_eta = q_jets[1].momentum.Eta();

      // Calculate qq pair vector
      TLorentzVector qjet_total = q_jets[0].momentum + q_jets[1].momentum;

      //
      // Compute azimuthal angle between qq pair and reco Higgs
      //

      phi_qq_higgs = abs(bjet_total.Phi()-qjet_total.Phi());

      //
      // Compute angles of vbf jets with respect to the qq pair rest frame
      //

      // qq pair rest frame
      
      std::vector<double> b_qjet = findBoostandRotation(qjet_total).b;
      TVector3 k_qjet= findBoostandRotation(qjet_total).k;
      double theta_qjet = findBoostandRotation(qjet_total).theta;
      
      TLorentzVector qjet1_boosted = q_jets[0].momentum;
      TLorentzVector qjet2_boosted = q_jets[1].momentum;

      qjet1_boosted.Boost(-b_qjet[0], -b_qjet[1], -b_qjet[2]);
      qjet1_boosted.Rotate(theta_qjet, k_qjet);

      qjet2_boosted.Boost(-b_qjet[0], -b_qjet[1], -b_qjet[2]);
      qjet2_boosted.Rotate(theta_qjet, k_qjet);

      double alpha_q1 = qjet1_boosted.Theta();
      double alpha_q2 = qjet2_boosted.Theta();

      alpha_qq = alpha_q1 < alpha_q2 ? alpha_q1 : alpha_q2;

      //
      // Compute the deltaR between higgs and the two vbf jets
      //

      dR_higgs_q1 = getDeltaR(bjet_total, q_jets[0].momentum);
      dR_higgs_q2 = getDeltaR(bjet_total, q_jets[1].momentum);

      //
      // Compute the average delta R of the b's
      //
      
      // Variable to hold the sum of all Delta R values
      double sumDeltaR = 0.0;
      int count_pairs = 0;

      for (size_t i = 0; i < bjet_total_vec.size(); ++i) {
        for (size_t j = i + 1; j < bjet_total_vec.size(); ++j) {
	  double deltaR = getDeltaR(bjet_total_vec[i], bjet_total_vec[j]);
	  sumDeltaR += deltaR;
	  ++count_pairs; // Count the number of pairs
        }
      }
      
      // Calculate the average Delta R
      avDeltaR = sumDeltaR / count_pairs;
      
      //
      // Calculate DeltaMmin
      //
      
      double unpairedMass;
      double minDeltaMold = getDeltaMmin(bjet_total_vec, unpairedMass);
      minDeltaM = getDeltaMminNew(bjet_total_vec);

      //
      // Calculate Mbbj
      //

      MassOfbbj = getMassOfbbj(non_q_jets, iEvt);

      // Find angle of a random b in bbb(b) ("Higgs") rest frame
      
      std::vector<double> b_bjet = findBoostandRotation(bjet_total).b;
      TVector3 k_bjet= findBoostandRotation(bjet_total).k;
      double theta_bjet = findBoostandRotation(bjet_total).theta;     
      
      // Generate a random number
      int random_number_jet;
      if (bjet_total_vec.size()>3) random_number_jet = rand_int(0, 3);
      else random_number_jet = rand_int(0, 2);
      
      TLorentzVector bjet_random = bjet_total_vec[random_number_jet];

      bjet_random.Boost(-b_bjet[0], -b_bjet[1], -b_bjet[2]);
      bjet_random.Rotate(theta_bjet, k_bjet);

      costheta0 = bjet_random.CosTheta();
    
      //
      // Find H_T
      //
      
      Float_t H_T_draft = 0.;
      for (const auto& p_jet : jets) H_T_draft += p_jet.Pt();
      H_T = H_T_draft;

      //
      // Find Sum(p->)_T / Sum(p_T)
      //
      
      TLorentzVector jet_total;
      for (size_t i = 0; i<jets.size(); i++) jet_total+=jets[i];

      H_Tvec = jet_total.Pt() / H_T;

      //
      // Find pz_q1 + pz_q2 + pz_b1 + pz_b2 + pz_b3 + (pz_b4)
      //
      
      Float_t H_z_draft = 0.;

      if(iEvt==1 || (iEvt>=4 && iEvt<=6)){
	for (size_t i = 0; i<3; i++) H_z_draft+=non_q_jets[i].momentum.Pz();
      }
      else{
	for (size_t i = 0; i<3; i++) H_z_draft+=non_q_jets[i].momentum.Pz();
      }
     
      for (size_t i = 0; i<2; i++) H_z_draft+=q_jets[i].momentum.Pz();

      H_z = H_z_draft;
      
      //
      // Find final jet multiplicity
      //
      
      N_jet = jets.size();
      N_bjet = b_jets.size();
      N_qjet = q_jets.size();
      N_untaggedjet = untagged_jets.size();

      int jet_multi_eta_cut(0);
      
      for (size_t i = 0; i<jets.size(); i++) {
	if (abs(jets[i].Eta()<2.4)) jet_multi_eta_cut++;
      }
      
      N_jet_eta_cut = jet_multi_eta_cut;

      //
      //find met pt and phi
      //
      
      MET_pt = ev.met_pt;
      MET_phi = ev.met_phi;

      //
      // Find p_T sum and energy sum of the untagged jets
      //

      double pt_rest_draft = 0;
      double E_rest_draft = 0;
      if(iEvt==5 || iEvt==6){
	for (size_t i = 3; i<non_q_jets.size(); i++)
	  {
	    pt_rest_draft += non_q_jets[i].momentum.Pt();
	    E_rest_draft += non_q_jets[i].momentum.E();
	  }
      }
      if(iEvt==3 || iEvt>=9){
	for (size_t i = 4; i<non_q_jets.size(); i++)
	  {
	    pt_rest_draft += non_q_jets[i].momentum.Pt();
	    E_rest_draft += non_q_jets[i].momentum.E();
	  }
      }
      
      pt_rest = pt_rest_draft;
      E_rest = E_rest_draft;

      //----------------------------------------------------------------------------//
      //-------------------------------Fill histogramms-----------------------------//
      //----------------------------------------------------------------------------//
      
      //
      // qq quantities
      //

      mon.fillHisto("deltaEta", tag + "qjet1_qjet2", DeltaEta, weight);
      mon.fillHisto("prodEta", tag + "qjet1_qjet2", ProductEta, weight);
      mon.fillHisto("Mqq", tag + "qjet1_qjet2", qq_m, weight);
      mon.fillHisto("deltaR", tag + "qjet1_qjet2", qq_dR, weight);
      mon.fillHisto("deltaPhi", tag + "qjet1_qjet2", qq_dphi, weight);

      if(hasvbfTrigger1) mon.fillHisto("Mqq", tag + "qjet1_qjet2_trigger1", qq_m, weight);
      if(hasvbfTrigger2) mon.fillHisto("Mqq", tag + "qjet1_qjet2_trigger2", qq_m, weight);

      //
      // jet kinematics
      //

      mon.fillHisto("pt", tag + "jet1", jet1_pt, weight);
      mon.fillHisto("pt", tag + "jet2", jet2_pt, weight);
      mon.fillHisto("pt", tag + "jet3", jet3_pt, weight);
      mon.fillHisto("pt", tag + "jet4", jet4_pt, weight);
      mon.fillHisto("pt", tag + "jet5", jet5_pt, weight);

      mon.fillHisto("eta", tag + "jet1", jet1_eta, weight);
      mon.fillHisto("eta", tag + "jet2", jet2_eta, weight);
      mon.fillHisto("eta", tag + "jet3", jet3_eta, weight);
      mon.fillHisto("eta", tag + "jet4", jet4_eta, weight);
      mon.fillHisto("eta", tag + "jet5", jet5_eta, weight);

      if (jets.size()>=6)
	{
	  mon.fillHisto("pt", tag + "jet6", jet6_pt, weight);
	  mon.fillHisto("eta", tag + "jet6", jet6_eta, weight);
	}

      //
      // bjet, qjet kinematics (pt,eta) and btag values
      //

      mon.fillHisto("pt", tag + "bjet1", bjet1_pt, weight);
      mon.fillHisto("pt", tag + "bjet2", bjet2_pt, weight);
      mon.fillHisto("pt", tag + "bjet3", bjet3_pt, weight);

      mon.fillHisto("eta", tag + "bjet1", bjet1_eta, weight);
      mon.fillHisto("eta", tag + "bjet2", bjet2_eta, weight);
      mon.fillHisto("eta", tag + "bjet3", bjet3_eta, weight);

      mon.fillHisto("btag", tag + "bjet1", bjet1_btag, weight);
      mon.fillHisto("btag", tag + "bjet2", bjet2_btag, weight);
      mon.fillHisto("btag", tag + "bjet3", bjet3_btag, weight);

      mon.fillHisto("deltaPhi", tag + "bjet1bjet2", deltaPhibjet1bjet2, weight);
      mon.fillHisto("deltaEta", tag + "bjet1bjet2", deltaEtabjet1bjet2, weight);
      mon.fillHisto("deltaR", tag + "bjet1bjet2", deltaRbjet1bjet2, weight);

      mon.fillHisto("pt", tag + "qjet1", qjet1_pt, weight);
      mon.fillHisto("pt", tag + "qjet2", qjet2_pt, weight);

      mon.fillHisto("eta", tag + "qjet1", qjet1_eta, weight);
      mon.fillHisto("eta", tag + "qjet2", qjet2_eta, weight);

      //
      // bbb(b) ("higgs") kinematic variables
      //

      if (b_jets.size()>=4)
	{
	  mon.fillHisto("M", tag + "bjet_total_4", bjet_total.M(), weight);
	}
  
      mon.fillHisto("M", tag + "bjet_total", higgs_m, weight);
      mon.fillHisto("eta", tag + "bjet_total", higgs_eta, weight);
      mon.fillHisto("pt", tag + "bjet_total", higgs_pt, weight);

      //
      // Angle between qq pair and reco Higgs
      //

      mon.fillHisto("deltaPhi", tag + "qq_higgs", phi_qq_higgs, weight);

      //
      // Angle of vbf jets with respect to the qq pair rest frame
      //
      
      mon.fillHisto("theta", tag + "alpha_qq", alpha_qq, weight);

      //
      // deltaR between higgs and the two vbf jets
      //

      mon.fillHisto("deltaR", tag + "H_qjet1", dR_higgs_q1, weight);
      mon.fillHisto("deltaR", tag + "H_qjet2", dR_higgs_q2, weight);

      //
      // Minimum bb pair azimuthal angle
      //
      mon.fillHisto("deltaPhi", tag + "bjet1_bjet2_min", minDeltaPhi, weight);
      
      //
      // Average delta R of the b's
      //
      
      mon.fillHisto("deltaR", tag + "av_bb", avDeltaR, weight);

      //
      // DeltaMmin
      //
      
      mon.fillHisto("M", tag + "Delta_pair_min", minDeltaMold, weight);
      mon.fillHisto("M", tag + "Delta_pair_min_new", minDeltaM, weight);

      if (bjet_total_vec.size()==3) mon.fillHisto("M", tag + "unpaired_b", unpairedMass, weight);

      //
      // Mbbj
      //

      mon.fillHisto("M", tag + "bbj", MassOfbbj, weight);	  

      //
      // costheta0
      //
      
      mon.fillHisto("costheta", tag + "0_bjet", costheta0, weight);

      //
      // H_T, H_z, HT_vec
      //
      
      mon.fillHisto("HT", tag, H_T, weight);
      mon.fillHisto("Hz", tag, H_z, weight);
      mon.fillHisto("HT_vec", tag, H_Tvec, weight);	    

      //
      // Multiplicities
      //
      
      mon.fillHisto("multi", tag + "jet_final", N_jet, weight);
      mon.fillHisto("multi", tag + "jet_final_eta_cuts", jet_multi_eta_cut, weight);
      mon.fillHisto("multi", tag + "bjet_final", N_bjet, weight);
      mon.fillHisto("multi", tag + "qjet_final", N_qjet, weight);
      mon.fillHisto("multi", tag + "untaggedjet_final", N_untaggedjet, weight);

      //
      // Met
      //
      
      mon.fillHisto("pt", tag + "met", MET_pt, weight);
      mon.fillHisto("phi", tag + "met", MET_phi, weight);

      //
      // Rest
      //
      
      mon.fillHisto("pt", tag + "rest", pt_rest, weight);
      mon.fillHisto("E", tag + "rest", E_rest, weight);

      //
      // Multiplicity in all regions
      //
      mon.fillHisto("event_cat", "weighted", iEvt-0.5, weight);
      mon.fillHisto("event_cat", "unweighted", iEvt-0.5, 1.);
      
      //#######################################################################################################################################//
      //--------------------------------------------------------------TMVA Reader--------------------------------------------------------------//
      //#######################################################################################################################################//

      float mvaBDT;
      mvaBDT = VBFhAnalysisReader.GenReMVAReader(higgs_m, higgs_pt, dR_higgs_q2, phi_qq_higgs, avDeltaR, minDeltaM, MassOfbbj, minDeltaPhi, qq_dphi, alpha_qq, jet3_pt, jet5_pt, H_T, DeltaEta, qq_m, pt_rest, E_rest, qjet2_pt, "VBFhAnalysisClass");
      mon.fillHisto("BDT", tag + "VBFhAnalysis", mvaBDT, weight);
  
      //#######################################################################################################################################//
      //--------------------------------------------------------------MVA Handler--------------------------------------------------------------//
      //#######################################################################################################################################//

      if (runMVA && (isSR != runQCD))
	{
	  myMVAHandler_.getEntry(higgs_m, higgs_pt, higgs_eta, costheta0, dR_higgs_q1, dR_higgs_q2, phi_qq_higgs, avDeltaR, minDeltaM, MassOfbbj, minDeltaPhi, DeltaEta, qq_m, ProductEta, qq_dR, qq_dphi, alpha_qq, N_jet, N_jet_eta_cut, jet1_pt, jet2_pt, jet3_pt, jet4_pt, jet5_pt, jet6_pt, jet1_eta, jet2_eta, jet3_eta, jet4_eta, jet5_eta, jet6_eta, H_T, H_z, H_Tvec, pt_rest, E_rest, N_bjet, bjet1_pt, bjet2_pt, bjet3_pt, bjet1_eta, bjet2_eta, bjet3_eta, bjet1_btag, bjet2_btag, bjet3_btag, N_qjet, qjet1_pt, qjet2_pt, qjet1_eta, qjet2_eta, N_untaggedjet, MET_pt, MET_phi, weight);
	  myMVAHandler_.fillTree();
	} // runMVA
      
    } // END loop over entries

  //#######################################################################################################################################//
  //-----------------------------------------------------------------SAVE------------------------------------------------------------------//
  //#######################################################################################################################################//

  // Write MVA files
  TString mvaout = TString ( runProcess.getParameter<std::string>("outdir") ) + "/mva_" + outFileUrl + ".root";
  if(runMVA){ myMVAHandler_.writeTree(mvaout); }

  outUrl += "/";
  outUrl += outFileUrl + ".root";
  printf("Results saved in %s\n", outUrl.Data());

  // Save plots to file
  int nTrial = 0;
  TFile *ofile=TFile::Open(outUrl, "recreate");
  while( !ofile->IsOpen() || ofile->IsZombie() )
    {
      if(nTrial > 3)
	{
	  printf("Output file open failed!");
	  if ( outTxtFile_final ) fclose(outTxtFile_final);
	  return -1;
	}
      nTrial++;
      usleep(1000000*nTrial);
      ofile=TFile::Open(outUrl, "update");
    }
  mon.Write();
  ofile->Close();
  if( outTxtFile_final ) fclose( outTxtFile_final );

  //#######################################################################################################################################//
  //---------------------------------------------------------------PRINTOUTS---------------------------------------------------------------//
  //#######################################################################################################################################//

  std::cout << "number of events: " << nevts << std::endl;
  cout << endl;
  cout << "weight is equal to: " << weight << endl;
  cout << endl;
  
} // END main function

//#######################################################################################################################################//
//---------------------------------------------------------------FUNCTIONS---------------------------------------------------------------//
//#######################################################################################################################################//

// Function to calculate the deltaR distance in the (eta, phi) plane between two four-vectors
double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2)
  
{
  double delta_phi;
  double delta_eta;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_eta = vec_1.Eta() - vec_2.Eta();

  return std::sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function to calculate the deltaR in the (psi, phi) plane between two four-vectors
double getDeltaRyy(TLorentzVector vec_1, TLorentzVector vec_2)
  
{
  double delta_phi;
  double delta_psi;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_psi = vec_1.Rapidity() - vec_2.Rapidity();

  return std::sqrt(delta_phi * delta_phi + delta_psi * delta_psi);
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function to find boosting direction and rotation direction and angle
triplet findBoostandRotation(TLorentzVector p_mom)
{
  //boosting (b) vector
  double b1 = p_mom.Px() / p_mom.E();
  double b2 = p_mom.Py() / p_mom.E();
  double b3 = p_mom.Pz() / p_mom.E();

  std::vector<double> b;
  b.push_back(b1);
  b.push_back(b2);
  b.push_back(b3);
  
  //rotation vector (k)
  double k1 = -p_mom.Py();
  double k2 = p_mom.Px();
  double k3 = 0;
  double k_mag = sqrt(pow(k1,2)+pow(k2,2)+pow(k3,2));
  TVector3 k(k1/k_mag, k2/k_mag, k3/k_mag);

  //rotation angle
  double p_mag = sqrt(pow(p_mom.Px(),2)+pow(p_mom.Py(),2)+pow(p_mom.Pz(),2));
  double theta = TMath::ACos(p_mom.Pz()/p_mag);

  triplet trip = {b, k, theta};
  return trip;
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function that returns true if a four-vector has an angular distance less than a threshold value with at least one four-vector within a vector
bool match_jets(TLorentzVector det_vec, std::vector<TLorentzVector> gen_vec, double threshold)
{
  for (size_t i = 0; i < gen_vec.size(); i++)
  {
    // Calculate deltaR between det_vec and the current element in gen_vec
    double deltaR = getDeltaR(gen_vec[i], det_vec);

    // If a match is found, return true immediately
    if (deltaR < threshold) 
      return true;
  }

  // If no match is found, return false
  return false;
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function that returns the minimum angular distance between a four-vector and a collection of four-vectors along with the position of the four-vector that minimizes the distance
std::pair<double, int> getDeltaRmin(TLorentzVector det_vec, std::vector<TLorentzVector> gen_vec)
{
  // Initialize dR and index of minimum dR
  double DeltaRmin(1000);
  int minIndex(-1);
      
  for (size_t i = 0; i < gen_vec.size(); i++)
    {
      // Calculate deltaR between det_vec and the current element in gen_vec
      double deltaR = getDeltaR(gen_vec[i], det_vec);
      if (deltaR<DeltaRmin)
	{
	  DeltaRmin=deltaR;
	  minIndex = i;
	}
    }
  return std::make_pair(DeltaRmin, minIndex);
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function to generate a random integer number between two values
int rand_int(int min, int max)
{
  // Random number generator
  std::random_device rd; 
  std::mt19937 gen(rd()); // Mersenne Twister engine seeded with rd()

  // Define a distribution that generates integers between 0 and 3 inclusive
  std::uniform_int_distribution<> dis(min, max);

  // Generate a random number
  int random_number = dis(gen);

  return random_number;
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function to calculate the minimum mass difference of b-jet pairs
// Input: 1) vec_bjets: a vector that contains all of the bjets momenta
// Output: 1) the minimum difference of the masses of the pairs (in the case of exactly three bjets pair2 is simply the mass of the unpaired bjet)
//         2) the mass of the unpaired 3rd bjet in the case of 3 bjets
double getDeltaMmin(const std::vector<TLorentzVector>& vec_bjets, double& unpairedMass) 
{
  double DeltaMmin = 1e9; // Initialize to a large value
  double DeltaM;

  if (vec_bjets.size() == 3) {
    // Case for exactly 3 b-jets
    TLorentzVector b1 = vec_bjets[0];
    TLorentzVector b2 = vec_bjets[1];
    TLorentzVector b3 = vec_bjets[2];
        
    // Calculate all possible pairs and track minimum DeltaM
    double m_pair1 = (b1 + b2).M();
    double m_pair2 = b3.M();
    DeltaM = fabs(m_pair1 - m_pair2);
    if (DeltaM < DeltaMmin) {
      DeltaMmin = DeltaM;
      unpairedMass = b3.M(); // Track the unpaired b-jet mass
    }

    m_pair1 = (b1 + b3).M();
    m_pair2 = b2.M();
    DeltaM = fabs(m_pair1 - m_pair2);
    if (DeltaM < DeltaMmin) {
      DeltaMmin = DeltaM;
      unpairedMass = b2.M(); // Track the unpaired b-jet mass
    }

    m_pair1 = (b2 + b3).M();
    m_pair2 = b1.M();
    DeltaM = fabs(m_pair1 - m_pair2);
    if (DeltaM < DeltaMmin) {
      DeltaMmin = DeltaM;
      unpairedMass = b1.M(); // Track the unpaired b-jet mass
    }
  } 
  else if (vec_bjets.size() >= 4) {
    // Case for 4 or more b-jets: no unpaired jets
    unpairedMass = -1; // Sentinel value to indicate no unpaired jet

    for (size_t i = 0; i < 4; i++) {
      for (size_t j = i + 1; j < 4; j++) {
	for (size_t k = 0; k < 4; k++) {
	  for (size_t l = k + 1; l < 4; l++) {
	    if (i != k && i != l && j != k && j != l) {
	      TLorentzVector b1 = vec_bjets[i];
	      TLorentzVector b2 = vec_bjets[j];
	      TLorentzVector b3 = vec_bjets[k];
	      TLorentzVector b4 = vec_bjets[l];
                            
	      double m_pair1 = (b1 + b2).M();
	      double m_pair2 = (b3 + b4).M();
	      DeltaM = fabs(m_pair1 - m_pair2);

	      // Track minimum DeltaM
	      DeltaMmin = std::min(DeltaMmin, DeltaM);
	    }
	  }
	}
      }
    }
  }
  return DeltaMmin;
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function to calculate the minimum mass difference of b-jet pairs
// Input: 1) vec_bjets: a vector that contains all of the bjets momenta as configured with the addition of an extra untagged jet if it exists and if the total number of bjets is less than four
// Output: 1) the minimum difference of the masses of the pairs (in the case of exactly three bjets all the possible combinations (including common jets in both pairs) are examined)
//         2) the mass of the unpaired 3rd bjet in the case of 3 bjets
double getDeltaMminNew(const std::vector<TLorentzVector>& vec_bjets) 
{
  double DeltaMmin = 1e9; // Initialize to a large value
  double DeltaM;

  if (vec_bjets.size() == 3) {
    // Case for exactly 3 b-jets
    TLorentzVector b1 = vec_bjets[0];
    TLorentzVector b2 = vec_bjets[1];
    TLorentzVector b3 = vec_bjets[2];
        
    // Calculate all possible pairs and track minimum DeltaM
    double m_pair1 = (b1 + b2).M();
    double m_pair2 = (b1 + b3).M();
    double m_pair3 = (b2 + b3).M();

    double DM1 = abs(m_pair1 - m_pair2);
    double DM2 = abs(m_pair2 - m_pair3);
    double DM3 = abs(m_pair3 - m_pair2);

    DeltaMmin = std::min(DM1, std::min(DM2, DM3));
  } 
  else if (vec_bjets.size() >= 4) {
    // Case for 4 or more b-jets: no unpaired jets

    for (size_t i = 0; i < 4; i++) {
      for (size_t j = i + 1; j < 4; j++) {
	for (size_t k = 0; k < 4; k++) {
	  for (size_t l = k + 1; l < 4; l++) {
	    if (i != k && i != l && j != k && j != l) {
	      TLorentzVector b1 = vec_bjets[i];
	      TLorentzVector b2 = vec_bjets[j];
	      TLorentzVector b3 = vec_bjets[k];
	      TLorentzVector b4 = vec_bjets[l];
                            
	      double m_pair1 = (b1 + b2).M();
	      double m_pair2 = (b3 + b4).M();
	      DeltaM = fabs(m_pair1 - m_pair2);

	      // Track minimum DeltaM
	      DeltaMmin = std::min(DeltaMmin, DeltaM);
	    }
	  }
	}
      }
    }
  }
  return DeltaMmin;
}
//#######################################################################################################################################//
//#######################################################################################################################################//
// Function to calculate the Mbbj variable.
// Inputs: 1) vec_nonqjets: vector that contains the btagged jets and the untagged jets (all the non outgoing jets) sorted in their btag value
//         2) iEvt: integer for the event category
// Output: Mbbj

double getMassOfbbj(const std::vector<jets_struct>& vec_nonqjets, int iEvt)
{
  double deltaRmin = 1e9;
  double Mbbj = -1;   // Initialize with an invalid mass value
  TLorentzVector bb;  // 4-Vector to store bb pair
  std::vector<jets_struct> unused_bjets; //Vector to store unused b pair

  //Case 1: 2b's or 3b's and no other non q's : simply sum the first three jets
  if(iEvt>=1 && iEvt<=4) Mbbj = (vec_nonqjets[0].momentum + vec_nonqjets[1].momentum + vec_nonqjets[2].momentum).M();
  //Case 2: 3b's and one or more non q's : find the closest pair among the three first and then add the fourth
  else if(iEvt<=7){
    for (size_t i = 0; i < 3; ++i)
      {
	for (size_t j = i + 1; j < 3; ++j)
	  {
	    TLorentzVector b1 = vec_nonqjets[i].momentum;
	    TLorentzVector b2 = vec_nonqjets[j].momentum;
	    double deltaR = getDeltaR(b1, b2);
	    if (deltaR < deltaRmin)
	      {
		deltaRmin = deltaR;
		bb = b1 + b2;
	      }
	  }
      }
    Mbbj = (bb + vec_nonqjets[3].momentum).M();
  }
  //Case 3: 4 b's : pick the closest pair and then the least b tagged among the other two
  else{
    for (size_t i = 0; i < 4; ++i)
      {
	for (size_t j = i + 1; j < 4; ++j)
	  {
	    TLorentzVector b1 = vec_nonqjets[i].momentum;
	    TLorentzVector b2 = vec_nonqjets[j].momentum;
	    double deltaR = getDeltaR(b1, b2);
	    if (deltaR < deltaRmin)
	      {
		deltaRmin = deltaR;
		bb = b1 + b2;
		unused_bjets.clear();
		for (size_t k = 0; k < 4; ++k)
		  {
		    if (k != i && k != j) unused_bjets.push_back(vec_nonqjets[k]);              
		  }
	      }
	  }
      }
    TLorentzVector unused_bjet = unused_bjets[0].btag < unused_bjets[1].btag ? unused_bjets[0].momentum : unused_bjets[1].momentum;
    Mbbj = (bb+unused_bjet).M();
  }
  return Mbbj;
}
//#######################################################################################################################################//
//#######################################################################################################################################//
int EvtCatIndex(const std::vector<jets_struct>& vec_bjets,
                const std::vector<jets_struct>& vec_nonqjets,
                const std::vector<jets_struct>& vec_qjets)
{
    // Count the number of jets
    int n_bjets = vec_bjets.size();
    int n_nonqjets = vec_nonqjets.size();
    int n_qjets = vec_qjets.size();

    // Minimum requirements: at least 2 q-jets, 2 b-jets, and 3 non-q jets.
    if (n_qjets < 2 || n_bjets < 2 || n_nonqjets < 3)
      return 0;

    // Determine the offset based on the number of non-q jets.
    int offset;
    if (n_nonqjets == 3) offset = 1;
    else if (n_nonqjets == 4) offset = 2;
    else offset = 3;
       
    // Cap the number of b-jets at 5.
    int effective_bjets = (n_bjets > 5) ? 5 : n_bjets;

    // Compute the category index.
    int iEvt = (effective_bjets - 2) * 3 + offset;

    return iEvt;
}
