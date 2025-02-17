#define YEAR_2017
#include <iostream>
#include <map>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakePyBind11ParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/MVAHandler_wh.h"
#include "UserCode/bsmhiggs_fwk/interface/TMVAReader_wh.h"
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

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

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
#include "TMath.h"
#include "TLatex.h"

#include <limits>
#include <cmath>
#include <utility>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <string>
#include <unistd.h>
using namespace std;

//######################### FUNCTIONS ###########################//
float getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2);
float getDeltaRmin(TLorentzVector vec_1, vector<TLorentzVector> vec_2);
float getAverageDeltaR(vector<TLorentzVector> vec);
float getdeltaMass(TLorentzVector b1, TLorentzVector b2, float targetMass);
float getDeltaMassMin(vector<TLorentzVector> vec, int n);
float getDeltaMassMinNew(vector<TLorentzVector> vec, float& unpairedmass);
float fjet_matched(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2);
float getmbbb(vector<TLorentzVector> vec);
float computeMT2(TLorentzVector vis1, TLorentzVector vis2, TLorentzVector invis);
float getDeltaPhiMin(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2);
float getMT(TLorentzVector vec_1, TLorentzVector vec_2);
float HTjets(std::vector<TLorentzVector> vec);
TLorentzVector getbbwithDRMin(vector<TLorentzVector> vec);

std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec);
bool sortlep(std::pair<TLorentzVector, TString>  vec_id_i, std::pair<TLorentzVector, TString> vec_id_j)
{
  return  vec_id_i.first.Pt() > vec_id_j.first.Pt();
}
std::pair<int, int> findFirstPair(std::vector<TLorentzVector> vec_b, float targetMass);

bool sortPair(std::pair<TLorentzVector, float> vec_1, std::pair<TLorentzVector, float> vec_2);

struct PairResult {
    TLorentzVector pair1;
    TLorentzVector pair2;
};
PairResult getPair(vector<TLorentzVector> vec);




int main(int argc, char* argv[])
{
  //##################################################################################//
  //##########################    GLOBAL INITIALIZATION     ##########################//
  //##################################################################################//
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

  bool isMC       = runProcess.getParameter<bool>("isMC");
  double xsec     = runProcess.getParameter<double>("xsec");
  double nevts    = runProcess.getParameter<double>("nevts");
  int mctruthmode = runProcess.getParameter<int>("mctruthmode");

  TString proc   = runProcess.getParameter<std::string>("proc");
  TString dtag   = runProcess.getParameter<std::string>("tag");
  TString suffix = runProcess.getParameter<std::string>("suffix");

  bool is2017data  = (!isMC && dtag.Contains("2017"));
  bool is2017MC    = (isMC && dtag.Contains("2017"));
  bool is2018data  = (!isMC && dtag.Contains("2018"));
  bool is2018MC    = (isMC && dtag.Contains("2018"));
  bool is2017_2018 = (is2017MC || is2017data || is2018MC || is2018data);

  bool verbose = runProcess.getParameter<bool>("verbose");

  TString url = runProcess.getParameter<std::string>("input");
  TString outFileUrl(dtag); 
  gSystem->BaseName( url );

  bool isMC_Wh   = isMC && (string(url.Data()).find("Wh_amass")  != string::npos); 
  bool isMC_Zh   = isMC && (string(url.Data()).find("Zh_amass")  != string::npos);
  bool isMC_VBFh = isMC && (string(url.Data()).find("VBFh_amass") !=string::npos);
  bool isSignal  = (isMC_Wh || isMC_Zh || isMC_VBFh);
  
  bool isMC_WJets       = isMC && (string(url.Data()).find("MC13TeV_WJets") != string::npos);
  bool isMC_WJets_HTbin = isMC_WJets && dtag.Contains("HT") ;
  bool isMC_ttbar       = isMC && (string(url.Data()).find("TTTo")  != string::npos);
  if(is2017MC || is2018MC) isMC_ttbar = isMC && (string(url.Data()).find("TTTo")  != string::npos);

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

  // Calculate xsec weight
  float xsec_weight = isMC ? xsec / nevts : 1;
  float weight      = 1.0;
  if(verbose)
    {
      cout << "Cross section   : " << xsec   << " [pb]" << endl;
      cout << "Weight per event: " << xsec_weight       << endl;
      cout << "Weight          : " << weight            << endl;
    }

  // B-tag Working Points
  float looseWP  = 0.0532;
  float mediumWP = 0.3040;
  float tightWP  = 0.7476;
  
 
  //##################################################################################//
  //##########################    INITIATING HISTOGRAMS     ##########################//
  //##################################################################################//
  SmartSelectionMonitor mon;

  // Event flow 
  TH1F *h = (TH1F*) mon.addHistogram ( new TH1F ("eventflow", ";;N_{Events}", 5, 0, 5 ) );
  h->GetXaxis()->SetBinLabel(1,"Raw");
  h->GetXaxis()->SetBinLabel(2,"1 lepton");
  h->GetXaxis()->SetBinLabel(3,"MET>25 & MTW>50");
  h->GetXaxis()->SetBinLabel(4,">=3 jets");
  h->GetXaxis()->SetBinLabel(5,">=3 b-tag jets");            //SR1
  //h->GetXaxis()->SetBinLabel(6,">=1 b-tag and >=1 f-jet");   //SR2

  // Multiplicity 
  mon.addHistogram ( new TH1F ("multi", ";N_{lepton} ;N_{Events}", 5, 0, 5 ) );
  mon.addHistogram ( new TH1F ("qgmulti", ";N_{qg} ;N_{Events}", 15, 0, 15 ) );
  mon.addHistogram ( new TH1F ("jetmulti", ";N_{jet} ;N_{Events}", 15, 0, 15 ) );

  // Particles/Objects kinematics
  mon.addHistogram ( new TH1F ("pt", ";p_{T} [GeV] ;N_{Events}", 100, 0, 500 ) );
  mon.addHistogram ( new TH1F ("eta", ";#eta ;N_{Events}", 100, -6, 6 ) );
  mon.addHistogram ( new TH1F ("phi", ";#phi ;N_{Events}", 100, -TMath::Pi(), TMath::Pi() ) );
  mon.addHistogram ( new TH1F ("M", ";M [GeV] ;N_{Events}", 100, 0, 800) );
  mon.addHistogram ( new TH2F("fjetpt_vs_bbpt", ";p_{T} ;p_{T}", 50, 0, 500, 50, 0, 500 ) );

  // quarks kinematics
  mon.addHistogram ( new TH1F ("max_pt", ";p_{T} [GeV] ;N_{Events}", 100, 0, 200 ) );
  mon.addHistogram ( new TH1F ("min_pt", ";p_{T} [GeV] ;N_{Events}", 100, 0, 200 ) );

  // Delta R
  mon.addHistogram ( new TH1F( "DR", ";#DeltaR ;N_{events}", 100, 0, 8 ) );

  // Delta phi
  mon.addHistogram ( new TH1F( "dphi", ";|#Delta#phi| ;N_{Events}", 100, 0, TMath::Pi() ) );

  // Delta eta
  mon.addHistogram ( new TH1F( "dEta", ";|#Delta#eta| ;N_{Events}", 100, 0, 10 ) );

  // HT
  mon.addHistogram ( new TH1F( "HTWjets", ";H_{T} ;N_{Events}", 100, 0, 5000 ) );

  // MT W
  mon.addHistogram ( new TH1F( "MT", ";M_{T} ;N_{Events}", 100, 0, 500) );

  // <DR> vs pt 
  mon.addHistogram ( new TProfile ( "DR_vs_pt", ";p_{T} [GeV] ;<#DeltaR>", 100, 0, 500 ) );
  mon.addHistogram ( new TH2F ( "DRbb_vs_ptbb", ";p_{T}^{bb} [GeV]; ;#DeltaR", 100, 0, 500, 100, 0, 8 ));

  // b-tag discriminator
  mon.addHistogram ( new TH1F( "btag", ";b-tag ;N_{Events}", 100, 0, 1 ) );
  
  // Filter efficiency for ttbar
  mon.addHistogram ( new TH1F( "filter", ";filter ;" , 1, 0, 1 ) );

  // Fat jets and soft-drop mass
  mon.addHistogram ( new TH1F ( "SDM", ";M_{SD} [GeV] ;N_{Events}", 100, 0, 200 ) );
  mon.addHistogram ( new TH2F ( "fjet_pt_vs_SDM", ";p_{T} [GeV] ;M_{SD} [Gev]", 100, 0, 500, 100, 0, 200 ) );
  mon.addHistogram ( new TH2F ( "fjet_eta_vs_SDM", ";#eta ;M_{SD} [Gev]", 100, -7, 7, 100, 0, 200 ) );
  mon.addHistogram ( new TH2F ( "fjet_subjet_vs_SDM", ";N_{subjet} ;M_{SD} [Gev]", 5, 0, 5, 100, 0, 200 ) );
  mon.addHistogram ( new TProfile ( "fjetpt_vs_SDM", ";p_{T} [GeV] ;<M_{SD}> [Gev]", 100, 0, 500 ) );
  mon.addHistogram ( new TProfile ( "fjeteta_vs_SDM", ";#eta ;<M_{SD}> [Gev]", 100, -7, 7 ) );
  mon.addHistogram ( new TProfile ( "fjetsubjet_vs_SDM", ";N_{subjet} ;<M_{SD}> [Gev]", 5, 0, 5 ) );

  mon.addHistogram( new TH2F ( "pt_bb_vs_met", ";p_{T} [GeV] ;E_{T}^{miss} [GeV]", 100, 0, 500, 100, 0, 500 ) );
  mon.addHistogram( new TProfile ( "pt_fjet_vs_met_prof", ";p_{T} [GeV] ;<E_{T}^{miss}> [GeV]", 100, 0, 500 ) );

  // Subjet count
  mon.addHistogram (new TH1F ( "count", ";N_{subjets} ;N_{Events}", 5, 0, 5 ) );

  // Discriminants
  mon.addHistogram ( new TH1F( "discriminant", ";discriminant ;N_{Events}", 100, 0, 1 ) );
  mon.addHistogram ( new TH2F( "fjet_Xbb_pt", ";p_{T} [GeV] ;Xbb", 100, 0, 500, 50, 0, 1 ) );
  mon.addHistogram ( new TH2F( "fjet_XbbXccXqq_pt", ";p_{T} [GeV] ;XbbXccXqq", 100, 0, 500, 50, 0, 1 ) );
  mon.addHistogram ( new TProfile( "Xbb_pt", ";p_{T} [GeV] ;<Xbb>", 100, 0, 500 ) );
  mon.addHistogram ( new TProfile( "XbbXccXqq_pt", ";p_{T} [GeV] ;<XbbXccXqq>", 100, 0, 500 ) );

  // Final step histograms for /submit 3.1
  mon.addHistogram ( new TH1F ( "BDT", ";BDT ;Events", 50, -1, 1 ) );
  mon.addHistogram ( new TH1F ( "H_pt", ";p_{T}^{H} [GeV] ;Events", 50, 0, 800 ) );
  mon.addHistogram ( new TH1F ( "W_pt", ";p_{T}^{W} [GeV] ;Events", 50, 0, 800 ) );
  mon.addHistogram ( new TH1F ( "l_pt", ";p_{T}^{lep} [GeV] ;Events", 50, 0, 500 ) );
  mon.addHistogram ( new TH1F ( "bjet1_pt", ";p_{T}^{bjet1} [GeV] ;Events", 50, 0, 800 ) );
  mon.addHistogram ( new TH1F ( "HT", ";H_{T} ;Events", 50, 0, 1000 ) );
  mon.addHistogram ( new TH1F ( "H_M", ";M_{H} [GeV] ;Events", 50, 0, 800) );
  mon.addHistogram ( new TH1F ( "bbj_M", ";M_{bbj} [GeV ;Events", 50, 0, 1000) );
  mon.addHistogram ( new TH1F ( "W_MT", ";M_{T}^{W} [GeV] ;Events", 50, 0, 500) );
  mon.addHistogram ( new TH1F ( "MET", ";E_{T}^{miss} [GeV] ;Events", 50, 0, 500 ) );
  mon.addHistogram ( new TH1F ( "WH_dphi", ";|#Delta#phi(WH)| ;Events", 50, 0, TMath::Pi() ) );
  mon.addHistogram ( new TH1F ( "METjet_dphi", ";min|#Delta#phi(E_{T}^{miss}-jet)| ;Events", 50, 0, TMath::Pi() ) );
  mon.addHistogram ( new TH1F ( "WH_deta", ";|#Delta#eta(WH)| ;Events", 50, 0, 7 ) );
  mon.addHistogram ( new TH1F ( "avbb_dR", ";#DeltaR(bb)^{aver} ;Events", 50, 0, 7 ) );
  mon.addHistogram ( new TH1F ( "jet-mn_dR", ";#DeltaR(jet-#mu) ;Events", 100, 0, 7 ) );
  mon.addHistogram ( new TH1F ( "jet-mn_cc_dR", ";#DeltaR(jet-#mu)_{cc} ;Events", 100, 0, 7 ) );
  mon.addHistogram ( new TH1F ( "jet-en_dR", ";#DeltaR(jet-e) ;Events", 100, 0, 7 ) );
  mon.addHistogram ( new TH1F ( "jet-en_cc_dR", ";#DeltaR(jet-#mu)_{cc} ;Events", 100, 0, 7 ) );
  mon.addHistogram ( new TH1F ( "minbb_dM", ";#DeltaM(bb)^{min} [GeV] ;Events", 50, 0, 400) );
  mon.addHistogram ( new TH1F ( "btag1_bef", ";b-tag1 ;Events", 100, mediumWP, 1 ) );
  mon.addHistogram ( new TH1F ( "btag2_bef", ";b-tag2 ;Events", 100, mediumWP, 1 ) );
  mon.addHistogram ( new TH1F ( "btag3_bef", ";b-tag3 ;Events", 100, mediumWP, 1 ) );
  mon.addHistogram ( new TH1F ( "btag1", ";b-tag1 ;Events", 50, tightWP, 1 ) );
  mon.addHistogram ( new TH1F ( "btag2", ";b-tag2 ;Events", 50, mediumWP, 1 ) );
  mon.addHistogram ( new TH1F ( "btag3", ";b-tag3 ;Events", 50, mediumWP, 1 ) );
  
  
  // Tree info
  int evStart     = runProcess.getParameter<int>("evStart");
  int evEnd       = runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  //##################################################################################//
  //###############         GET READY FOR THE EVENT LOOP           ###################//
  //##################################################################################//
  
  // Open the file and get events tree
  DataEvtSummaryHandler summaryHandler_;

  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0) {return -1; printf("file is 0");}

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


  //------Event Counters-------//
  int n_event_lepton = 0;
  int n_event_jet    = 0;
  int n_event_bjet   = 0;
  int n_event_metcut = 0;
  int n_event_fjet   = 0;

  int n_SR1 = 0;
  int n_SR2 = 0;
  
  // Variables for TTJets 
  int nTot   = 0; 
  int nFilt1 = 0;
  int nFilt4 = 0;
  int nFilt5 = 0;

 
  //####################################################################################################################//
  //###########################################           TMVAReader         ###########################################//
  //####################################################################################################################//
  
  std::string SR1_chpath = "WhAnalysis/SR1/";
  std::string SR2_chpath = "WhAnalysis/SR2/";

  TMVAReader SR1adaReader;
  SR1adaReader.InitTMVAReader();
  std::string SR1ada_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+SR1_chpath+"TMVA_BDT_ada.weights.xml";
  SR1adaReader.SetupMVAReader("SR1adaClass", SR1ada_xml_path);

  TMVAReader SR1gradReader;
  SR1gradReader.InitTMVAReader();
  std::string SR1grad_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+SR1_chpath+"TMVA_BDT_grad.weights.xml";
  SR1gradReader.SetupMVAReader("SR1gradClass", SR1grad_xml_path);

  TMVAReader SR2adaReader;
  SR2adaReader.InitTMVAReader();
  std::string SR2ada_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+SR2_chpath+"TMVA_BDT_ada.weights.xml";
  SR2adaReader.SetupMVAReader("SR2adaClass", SR2ada_xml_path);

  TMVAReader SR2gradReader;
  SR2gradReader.InitTMVAReader();
  std::string SR2grad_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+SR2_chpath+"TMVA_BDT_grad.weights.xml";
  SR2gradReader.SetupMVAReader("SR2gradClass", SR2grad_xml_path);
  
 
  //####################################################################################################################//
  //###########################################           EVENT LOOP         ###########################################//
  //####################################################################################################################//
  // Loop on all the events
  int treeStep = (evEnd-evStart)/50; 
  if(treeStep==0)treeStep=1;
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);

  for(int iev=evStart; iev<evEnd; iev++)
    {
      // Initializa weight
      weight = 1.0;

      // Step 1 : xsec weighting
      weight *= xsec_weight;
      
      // Load the event content from tree
      summaryHandler_.getEntry(iev);
      DataEvtSummary_t &ev=summaryHandler_.getEvent();

      mon.fillHisto("eventflow", "histo", 0, 1);
      mon.fillHisto("eventflow", "histo_weighted", 0, weight);

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 0, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 0, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 0, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 0, weight);
	}
            
      if(!isMC && duplicatesChecker.isDuplicate(ev.run, ev.lumi, ev.event))
	{
	  nDuplicates++;
	  cout << "nDuplicates: " << nDuplicates << endl;
	  continue;
	}

      if (isMC_WJets && !(isMC_WJets_HTbin) ) { if(ev.lheHt >= 70) continue; }
	

     
      //############################################################//
      //################## GENERATOR LEVEL #########################//
      //###########################################################//
      
      // Create new object vectors
      std::vector<TLorentzVector> vec_genlep;  // store the lepton
      std::vector<TLorentzVector> vec_genqg;   // store the quarks/gluons
      std::vector<TLorentzVector> vec_genb;    // store the b quarks
      std::vector<TLorentzVector> vec_H;       // store the Higgs
      std::vector<TLorentzVector> vec_A;       // store the A
      std::vector<TLorentzVector> vec_W;       // store the W
      std::vector<TLorentzVector> vec_bb1;     // store the bb pair from mother A1
      std::vector<TLorentzVector> vec_bb2;     // store the bb pair from mother A2

      TLorentzVector met_gen;                 // Missing transverse energy vector

      // LOOP OVER MC PARTICLES at GENERATOR LEVEL
      for (int imc=0; imc<ev.nmcparticles; imc++)
	{
	  if(verbose && iev<5) 
	    {
	      std::cout << " imcparticle " << imc << " : is a " << ev.mc_id[imc] <<",has a mom:"<<ev.mc_mom[imc]<< " , and has a mother at: " << ev.mc_momidx[imc]  <<"has status   "<<ev.mc_status[imc] << "  and has a 4-vector p = (" << ev.mc_en[imc] << ", " << ev.mc_px[imc] << ", " << ev.mc_py[imc] << ", " << ev.mc_pz[imc] << " ) " << std::endl;
	    }

	  if(ev.mc_id[imc]==25) 
	    { // Found the Higgs boson
	      TLorentzVector pH;
	      pH.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Make histograms of Higgs boson kinematics:
	      mon.fillHisto("pt", "gen_H", pH.Pt(), weight);
	      mon.fillHisto("eta", "gen_H", pH.Eta(), weight);
	      mon.fillHisto("phi", "gen_H", pH.Phi(), weight);
	      mon.fillHisto("M", "gen_H", pH.M(), weight);
	      
	      vec_H.push_back(pH);        
	    }

	  if(abs(ev.mc_id[imc])==24)
	    { // Found the W
	      TLorentzVector pW;
	      pW.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Make histograms of the W boson kinematics
	      mon.fillHisto("pt", "gen_W", pW.Pt(), weight);
	      mon.fillHisto("eta", "gen_W", pW.Eta(), weight);
	      mon.fillHisto("phi", "gen_W", pW.Phi(), weight);
	      mon.fillHisto("M", "gen_W", pW.M(), weight);
	      
	      vec_W.push_back(pW);
	    }

	  if(ev.mc_id[imc]==36)
            { // Found the A                                                              
              TLorentzVector pA;
              pA.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

              // Make histograms of the A boson kinematics
	      mon.fillHisto("pt", "gen_A", pA.Pt(), weight);
	      mon.fillHisto("eta", "gen_A", pA.Eta(), weight);
	      mon.fillHisto("phi", "gen_A", pA.Phi(), weight);
	      mon.fillHisto("M", "gen_A", pA.M(), weight);
	      
	      vec_A.push_back(pA);
            }

	  if(abs(ev.mc_id[imc])==5 && ev.mc_mom[imc]==36)
	    { // Found the b quarks
	      TLorentzVector pbb;
	      pbb.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Apply acceptance cuts
	      if(fabs(pbb.Eta()) > 2.4 || pbb.Pt() < 20) continue;

	      // b-quarks from mother1
	      if(ev.mc_momidx[imc] == 4) vec_bb1.push_back(pbb);

	      // b-quarks from mother2
	      else if(ev.mc_momidx[imc] == 5) vec_bb2.push_back(pbb);

	      vec_genb.push_back(pbb);
	    }
	  	
      
	  if(ev.mc_momidx[imc] == 0 && abs(ev.mc_mom[imc]) == 24)
	    {
	      TLorentzVector plep;
	      plep.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
		  
	      if(abs(ev.mc_id[imc]) == 11 || abs(ev.mc_id[imc]) == 13 || abs(ev.mc_id[imc]) == 15)
		{ // Found the lepton
		 
		  // Apply acceptance cuts
		  if(fabs(plep.Eta()) > 2.4 || plep.Pt() < 20) continue;

		  vec_genlep.push_back(plep);
		}

	      else if(abs(ev.mc_id[imc]) == 12 || abs(ev.mc_id[imc]) == 14 || abs(ev.mc_id[imc]) == 16)
		{ // Found the lepton neutrino
		  met_gen += plep;
		}  
	    }


	  if(ev.mc_status[imc] == 23 || ev.mc_status[imc] == 24)
	    {
	      if((abs(ev.mc_id[imc]) > 0 && abs(ev.mc_id[imc]) < 6) || ev.mc_id[imc] == 21)
		{ // Final state particles (quarks or gluons)
		  TLorentzVector pqg;
		  pqg.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
		  
		  // Apply acceptance cuts
		  if(pqg.Pt() < 20.0 || fabs(pqg.Eta()) > 2.4) continue;

		  vec_genqg.push_back(pqg);
		}
	    }
	} // END MC PARTICLE LOOP 

      

      // Fill the multiplicity histograms before cuts
      mon.fillHisto("multi", "step0_gen_lepton", vec_genlep.size(), weight);
      mon.fillHisto("qgmulti", "step0_gen", vec_genqg.size(), weight);

      // Sorting the quarks/gluons and b quarks based on pT in descending order
      vec_genqg = sort_vec_pt(vec_genqg);
      vec_genb  = sort_vec_pt(vec_genb);
      vec_bb1   = sort_vec_pt(vec_bb1);
      vec_bb2   = sort_vec_pt(vec_bb2);

      // Fill histogram with lheHT for the W+Jets background
      mon.fillHisto("HTWJets", "gen", ev.lheHt, weight); 
	
      if(vec_bb1.size()>1)
	{
	  float DRbb1 = getDeltaR(vec_bb1[0], vec_bb1[1]);
	  TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];

	  // Fill histograms with bb-pair kinematics
	  mon.fillHisto("pt", "step0_gen_bb1", pbb1.Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_bb1", pbb1.Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_bb1", pbb1.Phi(), weight);
	  mon.fillHisto("M", "step0_gen_bb1", pbb1.M(), weight);

	  // Fill histograms of angular distance DR
	  mon.fillHisto("DR", "step0_gen_bb1", DRbb1, weight);
	  mon.fillProfile("DR_vs_pt", "step0_gen_bb1", pbb1.Pt(), DRbb1, weight);

	  if(DRbb1 < 0.8)
	    {
	      mon.fillHisto("DRbb_vs_ptbb", "1_DR<0.8", (vec_bb1[0]+vec_bb1[1]).Pt(), DRbb1, weight);
	      mon.fillProfile("DR_vs_pt", "1_DR<0.8", (vec_bb1[0]+vec_bb1[1]).Pt(), DRbb1, weight);
	      mon.fillHisto("pt", "1_DR<0.8", (vec_bb1[0]+vec_bb1[1]).Pt(), weight);
	    }	  
	}
      
      if(vec_bb2.size()>1)
	{
	  float DRbb2 = getDeltaR(vec_bb2[0], vec_bb2[1]);
	  TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];

	  // Fill histograms with bb-pair kinematics
	  mon.fillHisto("pt", "step0_gen_bb2", pbb2.Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_bb2", pbb2.Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_bb2", pbb2.Phi(), weight);
	  mon.fillHisto("M", "step0_gen_bb2", pbb2.M(), weight);

	   // Fill histograms of angular distance DR
	  mon.fillHisto("DR", "step0_gen_bb2", DRbb2, weight);
	  mon.fillProfile("DR_vs_pt", "step0_gen_bb2", pbb2.Pt(), DRbb2, weight);

	  if(DRbb2 < 0.8)
	    {
	      mon.fillHisto("DRbb_vs_ptbb", "2_DR<0.8", (vec_bb2[0]+vec_bb2[1]).Pt(), DRbb2, weight);
	      mon.fillProfile("DR_vs_pt", "2_DR<0.8", (vec_bb2[0]+vec_bb2[1]).Pt(), DRbb2, weight);
	      mon.fillHisto("pt", "2_DR<0.8", (vec_bb2[0]+vec_bb2[1]).Pt(), weight);
	    }
	}

       if(vec_bb1.size()>1 && vec_bb2.size()>1)
	{
	  TLorentzVector vec_A1 = vec_bb1[0] + vec_bb1[1];
	  TLorentzVector vec_A2 = vec_bb2[0] + vec_bb2[1];
       
	  float DRAA = getDeltaR(vec_A1, vec_A2);
	  mon.fillHisto("DR", "step0_gen_AA", DRAA, weight);
	}

      // Fill the kinematics histograms for the quarks/gluons
      if(vec_genqg.size()>0)
	{
	  mon.fillHisto("pt", "step0_gen_qg1", vec_genqg[0].Pt(), weight);
	  mon.fillHisto("max_pt", "step0_gen_qg1", vec_genqg[0].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_qg1", vec_genqg[0].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_qg1", vec_genqg[0].Phi(), weight);
	}

      if(vec_genqg.size()>1)
	{
	  mon.fillHisto("pt", "step0_gen_qg2", vec_genqg[1].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_qg2", vec_genqg[1].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_qg2", vec_genqg[1].Phi(), weight);
	}

      if(vec_genqg.size()>2)
	{
	  mon.fillHisto("pt", "step0_gen_qg3", vec_genqg[2].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_qg3", vec_genqg[2].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_qg3", vec_genqg[2].Phi(), weight);
	}

      if(vec_genqg.size()>3)
	{
	  mon.fillHisto("pt", "step0_gen_qg4", vec_genqg[3].Pt(), weight);
	  mon.fillHisto("min_pt", "step0_gen_qg4", vec_genqg[3].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_qg4", vec_genqg[3].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_qg4", vec_genqg[3].Phi(), weight);
	}

      
      // Fill the kinematics histograms for the b-quarks
      if(vec_genb.size()>0)
	{
	  mon.fillHisto("pt", "step0_gen_b1", vec_genb[0].Pt(), weight);
	  mon.fillHisto("max_pt", "step0_gen_b1", vec_genb[0].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_b1", vec_genb[0].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_b1", vec_genb[0].Phi(), weight);
	}

      if(vec_genb.size()>1)
	{
	  mon.fillHisto("pt", "step0_gen_b2", vec_genb[1].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_b2", vec_genb[1].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_b2", vec_genb[1].Phi(), weight);
	}

      if(vec_genb.size()>2)
	{
	  mon.fillHisto("pt", "step0_gen_b3", vec_genb[2].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_b3", vec_genb[2].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_b3", vec_genb[2].Phi(), weight);
	}

      if(vec_genb.size()>3)
	{
	  mon.fillHisto("pt", "step0_gen_b4", vec_genb[3].Pt(), weight);
	  mon.fillHisto("min_pt", "step0_gen_b4", vec_genb[3].Pt(), weight);
	  mon.fillHisto("eta", "step0_gen_b4", vec_genb[3].Eta(), weight);
	  mon.fillHisto("phi", "step0_gen_b4", vec_genb[3].Phi(), weight);
	}

      

      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is one lepton       
      if(vec_genlep.size() >= 1)
	{
	  // Fill the multiplicity histograms
	  mon.fillHisto("multi", "step1_gen_lepton", vec_genlep.size(), weight);
	  mon.fillHisto("qgmulti", "step1_gen", vec_genqg.size(), weight);
			
	  // Fill histograms of angular distance DR
	  if(vec_bb1.size()>1)
	    {
	      float DRbb1 = getDeltaR(vec_bb1[0], vec_bb1[1]);
	      mon.fillHisto("DR", "step1_gen_bb1", DRbb1, weight);
	      
	      TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	      mon.fillProfile("DR_vs_pt", "step1_gen_bb1", pbb1.Pt(), DRbb1, weight);
	    }

	  if(vec_bb2.size()>1)
	    {
	      float DRbb2 = getDeltaR(vec_bb2[0], vec_bb2[1]);
	      mon.fillHisto("DR", "step1_gen_bb2", DRbb2, weight);

	      TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	      mon.fillProfile("DR_vs_pt", "step1_gen_bb2", pbb2.Pt(), DRbb2, weight);
	    }

	  if(vec_bb1.size()>1 && vec_bb2.size()>1)
	    {
	      TLorentzVector vec_A1 = vec_bb1[0] + vec_bb1[1];
	      TLorentzVector vec_A2 = vec_bb2[0] + vec_bb2[1];
       
	      float DRAA = getDeltaR(vec_A1, vec_A2);
	      mon.fillHisto("DR", "step1_gen_AA", DRAA, weight);
	    }
	}

	  
      //### STEP 2 ###//
      // If there is at least one lepton and at least three quarks/gluons                                                                                     
      if(vec_genlep.size()>0 && vec_genqg.size()>2)
	{
	  // Fill the histograms of qj kinematics
	  mon.fillHisto("pt", "step2_gen_qg1", vec_genqg[0].Pt(), weight);
	  mon.fillHisto("eta", "step2_gen_qg1", vec_genqg[0].Eta(), weight);
	  mon.fillHisto("phi", "step2_gen_qg1", vec_genqg[0].Phi(), weight);
	      
	  mon.fillHisto("pt", "step2_gen_qg2", vec_genqg[1].Pt(), weight);
	  mon.fillHisto("eta", "step2_gen_qg2", vec_genqg[1].Eta(), weight);
	  mon.fillHisto("phi", "step2_gen_qg2", vec_genqg[1].Phi(), weight);

	  mon.fillHisto("pt", "step2_gen_qg3", vec_genqg[2].Pt(), weight);
	  mon.fillHisto("eta", "step2_gen_qg3", vec_genqg[2].Eta(), weight);
	  mon.fillHisto("phi", "step2_gen_qg3", vec_genqg[2].Phi(), weight);
	    
	  if(vec_genqg.size() > 3)
	    {
	      mon.fillHisto("pt", "step2_gen_qg4", vec_genqg[3].Pt(), weight);
	      mon.fillHisto("eta", "step2_gen_qg4", vec_genqg[3].Eta(), weight);
	      mon.fillHisto("phi", "step2_gen_qg4", vec_genqg[3].Phi(), weight);
	    }

	  // Fill the lepton kinematics histograms
	  mon.fillHisto("pt", "gen_lepton", vec_genlep[0].Pt(), weight);
	  mon.fillHisto("eta", "gen_lepton", vec_genlep[0].Eta(), weight);
	  mon.fillHisto("phi", "gen_lepton", vec_genlep[0].Phi(), weight);
	  
	  // Fill the histogram of HT q/g
	  float HTqg = HTjets(vec_genqg);
	  mon.fillHisto("HT", "gen", HTqg, weight);
	  
	  // Fill the histogram of MET pT
	  mon.fillHisto("pt", "gen_MET", met_gen.Pt(), weight);

	  // Fill the histogram of W transverse mass
	  float MTW_gen = getMT(vec_genlep[0], met_gen);
	  mon.fillHisto("MT", "gen_recoW", MTW_gen, weight);

	  // Define the hadronic system
	  TLorentzVector phad_gen = vec_genqg[0] + vec_genqg[1] + vec_genqg[2];
	  if(vec_genqg.size()>3) phad_gen += vec_genqg[3];

	  // Fill the kinematics histograms of hadronic system
	  mon.fillHisto("pt", "gen_recoH", phad_gen.Pt(), weight);
	  mon.fillHisto("eta", "gen_recoH", phad_gen.Eta(), weight);
	  mon.fillHisto("phi", "gen_recoH", phad_gen.Phi(), weight);
	  mon.fillHisto("M", "gen_recoH", phad_gen.M(), weight);

	  // Define the leptonic system
	  TLorentzVector plep_gen = vec_genlep[0] + met_gen;
	  mon.fillHisto("pt", "gen_recoW", plep_gen.Pt(), weight);
	  mon.fillHisto("eta", "gen_recoW", plep_gen.Eta(), weight);
	  mon.fillHisto("phi", "gen_recoW", plep_gen.Phi(), weight);

	  // Fill Delta phi WH histogram
	  float dphiWH_gen = fabs(plep_gen.DeltaPhi(phad_gen));
	  mon.fillHisto("dphi", "gen_WH", dphiWH_gen, weight);

	  // Fill Delta phi of MET and lepton
	  float dphiMETlep_gen = fabs(met_gen.DeltaPhi(vec_genlep[0]));
	  mon.fillHisto("dphi", "gen_met-lepton", dphiMETlep_gen, weight);

	  
	} // END OF GENERATOR LEVEL ANALYSIS 



      
	    
      //######################################################################//
      //##########################  DETECTOR LEVEL ##########################//
      //####################################################################//

      // Create new object vectors
      std::vector<TLorentzVector> vec_en;                                     // store the electron
      std::vector<TLorentzVector> vec_mn;                                     // store the muon
      std::vector<TLorentzVector> vec_lepton;                                 // store the leptons
      std::vector<TLorentzVector> vec_jets;                                   // store all the jets based on their b-tag
      std::vector<TLorentzVector> vec_fjets_raw;                              // store all the fat jets with pT>30 GeV and at least 2 subjets
      std::vector<TLorentzVector> vec_fjets;                                  // store all the fat jets with pT>150 GeV
      std::vector<TLorentzVector> vec_bjets;                                  // store the b-tagged jets based on their pt
      std::vector<TLorentzVector> vec_btag;                                   // store the b-tagged jets based on their b-tag
      std::vector<TLorentzVector> vec_untag;                                  // store the untagged jets based on their b-tag
      std::vector<TLorentzVector> vec_jets_cc;                                // store the jets that dont overlap with fat jets
      std::vector<TLorentzVector> vec_bjets_cc;                               // store the b-jets that dont overlap with fat jets
      std::vector<std::pair<TLorentzVector, int>>   vec_lepton_id;            // store the lepton with its id index
      std::vector<std::pair<TLorentzVector, float>> vec_untagbtag;            // store the untagged jets with their b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_jetsbtag;             // store the jets with their b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_bjetsbtag;            // store the b tagged jets with their b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_jetsbtag_cc;          // store the jets with their b-tag index that dont overlap with fat jets
      std::vector<std::pair<TLorentzVector, float>> vec_bjetsbtag_cc;         // store the b-jets with their b-tag index that dont overlap with fat jets
      std::vector<std::pair<TLorentzVector, int>>   vec_fjets_index;          // store the fat jet with its index when pT>150 GeV
      std::vector<std::pair<TLorentzVector, int>>   vec_fjets_index_raw;      // store the fat jet with its index when pT>30 GeV and at least 2 subjets

      
      // Define MET
      TLorentzVector met;
      met.SetPtEtaPhiM(ev.met_pt, 0, ev.met_phi, 0);

      // Define the bb pair from the same mother
      TLorentzVector bb;

      // ============================ //
      // LOOP OVER DETECTOR PARTICLES //
      // ============================ //

      // Lepton configuration
      
      // Muon 
      for(int i=0; i<ev.mn; i++)
	{
	  TLorentzVector pmn;
	  pmn.SetPxPyPzE(ev.mn_px[i], ev.mn_py[i], ev.mn_pz[i], ev.mn_en[i]);

	  // Apply acceptance cuts
	  if(pmn.Pt()<20 || abs(pmn.Eta())>2.5) continue;

	  // Apply ID and Isolation
	  if(ev.mn_passId[i] && ev.mn_passIso[i])
	    {
	      vec_mn.push_back(pmn);
	      vec_lepton_id.push_back(std::make_pair(pmn, ev.mn_id[i]));
	    }
	}
      
      // Electron
      for(int i=0; i<ev.en; i++)
	{
	  TLorentzVector pen;
	  pen.SetPxPyPzE(ev.en_px[i], ev.en_py[i], ev.en_pz[i], ev.en_en[i]);

	  // Apply acceptance cuts
	  if(pen.Pt()<20 || fabs(pen.Eta())>2.4) continue;

	  // Apply ID and Isolation
	  if(ev.en_passId[i] && ev.en_passIso[i])
	    {
	      vec_en.push_back(pen);
	      vec_lepton_id.push_back(std::make_pair(pen, ev.en_id[i]));
	    }
	}

      // Sort the lepton vectors based on their pT
      vec_mn = sort_vec_pt(vec_mn);
      vec_en = sort_vec_pt(vec_en);
      std::sort(vec_lepton_id.begin(), vec_lepton_id.end(), sortlep);
     

      //-------------------------//
      // Jets and cross-cleaning //
      //-------------------------//
      
      int nHFb = 0;
      int nHFc = 0;

      float SR1_btag1_WP = 0.0;
      float SR1_btag2_WP = 0.0;
      float SR1_btag3_WP = 0.0;
      
      for(int i=0; i<ev.jet; i++)
	{
	  bool overlap = false;

	  TLorentzVector pjet;
	  pjet.SetPxPyPzE(ev.jet_px[i], ev.jet_py[i], ev.jet_pz[i], ev.jet_en[i]);

	  // Apply acceptance cuts
	  if(pjet.Pt()<20 || fabs(pjet.Eta())>2.5) continue;

	  // Jets id 
	  if(!ev.jet_PFTight[i]) continue;

	  for(int mn_count=0; mn_count<vec_mn.size(); mn_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_mn[mn_count]);

	      // Fill the histogram of DR(jet-muon) before cross cleaning
	      mon.fillHisto("jet-mn_dR", "step0", DRjet, weight);

	      if(vec_mn.size()>1) mon.fillHisto("jet-mn_dR", "step1", DRjet, weight);
		
	      if(DRjet<0.4) overlap = true;
		
	      else
		{
		  mon.fillHisto("jet-mn_cc_dR", "step0", DRjet, weight);

		  if(vec_mn.size()>1) mon.fillHisto("jet-mn_cc_dR", "step1", DRjet, weight);
		}
	    }
      
	  
	  for(int en_count=0; en_count<vec_en.size(); en_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_en[en_count]);

	      // Fill the histogram of DR(jet-electron) before cross cleaning
	      mon.fillHisto("jet-en_dR", "step0", DRjet, weight);

	      if(vec_en.size()>1) mon.fillHisto("jet-en_dR", "step1", DRjet, weight);

	      if(DRjet<0.4) overlap = true;
	
	      else
		{
		  mon.fillHisto("jet-en_cc_dR", "step0", DRjet, weight);

		  if(vec_en.size()>1) mon.fillHisto("jet-en_cc_dR", "step1", DRjet, weight);
		}
	    }
	  
	  if(overlap) continue;
	    
	  vec_jetsbtag.push_back(std::make_pair(pjet,ev.jet_btag1[i]));
	  std::sort(vec_jetsbtag.begin() , vec_jetsbtag.end(),  sortPair);

	  // b-tag discriminator for 3 most energetic b-jets before WP
	  SR1_btag1_WP = vec_jetsbtag[0].second;
	  SR1_btag2_WP = vec_jetsbtag[1].second;
	  SR1_btag3_WP = vec_jetsbtag[2].second;

	  // Keep the objects with DeepJetMediumWP
	  if(ev.jet_btag1[i] >= mediumWP)
	    { // Found the b-jets and store them in a vector with their corresponding b-tag index
	      vec_bjetsbtag.push_back(std::make_pair(pjet, ev.jet_btag1[i]));
	      std::sort(vec_bjetsbtag.begin(), vec_bjetsbtag.end(), sortPair);
	    }

	  else
	    { // Store the untagged jets with their corresponding b-tag index
	      vec_untagbtag.push_back(std::make_pair(pjet,ev.jet_btag1[i]));
	      std::sort(vec_untagbtag.begin(), vec_untagbtag.end(), sortPair);
	    }
	  
	  if(isMC_ttbar)
	    {
	      bool isMatched = false;

	      for(int imc=0; imc<ev.nmcparticles; imc++)
	  	{
	  	  if(fabs(ev.mc_id[imc])!=5) continue;
		  
	  	  TLorentzVector pgen;
	  	  pgen.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	  	  float DR = getDeltaR(pjet, pgen);
	  	  if(DR<0.4)
	  	    {
	  	      isMatched = true;
	  	      break;
	  	    }
	  	}

	      if(abs(ev.jet_partonFlavour[i])==5)
	  	{
	  	  if(!isMatched) nHFb++;
	  	}
	      else if(abs(ev.jet_partonFlavour[i])==4)
	  	{
	  	  nHFc++;
	  	}	
	    }   
	} //End of jet loop
      
      //Split inclusive TTJets POWHEG sample into tt+bb, tt+cc, tt+light
      if(isMC_ttbar)
	{	  
	  nTot++;

	  if(mctruthmode==5)
	    {
	      if(!(nHFb>0)) {continue;}
	      nFilt5++;
	    }
	  
	  if(mctruthmode==4)
	    {
	      if(!(nHFc>0 && nHFb==0)) {continue;}
	      nFilt4++;
	    }
	  
	  if(mctruthmode==1)
	    {
	      if(nHFb>0 || (nHFc>0 && nHFb == 0)) continue;
	      nFilt1++;
	    }
	}

      for (const auto& pair : vec_jetsbtag)
	{
	  vec_jets.push_back(pair.first);
	}

      for (const auto& pair : vec_bjetsbtag)
	{
	  vec_btag.push_back(pair.first);
	}

      for (const auto& pair : vec_untagbtag)
	{
	  vec_untag.push_back(pair.first);
	}

      //================//
      // F A T  J E T S //
      //================//

      // Lepton-fat jet cross cleaning
      for(int i=0; i<ev.fjet; i++)
	{
	  bool overlap = false;

	  TLorentzVector pfjet;
	  pfjet.SetPxPyPzE(ev.fjet_px[i], ev.fjet_py[i], ev.fjet_pz[i], ev.fjet_en[i]);

	  if(ev.fjet_subjet_count[i]<2)continue;

	  // Apply acceptance cuts
	  if(pfjet.Pt()<30 || fabs(pfjet.Eta())>2.5) continue;
	  
	  for(int mn_count=0; mn_count<vec_mn.size(); mn_count++)
	    {
	      float DRjet = getDeltaR(pfjet, vec_mn[mn_count]);

	      // Fill the histogram of DR(fjet-muon) before cross cleaning
	      mon.fillHisto("DR", "step0_fjet-mn", DRjet, weight);

	      if(vec_mn.size()>1) mon.fillHisto("DR", "step1_fjet-mn", DRjet, weight);
		
	      if(DRjet<0.8) overlap = true;
		
	      else
		{
		  mon.fillHisto("DR", "step0_fjet-mn_cc", DRjet, weight);

		  if(vec_mn.size()>1) mon.fillHisto("DR", "step1_fjet-mn_cc", DRjet, weight);
		}
	    }

	  for(int en_count=0; en_count<vec_en.size(); en_count++)
	    {
	      float DRjet = getDeltaR(pfjet, vec_en[en_count]);

	      // Fill the histogram of DR(fjet-electron) before cross cleaning
	      mon.fillHisto("DR", "step0_fjet-en", DRjet, weight);

	      if(vec_en.size()>1) mon.fillHisto("DR", "step1_fjet-en", DRjet, weight);
		
	      if(DRjet<0.8) overlap = true;
		
	      else
		{
		  mon.fillHisto("DR", "step0_fjet-en_cc", DRjet, weight);

		  if(vec_en.size()>1) mon.fillHisto("DR", "step1_fjet-en_cc", DRjet, weight);
		}
	    }

	  if (!overlap)
	    {
	      vec_fjets_raw.push_back(pfjet);
	      vec_fjets_index_raw.push_back(make_pair(pfjet, i));

	      // pT cut for the fat jets
	      if(pfjet.Pt() > 150)
		{
		  if( (ev.fjet_btag10[i]+ev.fjet_btag11[i]+ev.fjet_btag12[i]) / (ev.fjet_btag10[i]+ev.fjet_btag11[i]+ev.fjet_btag12[i]+ev.fjet_btag13[i]+ev.fjet_btag14[i]+ev.fjet_btag15[i]+ev.fjet_btag16[i]+ev.fjet_btag17[i]) > 0.5 )
		    {  
		      vec_fjets.push_back(pfjet);
		      vec_fjets_index.push_back(make_pair(pfjet, i));
		    }
		}

	    } //end overlap
	} // end fjets loop
      

      // Sort fat jets based on their pt
      vec_fjets_raw = sort_vec_pt(vec_fjets_raw);
      vec_fjets     = sort_vec_pt(vec_fjets);

      // Jet - fat-jet cross cleaning
      for (const auto& pair: vec_jetsbtag)
	{ 
	  if(getDeltaRmin(pair.first, vec_fjets) > 0.8)
	    {
	      vec_jetsbtag_cc.push_back(make_pair(pair.first, pair.second));
	      std::sort(vec_jetsbtag_cc.begin() , vec_jetsbtag_cc.end(),  sortPair);
	    }
	}

      // b-jet - fat-jet cross-cleaning
      for (const auto& pair: vec_bjetsbtag)
	{ 
	  if(getDeltaRmin(pair.first, vec_fjets) > 0.8)
	    {
	      vec_bjetsbtag_cc.push_back(make_pair(pair.first, pair.second));
	      std::sort(vec_bjetsbtag_cc.begin() , vec_bjetsbtag_cc.end(),  sortPair);
	    }
	}

      for (const auto& pair : vec_jetsbtag_cc)
	{
	  vec_jets_cc.push_back(pair.first);
	}

      for (const auto& pair : vec_bjetsbtag_cc)
	{
	  vec_bjets_cc.push_back(pair.first);
	}

      /////////////////////////////////////////////////////
      // E N D  O B J E C T S  C O N F I G U R A T I O N //
      /////////////////////////////////////////////////////

      // STEP 0  :  R A W  E V E N T S

      // Fill the multiplicity histograms for lepton
      mon.fillHisto("multi", "step0_en", vec_en.size(), weight);
      mon.fillHisto("multi", "step0_mn", vec_mn.size(), weight);
      mon.fillHisto("multi", "step0_lepton", vec_lepton_id.size(), weight);

      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "step0", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "step0_b", vec_btag.size(), weight);
      mon.fillHisto("jetmulti", "step0_untag", vec_untag.size(), weight);
      mon.fillHisto("jetmulti", "step0_f", vec_fjets.size(), weight);

      // Scalar sum of jets
      float HTWJet = HTjets(vec_jets);
      mon.fillHisto("HTWJets", "det", HTWJet, weight);

      
      if(vec_fjets_raw.size()>0)
	{
	  mon.fillHisto("pt", "step0_AK8_1", vec_fjets_raw[0].Pt(), weight);
	  mon.fillHisto("eta", "step0_AK8_1", vec_fjets_raw[0].Eta(), weight);
	  mon.fillHisto("phi", "step0_AK8_1", vec_fjets_raw[0].Phi(), weight);

	  mon.fillHisto("SDM", "step0_1", ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	  mon.fillHisto("fjet_pt_vs_SDM", "step0_1", vec_fjets_raw[0].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	  mon.fillProfile("fjetpt_vs_SDM", "step0_1", vec_fjets_raw[0].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	  mon.fillHisto("fjet_eta_vs_SDM", "step0_1", vec_fjets_raw[0].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	  mon.fillProfile("fjeteta_vs_SDM", "step0_1", vec_fjets_raw[0].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	  mon.fillHisto("fjet_subjet_vs_SDM", "step0_1", ev.fjet_subjet_count[0], ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	  mon.fillHisto("fjetsubjet_vs_SDM", "step0_1", ev.fjet_subjet_count[0], ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
      
	  bool matched(false);
	  bool matched1(false);
	  bool matched2(false);

	  // Discriminants
	  float denominator = ev.fjet_btag10[vec_fjets_index_raw[0].second] + ev.fjet_btag11[vec_fjets_index_raw[0].second] + ev.fjet_btag12[vec_fjets_index_raw[0].second] + ev.fjet_btag13[vec_fjets_index_raw[0].second] + ev.fjet_btag14[vec_fjets_index_raw[0].second] + ev.fjet_btag15[vec_fjets_index_raw[0].second] + ev.fjet_btag16[vec_fjets_index_raw[0].second] + ev.fjet_btag17[vec_fjets_index_raw[0].second];
     
	  float xbbccqq1 = (ev.fjet_btag10[vec_fjets_index_raw[0].second] + ev.fjet_btag11[vec_fjets_index_raw[0].second] + ev.fjet_btag12[vec_fjets_index_raw[0].second]) / denominator;
	  float xbb1     = ev.fjet_btag10[vec_fjets_index_raw[0].second] / (denominator - ev.fjet_btag11[vec_fjets_index_raw[0].second] - ev.fjet_btag12[vec_fjets_index_raw[0].second]);

	  mon.fillHisto("fjet_XbbXccXqq_pt", "1", vec_fjets_raw[0].Pt(), xbbccqq1, weight);
	  mon.fillProfile("XbbXccXqq_pt", "1", vec_fjets_raw[0].Pt(), xbbccqq1, weight);
	  mon.fillHisto("fjet_Xbb_pt", "1", vec_fjets_raw[0].Pt(), xbb1, weight);
	  mon.fillProfile("Xbb_pt", "1", vec_fjets_raw[0].Pt(), xbb1, weight);
	  mon.fillHisto("discriminant", "xbbccqq1", xbbccqq1, weight);
	  mon.fillHisto("discriminant", "xbb1", xbb1, weight);

	  if(isSignal)
	    {
	      if(fjet_matched(vec_fjets_raw[0], vec_bb1)>0)
		{
		  mon.fillHisto("fjetpt_vs_bbpt", "step0_1_matched", vec_fjets_raw[0].Pt(), (vec_bb1[0]+vec_bb1[1]).Pt(), weight);
		  matched1 = true;
		} 
	      else if(fjet_matched(vec_fjets_raw[0], vec_bb2)>0)
		{
		  mon.fillHisto("fjetpt_vs_bbpt", "step0_1_matched", vec_fjets_raw[0].Pt(), (vec_bb2[0]+vec_bb2[1]).Pt(), weight);
		  matched2 = true;
		}

	      if(matched1||matched2) matched = true;
	    }

	  if(matched)
	    {
	      mon.fillHisto("pt", "step0_AK8_1_matched", vec_fjets_raw[0].Pt(), weight);
	      mon.fillHisto("SDM", "step0_1_matched", ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillHisto("fjet_pt_vs_SDM", "step0_1_matched", vec_fjets_raw[0].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillProfile("fjetpt_vs_SDM", "step0_1_matched", vec_fjets_raw[0].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillHisto("fjet_eta_vs_SDM", "step0_1_matched", vec_fjets_raw[0].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillProfile("fjeteta_vs_SDM", "step0_1_matched", vec_fjets_raw[0].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	    }
	}

      if(vec_fjets_raw.size()>1)
	{
	  mon.fillHisto("pt", "step0_AK8_2", vec_fjets_raw[1].Pt(), weight);
	  mon.fillHisto("eta", "step0_AK8_2", vec_fjets_raw[1].Eta(), weight);
	  mon.fillHisto("phi", "step0_AK8_2", vec_fjets_raw[1].Phi(), weight);

	  mon.fillHisto("SDM", "step0_2", ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	  mon.fillHisto("fjet_pt_vs_SDM", "2", vec_fjets_raw[1].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	  mon.fillProfile("fjetpt_vs_SDM", "step0_2", vec_fjets_raw[1].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	  mon.fillHisto("fjet_eta_vs_SDM", "step0_2", vec_fjets_raw[1].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	  mon.fillProfile("fjeteta_vs_SDM", "step0_2", vec_fjets_raw[1].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	  mon.fillHisto("fjet_subjet_vs_SDM", "step0_2", ev.fjet_subjet_count[1], ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	  mon.fillHisto("fjetsubjet_vs_SDM", "step0_2", ev.fjet_subjet_count[1], ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
      
      
	  bool matched(false);
	  bool matched1(false);
	  bool matched2(false);

	  // Discriminants
	  float denominator = ev.fjet_btag10[vec_fjets_index_raw[1].second] + ev.fjet_btag11[vec_fjets_index_raw[1].second] + ev.fjet_btag12[vec_fjets_index_raw[1].second] + ev.fjet_btag13[vec_fjets_index_raw[1].second] + ev.fjet_btag14[vec_fjets_index_raw[1].second] + ev.fjet_btag15[vec_fjets_index_raw[1].second] + ev.fjet_btag16[vec_fjets_index_raw[1].second] + ev.fjet_btag17[vec_fjets_index_raw[1].second];
     
	  float xbbccqq2 = (ev.fjet_btag10[vec_fjets_index_raw[1].second] + ev.fjet_btag11[vec_fjets_index_raw[1].second] + ev.fjet_btag12[vec_fjets_index_raw[1].second]) / denominator;
	  float xbb2     = ev.fjet_btag10[vec_fjets_index_raw[1].second] / (denominator - ev.fjet_btag11[vec_fjets_index_raw[1].second] - ev.fjet_btag12[vec_fjets_index_raw[1].second]);

	  mon.fillHisto("fjet_XbbXccXqq_pt", "2", vec_fjets_raw[1].Pt(), xbbccqq2, weight);
	  mon.fillProfile("XbbXccXqq_pt", "2", vec_fjets_raw[1].Pt(), xbbccqq2, weight);
	  mon.fillHisto("fjet_Xbb_pt", "2", vec_fjets_raw[1].Pt(), xbb2, weight);
	  mon.fillProfile("Xbb_pt", "2", vec_fjets_raw[1].Pt(), xbb2, weight);
	  mon.fillHisto("discriminant", "xbbccqq2", xbbccqq2, weight);
	  mon.fillHisto("discriminant", "xbb2", xbb2, weight);

	  if(isSignal)
	    {
	      if(fjet_matched(vec_fjets_raw[1], vec_bb1)>0)
		{
		  mon.fillHisto("fjetpt_vs_bbpt", "step0_2_matched", vec_fjets_raw[1].Pt(), (vec_bb1[0]+vec_bb1[1]).Pt(), weight);
		  matched1 = true;
		} 
	      else if(fjet_matched(vec_fjets_raw[1], vec_bb2)>0)
		{
		  mon.fillHisto("fjetpt_vs_bbpt", "step0_2_matched", vec_fjets_raw[1].Pt(), (vec_bb2[0]+vec_bb2[1]).Pt(), weight);
		  matched2 = true;
		}
	      
	      if(matched1||matched2) matched = true;
	    }

	  if(matched)
	    {
	      mon.fillHisto("pt", "step0_AK8_2_matched", vec_fjets_raw[1].Pt(), weight);
	      mon.fillHisto("SDM", "step0_2_matched", ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	      mon.fillHisto("fjet_pt_vs_SDM", "step0_2_matched", vec_fjets_raw[1].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	      mon.fillProfile("fjetpt_vs_SDM", "step0_2_matched", vec_fjets_raw[1].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	      mon.fillHisto("fjet_eta_vs_SDM", "step0_2_matched", vec_fjets_raw[1].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	      mon.fillProfile("fjeteta_vs_SDM", "step0_2_matched", vec_fjets_raw[1].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[1].second], weight);
	    }
	}

      std::vector<TLorentzVector> vec_btag_pt;
      for (const auto& bjet: vec_btag)
	{
	  vec_btag_pt.push_back(bjet);
	}
      vec_btag_pt = sort_vec_pt(vec_btag_pt);

      if(vec_btag.size()>0)
       {
	 mon.fillHisto("pt", "AK4_1", vec_btag_pt[0].Pt(), weight);
	 mon.fillHisto("eta", "AK4_1", vec_btag_pt[0].Eta(), weight);
	 mon.fillHisto("phi", "AK4_1", vec_btag_pt[0].Phi(), weight);
	 bool matched(false);
	 bool matched1(false);bool matched2(false);
	 
	 if(isSignal)
	   {
	     if(fjet_matched(vec_btag_pt[0], vec_bb1)>0) matched1=true;
	       
	     else if(fjet_matched(vec_btag_pt[0], vec_bb2)>0) matched2=true;
	      
	     if(matched1 || matched2) matched=true;
	   }
	 
	 if(matched) mon.fillHisto("pt", "step0_AK4_matched", vec_btag_pt[0].Pt(), weight);	 
       }
      
      if(vec_btag.size()>1)
	{
	  mon.fillHisto("pt", "AK4_2", vec_btag_pt[1].Pt(), weight);
	  mon.fillHisto("eta", "AK4_2", vec_btag_pt[1].Eta(), weight);
	  mon.fillHisto("phi", "AK4_2", vec_btag_pt[1].Phi(), weight);
	}
      if(vec_btag.size()>2)
	{
	  mon.fillHisto("pt", "AK4_3", vec_btag_pt[2].Pt(), weight);
	  mon.fillHisto("eta", "AK4_3", vec_btag_pt[2].Eta(), weight);
	  mon.fillHisto("phi", "AK4_3", vec_btag_pt[2].Phi(), weight);
	}
       if(vec_btag.size()>3)
       {
	 mon.fillHisto("pt", "AK4_4", vec_btag_pt[3].Pt(), weight);
	 mon.fillHisto("eta", "AK4_4", vec_btag_pt[3].Eta(), weight);
	 mon.fillHisto("phi", "AK4_4", vec_btag_pt[3].Phi(), weight);
       }
      
       
      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//
      
      // Lepton pT cuts
      if (!vec_lepton_id.empty())
	{
	  if((abs(vec_lepton_id[0].second) == 11 && vec_lepton_id[0].first.Pt() >= 35) || (abs(vec_lepton_id[0].second) == 13 && vec_lepton_id[0].first.Pt() >= 30))
	    {
	      for (std::vector<TLorentzVector>::size_type i=0; i<vec_lepton_id.size(); i++)
		{
		  vec_lepton.push_back(vec_lepton_id[i].first);
		}
	    }
	}

      //###############//
      //### STEP 1 ###//
      //##############//
      
      // If there is at least one lepton
      if(vec_lepton.size() < 1) continue;

      n_event_lepton++;

      mon.fillHisto("eventflow", "histo", 1, 1);
      mon.fillHisto("eventflow", "histo_weighted", 1, weight);

     if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 1, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 1, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 1, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 1, weight);
	}

      // Fill the multiplicity histograms
      mon.fillHisto("multi", "step1_lepton", vec_lepton.size(), weight);
      mon.fillHisto("multi", "step1_en", vec_en.size(), weight);
      mon.fillHisto("multi", "step1_mn", vec_mn.size(), weight);
      
      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "step1", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "step1_b", vec_btag.size(), weight);
      mon.fillHisto("jetmulti", "step1_untag", vec_untag.size(), weight);
      mon.fillHisto("jetmulti", "step1_cc_b", vec_bjets_cc.size(), weight);
      mon.fillHisto("jetmulti", "step1_f", vec_fjets.size(), weight);
      mon.fillHisto("jetmulti", "step1_cc", vec_jets_cc.size(), weight);

      // MET vs pt(bb)
      if(vec_fjets_raw.size()>0)
	{
	  mon.fillHisto("pt_bb_vs_met","1", met.Pt(), vec_fjets_raw[0].Pt(), weight);
	  mon.fillProfile("pt_fjet_vs_met_prof","1", met.Pt(), vec_fjets_raw[0].Pt(), weight);
	}                                                                                                          
      if(vec_fjets_raw.size()>1)
	{
	  mon.fillHisto("pt_bb_vs_met","2", met.Pt(), vec_fjets_raw[1].Pt(), weight);
	  mon.fillProfile("pt_fjet_vs_met_prof", "2", met.Pt(), vec_fjets_raw[1].Pt(), weight);
	}
    
      mon.fillHisto("pt", "step1_MET", met.Pt(), weight);

      //###############//
      //### STEP 2 ###//
      //##############//
      
      if(met.Pt() < 25 || getMT(vec_lepton[0], met) < 50) continue;

      n_event_metcut++;

      mon.fillHisto("eventflow", "histo", 2, 1);
      mon.fillHisto("eventflow", "histo_weighted", 2, weight);

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 2, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 2, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 2, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 2, weight);
	}

      
      //###############//
      //### STEP 3 ###//
      //##############//
      // We define two regions, SR1 for the resolved regime and SR2 for the merged regime

      bool isSR1(false);
      bool isSR2(false);

      if(vec_btag.size()>2) isSR1=true;
      if(vec_fjets.size()>0 && vec_bjets_cc.size()>0) isSR2=true;

      // Define the variables for mva
      Float_t mvaBDTada_SR1(-10.0), mvaBDTgrad_SR1(-10.0), mvaBDTada_SR2(-10.0), mvaBDTgrad_SR2(-10.0);
      Float_t SR1_Njets(0.0), SR2_Njets(0.0), SR1_HT(0.0), SR2_HT(0.0), SR1_H_pt(0.0), SR2_H_pt(0.0), SR1_W_pt(0.0), SR2_W_pt(0.0), SR1_lepton_pt(0.0), SR2_lepton_pt(0.0), SR1_bjet1_pt(0.0), SR2_bjet1_pt(0.0), SR2_fjet1_pt(0.0), SR1_MET(0.0), SR2_MET(0.0), SR1_H_M(0.0), SR2_H_M(0.0), SR1_MTW(0.0), SR2_MTW(0.0), SR1_dphiWH(0.0), SR2_dphiWH(0.0), SR1_detaWH(0.0), SR2_detaWH(0.0), SR1_dphiMetJetMin(0.0), SR2_dphiMetJetMin(0.0), SR1_dphiHMet(0.0), SR2_dphiHMet(0.0), SR1_DRbbav(0.0), SR2_DRjj(0.0), SR1_DMbbMin(0.0), SR1_mbbj(0.0), SR1_btag1(0.0), SR2_btag1(0.0), SR1_btag2(0.0), SR1_btag3(0.0), SR2_sdMass1(0.0), SR2_xbbccqq(0.0), SR2_xbb(0.0) ;


      //########################################//
      //########## SIGNAL REGION SR1 ###########//
      //########################################//
      
      if(isSR1)
	{
	  if(vec_jets.size()<3) continue;

	  n_event_jet++;

	  mon.fillHisto("eventflow", "histo", 3, 1);
	  mon.fillHisto("eventflow", "histo_weighted", 3, weight);

	  if (isSignal)
	    {
	      mon.fillHisto("eventflow", "Signal", 3, weight);

	      if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 3, weight);  
	      if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 3, weight);   
	      if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 3, weight);
	    }

	  // Fill the multiplicity histograms of jets
	  mon.fillHisto("jetmulti", "step3", vec_jets.size(), weight);
	  mon.fillHisto("jetmulti", "step3_b", vec_btag.size(), weight);
	  mon.fillHisto("jetmulti", "step3_untag", vec_untag.size(), weight);
	  mon.fillHisto("jetmulti", "step3_cc_b", vec_bjets_cc.size(), weight);
	  mon.fillHisto("jetmulti", "step3_f", vec_fjets.size(), weight);
	  mon.fillHisto("jetmulti", "step3_cc", vec_jets_cc.size(), weight);
	  
	  // b tagged jets working points
	  // Tight WP: 0.7476, Medium WP: 0.3040, Loose WP: 0.0532
	  if(vec_bjetsbtag[0].second < tightWP) continue;
     
	  for (const auto& pair : vec_bjetsbtag)
	    {
	      vec_bjets.push_back(pair.first);
	    }
	  
	  // Sort the b-jets based on pT in descending order
	  vec_bjets = sort_vec_pt(vec_bjets);

	  n_event_bjet++;

	  mon.fillHisto("eventflow", "histo", 4, 1);
	  mon.fillHisto("eventflow", "histo_weighted", 4, weight);

	  if (isSignal)
	    {
	      mon.fillHisto("eventflow", "Signal", 4, weight);

	      if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 4, weight);  
	      if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 4, weight);   
	      if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 4, weight);
	    }

	  // Define the variables

	  // b-jets kinematics
	  SR1_bjet1_pt    = vec_bjets[0].Pt();
	  float bjet1_eta = vec_bjets[0].Eta();
	  float bjet1_phi = vec_bjets[0].Phi();

	  float bjet2_pt  = vec_bjets[1].Pt();
	  float bjet2_eta = vec_bjets[1].Eta();
	  float bjet2_phi = vec_bjets[1].Phi();

	  float bjet3_pt  = vec_bjets[2].Pt();
	  float bjet3_eta = vec_bjets[2].Eta();
	  float bjet3_phi = vec_bjets[2].Phi();

	  float bjet4_pt  = vec_bjets[3].Pt();
	  float bjet4_eta = vec_bjets[3].Eta();
	  float bjet4_phi = vec_bjets[3].Phi();

	  // Scalar sum of jets
	  SR1_HT = HTjets(vec_jets);

	  // Missing transverse energy
	  SR1_MET = met.Pt();

	  // W transverse mass
	  SR1_MTW = getMT(vec_lepton[0], met);
        
	  // Lepton kinematics
	  SR1_lepton_pt    = vec_lepton[0].Pt();
	  float lepton_eta = vec_lepton[0].Eta();
	  float lepton_phi = vec_lepton[0].Phi();
     
	  // Define the leptonic system
	  TLorentzVector plepsys = vec_lepton[0] + met;
	  SR1_W_pt = plepsys.Pt();

	  // Define the hadronic system
	  TLorentzVector phadb;

	  if(vec_btag.size() > 3) phadb = vec_btag[0] + vec_btag[1] + vec_btag[2] + vec_btag[3];
	
	  else if(vec_btag.size() == 3 && vec_untag.size() > 0)  phadb = vec_btag[0] + vec_btag[1] + vec_btag[2] + vec_untag[0];

	  else phadb = vec_btag[0] + vec_btag[1] + vec_btag[2];

	  // Higgs kinematics
	  SR1_H_M   = phadb.M();
	  SR1_H_pt  = phadb.Pt();
	  float H_eta = phadb.Eta(); 
         
	  // WH system variables
	  SR1_dphiWH   = fabs(vec_lepton[0].DeltaPhi(phadb));
	  SR1_detaWH   = fabs(plepsys.Eta()-phadb.Eta());
	  SR1_dphiHMet = fabs(met.DeltaPhi(phadb));

	  // Delta phi of MET and lepton
	  float dphiMetLepton = fabs(met.DeltaPhi(vec_lepton[0]));
     
	  // Minimum Delta phi of MET and jet
	  SR1_dphiMetJetMin = getDeltaPhiMin(met, vec_jets);

	  // DeltaR bb average
	  SR1_DRbbav = getAverageDeltaR(vec_btag);;

	  // Minimum DeltaMass of bb pairs
	  if (vec_btag.size() > 3)
	    {
	      SR1_DMbbMin = getDeltaMassMin(vec_btag, 4);
	    }
	  
	  else if (vec_btag.size() == 3) 
	    {
	      if (vec_untag.size() > 0) 
		{
		  vector<TLorentzVector> vec_b = vec_btag;
		  vec_b.push_back(vec_untag[0]);
		  SR1_DMbbMin = getDeltaMassMin(vec_b, 4);
		} 
	      else 
		{
		  SR1_DMbbMin = getDeltaMassMin(vec_btag, 3);
		}
	    }
      
	  // Mass of bb pair with minimum DR and one jet
	  if(vec_untag.size() > 0)
	    {
	      SR1_mbbj = (getbbwithDRMin(vec_bjets) + vec_untag[0]).M();
	      mon.fillHisto("M", "bbj", SR1_mbbj, weight);
	    }
	  else
	    {
	      SR1_mbbj = getmbbb(vec_btag);
	      mon.fillHisto("M", "bbb", SR1_mbbj, weight);
	    }
	  
	  // Multiplicity of jets
	  SR1_Njets = vec_jets.size();

	  // b-tag discriminator for 3 most energetic b-jets
	  SR1_btag1 = vec_bjetsbtag[0].second;
	  SR1_btag2 = vec_bjetsbtag[1].second;
	  SR1_btag3 = vec_bjetsbtag[2].second;

	  // Invariant mass of bb pair from the same mother
	  PairResult result = getPair(vec_btag);
	  float bb_M;

	  if(vec_btag.size() == 3) bb_M = result.pair1.M() ;
	  else bb_M = (result.pair1.M() + result.pair2.M())/2.0;


	  // Fill the histograms
	  //mon.fillHisto("pt"         , "SR1_H"          , SR1_H_pt          , weight);
	  mon.fillHisto("H_pt"       , "SR1"            , SR1_H_pt          , weight);  // final
	  //mon.fillHisto("pt"         , "SR1_W"          , SR1_W_pt          , weight);
	  mon.fillHisto("W_pt"       , "SR1"            , SR1_W_pt          , weight);  // final
	  //mon.fillHisto("pt"         , "SR1_MET"        , SR1_MET           , weight);
	  mon.fillHisto("MET"        , "SR1"            , SR1_MET           , weight);  // final
	  mon.fillHisto("l_pt"       , "SR1"            , SR1_lepton_pt     , weight);  // final
	  //mon.fillHisto("pt"         , "SR1_lepton"     , SR1_lepton_pt     , weight);
	  //mon.fillHisto("pt"         , "SR1_bjet1"      , SR1_bjet1_pt      , weight);
	  mon.fillHisto("bjet1_pt"   , "SR1"            , SR1_bjet1_pt      , weight);  // final
	  mon.fillHisto("pt"         , "SR1_bjet2"      , bjet2_pt          , weight);
	  mon.fillHisto("pt"         , "SR1_bjet3"      , bjet3_pt          , weight);
	  mon.fillHisto("pt"         , "SR1_bjet4"      , bjet4_pt          , weight);
	  mon.fillHisto("HT"         , "SR1_jets"       , SR1_HT            , weight);
	  mon.fillHisto("eta"        , "SR1_H"          , H_eta             , weight);
	  //mon.fillHisto("deta"       , "SR1_WH"         , SR1_detaWH        , weight);
	  mon.fillHisto("WH_deta"    , "SR1"            , SR1_detaWH        , weight);  // final
	  mon.fillHisto("eta"        , "SR1_lepton"     , lepton_eta        , weight);
	  mon.fillHisto("eta"        , "SR1_bjet1"      , bjet1_eta         , weight);
	  mon.fillHisto("eta"        , "SR1_bjet2"      , bjet2_eta         , weight);
	  mon.fillHisto("eta"        , "SR1_bjet3"      , bjet3_eta         , weight);
	  mon.fillHisto("eta"        , "SR1_bjet4"      , bjet4_eta         , weight);
	  //mon.fillHisto("MT"         , "SR1_W"          , SR1_MTW           , weight);
	  mon.fillHisto("W_MT"       , "SR1"            , SR1_MTW           , weight);  // final
	  //mon.fillHisto("M"          , "SR1_H"          , SR1_H_M           , weight);
	  mon.fillHisto("H_M"        , "SR1"            , SR1_H_M           , weight);  // final
	  //mon.fillHisto("M"          , "SR1_dMbbmin"    , SR1_DMbbMin       , weight);
	  mon.fillHisto("minbb_dM"   , "SR1"            , SR1_DMbbMin       , weight);  // final
	  //mon.fillHisto("M"          , "SR1_bbj"        , SR1_mbbj          , weight);
	  mon.fillHisto("bbj_M"      , "SR1"            , SR1_mbbj          , weight);  // final
	  mon.fillHisto("M"          , "SR1_bb"         , bb_M              , weight);
	  //mon.fillHisto("dphi"       , "SR1_WH"         , SR1_dphiWH        , weight);
	  mon.fillHisto("WH_dphi"    , "SR1"            , SR1_dphiWH        , weight);  // final
	  mon.fillHisto("dphi"       , "SR1_MET-lepton" , dphiMetLepton     , weight);
	  //mon.fillHisto("dphi"       , "SR1_METjet-min" , SR1_dphiMetJetMin , weight);
	  mon.fillHisto("METjet_dphi", "SR1"            , SR1_dphiMetJetMin , weight);  // final
	  //mon.fillHisto("DR"         , "SR1_bbaver"     , SR1_DRbbav        , weight);
	  mon.fillHisto("avbb_dR"    , "SR1"            , SR1_DRbbav        , weight);  // final
	  mon.fillHisto("jetmulti"   , "SR1_step3"      , vec_jets.size()   , weight);
	  mon.fillHisto("jetmulti"   , "SR1_step3_b"    , vec_bjets.size()  , weight);
	  mon.fillHisto("jetmulti"   , "SR1_step3_untag", vec_untag.size()  , weight);
	  mon.fillHisto("btag1_bef"  , "SR1"            , SR1_btag1_WP      , weight);  // final
	  mon.fillHisto("btag2_bef"  , "SR1"            , SR1_btag2_WP      , weight);  // final
	  mon.fillHisto("btag3_bef"  , "SR1"            , SR1_btag3_WP      , weight);  // final
	  mon.fillHisto("btag1"      , "SR1"            , SR1_btag1         , weight);  // final
	  mon.fillHisto("btag2"      , "SR1"            , SR1_btag2         , weight);  // final
	  mon.fillHisto("btag3"      , "SR1"            , SR1_btag3         , weight);  // final

	  //#####################################//
	  //############ TMVA Reader ############//
	  //#####################################//
	  
	  mvaBDTada_SR1  = SR1adaReader.GenReMVAReaderSR1(SR1_Njets, SR1_HT, SR1_H_pt, SR1_W_pt, SR1_lepton_pt, SR1_bjet1_pt, SR1_MET, SR1_H_M, SR1_MTW, SR1_dphiWH, SR1_dphiMetJetMin, SR1_DRbbav, SR1_DMbbMin, SR1_mbbj, SR1_btag1, SR1_btag2, SR1_btag3, "SR1adaClass");
	  mvaBDTgrad_SR1 = SR1gradReader.GenReMVAReaderSR1(SR1_Njets, SR1_HT, SR1_H_pt, SR1_W_pt, SR1_lepton_pt, SR1_bjet1_pt, SR1_MET, SR1_H_M, SR1_MTW, SR1_dphiWH, SR1_dphiMetJetMin, SR1_DRbbav, SR1_DMbbMin, SR1_mbbj, SR1_btag1, SR1_btag2, SR1_btag3, "SR1gradClass");
	  mon.fillHisto("BDT", "SR1_ada", mvaBDTada_SR1, weight);
	  mon.fillHisto("BDT", "SR1_grad", mvaBDTgrad_SR1, weight);
	  
	}


      //########################################//
      //########## SIGNAL REGION SR2 ###########//
      //########################################//

      if(isSR2)
	{
	  n_SR2++;

	  // mon.fillHisto("eventflow", "histo", 5, 1);
	  // mon.fillHisto("eventflow", "histo_weighted", 5, weight);

	  // if (isSignal)
	  //   {
	  //     mon.fillHisto("eventflow", "Signal", 5, weight);

	  //     if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 5, weight);  
	  //     if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 5, weight);   
	  //     if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 5, weight);
	  //   }

	  mon.fillHisto("count", "step3_subjet", ev.fjet_subjet_count[vec_fjets_index[0].second], weight);

	  if(vec_fjets_raw.size()>0)
	    {
	      mon.fillHisto("pt", "step3_AK8_1", vec_fjets_raw[0].Pt(), weight);
	      mon.fillHisto("eta", "step3_AK8_1", vec_fjets_raw[0].Eta(), weight);
	      mon.fillHisto("phi", "step3_AK8_1", vec_fjets_raw[0].Phi(), weight);

	      mon.fillHisto("SDM", "step3_1", ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillHisto("fjet_pt_vs_SDM", "step3_1", vec_fjets_raw[0].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillProfile("fjetpt_vs_SDM", "step3_1", vec_fjets_raw[0].Pt(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillHisto("fjet_eta_vs_SDM", "step3_1", vec_fjets_raw[0].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	      mon.fillProfile("fje_eta_vs_SDM", "step3_1", vec_fjets_raw[0].Eta(), ev.fjet_softdropM[vec_fjets_index_raw[0].second], weight);
	    }

	  // Fill the multiplicity histograms of jets
	  mon.fillHisto("jetmulti", "step3_cc_b", vec_bjets_cc.size(), weight);
	  mon.fillHisto("jetmulti", "step3_f", vec_fjets.size(), weight);
	  mon.fillHisto("jetmulti", "step3_cc", vec_jets_cc.size(), weight);

	  // b-jets kinematics
	  SR2_bjet1_pt    = vec_bjets_cc[0].Pt();
	  float bjet1_eta = vec_bjets_cc[0].Eta();
	  float bjet1_phi = vec_bjets_cc[0].Phi();

	  float bjet2_pt  = vec_bjets_cc[1].Pt();
	  float bjet2_eta = vec_bjets_cc[1].Eta();
	  float bjet2_phi = vec_bjets_cc[1].Phi();

	   // Multiplicity of jets
	  SR2_Njets = vec_jets_cc.size();
	  int Nsubjet = ev.fjet_subjet_count[vec_fjets_index[0].second];

	  // Scalar sum of jets
	  if(vec_jets_cc.size()>0) { for(int i=0; i<vec_jets_cc.size(); i++) SR2_HT += vec_jets_cc[i].Pt(); }
	  if(vec_fjets.size()>0)   { for(int j=0; j<vec_fjets.size(); j++)   SR2_HT += vec_fjets[j].Pt(); }

	  // Missing transverse energy
	  SR2_MET = met.Pt();

	  // W transverse mass
	  SR2_MTW = getMT(vec_lepton[0], met);
        
	  // Lepton kinematics
	  SR2_lepton_pt    = vec_lepton[0].Pt();
	  float lepton_eta = vec_lepton[0].Eta();
	  float lepton_phi = vec_lepton[0].Phi();
     
	  // Define the leptonic system
	  TLorentzVector plepsys = vec_lepton[0] + met;
	  SR2_W_pt = plepsys.Pt();

	  // Delta phi of MET and lepton
	  float dphiMetLepton = fabs(met.DeltaPhi(vec_lepton[0]));

	  //## FAT jet Analysis ##//
	  //----------------------//

	  // Fat jet kinematics
	  SR2_fjet1_pt    = vec_fjets[0].Pt();
	  float fjet1_eta = vec_fjets[0].Eta();
	  float fjet1_phi = vec_fjets[0].Phi();

	  float fjet2_pt  = vec_fjets[1].Pt();
	  float fjet2_eta = vec_fjets[1].Eta();
	  float fjet2_phi = vec_fjets[1].Phi();

	  // Soft drop mass
	  SR2_sdMass1 = ev.fjet_softdropM[vec_fjets_index[0].second];
	  
	  // Discriminators
	  SR2_btag1 = vec_bjetsbtag_cc[0].second;
	  
	  float denominator = ev.fjet_btag10[vec_fjets_index[0].second] + ev.fjet_btag11[vec_fjets_index[0].second] + ev.fjet_btag12[vec_fjets_index[0].second] + ev.fjet_btag13[vec_fjets_index[0].second] + ev.fjet_btag14[vec_fjets_index[0].second] + ev.fjet_btag15[vec_fjets_index[0].second] + ev.fjet_btag16[vec_fjets_index[0].second] + ev.fjet_btag17[vec_fjets_index[0].second];
     
	  SR2_xbbccqq = (ev.fjet_btag10[vec_fjets_index[0].second] + ev.fjet_btag11[vec_fjets_index[0].second] + ev.fjet_btag12[vec_fjets_index[0].second]) / denominator;
	  SR2_xbb     = ev.fjet_btag10[vec_fjets_index[0].second] / (denominator - ev.fjet_btag11[vec_fjets_index[0].second] - ev.fjet_btag12[vec_fjets_index[0].second]);

	  // DR b-jet - fat-jet
	  SR2_DRjj = getDeltaR(vec_fjets[0], vec_bjets_cc[0]);
	  
	   // Minimum Delta phi of MET and jet
	  SR2_dphiMetJetMin = std::min(fabs(getDeltaPhiMin(met, vec_fjets)), fabs(getDeltaPhiMin(met, vec_bjets_cc)));

	  // Define the hadronic system
	  TLorentzVector phad;

	  if(vec_bjets_cc.size() == 1)
	    {
	      if(vec_jets_cc.size() == 1) phad = vec_fjets[0] + vec_bjets_cc[0];
	
	      else if(vec_jets_cc.size() > 1)  phad = vec_fjets[0] + vec_jets_cc[0] + vec_jets_cc[1];
	    }

	  else if(vec_bjets_cc.size() > 1) phad = vec_fjets[0] + vec_bjets_cc[0] + vec_bjets_cc[1];
	      
	  // Higgs kinematics
	  SR2_H_M   = phad.M();
	  SR2_H_pt  = phad.Pt();
	  float H_eta = phad.Eta(); 

	  // WH system variables
	  if(fabs(vec_lepton[0].DeltaPhi(phad)) < TMath::Pi()) SR2_dphiWH = fabs(vec_lepton[0].DeltaPhi(phad));
	  else  SR2_dphiWH   = 2*TMath::Pi() - fabs(vec_lepton[0].DeltaPhi(phad));
	  SR2_detaWH   = fabs(plepsys.Eta()-phad.Eta());
	  SR2_dphiHMet = fabs(met.DeltaPhi(phad));

	  
	  // Fill the histograms
	  mon.fillHisto("pt"               , "SR2_H"          , SR2_H_pt                  , weight);
	  mon.fillHisto("pt"               , "SR2_W"          , SR2_W_pt                  , weight);
	  mon.fillHisto("pt"               , "SR2_MET"        , SR2_MET                   , weight);
	  mon.fillHisto("pt"               , "SR2_lepton"     , SR2_lepton_pt             , weight);
	  mon.fillHisto("pt"               , "SR2_bjet1"      , SR2_bjet1_pt              , weight);
	  mon.fillHisto("pt"               , "SR2_bjet2"      , bjet2_pt                  , weight);
	  mon.fillHisto("pt"               , "SR2_fjet1"      , SR2_fjet1_pt              , weight);
	  mon.fillHisto("pt"               , "SR2_fjet2"      , fjet2_pt                  , weight);
	  mon.fillHisto("HT"               , "SR2_jets"       , SR2_HT                    , weight);
	  mon.fillHisto("eta"              , "SR2_H"          , H_eta                     , weight);
	  mon.fillHisto("eta"              , "SR2_WH"         , SR2_detaWH                , weight);
	  mon.fillHisto("eta"              , "SR2_lepton"     , lepton_eta                , weight);
	  mon.fillHisto("eta"              , "SR2_bjet1"      , bjet1_eta                 , weight);
	  mon.fillHisto("eta"              , "SR2_bjet2"      , bjet2_eta                 , weight);
	  mon.fillHisto("eta"              , "SR2_fjet1"      , fjet1_eta                 , weight);
	  mon.fillHisto("eta"              , "SR2_fjet2"      , fjet2_eta                 , weight);
	  mon.fillHisto("MT"               , "SR2_W"          , SR2_MTW                   , weight);
	  mon.fillHisto("M"                , "SR2_H"          , SR2_H_M                   , weight);
	  mon.fillHisto("dEta"             , "SR2_WH"         , SR2_detaWH                , weight);
	  mon.fillHisto("dphi"             , "SR2_WH"         , SR2_dphiWH                , weight);
	  mon.fillHisto("dphi"             , "SR2_MET-lepton" , dphiMetLepton             , weight);
	  mon.fillHisto("dphi"             , "SR2_METjet-min" , SR2_dphiMetJetMin         , weight);
	  mon.fillHisto("dphi"             , "SR2_H-MET"      , SR2_dphiHMet              , weight);
	  mon.fillHisto("DR"               , "SR2_jj"         , SR2_DRjj                  , weight);
	  mon.fillHisto("jetmulti"         , "SR2_step3"      , vec_jets_cc.size()        , weight);
	  mon.fillHisto("jetmulti"         , "SR2_step3_b"    , vec_bjets_cc.size()       , weight);
	  mon.fillHisto("jetmulti"         , "SR2_step3_f"    , vec_fjets.size()          , weight);
	  mon.fillHisto("btag"             , "SR2_jet1"       , SR2_btag1                 , weight);
	  mon.fillHisto("fjet_XbbXccXqq_pt", "SR2_step3"      , SR2_fjet1_pt, SR2_xbbccqq , weight);
	  mon.fillHisto("fjet_Xbb_pt"      , "SR2_step3"      , SR2_fjet1_pt, SR2_xbb     , weight);
	  mon.fillHisto("discriminant"     , "SR2_xbbccqq1"   , SR2_xbbccqq               , weight);
	  mon.fillHisto("discriminant"     , "SR2_xbb1"       , SR2_xbb                   , weight);
	  mon.fillHisto("count"            , "SR2_subjet"     , Nsubjet                   , weight);
	  mon.fillHisto("SDM"              , "SR2_step3"      , SR2_sdMass1               , weight);
	  mon.fillHisto("fjet_pt_vs_SDM"   , "SR2_step3"      , SR2_fjet1_pt, SR2_sdMass1 , weight);
	  mon.fillProfile("fjetpt_vs_SDM"  , "SR2_step3"      , SR2_fjet1_pt, SR2_sdMass1 , weight);

	  //#####################################//
	  //############ TMVA Reader ############//
	  //#####################################//
	  
	  mvaBDTada_SR2  = SR2adaReader.GenReMVAReaderSR2(SR2_Njets, SR2_HT, SR2_H_pt, SR2_W_pt, SR2_lepton_pt, SR2_bjet1_pt, SR2_fjet1_pt, SR2_MET, SR2_H_M, SR2_MTW, SR2_dphiWH, SR2_dphiMetJetMin, SR2_DRjj, SR2_btag1, SR2_sdMass1, SR2_xbbccqq, SR2_xbb,  "SR2adaClass");
	  mvaBDTgrad_SR2 = SR2gradReader.GenReMVAReaderSR2(SR2_Njets, SR2_HT, SR2_H_pt, SR2_W_pt, SR2_lepton_pt, SR2_bjet1_pt, SR2_fjet1_pt, SR2_MET, SR2_H_M, SR2_MTW, SR2_dphiWH, SR2_dphiMetJetMin, SR2_DRjj, SR2_btag1, SR2_sdMass1, SR2_xbbccqq, SR2_xbb, "SR2gradClass");
	  mon.fillHisto("BDT", "SR2_ada", mvaBDTada_SR2, weight);
	  mon.fillHisto("BDT", "SR2_grad", mvaBDTgrad_SR2, weight);
	  
	}

	

      //#####################################//
      //############ MVA Handler ############//
      //#####################################//
   
      if(runMVA)
	{
	  float mvaweight = 1.0; 
	  myMVAHandler_.getEntry(isSR1, isSR2, SR1_Njets, SR2_Njets, SR1_HT, SR2_HT, SR1_H_pt, SR2_H_pt, SR1_W_pt, SR2_W_pt, SR1_lepton_pt, SR2_lepton_pt, SR1_bjet1_pt, SR2_bjet1_pt, SR2_fjet1_pt, SR1_MET, SR2_MET, SR1_H_M, SR2_H_M, SR1_MTW, SR2_MTW, SR1_dphiWH, SR2_dphiWH, SR1_detaWH, SR2_detaWH, SR1_dphiMetJetMin, SR2_dphiMetJetMin, SR1_dphiHMet, SR2_dphiHMet, SR1_DRbbav, SR2_DRjj, SR1_DMbbMin, SR1_mbbj, SR1_btag1, SR2_btag1, SR1_btag2, SR1_btag3, SR2_sdMass1, SR2_xbbccqq, SR2_xbb, weight);
	  myMVAHandler_.fillTree();
	    
	}  
      
    } // END OF EVENT LOOP

  // Write MVA files
  TString mvaout = TString ( runProcess.getParameter<std::string>("outdir") ) + "/mva_" + outFileUrl + ".root";
  if(runMVA){ myMVAHandler_.writeTree(mvaout); }

	
  //###############################################//
  //########     SAVING HISTOS TO FILE     ########//
  //###############################################//
  outUrl += "/";
  outUrl += outFileUrl + ".root";
  printf("Results saved in %s\n", outUrl.Data());
  
  if(isMC_ttbar)
    {
      float cFilt1=(float(nFilt1)/float(nTot));
      float cFilt4=(float(nFilt4)/float(nTot));
      float cFilt5=(float(nFilt5)/float(nTot));

      mon.fillHisto("filter", "1", cFilt1, weight);
      mon.fillHisto("filter", "4", cFilt4, weight);
      mon.fillHisto("filter", "5", cFilt5, weight);

      if(verbose) printf("From total = %i TTbar events ,  Found %.2f (filt1) , %.2f (filt4) , %.2f (filt5) \n\n", nTot, cFilt1, cFilt4, cFilt5);
    }
  

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
  
} // END MAIN

//#################//
//### FUNCTIONS ###//
//#################//

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
// Angular distance DR
float getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2)
{
  float delta_phi;
  float delta_eta;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_eta = vec_1.Eta() - vec_2.Eta();

  return std::sqrt(pow(delta_phi, 2) + pow(delta_eta, 2));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Get DR minimum
float getDeltaRmin(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2)
{
  float dR[vec_2.size()];
  float dRmin  = 0.0;
  int minInd   = 0;

  for(int i=0; i<vec_2.size(); i++)
    {
      dR[i] = getDeltaR(vec_1, vec_2[i]);

      if(dR[i] < dR[minInd]) minInd = i;
    }

  dRmin = dR[minInd];

  return dRmin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Average angular distance DR
float getAverageDeltaR(vector<TLorentzVector> vec)
{
  float sum = 0.0;
  int npairs = 0;

  for(size_t i=0; i<vec.size(); i++)
    {
      for(size_t j=i+1; j<vec.size(); j++)
	{
	  sum += getDeltaR(vec[i], vec[j]);
	  npairs ++ ;
	}
    }

  return sum/npairs;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Minimum Delta mass
float getDeltaMassMin(vector<TLorentzVector> vec, int n)
{
    float minDeltaMass = std::numeric_limits<float>::max();

    if (n == 4)
      {
	for (size_t i = 0; i < 4; ++i)
	  {
	    for (size_t j = i + 1; j < 4; ++j)
	      {
		TLorentzVector pair1 = vec[i] + vec[j];
		TLorentzVector pair2;

		for (size_t k = 0; k < 4; ++k)
		  {
		    if (k != i && k != j)
		      {
			pair2 += vec[k];
		      }
		  }
		
		float mass1 = pair1.M();
		float mass2 = pair2.M();
		float deltaMass = std::fabs(mass1 - mass2);

		if (deltaMass < minDeltaMass)
		  {
		    minDeltaMass = deltaMass;
		  }   
	      }
	  }
      }
    else if (n == 3)
      {
        for (size_t i = 0; i < 3; ++i)
	  {
	    size_t j = (i + 1) % 3;
	    size_t k = (i + 2) % 3;

	    TLorentzVector pair1 = vec[i] + vec[j];
	    TLorentzVector pair2 = vec[k] + vec[i];

	    float mass1 = pair1.M();
	    float mass2 = pair2.M();
	    float deltaMass = std::fabs(mass1 - mass2);

	    if (deltaMass < minDeltaMass)
	      {
		minDeltaMass = deltaMass;
	      }

	    pair1 = vec[i] + vec[k];
	    pair2 = vec[j] + vec[k];

	    mass1 = pair1.M();
	    mass2 = pair2.M();
	    deltaMass = std::fabs(mass1 - mass2);

	    if (deltaMass < minDeltaMass)
	      {
		minDeltaMass = deltaMass;
	      }
	  }
      }
    
    return minDeltaMass;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Another way to compute Delta mass bb min, for Nb=3
float getDeltaMassMinNew(vector<TLorentzVector> vec, float& unpairedmass)
{
  float minDeltaMass = std::numeric_limits<float>::max();

  for (size_t i = 0; i < 3; ++i)
    {
      size_t j = (i + 1) % 3;
      size_t k = (i + 2) % 3;

      TLorentzVector pair1 = vec[i] + vec[j];
      TLorentzVector pair2 = vec[k];

      float mass1 = pair1.M();
      float mass2 = pair2.M();
      float deltaMass = std::fabs(mass1 - mass2);

      if (deltaMass < minDeltaMass)
	{
	  minDeltaMass = deltaMass;
	  unpairedmass = vec[k].M();
	}

      pair1 = vec[i] + vec[k];
      pair2 = vec[j];

      mass1 = pair1.M();
      mass2 = pair2.M();
      deltaMass = std::fabs(mass1 - mass2);

      if (deltaMass < minDeltaMass)
	{
	  minDeltaMass = deltaMass;
	  unpairedmass = vec[j].M();
	}
    }

  return minDeltaMass;
}

 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Delta R minimum
TLorentzVector getbbwithDRMin(vector<TLorentzVector> vec)
{
  float mindR = std::numeric_limits<double>::max();
  TLorentzVector bjet1, bjet2;

  for(size_t i=0; i<vec.size(); i++)
    {
      for(size_t j=i+1; j<vec.size(); j++)
	{
	  float dR = vec[i].DeltaR(vec[j]);
	  if(dR<mindR)
	    {
	      mindR = dR;
	      bjet1 = vec[i];
	      bjet2 = vec[j];
	    }
	}
    }

  TLorentzVector pair = bjet1+bjet2;

  return pair;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float fjet_matched(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2)
{
  float qg_pt=-999.0;

  std::vector<TLorentzVector> vec_min;
  TLorentzVector vec_qgpt;
  
  bool ismatched(false);

  for (int i=0; i<vec_2.size(); i++)
    {
      ismatched=false;
      
      float dR = getDeltaR(vec_1, vec_2[i]);
      
      if(dR<0.8)
	{
	  ismatched=true;
	  vec_min.push_back(vec_2[i]);
	}
    }
  
  for (int i=0; i<vec_min.size(); i++)
    {
      vec_qgpt += vec_min[i];
    }
  
  if(ismatched) qg_pt = vec_qgpt.Pt();
  
  return qg_pt;  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float getmbbb(vector<TLorentzVector> vec)
{
  size_t n = vec.size();
  
  float mindR = std::numeric_limits<float>::max();
  TLorentzVector bjet1, bjet2;
  size_t index1 = 0, index2 = 0;

  for (size_t i = 0; i < n; ++i)
    {
      for (size_t j = i + 1; j < n; ++j)
        {
	  float dR = vec[i].DeltaR(vec[j]);
	  
	  if (dR < mindR || (dR == mindR && (i < index1 || (i == index1 && j < index2))))
            {
	      mindR = dR;
	      bjet1 = vec[i];
	      bjet2 = vec[j];
	      index1 = i;
	      index2 = j;
            }
        }
    }

  TLorentzVector remainingBJet;
  size_t remainingIndex = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < n; ++i)
    {
      if (i != index1 && i != index2)
        {
	  if (i < remainingIndex)
            {
	      remainingBJet = vec[i];
	      remainingIndex = i;
            }
        }
    }

  TLorentzVector combinedPair = bjet1 + bjet2 + remainingBJet;
  return combinedPair.M();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Scalar sum of pTs
float HTjets(std::vector<TLorentzVector> vec)
{
  float ht = 0.0;

  if(vec.size()>0)
    {
      for(int i=0; i<vec.size(); i++) ht += vec[i].Pt();
    }
  return ht;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Transverse mass
float getMT(TLorentzVector vec_1, TLorentzVector vec_2)
{
  float dphi = vec_1.DeltaPhi(vec_2);

  return std::sqrt(2 * vec_1.Pt() * vec_2.Pt() * (1 - TMath::Cos(dphi)));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculation of the transverse mass for mT2
float computeMT2(TLorentzVector vis1, TLorentzVector vis2, TLorentzVector invis)
{
  float minMT2 = 1e9;

  // Iterate over possible fractions of MET for the two hemispheres
  for (float alpha = 0; alpha<=1.0; alpha+=0.01)
    {
      // Split MET between the two hemispheres
      TLorentzVector met1 = alpha*invis;
      TLorentzVector met2 = (1-alpha)*invis;

      // Compute tranverse mass for the two hemispheres
      float mt1 = getMT(vis1, met1);
      float mt2 = getMT(vis2, met2);

      float maxMT = std::max(mt1, mt2);

      // Minimize the maximum transverse mass for all splits
      if(maxMT<minMT2)
	{
	  minMT2 = maxMT;
	}
    }

  return minMT2;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Find the bb pair from the same mother 
PairResult getPair(vector<TLorentzVector> vec)
{
  size_t n = vec.size();

  PairResult bestResult;
  
  if (n == 3)
    {
      float minDeltaR = std::numeric_limits<float>::max();

      for (size_t i = 0; i < n; ++i)
	{
	  for (size_t j = i + 1; j < n; ++j)
	    {
	      float deltaR = vec[i].DeltaR(vec[j]);

	      if (deltaR < minDeltaR)
		{
		  minDeltaR = deltaR;
		  bestResult.pair1 = vec[i] + vec[j];
		}
	    }
	}
    }
  
  else if(n > 3)
    {
      float minDeltaMass = std::numeric_limits<float>::max();

      for (size_t i = 0; i < 4; ++i)
	{
	  for (size_t j = i + 1; j < 4; ++j)
	    {
	      TLorentzVector pair1 = vec[i] + vec[j];
	      TLorentzVector pair2;

	      for (size_t k = 0; k < 4; ++k)
		{
		  if (k != i && k != j)
		    {
		      pair2 += vec[k];
		    }
		}

	      float mass1 = pair1.M();
	      float mass2 = pair2.M();
	      float deltaMass = std::fabs(mass1 - mass2);

	      if (deltaMass < minDeltaMass)
		{
		  minDeltaMass = deltaMass;
		  bestResult.pair1 = pair1;
		  bestResult.pair2 = pair2;
		}
	    }
	}
    }
  return bestResult;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Sorting based on pT
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec)
{
  std::sort(vec.begin(), vec.end(), [](const TLorentzVector& a, const TLorentzVector& b) {
    return a.Pt() > b.Pt();
  });

  return vec;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool sortPair(std::pair<TLorentzVector, float> vec_1, std::pair<TLorentzVector, float> vec_2)
{
  return vec_1.second > vec_2.second;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float getDeltaPhiMin(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2)
{
  float dphiMin = std::numeric_limits<float>::max();
  
  for (const auto& vec : vec_2)
    {
      float dphi = fabs(vec_1.DeltaPhi(vec));
      if (dphi < dphiMin)
	{
	  dphiMin = dphi;
	}
    }

  return dphiMin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
