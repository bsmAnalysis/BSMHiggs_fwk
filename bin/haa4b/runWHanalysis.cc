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
float getAverageDeltaR(vector<TLorentzVector> vec);
float getdeltaMass(TLorentzVector b1, TLorentzVector b2, float targetMass);
float getDeltaMassMin(vector<TLorentzVector> vec, int n);
float getDeltaMassMinNew(vector<TLorentzVector> vec, float& unpairedmass);
float getmbbb(vector<TLorentzVector> vec);
float computeMT2(TLorentzVector vis1, TLorentzVector vis2, TLorentzVector invis);
float getDeltaPhiMin(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2);
float getMT(TLorentzVector vec_1, TLorentzVector vec_2);
float HTjets(std::vector<TLorentzVector> vec);
TLorentzVector getbbwithDRMin(vector<TLorentzVector> vec);

std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec);
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
  float Lint   = 41.5;
  float Nexp   = xsec * 1000 * Lint;
  float weight = Nexp / nevts;
  if(verbose)
    {
      cout << "Cross section : " << xsec   << " [pb]" << endl;
      cout << "Weight        : " << weight            << endl;
      cout << "N expected    : " << Nexp              << endl; 
    }
 
  //##################################################################################//
  //##########################    INITIATING HISTOGRAMS     ##########################//
  //##################################################################################//
  SmartSelectionMonitor mon;

  // Count the events with 3 or more than 3 jets
  TH1F *hj = (TH1F*) mon.addHistogram ( new TH1F ("Events", ";;N_{Events}", 2,0,2) );
  hj->GetXaxis()->SetBinLabel(1,"3 jets");
  hj->GetXaxis()->SetBinLabel(2,">3 jets");
  
  // Event flow 
  TH1F *h = (TH1F*) mon.addHistogram ( new TH1F ("eventflow", ";;N_{Events}", 5,0,5) );
  h->GetXaxis()->SetBinLabel(1,"Raw");
  h->GetXaxis()->SetBinLabel(2,"1 lepton");
  h->GetXaxis()->SetBinLabel(3,">=3 jets");
  h->GetXaxis()->SetBinLabel(4,">=3 b-tags");
  h->GetXaxis()->SetBinLabel(5,"MET>25 & MTW>50");

  // Multiplicity 
  mon.addHistogram ( new TH1F ("multi", ";N_{lepton} ;N_{Events}", 5, 0, 5) );
  mon.addHistogram ( new TH1F ("qgmulti", ";N_{qg} ;N_{Events}", 15, 0, 15) );
  mon.addHistogram ( new TH1F ("jetmulti", ";N_{jet} ;N_{Events}", 15, 0, 15) );

  // Particles/Objects kinematics
  mon.addHistogram ( new TH1F ("pt", ";p_{T} [GeV] ;N_{Events}", 100, 0, 500) );
  mon.addHistogram ( new TH1F ("eta", ";#eta ;N_{Events}", 100, -6, 6));
  mon.addHistogram ( new TH1F ("phi", ";#phi ;N_{Events}", 100, -TMath::Pi(), TMath::Pi()) );
  mon.addHistogram ( new TH1F ("M", ";M [GeV] ;N_{Events}", 100, 0, 800) );
  mon.addHistogram ( new TH1F ("m", ";M [GeV] ;N_{Events}", 100, 0, 40 ) );

  // Delta R
  mon.addHistogram( new TH1F( "DR", ";#DeltaR ;N_{events}", 100, 0, 7) );

  // Delta phi
  mon.addHistogram( new TH1F( "dphi", ";|#Delta#phi| ;N_{Events}", 100, 0, TMath::Pi()) );

  // Delta eta
  mon.addHistogram( new TH1F( "dEta", ";|#Delta#eta| ;N_{Events}", 100, 0, 10 ) );

  // HT
  mon.addHistogram( new TH1F( "HT", ";H_{T} ;N_{Events}", 100, 0, 1000) );
  mon.addHistogram( new TH1F( "HTWjets", ";H_{T} ;N_{Events}", 100, 0, 5000) );

  // MT W
  mon.addHistogram( new TH1F( "MT", ";M_{T} ;N_{Events}", 100, 0, 500) );

  // <DR> vs pt 
  mon.addHistogram( new TProfile ("DR_vs_pt", ";p_{T} [GeV] ;<#DeltaR>", 100, 0, 500) );

  // b-tag discriminator
  mon.addHistogram( new TH1F( "btag", ";b-tag ;N_{Events}", 100, 0, 1) );

  // MVA BDT
  mon.addHistogram( new TH1F( "BDT", ";BDT ;N_{Events}", 100, -1, 1) );

  
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
  int n_event_SR     = 0;
  int n_event_njet   = 0;
  
  // Variables for TTJets 
  int nTot   = 0; 
  int nFilt1 = 0;
  int nFilt4 = 0;
  int nFilt5 = 0;
  
 
  //####################################################################################################################//
  //###########################################           TMVAReader         ###########################################//
  //####################################################################################################################//
  std::string chpath = "WhAnalysis/";

  TMVAReader SR0adaReader;
  SR0adaReader.InitTMVAReader();
  std::string SR0ada_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA_BDT_ada.weights.xml";
  SR0adaReader.SetupMVAReader("SR0adaClass", SR0ada_xml_path);

  TMVAReader SR0gradReader;
  SR0gradReader.InitTMVAReader();
  std::string SR0grad_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA_BDT_grad.weights.xml";
  SR0gradReader.SetupMVAReader("SR0gradClass", SR0grad_xml_path);

 
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

      float lheHT = ev.lheHt;
      if(isMC_WJets && !(isMC_WJets_HTbin))
	{
	  if (lheHT >= 70.0) continue; // reject the event
	}

     
      //############################################################//
      //################## GENERATOR LEVEL #########################//
      //###########################################################//
      
      // Create new object vectors
      std::vector<TLorentzVector> vec_genlep; // store the lepton
      std::vector<TLorentzVector> vec_genqg;  // store the quarks/gluons
      std::vector<TLorentzVector> vec_genb;   // store the b quarks
      std::vector<TLorentzVector> vec_H;      // store the Higgs
      std::vector<TLorentzVector> vec_A;      // store the A
      std::vector<TLorentzVector> vec_W;      // store the W
      std::vector<TLorentzVector> vec_bb1;    // store the bb pair from mother A1
      std::vector<TLorentzVector> vec_bb2;    // store the bb pair from mother A2

      TLorentzVector met_gen;                 // Missing transverse energy vector

      // LOOP OVER MC PARTICLES at GENERATOR LEVEL
      for (int imc=0; imc<ev.nmcparticles; imc++)
	{
	  /* if(verbose && iev <3) 
	    {
	      std::cout << " imcparticle " << imc << " : is a " << ev.mc_id[imc] <<",has a mom:"<<ev.mc_mom[imc]<< " , and has a mother at: " << ev.mc_momidx[imc]  <<"has status   "<<ev.mc_status[imc] << "  and has a 4-vector p = (" << ev.mc_en[imc] << ", " << ev.mc_px[imc] << ", " << ev.mc_py[imc] << ", " << ev.mc_pz[imc] << " ) " << std::endl;
	      }*/

	  if(ev.mc_id[imc]==25) 
	    { // Found the Higgs boson
	      TLorentzVector pH;
	      pH.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Make histograms of Higgs boson kinematics:
	      mon.fillHisto("pt", "gen_H", pH.Pt(), 1.0);
	      mon.fillHisto("eta", "gen_H", pH.Eta(), 1.0);
	      mon.fillHisto("phi", "gen_H", pH.Phi(), 1.0);
	      mon.fillHisto("M", "gen_H", pH.M(), 1.0);
	      
	      vec_H.push_back(pH);        
	    }

	  if(abs(ev.mc_id[imc])==24)
	    { // Found the W
	      TLorentzVector pW;
	      pW.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Make histograms of the W boson kinematics
	      mon.fillHisto("pt", "gen_W", pW.Pt(), 1.0);
	      mon.fillHisto("eta", "gen_W", pW.Eta(), 1.0);
	      mon.fillHisto("phi", "gen_W", pW.Phi(), 1.0);
	      mon.fillHisto("M", "gen_W", pW.M(), 1.0);
	      
	      vec_W.push_back(pW);
	    }

	  if(ev.mc_id[imc]==36)
            { // Found the A                                                              
              TLorentzVector pA;
              pA.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

              // Make histograms of the A boson kinematics
	      mon.fillHisto("pt", "gen_A", pA.Pt(), 1.0);
	      mon.fillHisto("eta", "gen_A", pA.Eta(), 1.0);
	      mon.fillHisto("phi", "gen_A", pA.Phi(), 1.0);
	      mon.fillHisto("M", "gen_A", pA.M(), 1.0);
	      
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
      mon.fillHisto("multi", "step0_gen_lepton", vec_genlep.size(), 1.0);
      mon.fillHisto("qgmulti", "step0_gen", vec_genqg.size(), 1.0);

      // Sorting the quarks/gluons and b quarks based on pT in descending order
      vec_genqg = sort_vec_pt(vec_genqg);
      vec_genb  = sort_vec_pt(vec_genb);
      vec_bb1   = sort_vec_pt(vec_bb1);
      vec_bb2   = sort_vec_pt(vec_bb2);

      // Fill histogram with lheHT for the W+Jets background
      mon.fillHisto("HTWJets", "gen", lheHT, weight); 
	
      // Fill histograms with bb-pair
      if(vec_bb1.size()>1)
	{
	  TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	 
	  mon.fillHisto("pt", "step0_gen_bb1", pbb1.Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_bb1", pbb1.Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_bb1", pbb1.Phi(), 1.0);
	  mon.fillHisto("M", "step0_gen_bb1", pbb1.M(), 1.0);
	}
      
      if(vec_bb2.size()>1)
	{
	  TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	  
	  mon.fillHisto("pt", "step0_gen_bb2", pbb2.Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_bb2", pbb2.Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_bb2", pbb2.Phi(), 1.0);
	  mon.fillHisto("M", "step0_gen_bb2", pbb2.M(), 1.0);
	}

      if(vec_genqg.size()>0)
	{
	  mon.fillHisto("pt", "step0_gen_qg1", vec_genqg[0].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_qg1", vec_genqg[0].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_qg1", vec_genqg[0].Phi(), 1.0);
	}

      if(vec_genqg.size()>1)
	{
	  mon.fillHisto("pt", "step0_gen_qg2", vec_genqg[1].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_qg2", vec_genqg[1].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_qg2", vec_genqg[1].Phi(), 1.0);
	}

      if(vec_genqg.size()>2)
	{
	  mon.fillHisto("pt", "step0_gen_qg3", vec_genqg[2].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_qg3", vec_genqg[2].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_qg3", vec_genqg[2].Phi(), 1.0);
	}

      if(vec_genqg.size()>3)
	{
	  mon.fillHisto("pt", "step0_gen_qg4", vec_genqg[3].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_qg4", vec_genqg[3].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_qg4", vec_genqg[3].Phi(), 1.0);
	}

      if(vec_genb.size()>0)
	{
	  mon.fillHisto("pt", "step0_gen_b1", vec_genb[0].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_b1", vec_genb[0].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_b1", vec_genb[0].Phi(), 1.0);
	}

      if(vec_genb.size()>1)
	{
	  mon.fillHisto("pt", "step0_gen_b2", vec_genb[1].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_b2", vec_genb[1].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_b2", vec_genb[1].Phi(), 1.0);
	}

      if(vec_genb.size()>2)
	{
	  mon.fillHisto("pt", "step0_gen_b3", vec_genb[2].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_b3", vec_genb[2].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_b3", vec_genb[2].Phi(), 1.0);
	}

      if(vec_genb.size()>3)
	{
	  mon.fillHisto("pt", "step0_gen_b4", vec_genb[3].Pt(), 1.0);
	  mon.fillHisto("eta", "step0_gen_b4", vec_genb[3].Eta(), 1.0);
	  mon.fillHisto("phi", "step0_gen_b4", vec_genb[3].Phi(), 1.0);
	}

      // Fill histograms of angular distance DR
      if(vec_bb1.size()>1)
	{
	  float DRbb1 = getDeltaR(vec_bb1[0], vec_bb1[1]);
	  mon.fillHisto("DR", "step0_gen_bb1", DRbb1, 1.0);
	      
	  TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	  mon.fillProfile("DR_vs_pt", "step0_gen_bb1", pbb1.Pt(), DRbb1, 1.0);
	}

      if(vec_bb2.size()>1)
	{
	  float DRbb2 = getDeltaR(vec_bb2[0], vec_bb2[1]);
	  mon.fillHisto("DR", "step0_gen_bb2", DRbb2, 1.0);

	  TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	  mon.fillProfile("DR_vs_pt", "step0_gen_bb2", pbb2.Pt(), DRbb2, 1.0);
	}

      if(vec_bb1.size()>1 && vec_bb2.size()>1)
	{
	  TLorentzVector vec_A1 = vec_bb1[0] + vec_bb1[1];
	  TLorentzVector vec_A2 = vec_bb2[0] + vec_bb2[1];
       
	  float DRAA = getDeltaR(vec_A1, vec_A2);
	  mon.fillHisto("DR", "step0_gen_AA", DRAA, 1.0);
	}
      

      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is one lepton       
      if(vec_genlep.size() >= 1)
	{
	  // Fill the multiplicity histograms
	  mon.fillHisto("multi", "step1_gen_lepton", vec_genlep.size(), 1.0);
	  mon.fillHisto("qgmulti", "step1_gen", vec_genqg.size(), 1.0);
			
	  // Fill histograms of angular distance DR
	  if(vec_bb1.size()>1)
	    {
	      float DRbb1 = getDeltaR(vec_bb1[0], vec_bb1[1]);
	      mon.fillHisto("DR", "step1_gen_bb1", DRbb1, 1.0);
	      
	      TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	      mon.fillProfile("DR_vs_pt", "step1_gen_bb1", pbb1.Pt(), DRbb1, 1.0);
	    }

	  if(vec_bb2.size()>1)
	    {
	      float DRbb2 = getDeltaR(vec_bb2[0], vec_bb2[1]);
	      mon.fillHisto("DR", "step1_gen_bb2", DRbb2, 1.0);

	      TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	      mon.fillProfile("DR_vs_pt", "step1_gen_bb2", pbb2.Pt(), DRbb2, 1.0);
	    }

	  if(vec_bb1.size()>1 && vec_bb2.size()>1)
	    {
	      TLorentzVector vec_A1 = vec_bb1[0] + vec_bb1[1];
	      TLorentzVector vec_A2 = vec_bb2[0] + vec_bb2[1];
       
	      float DRAA = getDeltaR(vec_A1, vec_A2);
	      mon.fillHisto("DR", "step1_gen_AA", DRAA, 1.0);
	    }
	}

	  
      //### STEP 2 ###//
      // If there is at least one lepton and at least three quarks/gluons                                                                                     
      if(vec_genlep.size()>0 && vec_genqg.size()>2)
	{
	  // Fill the histograms of qj kinematics
	  mon.fillHisto("pt", "step2_gen_qg1", vec_genqg[0].Pt(), 1.0);
	  mon.fillHisto("eta", "step2_gen_qg1", vec_genqg[0].Eta(), 1.0);
	  mon.fillHisto("phi", "step2_gen_qg1", vec_genqg[0].Phi(), 1.0);
	      
	  mon.fillHisto("pt", "step2_gen_qg2", vec_genqg[1].Pt(), 1.0);
	  mon.fillHisto("eta", "step2_gen_qg2", vec_genqg[1].Eta(), 1.0);
	  mon.fillHisto("phi", "step2_gen_qg2", vec_genqg[1].Phi(), 1.0);

	  mon.fillHisto("pt", "step2_gen_qg3", vec_genqg[2].Pt(), 1.0);
	  mon.fillHisto("eta", "step2_gen_qg3", vec_genqg[2].Eta(), 1.0);
	  mon.fillHisto("phi", "step2_gen_qg3", vec_genqg[2].Phi(), 1.0);
	    
	  if(vec_genqg.size() > 3)
	    {
	      mon.fillHisto("pt", "step2_gen_qg4", vec_genqg[3].Pt(), 1.0);
	      mon.fillHisto("eta", "step2_gen_qg4", vec_genqg[3].Eta(), 1.0);
	      mon.fillHisto("phi", "step2_gen_qg4", vec_genqg[3].Phi(), 1.0);
	    }

	  // Fill the lepton kinematics histograms
	  mon.fillHisto("pt", "gen_lepton", vec_genlep[0].Pt(), 1.0);
	  mon.fillHisto("eta", "gen_lepton", vec_genlep[0].Eta(), 1.0);
	  mon.fillHisto("phi", "gen_lepton", vec_genlep[0].Phi(), 1.0);
	  
	  // Fill the histogram of HT q/g
	  float HTqg = HTjets(vec_genqg);
	  mon.fillHisto("HT", "gen", HTqg, 1.0);
	  
	  // Fill the histogram of MET pT
	  mon.fillHisto("pt", "gen_MET", met_gen.Pt(), 1.0);

	  // Fill the histogram of W transverse mass
	  float MTW_gen = getMT(vec_genlep[0], met_gen);
	  mon.fillHisto("MT", "gen_recoW", MTW_gen, 1.0);

	  // Define the hadronic system
	  TLorentzVector phad_gen = vec_genqg[0] + vec_genqg[1] + vec_genqg[2];
	  if(vec_genqg.size()>3) phad_gen += vec_genqg[3];

	  // Fill the kinematics histograms of hadronic system
	  mon.fillHisto("pt", "gen_recoH", phad_gen.Pt(), 1.0);
	  mon.fillHisto("eta", "gen_recoH", phad_gen.Eta(), 1.0);
	  mon.fillHisto("phi", "gen_recoH", phad_gen.Phi(), 1.0);
	  mon.fillHisto("M", "gen_recoH", phad_gen.M(), 1.0);

	  // Define the leptonic system
	  TLorentzVector plep_gen = vec_genlep[0] + met_gen;
	  mon.fillHisto("pt", "gen_recoW", plep_gen.Pt(), 1.0);
	  mon.fillHisto("eta", "gen_recoW", plep_gen.Eta(), 1.0);
	  mon.fillHisto("phi", "gen_recoW", plep_gen.Phi(), 1.0);

	  // Fill Delta phi WH histogram
	  float dphiWH_gen = fabs(plep_gen.DeltaPhi(phad_gen));
	  mon.fillHisto("dphi", "gen_WH", dphiWH_gen, 1.0);

	  // Fill Delta phi of MET and lepton
	  float dphiMETlep_gen = fabs(met_gen.DeltaPhi(vec_genlep[0]));
	  mon.fillHisto("dphi", "gen_met-lepton", dphiMETlep_gen, 1.0);

	  
	} // END OF GENERATOR LEVEL ANALYSIS 



      
	    
      //######################################################################//
      //##########################  DETECTOR LEVEL ##########################//
      //####################################################################//

      // Create new object vectors
      std::vector<TLorentzVector> vec_en;                           // store the electron
      std::vector<TLorentzVector> vec_mn;                           // store the muon
      std::vector<TLorentzVector> vec_lepton;                       // store the lepton
      std::vector<TLorentzVector> vec_jets;                         // store all the jets based on their b-tag
      std::vector<TLorentzVector> vec_bjets;                        // store the b-tagged jets based on their pt
      std::vector<TLorentzVector> vec_btag;                         // store the b-tagged jets based on their b-tag
      std::vector<TLorentzVector> vec_untag;                        // store the untagged jets based on their b-tag
      std::vector<std::pair<TLorentzVector, float>> vec_untagbtag;  // store the untagged jets with their b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_jetsbtag;   // store the jets with its b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_bjetsbtag;  // store the b tagged jets with their b-tag index

      
      // Define MET
      TLorentzVector met;
      met.SetPtEtaPhiM(ev.met_pt, 0, ev.met_phi, 0);

      // Define the bb pair from the same mother
      TLorentzVector bb;


      // LOOP OVER DETECTOR PARTICLES

      // Muon
      for(int i=0; i<ev.mn; i++)
	{
	  TLorentzVector pmn;
	  pmn.SetPxPyPzE(ev.mn_px[i], ev.mn_py[i], ev.mn_pz[i], ev.mn_en[i]);

	  // Apply acceptance cuts
	  if(pmn.Pt()<20 || abs(pmn.Eta())>2.4) continue;

	  // Apply ID and Isolation
	  if(ev.mn_passId[i] && ev.mn_passIso[i])
	    {
	      vec_mn.push_back(pmn);
	      vec_lepton.push_back(pmn);
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
	      vec_lepton.push_back(pen);
	    }
	}


      // Fill the multiplicity histograms for lepton
      mon.fillHisto("multi", "step0_lepton", vec_lepton.size(), weight);
      mon.fillHisto("multi", "step0_en", vec_en.size(), weight);
      mon.fillHisto("multi", "step0_mn", vec_mn.size(), weight);


      // Jets and cross-cleaning
      int nHF  = 0;
      int nHFc = 0;
      
      for(int i=0; i<ev.jet; i++)
	{
	  bool overlap = false;

	  TLorentzVector pjet;
	  pjet.SetPxPyPzE(ev.jet_px[i], ev.jet_py[i], ev.jet_pz[i], ev.jet_en[i]);

	  // Apply acceptance cuts
	  if(pjet.Pt()<20 || fabs(pjet.Eta())>2.4) continue;

	  // Jets id 
	  if (!ev.jet_PFTight[i]) continue;

	  for(int mn_count=0; mn_count<vec_mn.size(); mn_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_mn[mn_count]);

	      // Fill the histogram of DR(jet-muon) before cross cleaning
	      mon.fillHisto("DR", "step0_jet-mn", DRjet, weight);

	      if(vec_lepton.size()>1) mon.fillHisto("DR", "step1_jet-mn", DRjet, weight);
		
	      if(DRjet<0.4) overlap = true;
		
	      else
		{
		  mon.fillHisto("DR", "step0_jet-mn_cc", DRjet, weight);

		  if(vec_lepton.size()>1) mon.fillHisto("DR", "step1_jet-mn_cc", DRjet, weight);
		}
	    }
      
	  
	  for(int en_count=0; en_count<vec_en.size(); en_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_en[en_count]);

	      // Fill the histogram of DR(jet-electron) before cross cleaning
	      mon.fillHisto("DR", "step0_jet-en", DRjet, weight);

	      if(vec_lepton.size()>1) mon.fillHisto("DR", "step1_jet-en", DRjet, weight);

	      if(DRjet<0.4) overlap = true;
	
	      else
		{
		  mon.fillHisto("DR", "step0_jet-en_cc", DRjet, weight);

		  if(vec_lepton.size()>1) mon.fillHisto("DR", "step1_jet-en_cc", DRjet, weight);
		}
	    }
	  
	  if(overlap) continue;
	    
	  vec_jetsbtag.push_back(std::make_pair(pjet,ev.jet_btag1[i]));
	  // Loose working point 
	  if(ev.jet_btag1[i]>0.0532)
	    { // Found the b-jets and store them in a vector with their corresponding b-tag index
	      vec_bjetsbtag.push_back(std::make_pair(pjet, ev.jet_btag1[i]));
	    }

	  else
	    { // Store the untagged jets with their corresponding b-tag index
	      vec_untagbtag.push_back(std::make_pair(pjet,ev.jet_btag1[i]));
	    }

	  if(isMC_ttbar)
	    {
	      bool isMatched = false;

	      for(int imc=0; imc<ev.nmcparticles; imc++)
		{
		  if(ev.mc_status[imc]==62) continue;     // 6*:particles produced by beam-remnant treatment/62:outgoing subprocess particle with primordial kT included
		  
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
		  if(!isMatched) nHF++;
		}
	      else if(abs(ev.jet_partonFlavour[i])==4)
		{
		  nHFc++;
		}	
	    }   
	} // End of jet loop

      // Split inclusive TTJets POWHEG sample into tt+bb, tt+cc, tt+light
      if(isMC_ttbar)
	{	  
	  nTot++;

	  if(mctruthmode==5)
	    {
	      if(!(nHF>0)) {continue;}
	      printf("tt+bb event detected\n");
	      nFilt5++;
	    }
	  
	  if(mctruthmode==4)
	    {
	      if(!(nHFc>0 && nHF==0)) {continue;}
	      printf("tt+cc event detected\n");
	      nFilt4++;
	    }
	  
	  if(mctruthmode==1)
	    {
	      if(nHF>0 || (nHFc>0 && nHF == 0)) continue;
	      printf("tt+light event detected\n");
	      nFilt1++;
	    }
	}

      // Sort the jets based on their b-tag
      std::sort(vec_jetsbtag.begin() , vec_jetsbtag.end(),  sortPair);
      std::sort(vec_bjetsbtag.begin(), vec_bjetsbtag.end(), sortPair);
      std::sort(vec_untagbtag.begin(), vec_untagbtag.end(), sortPair);

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
      
      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "step0", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "step0_b", vec_btag.size(), weight);
      mon.fillHisto("jetmulti", "step0_untag", vec_untag.size(), weight);

      // Scalar sum of jets
      float HT = HTjets(vec_jets);
      mon.fillHisto("HTWJets", "det", HT, weight); 

       
      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is one lepton
      if(vec_lepton.size()<1) continue;

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
      mon.fillHisto("multi", "step1_lepton_", vec_lepton.size(), weight);
      mon.fillHisto("multi", "step1_en", vec_en.size(), weight);
      mon.fillHisto("multi", "step1_mn", vec_mn.size(), weight);
      
      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "step1", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "step1_b", vec_btag.size(), weight);
      mon.fillHisto("jetmulti", "step1_untag", vec_untag.size(), weight);
      

           
      //### STEP 2 ###//
     
      // If there is at least three jets
      if(vec_jets.size()<3) continue;

      n_event_jet++;

      mon.fillHisto("eventflow", "histo", 2, 1);
      mon.fillHisto("eventflow", "histo_weighted", 2, weight);

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 2, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 2, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 2, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 2, weight);
	}
       
            
      //### STEP 3 ###//
      
      // If there is at least three b-jets
      if(vec_btag.size() < 3) continue;

      // b tagged jets working point : first b-jet tight WP, second b-jet medium WP, third (and fourth) b-jet loose WP 
      if( vec_bjetsbtag[0].second < 0.7476 || vec_bjetsbtag[1].second < 0.3040) continue;

      for (const auto& pair : vec_bjetsbtag)
	{
	  vec_bjets.push_back(pair.first);
	}

      // Sort the b-jets based on pT in descending order
      vec_bjets = sort_vec_pt(vec_bjets);

      n_event_bjet++;

      mon.fillHisto("eventflow", "histo", 3, 1);
      mon.fillHisto("eventflow", "histo_weighted", 3, weight);

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 3, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 3, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 3, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 3, weight);
	}

      //### STEP 4 ###//
      if(met.Pt() < 25 || getMT(vec_lepton[0], met) < 50) continue;

      n_event_metcut++;

      mon.fillHisto("eventflow", "histo", 4, 1);
      mon.fillHisto("eventflow", "histo_weighted", 4, weight);

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 4, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 4, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 4, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 4, weight);
	}

      
      // Missing transverse energy
      float MET = met.Pt();

      // W transverse mass
      float MTW = getMT(vec_lepton[0], met);
        
      // Lepton kinematics
      float lepton_pt = vec_lepton[0].Pt();
      float lepton_eta = vec_lepton[0].Eta();
      float lepton_phi = vec_lepton[0].Phi();
     
      // b-jet1 and b-jet2 kinematics
      float bjet1_pt  = vec_bjets[0].Pt();
      float bjet1_eta = vec_bjets[0].Eta();
      float bjet1_phi = vec_bjets[0].Phi();
      float bjet1_M   = vec_bjets[0].M();

      float bjet2_pt  = vec_bjets[1].Pt();
      float bjet2_eta = vec_bjets[1].Eta();
      float bjet2_phi = vec_bjets[1].Phi();
      float bjet2_M   = vec_bjets[1].M();

      float bjet3_pt  = vec_bjets[2].Pt();
      float bjet3_eta = vec_bjets[2].Eta();
      float bjet3_phi = vec_bjets[2].Phi();
      float bjet3_M   = vec_bjets[2].M();

      float bjet4_pt  = vec_bjets[3].Pt();
      float bjet4_eta = vec_bjets[3].Eta();
      float bjet4_phi = vec_bjets[3].Phi();
      float bjet4_M   = vec_bjets[3].M();
     

      // Define the leptonic system
      TLorentzVector plepsys = vec_lepton[0] + met;
      float W_pt = plepsys.Pt();

      // Define the hadronic system
      TLorentzVector phadb;

      if(vec_btag.size() > 3) phadb = vec_btag[0] + vec_btag[1] + vec_btag[2] + vec_btag[3];
	
      else if(vec_btag.size() == 3 && vec_untag.size() > 0)  phadb = vec_btag[0] + vec_btag[1] + vec_btag[2] + vec_untag[0];

      else phadb = vec_btag[0] + vec_btag[1] + vec_btag[2];

      // Higgs kinematics
      float H_M, H_pt, H_eta ;  
      H_M  = phadb.M();
      H_pt = phadb.Pt();
      H_eta = phadb.Eta(); 
         
      // WH system variables
      float dphiWH   = fabs(vec_lepton[0].DeltaPhi(phadb));
      float detaWH   = fabs(plepsys.Eta()-phadb.Eta());
      float DRWH     = getDeltaR(plepsys, phadb);
      float WH_M     = (plepsys + phadb).M();
      float dphiHMet = fabs(met.DeltaPhi(phadb));

      // Delta phi of MET and lepton
      float dphiMetLepton = fabs(met.DeltaPhi(vec_lepton[0]));
     
      // Minimum Delta phi of MET and jet
      float dphiMetJetMin = getDeltaPhiMin(met, vec_jets);

      // DeltaR bb average
      float DRbbav = getAverageDeltaR(vec_btag);;

      // Minimum DeltaMass of bb pairs
      float DeltaMassbbMin;
      if (vec_btag.size() > 3)
	{
	  DeltaMassbbMin = getDeltaMassMin(vec_btag, 4);
	}
	  
      else if (vec_btag.size() == 3) 
	{
	  if (vec_untag.size() > 0) 
	    {
	      vector<TLorentzVector> vec_b = vec_btag;
	      vec_b.push_back(vec_untag[0]);
	      DeltaMassbbMin = getDeltaMassMin(vec_b, 4);
	    } 
	  else 
	    {
	      DeltaMassbbMin = getDeltaMassMin(vec_btag, 3);
	    }
	}
      
      // Mass of bb pair with minimum DR and one jet
      float mbbuntag;
      if(vec_untag.size() > 0)
	{
	  mbbuntag = (getbbwithDRMin(vec_bjets) + vec_untag[0]).M();
	  mon.fillHisto("M", "bbj", mbbuntag, weight);
	}
      else
	{
	  mbbuntag = getmbbb(vec_btag);
	  mon.fillHisto("M", "bbb", mbbuntag, weight);
	}
	  
      // Multiplicity of jets
      float Njets = vec_jets.size();

      // b-tag discriminator for 3 most energetic b-jets
      float btag1 = vec_bjetsbtag[0].second;
      float btag2 = vec_bjetsbtag[1].second;
      float btag3 = vec_bjetsbtag[2].second;

      // Invariant mass of bb pair from the same mother
      PairResult result = getPair(vec_btag);
      float bb_M;

      if(vec_btag.size() == 3) bb_M = result.pair1.M() ;
      else bb_M = (result.pair1.M() + result.pair2.M())/2.0;


      // Count the events with exactly 3 jets and more than 3 jets
      if(vec_jets.size() == 3)
	{
	  n_event_njet++ ;
	  mon.fillHisto("Events", "jets", 0, weight); 
	}
      else if(vec_jets.size() > 3)
	{
	  n_event_njet++ ;
	  mon.fillHisto("Events", "jets", 1, weight); 
	}
	
      
      // Fill the histograms
      mon.fillHisto("pt"      , "H"          , H_pt               , weight);
      mon.fillHisto("pt"      , "W"          , W_pt               , weight);
      mon.fillHisto("pt"      , "MET"        , MET                , weight);
      mon.fillHisto("pt"      , "lepton"     , lepton_pt          , weight);
      mon.fillHisto("pt"      , "bjet1"      , bjet1_pt           , weight);
      mon.fillHisto("pt"      , "bjet2"      , bjet2_pt           , weight);
      mon.fillHisto("pt"      , "bjet3"      , bjet3_pt           , weight);
      mon.fillHisto("pt"      , "bjet4"      , bjet4_pt           , weight);
      mon.fillHisto("pt"      , "lepton"     , lepton_pt          , weight);
      mon.fillHisto("HT"      , "jets"       , HT                 , weight);
      mon.fillHisto("eta"     , "H"          , H_eta              , weight);
      mon.fillHisto("eta"     , "WH"         , detaWH             , weight);
      mon.fillHisto("eta"     , "lepton"     , lepton_eta         , weight);
      mon.fillHisto("eta"     , "bjet1"      , bjet1_eta          , weight);
      mon.fillHisto("eta"     , "bjet2"      , bjet2_eta          , weight);
      mon.fillHisto("eta"     , "bjet3"      , bjet3_eta          , weight);
      mon.fillHisto("eta"     , "bjet4"      , bjet4_eta          , weight);
      mon.fillHisto("MT"      , "W"          , MTW                , weight);
      mon.fillHisto("M"       , "H"          , H_M                , weight);
      mon.fillHisto("M"       , "WH"         , WH_M               , weight);
      mon.fillHisto("M"       , "dMbbmin"    , DeltaMassbbMin     , weight);
      mon.fillHisto("M"       , "bbuntag"    , mbbuntag           , weight);
      mon.fillHisto("M"       , "bb"         , bb_M               , weight);
      mon.fillHisto("dEta"    , "WH"         , detaWH             , weight);
      mon.fillHisto("dphi"    , "WH"         , dphiWH             , weight);
      mon.fillHisto("dphi"    , "MET-lepton" , dphiMetLepton      , weight);
      mon.fillHisto("dphi"    , "METjet-min" , dphiMetJetMin      , weight);
      mon.fillHisto("dphi"    , "H-MET"      , dphiHMet           , weight);
      mon.fillHisto("DR"      , "WH"         , DRWH               , weight);
      mon.fillHisto("DR"      , "bbaver"     , DRbbav             , weight);
      mon.fillHisto("jetmulti", "step4"      , vec_jets.size()    , weight);
      mon.fillHisto("jetmulti", "step4_b"    , vec_bjets.size()   , weight);
      mon.fillHisto("jetmulti", "step4_untag", vec_untag.size()   , weight);
      mon.fillHisto("btag"    , "jet1"       , btag1              , weight);
      mon.fillHisto("btag"    , "jet2"       , btag2              , weight);
      mon.fillHisto("btag"    , "jet3"       , btag3              , weight);
      
      
      


      
      //######################################//
      //########## SIGNAL REGIONS ###########//
      //#####################################//

      bool isSR1=false;      
      bool isSR2=false;     

      // for later .......


      //#####################################//
      //############ MVA Handler ############//
      //#####################################//
      
      if(runMVA)
	{
	  myMVAHandler_.getEntry(Njets, HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, WH_M, bb_M, MTW, dphiWH, detaWH, dphiMetLepton, dphiMetJetMin, dphiHMet, DRbbav, DRWH, DeltaMassbbMin, mbbuntag, btag1, btag2, btag3, weight);
	  myMVAHandler_.fillTree();
	}  

    
      //#####################################//
      //############ TMVA Reader ############//
      //#####################################//
      float mvaBDTada  = -10.0;
      float mvaBTDgrad = -10.0;

      mvaBDTada = SR0adaReader.GenReMVAReader(Njets, HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, detaWH, dphiMetJetMin, DRbbav, DeltaMassbbMin, mbbuntag, btag1, btag2, btag3, "SR0adaClass");
      mvaBTDgrad = SR0gradReader.GenReMVAReader(Njets, HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, detaWH, dphiMetJetMin, DRbbav, DeltaMassbbMin, mbbuntag, btag1, btag2, btag3, "SR0gradClass");
      mon.fillHisto("BDT", "SR0_ada", mvaBDTada, weight);
      mon.fillHisto("BDT", "SR0_grad", mvaBTDgrad, weight);
      
    
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

      //printf("From total = %i TTbar events ,  Found %.2f (filt1) , %.2f (filt4) , %.2f (filt5) \n\n", nTot, cFilt1, cFilt4, cFilt5);
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

// Sorting
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
    float dphiMin1 = fabs(vec_1.DeltaPhi(vec_2[0]));
    float dphiMin2 = fabs(vec_1.DeltaPhi(vec_2[1]));
    float dphiMin3 = fabs(vec_1.DeltaPhi(vec_2[2]));
   
    std::vector<float> dphivalues = {dphiMin1, dphiMin2, dphiMin3};
    
    auto minElem = std::min_element(dphivalues.begin(), dphivalues.end());
    
    return *minElem;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


