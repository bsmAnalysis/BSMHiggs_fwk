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
#include <utility>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <string>
#include <unistd.h>
using namespace std;

float getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2);
float getAverageDeltaR(vector<TLorentzVector> vec);
float getDeltaMassMin(vector<TLorentzVector> vec);
TLorentzVector getbbwithDRMin(vector<TLorentzVector> vec);
float getDeltaPhiMin(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2);
float getMT(TLorentzVector vec_1, TLorentzVector vec_2);
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec);
float HTjets(std::vector<TLorentzVector> vec);
bool sortPair(std::pair<TLorentzVector, float> vec_1, std::pair<TLorentzVector, float> vec_2); 

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

  bool isMC    = runProcess.getParameter<bool>("isMC");
  double xsec  = runProcess.getParameter<double>("xsec");
  double nevts = runProcess.getParameter<double>("nevts");

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

  bool isMC_ttbar_had   = isMC && (string(url.Data()).find("TTToHadronic") !=string::npos);
  bool isMC_ttbar_semi  = isMC && (string(url.Data()).find("TTToSemiLeptonic") !=string::npos);
  bool isMC_ttbar_2l2nu = isMC && (string(url.Data()).find("TTTo2L2Nu") !=string::npos);
  bool isMC_WJets_HTbin = isMC && (string(url.Data()).find("WJets") !=string::npos);
  bool isMC_QCD_MuEnr   = isMC && (string(url.Data()).find("MuEnr") !=string::npos);
  bool isMC_QCD_EMEnr   = isMC && (string(url.Data()).find("EMEnr") !=string::npos);
  bool isMC_QCD_bcToE   = isMC && (string(url.Data()).find("bcToE") !=string::npos);
  bool isBackground     = (isMC_ttbar_had || isMC_ttbar_semi || isMC_ttbar_2l2nu || isMC_WJets_HTbin || isMC_QCD_MuEnr || isMC_QCD_EMEnr || isMC_QCD_bcToE);

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
  double Lint = 41.5;
  double Nexp = xsec * 1000 * Lint;
  double weight = Nexp / nevts;
  if(verbose)
    {
      cout << "Cross section:" << xsec << " [pb]" << endl;
      cout << "Weight: " << weight << endl;
      cout << "N expected:" << Nexp << endl; 
    }
 
  //##################################################################################//
  //##########################    INITIATING HISTOGRAMS     ##########################//
  //##################################################################################//
  SmartSelectionMonitor mon;

  
  // Event flow 
  TH1F *h = (TH1F*) mon.addHistogram ( new TH1F ("eventflow", ";;N_{Events}", 7,0,7) );
  h->GetXaxis()->SetBinLabel(1,"Raw");
  h->GetXaxis()->SetBinLabel(2,"1 lepton");
  h->GetXaxis()->SetBinLabel(3,">=3 jets");
  h->GetXaxis()->SetBinLabel(4,">=3 b-tags");
  h->GetXaxis()->SetBinLabel(5,"MET>25 & MTW>50");
  h->GetXaxis()->SetBinLabel(6,"SR1");
  h->GetXaxis()->SetBinLabel(7,"SR2");
  
  // Multiplicity 
  mon.addHistogram ( new TH1F ("multi", ";N_{lepton} ;N_{Events}", 5, 0, 5) );
  mon.addHistogram ( new TH1F ("qgmulti", ";N_{qg} ;N_{Events}", 15, 0, 15) );
  mon.addHistogram ( new TH1F ("jetmulti", ";N_{jet} ;N_{Events}", 15, 0, 15) );

  // Particles/Objects kinematics
  mon.addHistogram ( new TH1F ("pt", ";p_{T} [GeV] ;N_{Events}", 100, 0, 500) );
  mon.addHistogram ( new TH1F ("eta", ";#eta ;N_{Events}", 100, -6, 6));
  mon.addHistogram ( new TH1F ("phi", ";#phi ;N_{Events}", 100, -TMath::Pi(), TMath::Pi()) );
  mon.addHistogram ( new TH1F ("M", ";M [GeV] ;N_{Events}", 100, 0, 800) );

  // Delta R
  mon.addHistogram( new TH1F( "DR", ";#DeltaR ;N_{events}", 100, 0, 7) );

  // Delta phi
  mon.addHistogram( new TH1F( "dphi", ";|#Delta#phi| ;N_{Events}", 100, 0, TMath::Pi()) );

  // HT
  mon.addHistogram( new TH1F( "HT", ";H_{T} ;N_{Events}", 100, 0, 1000) );

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


  //------Event Counters-------//
  int n_event_lepton = 0;
  int n_event_jet    = 0;
  int n_event_bjet   = 0;
  int n_event_metcut = 0;
  int n_event_SR     = 0;

  //####################################################################################################################//
  //###########################################           TMVAReader         ###########################################//
  //####################################################################################################################//
  std::string chpath = "WhAnalysis/old/";

  TMVAReader SR0Reader;
  SR0Reader.InitTMVAReader();
  std::string SR0_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA0_BDT.weights.xml";
  SR0Reader.SetupMVAReader("SR0Class", SR0_xml_path);

  TMVAReader SR1Reader;
  SR1Reader.InitTMVAReader();
  std::string SR1_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA1_BDT.weights.xml";
  SR1Reader.SetupMVAReader("SR1Class", SR1_xml_path);

  TMVAReader SR2Reader;
  SR2Reader.InitTMVAReader();
  std::string SR2_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA2_BDT.weights.xml";
  SR2Reader.SetupMVAReader("SR2Class", SR2_xml_path);

 
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

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 0, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 0, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 0, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 0, weight);
	}

      else if (isBackground)
	{
	  mon.fillHisto("eventflow", "Background", 0, weight);

	  if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 0, weight);
	  if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 0, weight);
	  if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 0, weight);
	  if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 0, weight);
	  if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 0, weight);
	  if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 0, weight);
	  if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 0, weight);
	}
            
      if(!isMC && duplicatesChecker.isDuplicate(ev.run, ev.lumi, ev.event))
	{
	  nDuplicates++;
	  cout << "nDuplicates: " << nDuplicates << endl;
	  continue;
	}

     
      //############################################################//
      //################## GENERATOR LEVEL #########################//
      //###########################################################//
      
      // Create new object vectors
      std::vector<TLorentzVector> vec_genlep; // store the lepton
      std::vector<TLorentzVector> vec_genqg;  // store the quarks/gluons
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
      mon.fillHisto("multi", "gen_lepton", vec_genlep.size(), weight);
      mon.fillHisto("qgmulti", "gen", vec_genqg.size(), weight);

      // Sorting the quarks/gluons based on pT in descending order
      vec_genqg = sort_vec_pt(vec_genqg);

      // Fill histograms with bb-pair
      if(vec_bb1.size()>1)
	{
	  TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	 
	  mon.fillHisto("pt", "gen_bb1", pbb1.Pt(), weight);
	  mon.fillHisto("eta", "gen_bb1", pbb1.Eta(), weight);
	  mon.fillHisto("phi", "gen_bb1", pbb1.Phi(), weight);
	  mon.fillHisto("M", "gen_bb1", pbb1.M(), weight);
	}
      

      if(vec_bb2.size()>1)
	{
	  TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	  
	  mon.fillHisto("pt", "gen_bb2", pbb2.Pt(), weight);
	  mon.fillHisto("eta", "gen_bb2", pbb2.Eta(), weight);
	  mon.fillHisto("phi", "gen_bb2", pbb2.Phi(), weight);
	  mon.fillHisto("M", "gen_bb2", pbb2.M(), weight);
	}

      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is at least one lepton       
      if(vec_genlep.size()>0)
	{
	  // Fill the multiplicity histograms
	  mon.fillHisto("multi", "step1_gen_lepton", vec_genlep.size(), weight);
	  mon.fillHisto("qgmulti", "step1_gen", vec_genqg.size(), weight);
			
	  // Fill histograms of angular distance DR
	  if(vec_bb1.size()>1)
	    {
	      float DRbb1 = getDeltaR(vec_bb1[0], vec_bb1[1]);
	      mon.fillHisto("DR", "gen_bb1", DRbb1, weight);
	      
	      TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	      mon.fillProfile("DR_vs_pt", "gen_bb1", pbb1.Pt(), DRbb1, weight);
	    }

	  if(vec_bb2.size()>1)
	    {
	      float DRbb2 = getDeltaR(vec_bb2[0], vec_bb2[1]);
	      mon.fillHisto("DR", "gen_bb2", DRbb2, weight);

	      TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	      mon.fillProfile("DR_vs_pt", "gen_bb2", pbb2.Pt(), DRbb2, weight);
	    }

	  if(vec_bb1.size()>1 && vec_bb2.size()>1)
	    {
	      TLorentzVector vec_A1 = vec_bb1[0] + vec_bb1[1];
	      TLorentzVector vec_A2 = vec_bb2[0] + vec_bb2[1];
       
	      float DRAA = getDeltaR(vec_A1, vec_A2);
	      mon.fillHisto("DR", "gen_AA", DRAA, weight);
	    }
	}

	  
      //### STEP 2 ###//
      // If there is at least one lepton and at least three quarks/gluons                                                                                     
      if(vec_genlep.size()>0 && vec_genqg.size()>2)
	{                                             
	  // Fill the histogram of HT q/g
	  float HTqg = HTjets(vec_genqg);
	  mon.fillHisto("HT", "gen", HTqg, weight);
	  
	  // Fill the histogram of MET pT
	  mon.fillHisto("pt", "gen_MET", met_gen.Pt(), weight);

	  // Fill the histogram of W transverse mass
	  float MTW_gen = getMT(vec_genlep[0], met_gen);
	  mon.fillHisto("MT", "gen_W", MTW_gen, weight);

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
	  mon.fillHisto("pt", "gen_lepton", plep_gen.Pt(), weight);

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
      std::vector<TLorentzVector> vec_en;                           // store the electron
      std::vector<TLorentzVector> vec_mn;                           // store the muon
      std::vector<TLorentzVector> vec_lepton;                       // store the lepton
      std::vector<TLorentzVector> vec_jets;                         // store all the jets based on their b-tag
      std::vector<TLorentzVector> vec_bjets;                        // store the b-tagged jets based on their pt
      std::vector<TLorentzVector> vec_btag;                         // store the b-tagged jets based on their b-tag
      std::vector<TLorentzVector> vec_untag;                        // store the untagged jets based on their b-tag
      std::vector<TLorentzVector> vec_fjets;                        // store the fat jets
      std::vector<std::pair<TLorentzVector, float>> vec_untagbtag;  // store the untagged jets with their b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_jetsbtag;   // store the jets with its b-tag index
      std::vector<std::pair<TLorentzVector, float>> vec_bjetsbtag;  // store the b tagged jets with their b-tag index

      
      // Define MET
      TLorentzVector met;
      met.SetPtEtaPhiM(ev.met_pt, 0, ev.met_phi, 0);


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

	      if(vec_lepton.size()>1)
		{
		  mon.fillHisto("DR", "step1_jet-mn", DRjet, weight);
		}
	      
	      if(DRjet<0.4)
		{
		  overlap = true;
		}
	      else
		{
		  mon.fillHisto("DR", "step0_jet-mn_cc", DRjet, weight);

		  if(vec_lepton.size()>1)
		    {
		      mon.fillHisto("DR", "step1_jet-mn_cc", DRjet, weight);
		    }
		}
	    }
      
	  
	  for(int en_count=0; en_count<vec_en.size(); en_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_en[en_count]);

	      // Fill the histogram of DR(jet-electron) before cross cleaning
	      mon.fillHisto("DR", "step0_jet-en", DRjet, weight);

	      if(vec_lepton.size()>1)
		{
		  mon.fillHisto("DR", "step1_jet-en", DRjet, weight);
		}

	      if(DRjet<0.4)
		{
		  overlap = true;
		}
	      else
		{
		  mon.fillHisto("DR", "step0_jet-en_cc", DRjet, weight);

		  if(vec_lepton.size()>1)
		    {
		      mon.fillHisto("DR", "step1_jet-en_cc", DRjet, weight);
		    }
		}
	    }
	  
	  if(!overlap)
	    {
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
	    } 
	} // End of jet loop

      //Fat-jets and cross-cleaning
      /* int index1=0;
      int index2=0;
      int count=0;
      
      for(int i=0; i<ev.fjet; i++)
	{
	  bool overlap = false;
	  
	  TLorentzVector pfjet;
	  pfjet.SetPxPyPzE(ev.fjet_px[i], ev.fjet_py[i], ev.fjet_pz[i], ev.fjet_en[i]);

	  if(ev.fjet_subjet_count[i]<2) continue;

	  // Apply acceptance cuts
	  if(pfjet.Pt()<20 || fabs(pfjet.Eta())>2.4) continue;

	  for(int mn_count=0; mn_count<vec_mn.size(); mn_count++)
	    {
	      double DRfjet = getDeltaR(pfjet, vec_mn[mn_count]);

	      // Fill the histogram of DR(fatjet-muon) before cross cleaning
	      mon.fillHisto("DR", "fjet-mn", DRfjet, weight);

	      if(DRfjet<0.8) overlap = true;

	      else mon.fillHisto("DR", "fjet-mn_cc", DRfjet, weight);
	    }

	  for(int en_count=0; en_count<vec_en.size(); i++)
	    {
	      double DRfjet = getDeltaR(pfjet, vec_en[en_count]);

	      // Fill the histogram of DR(fatjet-electron) before cross cleaning
	      mon.fillHisto("DR", "fjet-en", DRfjet, weight);

	      if(DRfjet<0.8) overlap = true;

	      else mon.fillHisto("DR", "fjet-en_cc", DRfjet, weight);
		  
	    }

	  if(!overlap)
	    {
	      count ++;
	      if(count==1) index1=i;
	      if(count==2) index2=i;
	  
	      vec_fjets.push_back(pfjet);
	    }
	    } // End of fat-jet loop */

     
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

       
      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is at least one lepton
      if(vec_lepton.size()<1) continue;

      n_event_lepton++;

     if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 1, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 1, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 1, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 1, weight);
	}

      else if (isBackground)
	{
	  mon.fillHisto("eventflow", "Background", 1, weight);

	  if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 1, weight);
	  if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 1, weight);
	  if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 1, weight);
	  if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 1, weight);
	  if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 1, weight);
	  if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 1, weight);
	  if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 1, weight);
	}

      // Fill the multiplicity histograms
      mon.fillHisto("multi", "step1_lepton_", vec_lepton.size(), weight);
      mon.fillHisto("multi", "step1_en", vec_en.size(), weight);
      mon.fillHisto("multi", "step1_mn", vec_mn.size(), weight);

      // Fill lepton kinematics
      mon.fillHisto("pt", "lepton", vec_lepton[0].Pt(), weight);
      mon.fillHisto("eta", "lepton", vec_lepton[0].Eta(), weight);
      mon.fillHisto("phi", "lepton", vec_lepton[0].Phi(), weight);
      
      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "step1", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "step1_b", vec_btag.size(), weight);
      mon.fillHisto("jetmulti", "step1_untag", vec_untag.size(), weight);
      

           
      //### STEP 2 ###//
     
      // If there is at least three jets
      if(vec_jets.size()<3) continue;

      n_event_jet++;

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 2, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 2, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 2, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 2, weight);
	}

      else if (isBackground)
	{
	  mon.fillHisto("eventflow", "Background", 2, weight);

	  if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 2, weight);
	  if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 2, weight);
	  if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 2, weight);
	  if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 2, weight);
	  if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 2, weight);
	  if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 2, weight);
	  if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 2, weight);
	}
       
            
      //### STEP 3 ###//
      
      // If there is at least three b-jets
      if(vec_btag.size() < 3) continue;

      // b tagged jets working point : first and second b-jet medium working point, third (and fourth) b-jet loose working point 
      if( vec_bjetsbtag[0].second < 0.3040 || vec_bjetsbtag[1].second < 0.3040) continue;

      for (const auto& pair : vec_bjetsbtag)
	{
	  vec_bjets.push_back(pair.first);
	}

      // Sort the b-jets based on pT in descending order
      vec_bjets = sort_vec_pt(vec_bjets);

      n_event_bjet++;

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 3, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 3, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 3, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 3, weight);
	}

      else if (isBackground)
	{
	  mon.fillHisto("eventflow", "Background", 3, weight);

	  if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 3, weight);
	  if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 3, weight);
	  if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 3, weight);
	  if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 3, weight);
	  if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 3, weight);
	  if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 3, weight);
	  if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 3, weight);
	}

 
      //### STEP 4 ###//
      if(met.Pt() < 25 || getMT(vec_lepton[0], met) < 50) continue;

      n_event_metcut++;

      if (isSignal)
	{
	  mon.fillHisto("eventflow", "Signal", 4, weight);

	  if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 4, weight);  
	  if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 4, weight);   
	  if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 4, weight);
	}

      else if (isBackground)
	{
	  mon.fillHisto("eventflow", "Background", 4, weight);

	  if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 4, weight);
	  if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 4, weight);
	  if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 4, weight);
	  if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 4, weight);
	  if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 4, weight);
	  if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 4, weight);
	  if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 4, weight);
	}

      
      // Missing transverse energy
      float MET = met.Pt();

      // W transverse mass
      float MTW = getMT(vec_lepton[0], met);
    
      // Define the hadronic system
      TLorentzVector phadb = vec_bjets[0] + vec_bjets[1] + vec_bjets[2];
      if(vec_bjets.size()>3) phadb += vec_bjets[3];

      float H_pt  = phadb.Pt();
      float H_M   = phadb.M();
      float H_eta = phadb.Eta();
      float H_phi = phadb.Phi();

      // Lepton pt
      float lepton_pt = vec_lepton[0].Pt();
     
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
     
      // Scalar sum of jets
      float HT = HTjets(vec_jets);
     
      // Delta phi WH
      float dphiWH = fabs(vec_lepton[0].DeltaPhi(phadb));

      // Define the leptonic system
      TLorentzVector plepsys = vec_lepton[0] + met;
      float W_pt = plepsys.Pt();

      // Delta phi of MET and lepton
      float dphiMetLepton = fabs(met.DeltaPhi(vec_lepton[0]));
     
      // Minimum Delta phi of MET and jet
      float dphiMetJetMin = getDeltaPhiMin(met, vec_jets);

      // DeltaR bb average
      float DRbbav = getAverageDeltaR(vec_bjets);;

      // Minimum DeltaMass of bb pair
      float DeltaMassbbMin = 0.0;
      
      // Mass of bb pair with minimum DR and one jet
      float mbbuntag = 0.0;
	  
      // Multiplicity of jets
      float Njets = vec_jets.size();

      // Multiplicity of lepton
      float Nlepton = vec_lepton.size();

      // b-tag discriminator for 3 most energetic b-jets
      float btag1 = vec_jetsbtag[0].second;
      float btag2 = vec_jetsbtag[1].second;
      float btag3 = vec_jetsbtag[2].second;
      
      // Fill the histograms
      mon.fillHisto("pt"      , "H", H_pt, weight);
      mon.fillHisto("pt"      , "W", W_pt, weight);
      mon.fillHisto("pt"      , "MET", MET, weight);
      mon.fillHisto("pt"      , "lepton", lepton_pt, weight);
      mon.fillHisto("pt"      , "bjet1", bjet1_pt, weight);
      mon.fillHisto("pt"      , "bjet2", bjet2_pt, weight);
      mon.fillHisto("pt"      , "bjet3", bjet3_pt, weight);
      mon.fillHisto("pt"      , "bjet4", bjet4_pt, weight);
      mon.fillHisto("HT"      , "jets", HT, weight);
      mon.fillHisto("eta"     , "H", H_eta, weight);
      mon.fillHisto("MT"      , "W", MTW, weight);
      mon.fillHisto("M"       , "H", H_M, weight);
      mon.fillHisto("M"       , "dMbbmin", DeltaMassbbMin, weight);  
      mon.fillHisto("dphi"    , "WH", dphiWH, weight);
      mon.fillHisto("dphi"    , "met-lepton", dphiMetLepton, weight);
      mon.fillHisto("dphi"    , "metjet-min", dphiMetJetMin, weight);
      mon.fillHisto("jetmulti", "step4", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "step4_b", vec_bjets.size(), weight);
      mon.fillHisto("jetmulti", "step4_untag", vec_untag.size(), weight);
      mon.fillHisto("btag"    , "jet1", btag1, weight);
      mon.fillHisto("btag"    , "jet2", btag2, weight);
      mon.fillHisto("btag"    , "jet3", btag3, weight);

 
      //######################################//
      //########## SIGNAL REGIONS ###########//
      //#####################################//

      bool isSR1=false;      // Nlepton>=1, Njets=3,  Nbjets=3
      bool isSR2=false;      // Nlepton>=1, Njets>=4

      if(vec_jets.size() == 3) isSR1=true;
      if(vec_jets.size() >= 4) isSR2=true;

      if (isSR1)
	{
	  if (isSignal)
	    {
	      n_event_SR++;
	      
	      mon.fillHisto("eventflow", "Signal", 5, weight);

	      if (isMC_Wh)   mon.fillHisto("eventflow", "Wh",   5, weight);  
	      if (isMC_Zh)   mon.fillHisto("eventflow", "Zh",   5, weight);   
	      if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 5, weight);
	    }

	  else if (isBackground)
	    {
	      mon.fillHisto("eventflow", "Background", 5, weight);

	      if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 5, weight);
	      if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 5, weight);
	      if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 5, weight);
	      if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 5, weight);
	      if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 5, weight);
	      if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 5, weight);
	      if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 5, weight);
	    }
	  
	  // Fill the histograms
	  mon.fillHisto("HT"      , "SR1_jets", HT, weight);
	  mon.fillHisto("pt"      , "SR1_H", H_pt, weight);
	  mon.fillHisto("pt"      , "SR1_W", W_pt, weight);
	  mon.fillHisto("pt"      , "SR1_lepton", lepton_pt, weight);
	  mon.fillHisto("pt"      , "SR1_bjet1", bjet1_pt, weight);
	  mon.fillHisto("pt"      , "SR1_bjet2", bjet2_pt, weight);
	  mon.fillHisto("pt"      , "SR1_bjet3", bjet3_pt, weight);
	  mon.fillHisto("pt"      , "SR1_bjet4", bjet4_pt, weight);
	  mon.fillHisto("pt"      , "SR1_MET", MET, weight);
	  mon.fillHisto("M"       , "SR1_H", H_M, weight);
	  mon.fillHisto("M"       , "SR1_dMbbmin", DeltaMassbbMin, weight);
	  mon.fillHisto("MT"      , "SR1_W", MTW, weight);
	  mon.fillHisto("eta"     , "SR1_H", H_eta, weight);
	  mon.fillHisto("dphi"    , "SR1_WH", dphiWH, weight);
	  mon.fillHisto("dphi"    , "SR1_met-lepton", dphiMetLepton, weight);
	  mon.fillHisto("dphi"    , "SR1_metjet-min", dphiMetJetMin, weight);
	  mon.fillHisto("DR"      , "SR1_bbav", DRbbav, weight);
	  mon.fillHisto("jetmulti", "SR1_step4", vec_jets.size(), weight);
	  mon.fillHisto("jetmulti", "SR1_step4_b", vec_bjets.size(), weight);
	  mon.fillHisto("jetmulti", "SR1_step4_untag", vec_untag.size(), weight);
	  mon.fillHisto("multi"   , "SR1_lepton_step4", vec_lepton.size(), weight);
	  mon.fillHisto("btag"    , "SR1_jet1", btag1, weight);
	  mon.fillHisto("btag"    , "SR1_jet2", btag2, weight);
	  mon.fillHisto("btag"    , "SR1_jet3", btag3, weight);
	} 
	
      else if(isSR2)
	{
	  n_event_SR++;

	  if (isSignal)
	    {
	      mon.fillHisto("eventflow", "Signal", 6, weight);

	      if (isMC_Wh)   mon.fillHisto("eventflow","Wh", 6, weight);  
	      if (isMC_Zh)   mon.fillHisto("eventflow","Zh", 6, weight);   
	      if (isMC_VBFh) mon.fillHisto("eventflow", "VBFh", 6, weight);
	    }

	  else if (isBackground)
	    {
	      mon.fillHisto("eventflow", "Background", 6, weight);

	      if (isMC_ttbar_had)    mon.fillHisto("eventflow", "TTToHadronic", 6, weight);
	      if (isMC_ttbar_semi)   mon.fillHisto("eventflow", "TTToSemiLeptonic", 6, weight);
	      if (isMC_ttbar_2l2nu)  mon.fillHisto("eventflow", "TTTo2L2Nu", 6, weight);
	      if (isMC_QCD_MuEnr)    mon.fillHisto("eventflow", "MuEnriched", 6, weight);
	      if (isMC_QCD_EMEnr)    mon.fillHisto("eventflow", "EMEnriched", 6, weight);
	      if (isMC_QCD_bcToE)    mon.fillHisto("eventflow", "bcToE", 6, weight);
	      if (isMC_WJets_HTbin)  mon.fillHisto("eventflow", "WJets_HTbin", 6, weight);
	    }

	  // Define the hadronic system
	  TLorentzVector phad;

	  // Higgs invariant mass and transverse momentum
	  if(vec_bjets.size() == 3)
	    {
	      phad  = vec_bjets[0] + vec_bjets[1] + vec_bjets[2] + vec_untag[0]; 
	      H_M   = phadb.M();
	      H_pt  = phadb.Pt();
	      H_eta = phadb.Eta(); 
	    }
	  
	  else if(vec_bjets.size() > 3)
	    {
	      phad = vec_bjets[0] + vec_bjets[1] + vec_bjets[2] + vec_bjets[3];
	      H_M  = phadb.M();
	      H_pt = phadb.Pt();
	      H_eta = phadb.Eta(); 
	    }	  

	  // Find the invariant mass of bb pair with minimum DeltaR and one untagged jet
	  if(vec_untag.size() > 0) mbbuntag = (getbbwithDRMin(vec_bjets) + vec_untag[0]).M();

	  // Fill the histograms
	  mon.fillHisto("HT"      , "SR2_jets", HT, weight);
	  mon.fillHisto("pt"      , "SR2_H", H_pt, weight);
	  mon.fillHisto("pt"      , "SR2_W", W_pt, weight);
	  mon.fillHisto("pt"      , "SR2_lepton", lepton_pt, weight);
	  mon.fillHisto("pt"      , "SR2_bjet1", bjet1_pt, weight);
	  mon.fillHisto("pt"      , "SR2_bjet2", bjet2_pt, weight);
	  mon.fillHisto("pt"      , "SR2_bjet3", bjet3_pt, weight);
	  mon.fillHisto("pt"      , "SR2_bjet4", bjet4_pt, weight);
	  mon.fillHisto("pt"      , "SR2_MET", MET, weight);
	  mon.fillHisto("eta"     , "SR2_H", H_eta, weight);
	  mon.fillHisto("M"       , "SR2_H", H_M, weight);
	  mon.fillHisto("M"       , "SR2_dMbbmin", DeltaMassbbMin, weight);
	  mon.fillHisto("M"       , "SR2_bbuntag", mbbuntag, weight);
	  mon.fillHisto("MT"      , "SR2_W", MTW, weight);
	  mon.fillHisto("dphi"    , "SR2_WH", dphiWH, weight);
	  mon.fillHisto("dphi"    , "SR2_met-lepton", dphiMetLepton, weight);
	  mon.fillHisto("dphi"    , "SR2_metjet-min", dphiMetJetMin, weight);
	  mon.fillHisto("DR"      , "SR2_bbav", DRbbav, weight);
	  mon.fillHisto("jetmulti", "SR2_step4", vec_jets.size(), weight);
	  mon.fillHisto("jetmulti", "SR2_step4_b", vec_bjets.size(), weight);
	  mon.fillHisto("jetmulti", "SR2_step4_untag", vec_untag.size(), weight);
	  mon.fillHisto("multi"   , "SR2_lepton_step4", vec_lepton.size(), weight);
	  mon.fillHisto("btag"    , "SR2_jet1", btag1, weight);
	  mon.fillHisto("btag"    , "SR2_jet2", btag2, weight);
	  mon.fillHisto("btag"    , "SR2_jet3", btag3, weight);
	}


      //#####################################//
      //############ MVA Handler ############//
      //#####################################//
      
      if(runMVA)
	{
	  myMVAHandler_.getEntry(isSR1, isSR2, Njets, Nlepton, HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, dphiMetLepton, dphiMetJetMin, DRbbav, DeltaMassbbMin, mbbuntag, btag1, btag2, btag3, weight);
	  myMVAHandler_.fillTree();
	}  

    
      //#####################################//
      //############ TMVA Reader ############//
      //#####################################//
      float mvaBDT0 = -10.0;
      float mvaBDT1 = -10.0;
      float mvaBDT2 = -10.0;

      mvaBDT0 = SR0Reader.GenReMVAReader(Njets, HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, dphiMetJetMin, DRbbav, btag1, btag2, btag3, "SR0Class");
      mon.fillHisto("BDT", "SR0", mvaBDT0, weight);
	
      if(isSR1)
	{
	  mvaBDT1 = SR1Reader.GenReMVAReader(HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, dphiMetJetMin, DRbbav, btag1, btag2, btag3, "SR1Class");
	  mon.fillHisto("BDT", "SR1", mvaBDT1, weight);
	}

      else if(isSR2)
	{
	  mvaBDT2 = SR2Reader.GenReMVAReader(Njets, HT, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, dphiMetJetMin, DRbbav, btag1, btag2, btag3, "SR2Class");
	  mon.fillHisto("BDT", "SR2", mvaBDT2, weight);
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
  
// Angular distance DR
float getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2)
{
  float delta_phi;
  float delta_eta;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_eta = vec_1.Eta() - vec_2.Eta();

  return std::sqrt(pow(delta_phi, 2) + pow(delta_eta, 2));
}

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

// Minimum Delta mass
float getDeltaMassMin(vector<TLorentzVector> vec)
{
  float minDeltaMass = std::numeric_limits<double>::max();

  for(size_t i=0; i<vec.size(); i++)
    {
      for(size_t j=i+1; j<vec.size(); j++)
	{
	  TLorentzVector pair1 = vec[i]+vec[j];
	  TLorentzVector pair2 = TLorentzVector();

	  for(size_t k=0; k<vec.size(); k++)
	    {
	      if(k!=i && k!=j)
		{
		  pair2 += vec[k];
		}
	    }
	  float mass1 = pair1.M();
	  float mass2 = pair2.M();
	  float deltaMass = std::fabs(mass1-mass2);

	  if(deltaMass<minDeltaMass)
	    {
	      minDeltaMass = deltaMass ;
	    }
	}
    }
  return minDeltaMass;
}

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

// Transverse mass
float getMT(TLorentzVector vec_1, TLorentzVector vec_2)
{
  float dphi = vec_1.DeltaPhi(vec_2);

  return std::sqrt(2 * vec_1.Pt() * vec_2.Pt() * (1 - TMath::Cos(dphi)));
}

// Sorting
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec)
{
  std::sort(vec.begin(), vec.end(), [](const TLorentzVector& a, const TLorentzVector& b) {
    return a.Pt() > b.Pt();
  });

  return vec;
}

bool sortPair(std::pair<TLorentzVector, float> vec_1, std::pair<TLorentzVector, float> vec_2)
{
  return vec_1.second > vec_2.second;
}

float getDeltaPhiMin(TLorentzVector vec_1, std::vector<TLorentzVector> vec_2)
{
    float dphiMin1 = fabs(vec_1.DeltaPhi(vec_2[0]));
    float dphiMin2 = fabs(vec_1.DeltaPhi(vec_2[1]));
    float dphiMin3 = fabs(vec_1.DeltaPhi(vec_2[2]));
   
    std::vector<float> dphivalues = {dphiMin1, dphiMin2, dphiMin3};
    
    auto minElem = std::min_element(dphivalues.begin(), dphivalues.end());
    
    return *minElem;
}
