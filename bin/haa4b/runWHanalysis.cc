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
float getMT(TLorentzVector vec_1, TLorentzVector vec_2);
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec);
float HT(std::vector<TLorentzVector> vec);
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

  bool isMC_Wh = isMC && (string(url.Data()).find("Wh_amass")  != string::npos); 
  bool isMC_Zh = isMC && (string(url.Data()).find("Zh_amass")  != string::npos);
  bool isMC_VBFh = isMC && (string(url.Data()).find("VBFh_amass") !=string::npos);
  bool isSignal = (isMC_Wh || isMC_Zh || isMC_VBFh);

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

  // Delta R
  mon.addHistogram( new TH1F( "DR", ";#DeltaR ;N_{events}", 100, 0, 7) );

  // Delta phi
  mon.addHistogram( new TH1F( "Dphi", ";|#Delta#phi| ;N_{Events}", 100, 0, TMath::Pi()) );

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

  //####################################################################################################################//
  //###########################################           TMVAReader         ###########################################//
  //####################################################################################################################//
  std::string chpath = "WhAnalysis/";
  TMVAReader WhAnalysisReader;
  WhAnalysisReader.InitTMVAReader();
  std::string WhAnalysis_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVAClassification_BDT.weights.xml";
  WhAnalysisReader.SetupMVAReader("WhAnalysisClass", WhAnalysis_xml_path);

 
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

      mon.fillHisto("eventflow","histo", 0, weight);

      if (isSignal)
	{
	  if (isMC_Wh)
	    {
	       mon.fillHisto("eventflow","Wh", 0, weight);
	    }
	  if (isMC_Zh)
	    {
	      mon.fillHisto("eventflow","Zh", 0, weight);
	    }
	  if (isMC_VBFh)
	    {
	      mon.fillHisto("eventflow", "VBFh", 0, weight);
	    }
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
	      // mon.fillHisto("pt", "Hgen", pH.Pt(), weight); etc 
	      vec_H.push_back(pH);        
	    }

	  if(abs(ev.mc_id[imc])==24)
	    { // Found the W
	      TLorentzVector pW;
	      pW.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Make histograms of the W boson kinematics
	      // mon.fillHisto("pt", "Wgen", pW.Pt(), weight); etc
	      vec_W.push_back(pW);
	    }

	  if(ev.mc_id[imc]==36)
            { // Found the A                                                              
              TLorentzVector pA;
              pA.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

              // Make histograms of the A boson kinematics
	      // mon.fillHisto("pt", "A", pA.Pt(), weight); etc
	      vec_A.push_back(pA);
            }

	  if(abs(ev.mc_id[imc])==5 && ev.mc_mom[imc]==36)
	    { // Found the b quarks
	      TLorentzVector pbb;
	      pbb.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Apply acceptance cuts
	      if(fabs(pbb.Eta()) > 2.4 || pbb.Pt() < 20) continue;

	      // b-quarks from mother1
	      if(ev.mc_momidx[imc]==4) vec_bb1.push_back(pbb);

	      // b-quarks from mother2
	      else if(ev.mc_momidx[imc]==5) vec_bb2.push_back(pbb);
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
		{
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
      mon.fillHisto("multi", "genlepton", vec_genlep.size(), weight);
      mon.fillHisto("qgmulti", "gen", vec_genqg.size(), weight);

      // Sorting the quarks/gluons based on pT in descending order
      vec_genqg = sort_vec_pt(vec_genqg);

      // Fill histograms with bb-pair
      if(vec_bb1.size()>1)
	{
	  TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	  float ptbb1 = pbb1.Pt();

	  mon.fillHisto("pt", "bb1", ptbb1, weight);
	}
      

      if(vec_bb2.size()>1)
	{
	  TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	  float ptbb2 = pbb2.Pt();

	  mon.fillHisto("pt", "bb2", ptbb2, weight);
	}

      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is at least one lepton       
      if(vec_genlep.size()>0)
	{
	  // Fill the multiplicity histograms
	  mon.fillHisto("multi", "genlepton_step1", vec_genlep.size(), weight);
	  mon.fillHisto("qgmulti", "gen_step1", vec_genqg.size(), weight);
			
	  // Fill histograms of angular distance DR
	  if(vec_bb1.size()>1)
	    {
	      float DRbb1 = getDeltaR(vec_bb1[0], vec_bb1[1]);
	      mon.fillHisto("DR", "bb1", DRbb1, weight);
	      
	      TLorentzVector pbb1 = vec_bb1[0] + vec_bb1[1];
	      mon.fillProfile("DR_vs_pt", "bb1", pbb1.Pt(), DRbb1, weight);
	    }

	  if(vec_bb2.size()>1)
	    {
	      float DRbb2 = getDeltaR(vec_bb2[0], vec_bb2[1]);
	      mon.fillHisto("DR", "bb2", DRbb2, weight);

	      TLorentzVector pbb2 = vec_bb2[0] + vec_bb2[1];
	      mon.fillProfile("DR_vs_pt", "bb2", pbb2.Pt(), DRbb2, weight);
	    }

	  if(vec_bb1.size()>1 && vec_bb2.size()>1)
	    {
	      TLorentzVector vec_A1 = vec_bb1[0] + vec_bb1[1];
	      TLorentzVector vec_A2 = vec_bb2[0] + vec_bb2[1];
       
	      float DRAA = getDeltaR(vec_A1, vec_A2);
	      mon.fillHisto("DR", "AA", DRAA, weight);
	    }
	}

	  
      //### STEP 2 ###//
      // If there is at least one lepton and at least three quarks/gluons                                                                                     
      if(vec_genlep.size()>0 && vec_genqg.size()>2)
	{                                             
	  // Fill the histogram of HT q/g
	  float HTqg = HT(vec_genqg);
	  mon.fillHisto("HT", "qg", HTqg, weight);
	  
	  // Fill the histogram of MET pT
	  mon.fillHisto("pt", "genmet", met_gen.Pt(), weight);

	  // Fill the histogram of W transverse mass
	  float MTW_gen = getMT(vec_genlep[0], met_gen);
	  mon.fillHisto("MT", "genW", MTW_gen, weight);

	  // Define the hadronic system
	  TLorentzVector phad_gen = vec_genqg[0] + vec_genqg[1] + vec_genqg[2];
	  if(vec_genqg.size()>3) phad_gen += vec_genqg[3];

	  // Fill the kinematics histograms of hadronic system
	  mon.fillHisto("pt", "hadgen", phad_gen.Pt(), weight);
	  mon.fillHisto("eta", "hadgen", phad_gen.Eta(), weight);
	  mon.fillHisto("phi", "hadgen", phad_gen.Phi(), weight);
	  mon.fillHisto("M", "hadgen", phad_gen.M(), weight);

	  // Define the leptonic system
	  TLorentzVector plep_gen = vec_genlep[0] + met_gen;
	  mon.fillHisto("pt", "genlep", plep_gen.Pt(), weight);

	  // Fill Delta phi WH histogram
	  float dphiWH_gen = fabs(plep_gen.DeltaPhi(phad_gen));
	  mon.fillHisto("Dphi", "gen_WH", dphiWH_gen, weight);

	  // Fill Delta phi of MET and lepton
	  float dphiMETlep_gen = fabs(met_gen.DeltaPhi(vec_genlep[0]));
	  mon.fillHisto("Dphi", "gen_met-lepton", dphiMETlep_gen, weight);	  
	
	} // END OF GENERATOR LEVEL ANALYSIS 
	  
	    
      //######################################################################//
      //##########################  DETECTOR LEVEL ##########################//
      //####################################################################//

      // Create new object vectors
      std::vector<TLorentzVector> vec_en;                           // store the electron
      std::vector<TLorentzVector> vec_mn;                           // store the muon
      std::vector<TLorentzVector> vec_lepton;                       // store the lepton
      std::vector<TLorentzVector> vec_jets;                         // store the jets
      std::vector<TLorentzVector> vec_bjets;                        // store the b-tags
      std::vector<TLorentzVector> vec_fjets;                        // store the fat jets
      std::vector<std::pair<TLorentzVector, float>> vec_jetbtag;    // store the jet with its b-tag index


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
      mon.fillHisto("multi", "lepton_step0", vec_lepton.size(), weight);
      mon.fillHisto("multi", "en_step0", vec_en.size(), weight);
      mon.fillHisto("multi", "mn_step0", vec_mn.size(), weight);


      // Jets and cross-cleaning
      for(int i=0; i<ev.jet; i++)
	{
	  bool overlap = false;

	  TLorentzVector pjet;
	  pjet.SetPxPyPzE(ev.jet_px[i], ev.jet_py[i], ev.jet_pz[i], ev.jet_en[i]);

	  // Apply acceptance cuts
	  if(pjet.Pt()<20 || fabs(pjet.Eta())>2.4) continue;

	  for(int mn_count=0; mn_count<vec_mn.size(); mn_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_mn[mn_count]);

	      // Fill the histogram of DR(jet-muon) before cross cleaning
	      mon.fillHisto("DR", "jet-mn_step0", DRjet, weight);

	      if(vec_lepton.size()>1)
		{
		  mon.fillHisto("DR", "jet-mn_step1", DRjet, weight);
		}
	      
	      if(DRjet<0.4)
		{
		  overlap = true;
		}
	      else
		{
		  mon.fillHisto("DR", "jet-mn_cc_step0", DRjet, weight);

		  if(vec_lepton.size()>1)
		    {
		      mon.fillHisto("DR", "jet-mn_cc_step1", DRjet, weight);
		    }
		}
	    }
      
	  
	  for(int en_count=0; en_count<vec_en.size(); en_count++)
	    {
	      float DRjet = getDeltaR(pjet, vec_en[en_count]);

	      // Fill the histogram of DR(jet-electron) before cross cleaning
	      mon.fillHisto("DR", "jet-en_step0", DRjet, weight);

	      if(vec_lepton.size()>1)
		{
		  mon.fillHisto("DR", "jet-en_step1", DRjet, weight);
		}

	      if(DRjet<0.4)
		{
		  overlap = true;
		}
	      else
		{
		  mon.fillHisto("DR", "jet-en_cc_step0", DRjet, weight);

		  if(vec_lepton.size()>1)
		    {
		      mon.fillHisto("DR", "jet-en_cc_step1", DRjet, weight);
		    }
		}
	    }
	  
	  if(!overlap)
	    {
	      vec_jets.push_back(pjet);
	      vec_jetbtag.push_back(std::make_pair(pjet,ev.jet_btag1[i]));

	      if(ev.jet_btag1[i]>0.4941)
	      	{ // Found the b-jets
	      	  vec_bjets.push_back(pjet);
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

      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "untag_step0", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "b_step0", vec_bjets.size(), weight);

       
      //##################################################//
      //########    EVENT SELECTION CRITERIA     ########//
      //#################################################//

      //### STEP 1 ###//
      
      // If there is at least one lepton
      if(vec_lepton.size()<1) continue;

      n_event_lepton++;

      mon.fillHisto("eventflow","histo", 1, weight);

      if (isSignal)
	{
	  if (isMC_Wh)
	    {
	      mon.fillHisto("eventflow","Wh", 1, weight);
	    }
	  if (isMC_Zh)
	    {
	      mon.fillHisto("eventflow","Zh", 1, weight);
	    }
	  if (isMC_VBFh)
	    {
	      mon.fillHisto("eventflow", "VBFh", 1, weight);
	    }
	}

      // Fill the multiplicity histograms
      mon.fillHisto("multi", "lepton_step1", vec_lepton.size(), weight);
      mon.fillHisto("multi", "en_step1", vec_en.size(), weight);
      mon.fillHisto("multi", "mn_step1", vec_mn.size(), weight);

      // Fill lepton kinematics
      mon.fillHisto("pt", "lepton", vec_lepton[0].Pt(), weight);
      mon.fillHisto("eta", "lepton", vec_lepton[0].Eta(), weight);
      mon.fillHisto("phi", "lepton", vec_lepton[0].Phi(), weight);
      
      // Fill the multiplicity histograms of jets
      mon.fillHisto("jetmulti", "untag_step1", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "b_step1", vec_bjets.size(), weight);
     
      // Sort the jets based on pT in descending order
      vec_bjets = sort_vec_pt(vec_bjets);
     
      // Sort the jets based on their b-tag
      std::sort(vec_jetbtag.begin(), vec_jetbtag.end(), sortPair);

           
      //### STEP 2 ###//
      
      // If there is at least three jets
      if(vec_jets.size()<3) continue;

      n_event_jet++;

      mon.fillHisto("eventflow","histo", 2, weight);

       if (isSignal)
	{
	  if (isMC_Wh)
	    {
	       mon.fillHisto("eventflow","Wh", 2, weight);
	    }
	  if (isMC_Zh)
	    {
	      mon.fillHisto("eventflow","Zh", 2, weight);
	    }
	  if (isMC_VBFh)
	    {
	      mon.fillHisto("eventflow", "VBFh", 2, weight);
	    }
	}

         
      //### STEP 3 ###//
      
      // If there is at least three b-jets
      if(vec_bjets.size()<3) continue;

      n_event_bjet++;

      mon.fillHisto("eventflow","histo", 3, weight);

      if (isSignal)
	{
	  if (isMC_Wh)
	    {
	      mon.fillHisto("eventflow","Wh", 3, weight);
	    }
	  if (isMC_Zh)
	    {
	      mon.fillHisto("eventflow","Zh", 3, weight);
	    }
	  if (isMC_VBFh)
	    {
	      mon.fillHisto("eventflow", "VBFh", 3, weight);
	    }
	}

      // Define MET
      TLorentzVector met;
      met.SetPtEtaPhiM(ev.met_pt, 0, ev.met_phi, 0);
      float MET = met.Pt();

      // Define W transverse mass
      float MTW = getMT(vec_lepton[0], met);
    
      //### STEP 4 ###//
      if(MET<25 || MTW<50) continue;

      n_event_metcut++;
      
      mon.fillHisto("eventflow","histo", 4, weight);

      if (isSignal)
	{
	  if (isMC_Wh)
	    {
	      mon.fillHisto("eventflow","Wh", 4, weight);
	    }
	  if (isMC_Zh)
	    {
	      mon.fillHisto("eventflow","Zh", 4, weight);
	    }
	  if (isMC_VBFh)
	    {
	      mon.fillHisto("eventflow", "VBFh", 4, weight);
	    }
	}

      // Define the hadronic system
      TLorentzVector phadb = vec_bjets[0] + vec_bjets[1] + vec_bjets[2];
      if(vec_bjets.size()>3) phadb += vec_bjets[3];

      float H_pt  = phadb.Pt();
      float H_eta = phadb.Eta();
      float H_phi = phadb.Phi();
      float H_M   = phadb.M();

      // Fill the hadronic kinematics histograms
      mon.fillHisto("pt", "H", H_pt, weight);
      mon.fillHisto("eta", "H", H_eta, weight);
      mon.fillHisto("phi", "H", H_phi, weight);
      mon.fillHisto("M", "H", H_M, weight);

      // Fill the lepton pt histogram
      float lepton_pt = vec_lepton[0].Pt();
      mon.fillHisto("pt", "lepton", lepton_pt, weight);

      // Fill the most energetic b-jet pt
      float bjet1_pt = vec_bjets[0].Pt();
      mon.fillHisto("pt", "bjet1", bjet1_pt, weight);

      // Fill the MET histogram
      mon.fillHisto("pt", "met", MET, weight);

      // Fill the histogram of W transverse mass
      mon.fillHisto("MT", "W", MTW, weight);
           
      // Fill the HT histogram of jets
      float HTjets = HT(vec_jets);
      mon.fillHisto("HT", "jets", HTjets, weight);

      // Fill the Delta phi WH histogram
      float dphiWH = fabs(vec_lepton[0].DeltaPhi(phadb));
      mon.fillHisto("Dphi", "WH", dphiWH, weight);

      // Define the leptonic system
      TLorentzVector plepsys = vec_lepton[0] + met;
      float W_pt = plepsys.Pt();
      
      // Fill the reconstructed W pt
      mon.fillHisto("pt", "W", W_pt, weight);

      // Fill the Delta phi of MET and lepton
      float dphiMetLepton = fabs(met.DeltaPhi(vec_lepton[0]));
      mon.fillHisto("Dphi", "met-lepton", dphiMetLepton, weight);

      // Find the minimum Delta phi of MET and jet
      float dphiMetJetMin1 = fabs(met.DeltaPhi(vec_jets[0]));
      float dphiMetJetMin2 = fabs(met.DeltaPhi(vec_jets[1]));
      float dphiMetJetMin3 = fabs(met.DeltaPhi(vec_jets[2]));
      // Store the values in a vector
      std::vector<float> dphivalues = {dphiMetJetMin1, dphiMetJetMin2, dphiMetJetMin3};
      // Find the minimum delta phi
      auto minElem = std::min_element(dphivalues.begin(), dphivalues.end());
      int minIndex = std::distance(dphivalues.begin(), minElem);
      float dphiMetJetMin = dphivalues[minIndex];
      // Fill the histogram
      mon.fillHisto("Dphi", "metjet-min", dphiMetJetMin, weight);
 
      // Fill the multiplicity histograms of jets
      float Njets = vec_jets.size();
      mon.fillHisto("jetmulti", "untag_step4", vec_jets.size(), weight);
      mon.fillHisto("jetmulti", "b_step4", vec_bjets.size(), weight);

      // b-tag discriminator for 3 most energetic b-jets
      float btag1 = vec_jetbtag[0].second;
      float btag2 = vec_jetbtag[1].second;
      float btag3 = vec_jetbtag[2].second;
      mon.fillHisto("btag", "jet1", btag1, weight);
      mon.fillHisto("btag", "jet2", btag2, weight);
      mon.fillHisto("btag", "jet3", btag3, weight);

     
      //#####################################//
      //############ TMVA Reader ############//
      //#####################################//
      float mvaBDT = -10.0;

      mvaBDT = WhAnalysisReader.GenReMVAReader(Njets, HTjets, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, dphiMetJetMin, btag1, btag2, btag3, "WhAnalysisClass");

      mon.fillHisto("BDT", "WhAnalysis", mvaBDT, weight);

     

      //#####################################//
      //############ MVA Handler ############//
      //#####################################//
      
      if(runMVA)
	{
	  myMVAHandler_.getEntry(Njets, HTjets, H_pt, W_pt, lepton_pt, bjet1_pt, MET, H_M, MTW, dphiWH, dphiMetLepton, dphiMetJetMin, btag1, btag2, btag3, weight);
	  myMVAHandler_.fillTree();
	} // runMVA      
      
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

// Scalar sum of pTs
float HT(std::vector<TLorentzVector> vec)
{
  float HT = 0.0;

  if(vec.size()>0)
    {
      for(int i=0; i<vec.size(); i++) HT += vec[i].Pt();
    }
  return HT;
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

