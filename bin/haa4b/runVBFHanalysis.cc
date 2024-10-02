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

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
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

// Initialize functions to use
double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2);
double getDeltaRyy(TLorentzVector vec_1, TLorentzVector vec_2);
bool match_jets(TLorentzVector det_vec, std::vector<TLorentzVector> gen_vec, double threshold);
std::pair<double, int> getDeltaRmin(TLorentzVector det_vec, std::vector<TLorentzVector> gen_vec);
triplet findBoostandRotation(TLorentzVector p_mom);
int rand_int(int min, int max);
double getDeltaMmin(const std::vector<TLorentzVector>& vec_bjets, double& unpairedMass);
double getMassOfbbj(const std::vector<jets_struct>& vec_bjets, const std::vector<jets_struct>& vec_untaggedjets);
double getDeltaMminNew(const std::vector<TLorentzVector>& vec_bjets);

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


  //#######################################################################################################################################//
  //--------------------------------------------------------INITIALIZE HISTOGRAMS----------------------------------------------------------//
  //#######################################################################################################################################//

  SmartSelectionMonitor mon;

  // Event flow 
  TH1F *h_event_flow = (TH1F*) mon.addHistogram ( new TH1F ("event_flow", "Cut Flow;Step;N_{Events}", 6,0,6) );
  h_event_flow->GetXaxis()->SetBinLabel(1,"0: Raw Events");
  h_event_flow->GetXaxis()->SetBinLabel(2,"1: N_leptons=0");
  h_event_flow->GetXaxis()->SetBinLabel(3,"2: N_jets>=5");
  h_event_flow->GetXaxis()->SetBinLabel(4,"3: N_bjets>=3 and N_qjets>=2");
  h_event_flow->GetXaxis()->SetBinLabel(5,"4: Tight btag for at least two jets");
  h_event_flow->GetXaxis()->SetBinLabel(6,"5: Dh>2.5 and Mqq>250");
  
  // Particles/Objects kinematics (single)
  mon.addHistogram ( new TH1F ("pt", "Transverse Momentum;p_{T} [GeV];N_{Events}", 150, 0, 500) );
  mon.addHistogram ( new TH1F ("eta", "Pseudorapidity;#eta [-];N_{Events}", 150, -6, 6) );
  mon.addHistogram ( new TH1F ("psi", "Rapidity;y [-];N_{Events}", 150, -6, 6) );
  mon.addHistogram ( new TH1F ("phi", "Azimuthal Angle;#phi [-];N_{Events}", 150, -TMath::Pi(), TMath::Pi()) );
  mon.addHistogram ( new TH1F ("M", "Mass;M [GeV];N_{Events}", 150, 0, 800) );
  mon.addHistogram ( new TH1F ("E", "Energy;E [GeV];N_{Events}", 150, 0, 500) );

  mon.addHistogram ( new TH1F ("theta", "Polar Angle;#theta [-];N_{Events}", 150, 0, TMath::Pi()) );
  mon.addHistogram ( new TH1F ("costheta", "Polar Angle;cos #theta [-];N_{Events}", 150, -1, 1) );

  // -//- (multi)
  mon.addHistogram ( new TH1F ("deltaR", "Angular Distance;#Delta R [-];N_{Events}", 150, 0, 16) );
  mon.addHistogram ( new TH1F ("deltaRyy", "Angular Distance;#Delta R [-];N_{Events}", 150, 0, 7) );
  mon.addHistogram ( new TH1F ("deltaEta", "Pseudorapidity Difference;#left|#Delta#eta#right| [-];N_{Events}", 150, 0, 10) );
  mon.addHistogram ( new TH1F ("deltaPhi", "Azimuthal Angle Difference;#left|#Delta#phi#right| [-];N_{Events}", 150, 0, 10) );
  mon.addHistogram ( new TH1F ("prodEta", "Pseudorapidity Product;#eta_{1} #times #eta_{2} [-];N_{Events}", 150, -25, 10) );
  mon.addHistogram ( new TH1F ("Mqq", "Mass;M [GeV];N_{Events}", 150, 0, 1400) );
  mon.addHistogram ( new TH1F ("avdeltaR", "Angular Distance;#bar{#Delta R}_{bb} [-];N_{Events}", 150, 0, 7) );

  // Multiplicity
  mon.addHistogram ( new TH1F ("multi", ";N [-] ;N_{Events}", 18, 0, 18) );

  // B-tag discriminator
  mon.addHistogram( new TH1F( "btag", "b-tag discriminator;score [-];N_{Events}", 150, 0, 1) );

  // HT, Hz
  mon.addHistogram( new TH1F( "HT", ";H_{T} [GeV];N_{Events}", 150, 0, 1400) );
  mon.addHistogram ( new TH1F ("Hz", "Longitudinal Momentum;#sum_{i}{p_{z,{b_i}}} + #sum_{j}{p_{z,{q_j}}} [GeV];N_{Events}", 300, -6000, 6000) );
  mon.addHistogram( new TH1F( "HT_vec", ";#frac{#sum_{i}{\vec{p_i}}_T}{#sum_{i}{p_{i,T}}} [-];N_{Events}", 150, 0, 2) );

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

  //#######################################################################################################################################//
  //----------------------------------------------------INITIALIZE KINEMATIC VARIABLES-----------------------------------------------------//
  //#######################################################################################################################################//

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

  // outgoing quark pair kinematics
  float DeltaEta;
  float qq_m;
  float ProductEta;
  float qq_dR;
  float qq_dphi;
  float alpha_qq;

  // jet kinematics
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
  float pt_rest;
  float E_rest;

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

  // untagged jets
  float N_untaggedjet;

  // met
  float MET_pt;
  float MET_phi;
  
  //#######################################################################################################################################//

  //----------------------------------------------------------Event Counters-----------------------------------------------------------------//
  Int_t count_leptons = 0;
  Int_t count_jets = 0;
  Int_t count_bjets = 0;
  Int_t count_btag = 0;
  Int_t count_jets_max = 0;
  Int_t count_final_events = 0;

  //-------------------------------------------------------Signal Region Counters------------------------------------------------------------//
  Int_t count_bjets_SR1 = 0;
  Int_t count_bjets_SR2 = 0;
  Int_t count_bjets_SR3 = 0;

  //#######################################################################################################################################//

  //---------------------------------------------------------------Get Tree info-------------------------------------------------------------f-//
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
            
      if(!isMC && duplicatesChecker.isDuplicate(ev.run, ev.lumi, ev.event))
	{
	  nDuplicates++;
	  cout << "nDuplicates: " << nDuplicates << endl;
	  continue;
	}

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

	  for (int imc=0; imc<ev.nmcparticles;imc++)
	    {
	      TLorentzVector p_particle;
	      p_particle.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

	      // Find how many categories exist in status and perform printouts

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

	  // Boost to A boson f.o.m

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

      // Leptons and jets vectors
      std::vector<TLorentzVector> leptons;
      std::vector<TLorentzVector> muons;
      std::vector<TLorentzVector> electrons;
      
      std::vector<TLorentzVector> jets;
      std::vector<jets_struct> b_jets;
      std::vector<jets_struct> non_b_jets;
      std::vector<jets_struct> q_jets;
      std::vector<jets_struct> untagged_jets;

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

      for (int ijet = 0; ijet < ev.jet; ijet++)
	{
	  TLorentzVector p_jet;
	  p_jet.SetPxPyPzE(ev.jet_px[ijet], ev.jet_py[ijet], ev.jet_pz[ijet], ev.jet_en[ijet]);
	  
	  // Acceptance cuts
	  if (p_jet.Pt()<20. || abs(p_jet.Eta())>4.7) continue;
	  
	  // Id
	  if (!ev.jet_PFTight[ijet]) continue;

	  // Cross cleaning with respect to electrons and muons
	  bool overlap = false;

	  for (int imn = 0; imn < muons.size(); imn++)
	    {
	      float dR_jet_mn = getDeltaR(p_jet, muons[imn]);
	      mon.fillHisto("deltaR", "jet_mn_before", dR_jet_mn, weight);
	      if (dR_jet_mn < 0.4) overlap = true;	     
	      else mon.fillHisto("deltaR", "jet_mn_after", dR_jet_mn, weight);	   
	    }
	  
	  for (int ien = 0; ien < electrons.size(); ien++)
	    {
	      float dR_jet_en = getDeltaR(p_jet, electrons[ien]);
	      mon.fillHisto("deltaR", "jet_en_before", dR_jet_en, weight);
	      if (dR_jet_en < 0.4) overlap = true;
	      else mon.fillHisto("deltaR", "jet_en_after", dR_jet_en, weight);		   
	    }
	  
	  if (!overlap)
	    {
	      jets.push_back(p_jet);
	      
	      // b tagging algorithm
	      double btag_value = ev.jet_btag1[ijet];

	      bool match_bjets = match_jets(p_jet, b_quarks, 0.3);
	      bool match_qjets = match_jets(p_jet, out_quarks, 0.5);
	      
	      // Configure b jets and non b jets
	      if (ev.jet_btag1[ijet] > 0.3040) //medium working point
		//if (jet_btag1[ijet] > 0.7476) //tight working point
		{
		  // Demand b jets to have eta less than 2.4
		  if (abs(p_jet.Eta())>2.4) continue;
		  
		  if (isSignal) // This analysis corresponds only to signal (it can be removed when the fine tuning is done) and it matches the tagged jets with the "true" ones
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
		  if (isSignal)
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
      
      //printout
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

      //sort jets on pt
      std::sort(jets.begin(), jets.end(), [](const TLorentzVector &a, const TLorentzVector &b){
	  return a.Pt() > b.Pt();
	});

      //sort non b jets on pt
      std::sort(non_b_jets.begin(), non_b_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.momentum.Pt() > b.momentum.Pt();
	});

      //----------configure qjets and untagged jets from non b jets----------//

      if (non_b_jets.size()>=2)
	{
	  q_jets.push_back(non_b_jets[0]);
	  q_jets.push_back(non_b_jets[1]);
	}
      
      if (non_b_jets.size()>=3)
	{
	  for (size_t i = 2; i < non_b_jets.size(); i++)
	    {
	      untagged_jets.push_back(non_b_jets[i]);
	    }
	}

      //sort q jets (unnecessary)
      std::sort(q_jets.begin(), q_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.momentum.Pt() > b.momentum.Pt();
	});

      //sort untagged jets (based on btag)
      std::sort(untagged_jets.begin(), untagged_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.btag > b.btag;
	});

      //----------------------------------------------------------------------------//
      //-------------------------------Fill histogramms-----------------------------//
      //----------------------------------------------------------------------------//

      // Lepton multiplicity
      
      mon.fillHisto("multi", "el", electrons.size(), weight);
      mon.fillHisto("multi", "mn", muons.size(), weight);
      mon.fillHisto("multi", "lep", leptons.size(), weight);	   

      // Compare with generator

      if (isSignal)
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
	} // END if isSignal

      //----------------------------------------------------------------------------//
      //---------------------------EVENT SELECTION CRITERIA-------------------------// 
      //----------------------------------------------------------------------------//

      mon.fillHisto("event_flow", "hist", 0.5, 1);

      // STEP1: Recquire zero leptons
      if (leptons.size()!=0) continue;
      
      count_leptons++;

      mon.fillHisto("event_flow", "hist", 1.5, 1);

      // Control plot: jet multiplicity
      mon.fillHisto("multi", "jet", jets.size(), weight);

      // STEP2: Recquire at least 5 jets
      if(jets.size()<5) continue;
      
      count_jets++;

      mon.fillHisto("event_flow", "hist", 2.5, 1);
 
      // Control plot: b and q jet multiplicity
      mon.fillHisto("multi", "bjet", b_jets.size(), weight);
      mon.fillHisto("multi", "qjet", q_jets.size(), weight);
      mon.fillHisto("multi", "untaggedjets", untagged_jets.size(), weight);

      // Compare with generator

      if (isSignal)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets","step2", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets","step2", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets","step2", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets","step2", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets","step2", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets","step2", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets","step2", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets","step2", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets","step2", 2.5, 1);
	    }
	} // END if isSignal

      //STEP3: Recquire at least 3 b tagged and 2 q tagged jets
      if(b_jets.size()<3 || q_jets.size()<2) continue;

      count_bjets++;

      mon.fillHisto("event_flow", "hist", 3.5, 1);

      DeltaEta = abs(q_jets[0].momentum.Eta()-q_jets[1].momentum.Eta());
      qq_m = (q_jets[0].momentum+q_jets[1].momentum).M();

      // Control plot: b jet, q jet and untagged jet multiplicity
      mon.fillHisto("multi", "bjet_step3", b_jets.size(), weight);
      mon.fillHisto("multi", "qjet_step3", q_jets.size(), weight);
      mon.fillHisto("multi", "untaggedjet_step3", untagged_jets.size(), weight);         

      // Control plot: qq mass, eta difference and eta product plots
      mon.fillHisto("deltaEta", "qjet1_qjet2_before", DeltaEta, weight);
      mon.fillHisto("prodEta", "qjet1_qjet2_before", q_jets[0].momentum.Eta()*q_jets[1].momentum.Eta(), weight);
      mon.fillHisto("Mqq", "qjet1_qjet2_before", qq_m, weight);

      // Compare with generator

      if (isSignal)
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
	} // END if isSignal

      //sort b jets on btag
      std::sort(b_jets.begin(), b_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.btag > b.btag;
	});

      //STEP4: Recquire tight working point for the btag algorithm for the first two b jets
      if (b_jets[0].btag < 0.7476 || b_jets[1].btag < 0.7476) continue;

      count_btag++;

      mon.fillHisto("event_flow", "hist", 4.5, 1);

      // Compare with generator

      if (isSignal)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets","step4", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets","step4", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets","step4", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets","step4", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets","step4", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets","step4", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets","step4", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets","step4", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets","step4", 2.5, 1);
	    }
	} // END if isSignal


      //sort b jets on pt
      std::sort(b_jets.begin(), b_jets.end(), [](const jets_struct &a, const jets_struct &b){
	  return a.momentum.Pt() > b.momentum.Pt();
	});

      //STEP5: Recquire Mqq more than 250 GeV and Dhqq more than 2.5
      
      if (DeltaEta<2.5 || qq_m<250) continue;

      mon.fillHisto("event_flow", "hist", 5.5, 1);

      count_final_events++;

      // Compare with generator

      if (isSignal)
	{
	  // Matching b tagged, q tagged and untagged jets with either b quarks, q quarks or none (categorical plots)

	  //b tagged
	  for (size_t i = 0; i < b_jets.size(); i++)
	    {
	      int category = b_jets[i].category;
	      if (category==1) mon.fillHisto("match_bjets","final", 0.5, 1);
	      if (category==2) mon.fillHisto("match_bjets","final", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_bjets","final", 2.5, 1);
	    }
	  
	  //q tagged
	  for (size_t i = 0; i < q_jets.size(); i++)
	    {
	      int category = q_jets[i].category;
	      if (category==1) mon.fillHisto("match_qjets","final", 0.5, 1);
	      if (category==2) mon.fillHisto("match_qjets","final", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_qjets","final", 2.5, 1);
	    }
	  
	  //untagged
	  for (size_t i = 0; i < untagged_jets.size(); i++)
	    {
	      int category = untagged_jets[i].category;
	      if (category==1) mon.fillHisto("match_untaggedjets","final", 0.5, 1);
	      if (category==2) mon.fillHisto("match_untaggedjets","final", 1.5, 1);
	      if (category==-1) mon.fillHisto("match_untaggedjets","final", 2.5, 1);
	    }

	  //untagged n.0 matchings 

	  if (untagged_jets.size()>0 && (abs(untagged_jets[0].momentum.Eta()) < 2.4 || abs(untagged_jets[0].momentum.Eta()) > 3.3)) 
	    {

	      if (untagged_jets[0].category == 1)
		{
		  mon.fillHisto(TString("match_eta"), "untagged0", untagged_jets[0].momentum.Eta(), 0.5, 1.);
		  mon.fillHisto(TString("match_pt"), "untagged0", untagged_jets[0].momentum.Pt(), 0.5, 1.);
		  mon.fillHisto(TString("match_untaggedjets"),"0_final", 0.5, 1);
		}
	      if (untagged_jets[0].category == 2)
		{
		  mon.fillHisto(TString("match_eta"), "untagged0", untagged_jets[0].momentum.Eta(), 1.5, 1.);
		  mon.fillHisto(TString("match_pt"), "untagged0", untagged_jets[0].momentum.Pt(), 1.5, 1.);
		  mon.fillHisto(TString("match_untaggedjets"),"0_final", 1.5, 1);
		}
	      if (untagged_jets[0].category == -1)
		{
		  mon.fillHisto(TString("match_eta"), "untagged0", untagged_jets[0].momentum.Eta(), 2.5, 1.);
		  mon.fillHisto(TString("match_pt"), "untagged0", untagged_jets[0].momentum.Pt(), 2.5, 1.);
		  mon.fillHisto(TString("match_untaggedjets"),"0_final", 2.5, 1);
		}		  
	    }

	  // Calculate the minimum angular distance between each jet and the collection of generated b's and q's (only for signal)

	  for (size_t i=0; i<jets.size(); i++)
	    {
	      TLorentzVector p_jet = jets[i];
	      
	      double dRmin_jet_b = getDeltaRmin(p_jet,b_quarks).first;
	      double dRmin_jet_q = getDeltaRmin(p_jet,out_quarks).first;

	      mon.fillHisto("deltaR", "min_reco_b", dRmin_jet_b, 1);
	      mon.fillHisto("deltaR", "min_reco_q", dRmin_jet_q, 1);
	    }

	  for (size_t i=0; i<b_jets.size(); i++)
	    {
	      TLorentzVector p_jet = b_jets[i].momentum;

	      double dRmin_bjet_b = getDeltaRmin(p_jet,b_quarks).first;
	      double dRmin_bjet_q = getDeltaRmin(p_jet,out_quarks).first;

	      int iter_b = getDeltaRmin(p_jet,b_quarks).second;

	      mon.fillHisto("deltaR", "min_bjet_b", dRmin_bjet_b, 1);
	      mon.fillHisto("deltaR", "min_bjet_q", dRmin_bjet_q, 1);
	      mon.fillHisto(TString("mindeltaR_vs_pt"), "bjet_b", b_quarks[iter_b].Pt(), dRmin_bjet_b, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), "bjet_b", b_quarks[iter_b].Eta(), dRmin_bjet_b, 1.);
	    }

	  for (size_t i=0; i<q_jets.size(); i++)
	    {
	      TLorentzVector p_jet = q_jets[i].momentum;
       
	      double dRmin_qjet_b =  getDeltaRmin(p_jet,b_quarks).first;
	      double dRmin_qjet_q = getDeltaRmin(p_jet,out_quarks).first;

	      int iter_q = getDeltaRmin(p_jet,out_quarks).second;

	      mon.fillHisto("deltaR", "min_qjet_b", dRmin_qjet_b, 1);
	      mon.fillHisto("deltaR", "min_qjet_q", dRmin_qjet_q, 1);
	      mon.fillHisto(TString("mindeltaR_vs_pt"), "qjet_q", out_quarks[iter_q].Pt(), dRmin_qjet_q, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), "qjet_q", out_quarks[iter_q].Eta(), dRmin_qjet_q, 1.);

	    }

	  if (untagged_jets.size()>0)
	    {		
	      TLorentzVector p_jet_0 = untagged_jets[0].momentum;
	      
	      double dRmin_0_b = getDeltaRmin(p_jet_0,b_quarks).first;
	      double dRmin_0_q = getDeltaRmin(p_jet_0,out_quarks).first;

	      int iter_q = getDeltaRmin(p_jet_0,out_quarks).second;
	      int iter_b = getDeltaRmin(p_jet_0,b_quarks).second;

	      mon.fillHisto(TString("mindeltaR_vs_pt"), "0_b", b_quarks[iter_b].Pt(), dRmin_0_b, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_pt"), "0_q", out_quarks[iter_q].Pt(), dRmin_0_q, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), "0_b", b_quarks[iter_b].Eta(), dRmin_0_b, 1.);
	      mon.fillHisto(TString("mindeltaR_vs_eta"), "0_q", out_quarks[iter_q].Eta(), dRmin_0_q, 1.);
	    
	      mon.fillHisto("deltaR", "min_0_b", dRmin_0_b, 1);
	      mon.fillHisto("deltaR", "min_0_q", dRmin_0_q, 1);
	    }
	  
	} // END if isSignal

      //
      // qq quantities
      //
      
      ProductEta = q_jets[0].momentum.Eta()*q_jets[1].momentum.Eta();
      qq_dR = getDeltaR(q_jets[0].momentum, q_jets[1].momentum);
      qq_dphi = abs(q_jets[0].momentum.Phi() - q_jets[1].momentum.Phi());

      mon.fillHisto("deltaEta","qjet1_qjet2", DeltaEta, weight);
      mon.fillHisto("prodEta","qjet1_qjet2", ProductEta, weight);
      mon.fillHisto("Mqq","qjet1_qjet2", qq_m, weight);
      mon.fillHisto("deltaR","qjet1_qjet2", qq_dR, weight);
      mon.fillHisto("deltaPhi","qjet1_qjet2", qq_dphi, weight);
	    
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

      mon.fillHisto("pt","jet1", jet1_pt, weight);
      mon.fillHisto("pt","jet2", jet2_pt, weight);
      mon.fillHisto("pt","jet3", jet3_pt, weight);
      mon.fillHisto("pt","jet4", jet4_pt, weight);
      mon.fillHisto("pt","jet5", jet5_pt, weight);

      mon.fillHisto("eta","jet1", jet1_eta, weight);
      mon.fillHisto("eta","jet2", jet2_eta, weight);
      mon.fillHisto("eta","jet3", jet3_eta, weight);
      mon.fillHisto("eta","jet4", jet4_eta, weight);
      mon.fillHisto("eta","jet5", jet5_eta, weight);

      if (jets.size()>=6)
	{
	  jet6_pt = jets[5].Pt();
	  jet6_eta = jets[5].Eta();
	  mon.fillHisto("pt","jet6", jet6_pt, weight);
	  mon.fillHisto("eta","jet6", jet6_eta, weight);
	}
      else
	{
	  jet6_pt = std::nan("");
	  jet6_eta = std::nan("");
	}

      //
      // bjet, qjet kinematics (pt,eta) and btag values
      //

      bjet1_pt = b_jets[0].momentum.Pt();
      bjet2_pt = b_jets[1].momentum.Pt();
      bjet3_pt = b_jets[2].momentum.Pt();
      
      bjet1_eta = b_jets[0].momentum.Eta();
      bjet2_eta = b_jets[1].momentum.Eta();
      bjet3_eta = b_jets[2].momentum.Eta();

      bjet1_btag = b_jets[0].btag;
      bjet2_btag = b_jets[1].btag;
      bjet3_btag = b_jets[2].btag;

      mon.fillHisto("pt","bjet1", bjet1_pt, weight);
      mon.fillHisto("pt","bjet2", bjet2_pt, weight);
      mon.fillHisto("pt","bjet3", bjet3_pt, weight);

      mon.fillHisto("eta","bjet1", bjet1_eta, weight);
      mon.fillHisto("eta","bjet2", bjet2_eta, weight);
      mon.fillHisto("eta","bjet3", bjet3_eta, weight);

      mon.fillHisto("btag","bjet1", bjet1_btag, weight);
      mon.fillHisto("btag","bjet2", bjet2_btag, weight);
      mon.fillHisto("btag","bjet3", bjet3_btag, weight);

      if (b_jets.size()>=4)
	{
	  bjet4_pt = b_jets[3].momentum.Pt();
	  bjet4_eta = b_jets[3].momentum.Eta();
	  bjet4_btag = b_jets[3].btag;

	  mon.fillHisto("pt","bjet4", bjet4_pt, weight);
	  mon.fillHisto("eta","bjet4", bjet4_eta, weight);
	  mon.fillHisto("btag","bjet4", bjet4_btag, weight);
	}
      else
	{
	  bjet4_pt = std::nan("");
	  bjet4_eta =  std::nan("");
	  bjet4_btag =  std::nan("");
	}

      qjet1_pt = q_jets[0].momentum.Pt();
      qjet1_eta = q_jets[0].momentum.Eta();
      qjet2_pt = q_jets[1].momentum.Pt();
      qjet2_eta = q_jets[1].momentum.Eta();

      mon.fillHisto("pt","qjet1", qjet1_pt, weight);
      mon.fillHisto("pt","qjet2", qjet2_pt, weight);

      mon.fillHisto("eta","qjet1", qjet1_eta, weight);
      mon.fillHisto("eta","qjet2", qjet2_eta, weight);

      // calculate qq pair vector
      TLorentzVector qjet_total = q_jets[0].momentum + q_jets[1].momentum;

      //
      // bbb(b) ("higgs") kinematic variables
      //

      // "Higgs" four vector and vector of four b's
      TLorentzVector bjet_total;
      std::vector<TLorentzVector> bjet_total_vec;
      
      for (int i = 0; i < 3; ++i)
	{
	  bjet_total += b_jets[i].momentum;
	  bjet_total_vec.push_back(b_jets[i].momentum);
	}
      if (b_jets.size()>=4)
	{
	  bjet_total += b_jets[3].momentum;
	  mon.fillHisto("M","bjet_total_4", bjet_total.M(), weight);
	  bjet_total_vec.push_back(b_jets[3].momentum);
	}
      /*
      //else if (b_jets.size() < 4 && untagged_jets.size() > 0 && abs(untagged_jets[0].momentum.Eta()) < 2.4)
      else if (b_jets.size() < 4 && untagged_jets.size() > 0 && (abs(untagged_jets[0].momentum.Eta()) < 2.4 || abs(untagged_jets[0].momentum.Eta()) > 3.3))
	//else if (b_jets.size() < 4 && untagged_jets.size() > 0 && abs(untagged_jets[0].momentum.Eta()) < 1)
	{
	  bjet_total += untagged_jets[0].momentum;
	  bjet_total_vec.push_back(untagged_jets[0].momentum);
	}
      */
      
      higgs_m = bjet_total.M();
      higgs_pt = bjet_total.Pt();
      higgs_eta = bjet_total.Eta();

      mon.fillHisto("M","bjet_total", higgs_m, weight);
      mon.fillHisto("eta","bjet_total", higgs_eta, weight);
      mon.fillHisto("pt","bjet_total", higgs_pt, weight);

      //
      // Compute azimuthal angle between qq pair and reco Higgs
      //

      phi_qq_higgs = abs(bjet_total.Phi()-qjet_total.Phi());

      mon.fillHisto("phi", "qq_higgs", phi_qq_higgs, weight);

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

      mon.fillHisto("theta", "alpha_qq", alpha_qq, weight);

      //
      // Compute the deltaR between higgs and the two vbf jets
      //

      dR_higgs_q1 = getDeltaR(bjet_total, q_jets[0].momentum);
      dR_higgs_q2 = getDeltaR(bjet_total, q_jets[1].momentum);

      mon.fillHisto("deltaR", "H_qjet1", dR_higgs_q1, weight);
      mon.fillHisto("deltaR", "H_qjet2", dR_higgs_q2, weight);

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
      mon.fillHisto("deltaR","av_bb", avDeltaR, weight);

      //
      // Calculate DeltaMmin
      //
      double unpairedMass;
      double minDeltaMold = getDeltaMmin(bjet_total_vec, unpairedMass);
      minDeltaM = getDeltaMminNew(bjet_total_vec);

      mon.fillHisto("M", "Delta_pair_min", minDeltaMold, weight);
      mon.fillHisto("M", "Delta_pair_min_new", minDeltaM, weight);

      if (bjet_total_vec.size()==3) mon.fillHisto("M", "unpaired_b", unpairedMass, weight);	  

      //
      // Calculate Mbbj
      //

      MassOfbbj = getMassOfbbj(b_jets, untagged_jets);
      mon.fillHisto("M", "bbj", MassOfbbj, weight);	  

      // Find angle of a random b in bbb(b) ("Higgs") rest frame
      
      std::vector<double> b_bjet = findBoostandRotation(bjet_total).b;
      TVector3 k_bjet= findBoostandRotation(bjet_total).k;
      double theta_bjet = findBoostandRotation(bjet_total).theta;     

      if (b_jets.size()>=4)
	{
	  // Generate a random number
	  int random_number_jet = rand_int(0, 3);
	  TLorentzVector bjet_random = b_jets[random_number_jet].momentum;

	  bjet_random.Boost(-b_bjet[0], -b_bjet[1], -b_bjet[2]);
	  bjet_random.Rotate(theta_bjet, k_bjet);

	  costheta0 = bjet_random.CosTheta();
	  mon.fillHisto("costheta","0_bjet", costheta0, weight);
	}
      /*
      //else if (b_jets.size() < 4 && untagged_jets.size() > 0 && abs(untagged_jets[0].momentum.Eta()) < 2.4)
      else if (b_jets.size() < 4 && untagged_jets.size() > 0 && (abs(untagged_jets[0].momentum.Eta()) < 2.4 || abs(untagged_jets[0].momentum.Eta()) > 3.3))
	//else if (b_jets.size() < 4 && untagged_jets.size() > 0 && abs(untagged_jets[0].momentum.Eta()) < 1)
	{
	  //create new "higgs" vector with 1st untagged as 4th b
	  std::vector<jets_struct> vec_bj;
	  for (size_t i = 0; i < b_jets.size(); i++) vec_bj.push_back(b_jets[i]);
	  vec_bj.push_back(untagged_jets[0]);

	  // Generate a random number
	  int random_number_jet = rand_int(0, 3);
	  TLorentzVector bjet_random = vec_bj[random_number_jet].momentum;

	  bjet_random.Boost(-b_bjet[0], -b_bjet[1], -b_bjet[2]);
	  bjet_random.Rotate(theta_bjet, k_bjet);
	  
	  costheta0 = bjet_random.CosTheta();
	  mon.fillHisto("costheta","0_bjet", costheta0, weight);
	}
      */
      else
	{
	  // Generate a random number
	  int random_number_jet = rand_int(0, 2);
	  TLorentzVector bjet_random = b_jets[random_number_jet].momentum;

	  bjet_random.Boost(-b_bjet[0], -b_bjet[1], -b_bjet[2]);
	  bjet_random.Rotate(theta_bjet, k_bjet);

	  costheta0 = bjet_random.CosTheta();
	  mon.fillHisto("costheta","0_bjet", costheta0, weight);
	}
      
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

      for (size_t i = 0; i<3; i++) H_z_draft+=b_jets[i].momentum.Pz();
      if (b_jets.size()>3) H_z_draft+=b_jets[3].momentum.Pz();
      for (size_t i = 0; i<2; i++) H_z_draft+=q_jets[i].momentum.Pz();

      H_z = H_z_draft;
	
      mon.fillHisto("HT","", H_T, weight);
      mon.fillHisto("Hz","", H_z, weight);
      mon.fillHisto("HT_vec","", H_Tvec, weight);	    

      //find final jet multiplicity
      N_jet = jets.size();
      N_bjet = b_jets.size();
      N_qjet = q_jets.size();
      N_untaggedjet = untagged_jets.size();

      int jet_multi_eta_cut(0);
      
      for (size_t i = 0; i<jets.size(); i++) {
	if (abs(jets[i].Eta()<2.4)) jet_multi_eta_cut++;
      }
      N_jet_eta_cut = jet_multi_eta_cut;

      mon.fillHisto("multi", "jet_final", N_jet, weight);
      mon.fillHisto("multi", "jet_final_eta_cuts", jet_multi_eta_cut, weight);
      mon.fillHisto("multi", "bjet_final", N_bjet, weight);
      mon.fillHisto("multi", "qjet_final", N_qjet, weight);
      mon.fillHisto("multi", "untaggedjet_final", N_untaggedjet, weight);

      //
      //find met pt and phi
      //
      
      MET_pt = ev.met_pt;
      MET_phi = ev.met_phi;

      mon.fillHisto("pt", "met", MET_pt, weight);
      mon.fillHisto("phi", "met", MET_phi, weight);

      //
      // Find p_T sum and energy sum of the untagged jets
      //

      double pt_rest_draft = 0;
      double E_rest_draft = 0;
      for (size_t i = 0; i<untagged_jets.size(); i++)
	{
	  pt_rest_draft += untagged_jets[i].momentum.Pt();
	  E_rest_draft += untagged_jets[i].momentum.E();
	}

      pt_rest = pt_rest_draft;
      E_rest = E_rest_draft;

      mon.fillHisto("pt", "rest", pt_rest, weight);
      mon.fillHisto("E", "rest", E_rest, weight);

      
      
      //--------------------------------discrete signal region analysis--------------------------------//
      bool isSR1 = false;
      bool isSR2 = false;
      bool isSR3 = false;

      if (b_jets.size()==3 && untagged_jets.size()==0) isSR1 = true;
      else if (b_jets.size()==3 && untagged_jets.size()>0) isSR2 = true;
      else isSR3 = true;
      
      if (isSR1) count_bjets_SR1++;
      else if (isSR2) count_bjets_SR2++;
      else count_bjets_SR3++;


      //#######################################################################################################################################//
      //--------------------------------------------------------------MVA Handler--------------------------------------------------------------//
      //#######################################################################################################################################//

      if(runMVA)
	{
	  myMVAHandler_.getEntry(higgs_m, higgs_pt, higgs_eta, costheta0, dR_higgs_q1, dR_higgs_q2, phi_qq_higgs, avDeltaR, minDeltaM, MassOfbbj, DeltaEta, qq_m, ProductEta, qq_dR, qq_dphi, alpha_qq, N_jet, N_jet_eta_cut, jet1_pt, jet2_pt, jet3_pt, jet4_pt, jet5_pt, jet6_pt, jet1_eta, jet2_eta, jet3_eta, jet4_eta, jet5_eta, jet6_eta, H_T, H_z, H_Tvec, pt_rest, E_rest, N_bjet, bjet1_pt, bjet2_pt, bjet3_pt, bjet4_pt, bjet1_eta, bjet2_eta, bjet3_eta, bjet4_eta, bjet1_btag, bjet2_btag, bjet3_btag, bjet4_btag, N_qjet, qjet1_pt, qjet2_pt, qjet1_eta, qjet2_eta, N_untaggedjet, MET_pt, MET_phi, weight);
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

// Function to calculate the deltaR in the (psi, phi) plane between two four-vectors
double getDeltaRyy(TLorentzVector vec_1, TLorentzVector vec_2)
  
{
  double delta_phi;
  double delta_psi;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_psi = vec_1.Rapidity() - vec_2.Rapidity();

  return std::sqrt(delta_phi * delta_phi + delta_psi * delta_psi);
}

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

// Function that returns the minimum angular distance between a four-vector and a collection of four-vectors along with the position of the four-vector that minimizes the distance
std::pair<double, int> getDeltaRmin(TLorentzVector det_vec, std::vector<TLorentzVector> gen_vec)
{
  //initialize dR and index of minimum dR
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

// Function to generate a random integer number between two values
int rand_int(int min, int max)
{
  //random number generator
  std::random_device rd; 
  std::mt19937 gen(rd()); // Mersenne Twister engine seeded with rd()

  // Define a distribution that generates integers between 0 and 3 inclusive
  std::uniform_int_distribution<> dis(min, max);

  // Generate a random number
  int random_number = dis(gen);

  return random_number;
}

// Function to calculate the minimum mass difference of b-jet pairs
// Input: 1) vec_bjets: a vector that contains all of the bjets momenta as configured with the addition of an extra untagged jet if it exists and if the total number of bjets is less than four
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

// Function to calculate the Mbbj variable.
// Inputs: 1) vec_bjets: vector that contains the btagged jets sorted in their pt
//         2) vec_untaggedjets: vector that contains the untagged jets sorted in their btag value
// Output: Mbbj

double getMassOfbbj(const std::vector<jets_struct>& vec_bjets, const std::vector<jets_struct>& vec_untaggedjets)
{
  double deltaRmin = 1e9;
  double Mbbj = -1; // Initialize with an invalid mass value
  TLorentzVector bb; 
  std::vector<jets_struct> unused_bjets; 

  if (vec_bjets.size() == 3 && vec_untaggedjets.size() == 0)
  {
    // Case 1: Sum all three b-jet momenta
    Mbbj = (vec_bjets[0].momentum + vec_bjets[1].momentum + vec_bjets[2].momentum).M();
  }
  else if (vec_bjets.size() == 3 && vec_untaggedjets.size()>0)
  {
    // Case 2: Find closest b-jet pair and add first untagged jet
    for (size_t i = 0; i < 3; ++i)
    {
      for (size_t j = i + 1; j < 3; ++j)
      {
        TLorentzVector b1 = vec_bjets[i].momentum;
        TLorentzVector b2 = vec_bjets[j].momentum;
        double deltaR = getDeltaR(b1, b2);
        if (deltaR < deltaRmin)
        {
          deltaRmin = deltaR;
          bb = b1 + b2;
        }
      }
    }
    Mbbj = (bb + vec_untaggedjets[0].momentum).M();
  }
  else
  {
    // Case 3: Find closest b-jet pair and select unused b-jet with appropriate btag
    for (size_t i = 0; i < 4; ++i)
    {
      for (size_t j = i + 1; j < 4; ++j)
      {
        TLorentzVector b1 = vec_bjets[i].momentum;
        TLorentzVector b2 = vec_bjets[j].momentum;
        double deltaR = getDeltaR(b1, b2);
        if (deltaR < deltaRmin)
        {
          deltaRmin = deltaR;
          bb = b1 + b2;
          unused_bjets.clear();
          for (size_t k = 0; k < 4; ++k)
          {
            if (k != i && k != j) unused_bjets.push_back(vec_bjets[k]);              
          }
        }
      }
    }
    TLorentzVector unused_bjet = unused_bjets[0].btag < unused_bjets[1].btag ? unused_bjets[0].momentum : unused_bjets[1].momentum;
    Mbbj = (bb+unused_bjet).M();
  }
  return Mbbj;
}

