#define YEAR_2017
#include <iostream>
#include <map>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakePyBind11ParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/TMVAReader.h"
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/bsmhiggs_fwk/interface/PDFInfo.h"
//#include "UserCode/bsmhiggs_fwk/interface/rochcor2016.h" 
#include "UserCode/bsmhiggs_fwk/interface/RoccoR.h" 
//#include "UserCode/bsmhiggs_fwk/interface/muresolution_run2.h" 
#include "UserCode/bsmhiggs_fwk/interface/LeptonEfficiencySF.h"
//#include "UserCode/bsmhiggs_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/bsmhiggs_fwk/interface/BtagUncertaintyComputer.h"
#include "UserCode/bsmhiggs_fwk/interface/METUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/XYMETCorrection.h"
//#include "UserCode/bsmhiggs_fwk/interface/BTagUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/EventCategory.h"
#include "UserCode/bsmhiggs_fwk/interface/statWgt.h"

#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
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

#include <limits>
#include <utility>
#include <vector>
#include <stdexcept>
#include <string>
#include <unistd.h>
using namespace std;

double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2);
bool jet_matched(TLorentzVector jet,std::vector<TLorentzVector> vecb);
double dR_j_fj_min(TLorentzVector jet,std::vector<TLorentzVector> vecfjet);
double dphi_met_j_min(float metphi,std::vector<TLorentzVector> jet,std::vector<TLorentzVector> fjet);
float fjet_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb) ;
double Ht(std::vector<TLorentzVector> vecjet,std::vector<TLorentzVector> vecfjet);
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec);
bool sortBtag( std::pair<TLorentzVector,float>  vec_jet_and_btagi,std::pair<TLorentzVector,float> vec_jet_and_btagj){
  return  vec_jet_and_btagi.second >  vec_jet_and_btagj.second ;
}

int main(int argc, char* argv[])
{
  //##################################################################################
  //##########################    GLOBAL INITIALIZATION     ##########################
  //##################################################################################
  gSystem->Load("libTree.so");
  //check arguments
  if(argc<2) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    exit(0);
  }
   
  //load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  // AutoLibraryLoader::enable();
  FWLiteEnabler::enable();
  MVAHandler myMVAHandler_;
  myMVAHandler_.initTree();
  //configure the process
  const edm::ParameterSet &runProcess = edm::cmspybind11::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

  bool run0lep = runProcess.getParameter<bool>("run0lep");
  
  // will produce the input root trees to BDT training (optional)
  bool runMVA = runProcess.getParameter<bool>("runMVA");
  
  bool isMC = runProcess.getParameter<bool>("isMC");
  double xsec = runProcess.getParameter<double>("xsec");
  double nevts = runProcess.getParameter<double>("nevts");


  TString proc=runProcess.getParameter<std::string>("proc");
  TString dtag=runProcess.getParameter<std::string>("tag");
  TString suffix=runProcess.getParameter<std::string>("suffix");

  bool is2017data = (!isMC && dtag.Contains("2017"));
  bool is2017MC = (isMC && dtag.Contains("2017"));
  bool is2018data = (!isMC && dtag.Contains("2018"));
  bool is2018MC = (isMC && dtag.Contains("2018"));
  bool is2017_2018 = (is2017MC || is2017data || is2018MC || is2018data);

  bool verbose = runProcess.getParameter<bool>("verbose");
  
  TString url = runProcess.getParameter<std::string>("input");
  TString outFileUrl( dtag ); 
  gSystem->BaseName(url);
  
  
  
  bool isMC_ZVV = isMC && ( (string(url.Data()).find("ZJetsToNuNu")  != string::npos) );
  bool isMC_DY_HTbin = isMC_ZVV && dtag.Contains("HT") ;

  bool isMC_ttbar = isMC && (string(url.Data()).find("TTTo")  != string::npos);
  if(is2017MC || is2018MC) isMC_ttbar = isMC && (string(url.Data()).find("TTTo")  != string::npos);
  

  bool isMC_Zh = isMC && (string(url.Data()).find("_Zh_amass")  != string::npos); 
  bool isSignal = isMC_Zh; 

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
  

 

  SmartSelectionMonitor mon;
  TH1F *drbf= new TH1F ("dr-min-bjet-fjet", ";dr;Events", 100,0,8);
  TH1F *drjf= new TH1F ("dr-min-jet-fjet", ";dr;Events", 100,0,8);
  TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 6,0,6) );
  h->GetXaxis()->SetBinLabel(1,"Raw");
  if (run0lep) {
  h->GetXaxis()->SetBinLabel(2,"0 lepton");
  h->GetXaxis()->SetBinLabel(3,"MET cut>170");  
  } else {
  h->GetXaxis()->SetBinLabel(2,"2 leptons");
  h->GetXaxis()->SetBinLabel(3,"Z-mass window"); 
   }
  
  h->GetXaxis()->SetBinLabel(4,"=2 AK4 b-jets");
  h->GetXaxis()->SetBinLabel(5,">=3 AK4 b-jets");
  h->GetXaxis()->SetBinLabel(6,">=1 f-jet and 1 AK4 b-jet");
 
  

  TH1F *hr= new TH1F ("raw_eventflow", ";;Events", 6,1,6);
  hr->GetXaxis()->SetBinLabel(1,"Raw");
  if (run0lep) {
    hr->GetXaxis()->SetBinLabel(2,"0 lepton");
    hr->GetXaxis()->SetBinLabel(3,"MET cut>170");
  } else {
    hr->GetXaxis()->SetBinLabel(2,"2 leptons");
    hr->GetXaxis()->SetBinLabel(3,"Z-mass window");
  }
  
  hr->GetXaxis()->SetBinLabel(4,"=2 AK4 b-jets");
  hr->GetXaxis()->SetBinLabel(5,">=3 AK4 b-jets");
  hr->GetXaxis()->SetBinLabel(6,">=1 f-jet and 1 AK4 b-jet");
 
  //MVA BDT
  TH1F *bdt= new TH1F ("bdt", ";BDT;Events", 100, -1,1);
  mon.addHistogram( new TH1F( "bdt", ";BDT;Events", 100, -1,1) );

  //reconstruction level
  //dR between leptons-jets/leptons-fjets/jets-fjets: before and after cross-cleaning
   // -- -- -- Delta R -- -- --
  //dR jets-leptons
  mon.addHistogram( new TH1F("dR",";#Delta R;Events", 100, 0, 8));
 
  //dR fjets-leptons 
 
 
  //dR jets-fjets
  mon.addHistogram( new TH1F("dR_j-fj", "jet-fjet", 100, 0, 8));
  mon.addHistogram( new TH1F("dR_j-fj_min", "jet-fjet_min", 100, 0,8));
  //2-lepton inv mass
  mon.addHistogram( new TH1F("inv_mass", "#it{m}_{ll} [GeV];Events", 80, 0, 200));
  //m bbb
  mon.addHistogram( new TH1F("m_bbb", ";#it{m}_{bbb} [GeV];Events", 50, 0, 500));
  mon.addHistogram( new TH1F("pt_bbb", ";p_{T}_{bbb} [GeV];Events", 60, 0, 600));
  mon.addHistogram( new TH1F("pt_bb", ";p_{T}_{bb} [GeV];Events", 50, 0, 500));
   // fjets
  //fjet mult
  mon.addHistogram( new TH1F("mult", ";multiplicity;Events", 10, 0, 10));

  //fjet kinematics
  mon.addHistogram( new TH1F("pt", ";p_{T} [GeV];Events", 50, 0, 500));
   mon.addHistogram( new TH1F("M", ";M [GeV];Events", 50, 0, 500));
  mon.addHistogram( new TH1F("eta", ";#eta ;Events", 50, -5, 5));
  mon.addHistogram( new TH1F("phi", ";#phi;Events", 50, -3*TMath::Pi(), 3*TMath::Pi()));
  //fjet matched pt vs met(vv)
  
  mon.addHistogram( new TH2F("pt_bb_vs_met", "gen/reco",50,0,500, 50, 0,500));
  mon.addHistogram(new TProfile("pt_fjet_vs_met_prof", "gen/reco",50,0,500));
  //fjet matched pt
 
  //fjet pt vs pt(bb)  
  mon.addHistogram( new TH2F("fj_pt_bb_pt", "fj1_pt_bb_pt",50,0,500, 50, 0,500));
 
  
  //fjet and fjet matched soft drop mass 
  mon.addHistogram(new TH1F("fj_sd_mass",";soft drop mass [GeV];Events",100,0,200));
  mon.addHistogram( new TH2F("fj_pt_sd_mass", ";P_{T} [GeV];soft drop mass [Gev]",50,0,500, 100, 0,200));
  //B TAG JET
   mon.addHistogram( new TH1F("btag1", ";b-tag1;Events", 50, 0, 1));
  //fjet discriminants
  mon.addHistogram( new TH1F("fjet_btagXbb", "fjet1_btagXbb", 50, 0, 1));
  mon.addHistogram( new TH1F("fjet_btagXbbXccXqq", "fjet1_btagXbbXccXqq", 50, 0, 1));
 
  
  mon.addHistogram( new TH2F("fj1_fj2_btagXbb", "fjet_btagXbb",50,0,1, 50, 0, 1));
  mon.addHistogram( new TH2F("fj1_fj2_btagXbbXccXqq", "fjet_btagXbbXccXqq",50,0,1,50, 0, 1));

  mon.addHistogram(new TH2F("fjet_Xbb_pt", "fjet1_btagXbb", 100,0,500,50,0, 1));
  mon.addHistogram( new TH2F("fjet_XbbXccXqq_pt", "fjet1_btagXbb", 100,0,500,50,0, 1));
  mon.addHistogram( new TH1F("dphi",";H-Z #Delta #phi ;Events", 50, -TMath::Pi(), TMath::Pi()));
 
  //fjet subjet count
 
 mon.addHistogram (new TH1F("subcount", "sub_mult1", 5, 0, 5));
  
					       
 
 
   //MET
  //total met
  mon.addHistogram(new TH1F("met",";MET [GeV];Events",50,0,500));
  mon.addHistogram(new TH1F("met_phi",";#phi;Events",100,-TMath::Pi(), TMath::Pi()));
 
  						   
 //Ht
  mon.addHistogram(new TH1F("Ht",";H_{T} [GeV];Events",150,0,1500));
						   
   
  						   
  
  //tree info
  int evStart     = runProcess.getParameter<int>("evStart");
  int evEnd       = runProcess.getParameter<int>("evEnd");
  TString dirname = runProcess.getParameter<std::string>("dirName");

  
  //##################################################################################
  //#############         GET READY FOR THE EVENT LOOP           #####################
  //##################################################################################
  //open the file and get events tree
  DataEvtSummaryHandler summaryHandler_;
  
  TFile *file = TFile::Open(url);
  printf("Looping on %s\n",url.Data());
  if(file==0) { return -1; printf("file is 0"); }
  
  if(file->IsZombie()) return -1;
  if( !summaryHandler_.attachToTree( (TTree *)file->Get(dirname) ) ) {
    file->Close();
    return -1;
  }
  //check run range to compute scale factor (if not all entries are used)
  const Int_t totalEntries= summaryHandler_.getEntries();
  float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
  if(evEnd<0 || evEnd>summaryHandler_.getEntries() ) evEnd=totalEntries;
  if(evStart > evEnd ) {
    file->Close();
    return -1;
  }
  
  //------Event Counters-------------
  Int_t n_event_lepton_test = 0;
  Int_t n_fjets=0;
  Int_t INVM=0;
  Int_t METCUT=0;
  Int_t sr1=0;
  Int_t sr2=0;
  //z->:
  Int_t zqq_light=0;
  Int_t zbb=0;
  Int_t zvv=0;
  Int_t zll=0;
  
  Int_t st4_qq=0;
  Int_t st4_bb=0;
  Int_t st4_vv=0;
  Int_t st4_ll=0;

  double totalNumberofEvents(0.);
  TH1F* nevtH = (TH1F *) file->Get("mainNtuplizer/nevents");
  totalNumberofEvents = nevtH->GetBinContent(1);
  hr->SetBinContent(1, totalNumberofEvents);
  float Lint=43.5*1000;
  double nev_exp=xsec*Lint;
  double weight=nev_exp/nevts;
   //####################################################################################################################
    //###########################################           MVAHandler         ###########################################
    //####################################################################################################################
    //construct MVA out put file name
  // TString mvaout = TString ( runProcess.getParameter<std::string>("outdir") ) + "/mva_" + outFileUrl + ".root";
  // MVAHandler myMVAHandler_;
  // if (runMVA) { myMVAHandler_.initTree(mvaout); }

  //####################################################################################################################
    //###########################################           TMVAReader         ###########################################
    //####################################################################################################################
    
  std::string chpath = "MVA_2lepton/";
   if (run0lep) chpath =  "MVA_0lepton/";

   TMVAReader SR1Reader;
   SR1Reader.InitTMVAReader();
   std::string SR1_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA1_BDT.weights.xml"; 
   SR1Reader.SetupMVAReader( "SR1Class", SR1_xml_path );

   TMVAReader SR2Reader;
   SR2Reader.InitTMVAReader();
   std::string SR2_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA2_BDT.weights.xml"; 
   SR2Reader.SetupMVAReader( "SR2Class", SR2_xml_path );

   TMVAReader SR3Reader;
   SR3Reader.InitTMVAReader();
   std::string SR3_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/"+chpath+"TMVA3_BDT.weights.xml"; 
   SR3Reader.SetupMVAReader( "SR3Class", SR3_xml_path );
    
  
  //####################################################################################################################
  //###########################################           EVENT LOOP         ###########################################
  //####################################################################################################################
  //loop on all the events
  int treeStep = (evEnd-evStart)/50;
  if(treeStep==0)treeStep=1;
  DuplicatesChecker duplicatesChecker;
  int nDuplicates(0);
  
  for( int iev=evStart; iev<evEnd; iev++) 
    {

      //if ( verbose ) printf("\n\n Event info %3d: \n",iev);
      //load the event content from tree
      summaryHandler_.getEntry(iev);
      DataEvtSummary_t &ev=summaryHandler_.getEvent();
      if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) {
	nDuplicates++;
	cout << "nDuplicates: " << nDuplicates << endl;
	continue;
      }
      //tree variables                                                                                                                                
      
      // Calculate xsec weight:                                                                                                                       
      

     // Create new object vectors after configuration
      // Jets
      std::vector<TLorentzVector> vec_jet;
      std::vector<TLorentzVector> vec_bjets;
      std::vector<TLorentzVector> vec_bjets_cc;
      std::vector<TLorentzVector> vec_fjet;
      std::vector<TLorentzVector> vec_jet_cc;
      // Leptons
      std::vector<TLorentzVector> vec_muons;
      std::vector<TLorentzVector> vec_ele;
      std::vector<TLorentzVector> vec_leptons;
    
      //generator level
      //before accepance cuts 
      std::vector<TLorentzVector> vec_genmu_bef;
      std::vector<TLorentzVector> vec_genele_bef;
      std::vector<TLorentzVector> vec_genlep_bef;
      std::vector<TLorentzVector> vec_genqg_bef; //quarks or gluons
      std::vector<TLorentzVector> vec_genb_bef; // only b's
    
      //after cuts
      std::vector<TLorentzVector> vec_genmu;
      std::vector<TLorentzVector> vec_genele;
      std::vector<TLorentzVector> vec_genlep;
      std::vector<TLorentzVector> vec_gena;
      std::vector<TLorentzVector> vec_h;
      std::vector<TLorentzVector> vec_z;
      std::vector<TLorentzVector> vec_genqg; //quarks or gluons
      std::vector<TLorentzVector> vec_genb; // only b's
      std::vector<TLorentzVector> vec_genbb1_truth;
      std::vector<TLorentzVector> vec_genbb2_truth;

      std::vector<TLorentzVector> vec_zll;
      std::vector<TLorentzVector> vec_zbb;
      std::vector<TLorentzVector> vec_zqq_light;
      std::vector<TLorentzVector> vec_zvv;
      std::vector<std::pair<TLorentzVector,int> > bjet_index;
      std::vector<std::pair<TLorentzVector,int> > fjet_index;
      std::vector<std::pair<TLorentzVector,int> > bjet_index_cc;
      std::vector<std::pair<TLorentzVector,float> > vec_jet_and_btag;
      //start generator level loop
      for (int imc=0; imc<ev.nmcparticles;imc++)	{
       if(verbose && iev <3) 
       {
	 std::cout << " imcparticle " << imc << " : is a " << ev.mc_id[imc] <<",has a mom:"<<ev.mc_mom[imc]<< " , and has a mother at: " << ev.mc_momidx[imc]  <<"has status   "<<ev.mc_status[imc] << "  and has a 4-vector p = (" << ev.mc_en[imc] << ", " << ev.mc_px[imc] << ", " << ev.mc_py[imc] << ", " << ev.mc_pz[imc] << " ) " << std::endl;
       }
       
        if(ev.mc_id[imc]==25) 
            { // found the Higgs boson
	      TLorentzVector pH;
	      pH.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);

              double ptH  = pH.Pt();
              double etaH = pH.Eta();
              double phiH = pH.Phi();
              double mH   = pH.M();

              // Make histograms of the pt, eta, phi, mass of the Higgs boson:
              mon.fillHisto("pt","H",pH.Pt(),weight); 
              mon.fillHisto("eta","H",etaH,weight); mon.fillHisto("phi","H",phiH,weight);
	      mon.fillHisto("M","Z",pH.M(),weight);
	      vec_h.push_back(pH);
              
            }
        if(ev.mc_id[imc]==23) 
            { // found the z boson

              TLorentzVector pz;
	            pz.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
              // Make histograms of the pt, eta, phi, mass of the z boson:
		    mon.fillHisto ("pt","Z",pz.Pt(),weight); 
		    mon.fillHisto("eta","Z",pz.Eta(),weight); mon.fillHisto("phi","Z",pz.Phi(),weight);
		    mon.fillHisto("M","Z",pz.M(),weight);
		    vec_z.push_back(pz);
            }
        if(ev.mc_id[imc]==36) 
            { // found the A
	      
              TLorentzVector pa;
	      pa.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
              // Make histograms of the pt, eta, phi, mass of the A:
	      mon.fillHisto("pt","A",pa.Pt(),weight); 
	      mon.fillHisto("eta","A",pa.Eta(),weight);mon.fillHisto("phi","A",pa.Phi(),weight);
	      mon.fillHisto("phi","A",pa.M(),weight);
	      vec_gena.push_back(pa);
            }

           
        if(abs(ev.mc_id[imc])==11) 
            { // electrons   
              TLorentzVector lep1;
              lep1.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
              vec_genele_bef.push_back(lep1);
              vec_genlep_bef.push_back(lep1);
              if (lep1.Pt()>20. && fabs(lep1.Eta())<2.4) 
                { // acceptance cuts
                  vec_genele.push_back(lep1);
                  vec_genlep.push_back(lep1);
                }
            }
        // electrons END
        if(abs(ev.mc_id[imc])==13) 
            { //muons
              TLorentzVector lep2;
              lep2.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
              vec_genmu_bef.push_back(lep2);
              vec_genlep_bef.push_back(lep2);
              if (lep2.Pt()>20. && fabs(lep2.Eta())<2.4) 
                {// acceptance cuts
                  vec_genmu.push_back(lep2);
                  vec_genlep.push_back(lep2);
                }
              
            } // muons END
        
        //quarks or gluons
       
        if( ((abs(ev.mc_id[imc])>=1 && abs(ev.mc_id[imc])<=5) || abs(ev.mc_id[imc])==21) &&  (ev.mc_status[imc]==23 || ev.mc_status[imc]==1))
	  {
            TLorentzVector pb;
            pb.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
            vec_genqg_bef.push_back(pb);
            if (abs(ev.mc_id[imc])==5) vec_genb_bef.push_back(pb);
	    // acceptance cuts :
            if (pb.Pt()>20. && fabs(pb.Eta())<2.4) 
            {
              vec_genqg.push_back(pb);
              if (abs(ev.mc_id[imc])==5)  vec_genb.push_back(pb);
	      //bb pair from same a boson 
              
              if (abs(ev.mc_id[imc])==5 && ev.mc_mom[imc]==36 && (ev.mc_momidx[imc]==4))  
                { 
                  vec_genbb1_truth.push_back(pb);
                }
              if (abs(ev.mc_id[imc])==5 && ev.mc_mom[imc]==36 && (ev.mc_momidx[imc]==5))  
                { 
                  vec_genbb2_truth.push_back(pb);
                }
		
	    } // end acceptance cuts
	 }// if q/g END    

         // z->qqlight/z->bb/z->vv/z->ll
         if ((abs(ev.mc_id[imc])>=1 && abs(ev.mc_id[imc])<=4)  && ev.mc_mom[imc]==23)
         {
            TLorentzVector q;
            q.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
            vec_zqq_light.push_back(q);
         } else if (abs(ev.mc_id[imc])==5  && ev.mc_mom[imc]==23)
         {
            TLorentzVector b;
            b.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
            vec_zbb.push_back(b);
         }
         if( (abs(ev.mc_id[imc])==12 || abs(ev.mc_id[imc])==14 || abs(ev.mc_id[imc])==16) && ev.mc_mom[imc]==23)
         {
            TLorentzVector vv;
            vv.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
            vec_zvv.push_back(vv);
            
         }
         if( (abs(ev.mc_id[imc])==11 || abs(ev.mc_id[imc])==13 || abs(ev.mc_id[imc])==15) && ev.mc_mom[imc]==23)
         {
            TLorentzVector ll;
            ll.SetPxPyPzE(ev.mc_px[imc],ev.mc_py[imc],ev.mc_pz[imc],ev.mc_en[imc]);
            vec_zll.push_back(ll);
         } 
      } // MC particle loop END
      
      // find dr between aas
      if (vec_gena.size()>1) mon.fillHisto("dR","AA",getDeltaR(vec_gena[0],vec_gena[1]),weight);
      
      //sort qluon-quark vector
      vec_genqg= sort_vec_pt( vec_genqg);
      //sort bquark vector
      vec_genb=sort_vec_pt(vec_genb);
      //check z->;
      bool is_vv(false);
      bool is_qq_light(false);
      bool is_bb(false);
      bool is_ll(false);
      if (vec_zll.size()>1)
	{ zll++;
	  is_ll=true;}
      if (vec_zqq_light.size()>1) 
	{ zqq_light++;
	  is_qq_light=true;}
      if (vec_zbb.size()>1) 
	{ zbb++;  
	  is_bb=true;}  
      if (vec_zvv.size()>1) 
	{ zvv++;  
	  is_vv=true;}  
      //mon.fillHisto("pt","a2",vec_gena[1].Pt(),weight);
      //mon.fillHisto("pt","a1",vec_gena[0].Pt(),weight);
      //plot dr and  vector sum pt of bbs of same mom and z->vv event
      if (vec_genbb1_truth.size()>1 ){
	mon.fillHisto ("dR","bb1",getDeltaR(vec_genbb1_truth[0],vec_genbb1_truth[1]),weight);
	mon.fillHisto("pt","bb1",(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
      }
      if (vec_genbb2_truth.size()>1 ){
	mon.fillHisto("dR","bb2",getDeltaR(vec_genbb2_truth[0],vec_genbb2_truth[1]),weight);
	mon.fillHisto ("pt","bb2",(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
      
       
      }
      
      //plot vec sum pt of vvs
      if(is_vv){
	mon.fillHisto("pt","gen_vv",(vec_zvv[0]+vec_zvv[1]).Pt(),weight);
	if(vec_genbb1_truth.size()>1){
          mon.fillHisto("pt_bb_vs_met","1_gen",(vec_zvv[0]+vec_zvv[1]).Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
	  mon.fillProfile("pt_fjet_vs_met_prof","1_gen",(vec_zvv[0]+vec_zvv[1]).Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
	  }
        else if(vec_genbb2_truth.size()>1){
	  mon.fillHisto("pt_bb_vs_met","2_gen",(vec_zvv[0]+vec_zvv[1]).Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
  mon.fillProfile("pt_fjet_vs_met_prof","2_gen",(vec_zvv[0]+vec_zvv[1]).Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
         }
       } 
   
    
    // Make Histograms: lepton mulciplicity, q/g multiplicity , b-quark multiplicity --> Check effect of the Acceptance cuts (pt/eta)
    mon.fillHisto("mult","en_gen_bef",vec_genele_bef.size(),weight);
    mon.fillHisto("mult","mn_gen_bef",vec_genmu_bef.size(),weight);
    mon.fillHisto("mult","leptgen_bef",vec_genlep_bef.size(),weight);
    mon.fillHisto("mult","qg_bef",vec_genqg_bef.size(),weight);
    mon.fillHisto("mult","bq_bef",vec_genb_bef.size(),weight);
    //after acceptance cuts
    mon.fillHisto("mult","en_gen",vec_genele.size(),weight);
    mon.fillHisto("mult","mn_gen",vec_genmu.size(),weight);
    mon.fillHisto("mult","leptgen",vec_genlep.size(),weight);
    mon.fillHisto("mult","qg",vec_genqg.size(),weight);
    mon.fillHisto("mult","bq",vec_genb.size(),weight);
    
    
    // END GENERATOR LEVEL ANALYSIS
    
    //start reco level analysis
     // RECONSTRUCTION LEVEL ANALYSIS:


    for (int i = 0; i <ev.mn; i++)
    {
      TLorentzVector p_muon;
      p_muon.SetPxPyPzE(ev.mn_px[i], ev.mn_py[i], ev.mn_pz[i], ev.mn_en[i]);

      if (p_muon.Pt() < 20. || std::fabs(p_muon.Eta()) > 2.4) continue;

      // id + Isolation
      if (ev.mn_passId[i] && ev.mn_passIso[i])
      {
        vec_muons.push_back(p_muon);
	     
        vec_leptons.push_back(p_muon);
      }
    }

    for (int i = 0; i <ev.en; i++)
    {
      TLorentzVector p_electron;
      p_electron.SetPxPyPzE(ev.en_px[i], ev.en_py[i], ev.en_pz[i], ev.en_en[i]);

      if (p_electron.Pt() < 20. || std::fabs(p_electron.Eta()) > 2.4) continue;

      if (ev.en_passId[i] && ev.en_passIso[i])
      {
        vec_ele.push_back(p_electron);
	vec_leptons.push_back(p_electron);
      }
    }

     // jets & cross cleaning
   
   
    for (int i = 0; i < ev.jet; i++)
    {
      bool overlap = false;
      double dR1 = 0;
      
      
      TLorentzVector p_jet;
      p_jet.SetPxPyPzE(ev.jet_px[i], ev.jet_py[i], ev.jet_pz[i], ev.jet_en[i]);

      if (p_jet.Pt() < 20. || std::fabs(p_jet.Eta()) > 2.5) continue;
      for (int mn_count = 0; mn_count < vec_muons.size(); mn_count++)
      {
        dR1 = getDeltaR(p_jet, vec_muons[mn_count]);
        mon.fillHisto ("dR","jet_mn_before",dR1,weight);
        if (dR1 < 0.4){
          overlap = true;
         }
         else {
	   mon.fillHisto("dR","jet_mn_after",dR1,weight);
         }
      }

      for (int en_count = 0; en_count < vec_ele.size(); en_count++)
      {
        dR1 = getDeltaR(p_jet, vec_ele[en_count]);
        mon.fillHisto("dR","jet_en_before",dR1,weight);
        if (dR1 < 0.4){
         overlap = true;
        }
        else {
	  mon.fillHisto("dR","jet_en_after",dR1,weight);
        }
      }

      //      if (!overlap) vec_jet.push_back(p_jet);
        
       //bjets
      if (! overlap ) {
	vec_jet.push_back(p_jet);
	vec_jet_and_btag.push_back(make_pair(p_jet,ev.jet_btag1[i]));

	if (ev.jet_btag1[i] > 0.4941)
	  {
	    vec_bjets.push_back(p_jet);
	    bjet_index.push_back(make_pair(p_jet,i));
	  }
      }

    } // end jet loop
    //sort the b-tagged jets by discriminator value   
    std::sort(vec_jet_and_btag.begin(), vec_jet_and_btag.end(), sortBtag);

 
    //fjets
    
    for (int i = 0; i < ev.fjet; i++)
      {  
         bool overlap = false;
         double dR1 = 0;
         TLorentzVector p_fjet;
         p_fjet.SetPxPyPzE(ev.fjet_px[i], ev.fjet_py[i], ev.fjet_pz[i], ev.fjet_en[i]);
         if(ev.fjet_subjet_count[i]<2)continue;
         if (p_fjet.Pt() <150 || std::fabs(p_fjet.Eta()) > 2.5 ) continue;
         if((ev.fjet_btag10[i]+ev.fjet_btag11[i]+ev.fjet_btag12[i])/(ev.fjet_btag10[i]+ev.fjet_btag11[i]+ev.fjet_btag12[i]+ev.fjet_btag13[i]+ev.fjet_btag14[i]+ev.fjet_btag15[i]+ev.fjet_btag16[i]+ev.fjet_btag17[i])<0.5)continue;
         for (int mn_count = 0; mn_count < vec_muons.size(); mn_count++)
         {
         dR1 = getDeltaR(p_fjet, vec_muons[mn_count]);
         mon.fillHisto("dR","fjet_mn_before",dR1,weight);
         if (dR1 < 0.8){
            overlap = true;
            }
            else {
	      mon.fillHisto("dR","fjet_mn_after",dR1,weight);
            }
         }

         for (int en_count = 0; en_count < vec_ele.size(); en_count++)
         {
         dR1 = getDeltaR(p_fjet, vec_ele[en_count]);
         mon.fillHisto("dR","fjet_en_before",dR1,weight);
         if (dR1 < 0.8){
            overlap = true;
         }
         else {
	   mon.fillHisto("dR","fjet_en_after",dR1,weight);
         }
         }

         if (!overlap)
         { 
          vec_fjet.push_back(p_fjet);
	  fjet_index.push_back(make_pair(p_fjet,i));
         }
      }
 
    
    // do it for vec_jet as well
    for (int j = 0; j < vec_jet.size(); j++)  {
      drjf->Fill(dR_j_fj_min(vec_jet[j],vec_fjet));
      if(dR_j_fj_min(vec_jet[j],vec_fjet)>0.8) {
        vec_jet_cc.push_back(vec_jet[j]);
      }
    }
    
      //jet f jet cross-cleaning
   
    for (const auto& bjet: bjet_index) {
      drbf->Fill(dR_j_fj_min(bjet.first,vec_fjet));
      if(dR_j_fj_min(bjet.first,vec_fjet)>0.8) {
	bjet_index_cc.push_back(make_pair(bjet.first,bjet.second));
	vec_bjets_cc.push_back(bjet.first);
      }
    }
    //-----End config ------------
    //raw events 
    
     mon.fillHisto("eventflow","histo",0,weight);
    
       // lepton cuts---------------------------------------------------

    if(run0lep) {
     if (vec_leptons.size() !=0) continue;
     n_event_lepton_test++;
     //mon.fillHisto("eventflow","histo",1,weight);
     hr->Fill(2);
    }
    else {
      if (vec_leptons.size() < 2) continue;
      n_event_lepton_test++;
      mon.fillHisto("eventflow","histo",1,weight);
      //mon.fillHisto("raw_eventflow","histo",1);
      hr->Fill(2);
    }
   
    // pt(bb)vs MET                                                                                                                                                             
    //(fjet1_pt vs met                                                                                                                                                  
    if(vec_fjet.size()>0){
      mon.fillHisto("pt_bb_vs_met","1_reco",ev.met_pt,vec_fjet[0].Pt(),weight);
      mon.fillProfile("pt_fjet_vs_met_prof","1_reco",ev.met_pt,vec_fjet[0].Pt(),weight);
    }
    //(fjet2_pt vs met                                                                                                                
    if(vec_fjet.size()>1){
      mon.fillHisto("pt_bb_vs_met","2_reco",ev.met_pt,vec_fjet[1].Pt(),weight);
      mon.fillProfile("pt_fjet_vs_met_prof","2_reco",ev.met_pt,vec_fjet[1].Pt(),weight);
    }
    //met
    mon.fillHisto("met","total",ev.met_pt,weight);
    mon.fillHisto("met_phi","total",ev.met_phi,weight);
    if(is_qq_light) mon.fillHisto("met","qq",ev.met_pt,weight);
    if(is_bb)mon.fillHisto("met","bb",ev.met_pt,weight);
    if(is_vv) mon.fillHisto("met","vv",ev.met_pt,weight);
    
    //trigger 
    bool hasMETtrigger1=(ev.triggerType>>11)&0x1;
    bool hasMETtrigger2=(ev.triggerType>>12)&0x1;
    

    if(hasMETtrigger1||hasMETtrigger2){
      mon.fillHisto("met","trig",ev.met_pt,weight);
      mon.fillHisto("met_phi","trig",ev.met_phi,weight);
      if(is_vv)mon.fillHisto("met","vv_trig",ev.met_pt,weight); 
    }
    
    // control plots mult
    mon.fillHisto("mult","en",vec_ele.size(),weight);
    mon.fillHisto("mult","mn",vec_muons.size(),weight);
    mon.fillHisto("mult","lept",vec_leptons.size(),weight);
 
    mon.fillHisto("mult","bjet_cc_1",vec_bjets_cc.size(),weight);
    mon.fillHisto("mult","fjet_1",vec_fjet.size(),weight);
    mon.fillHisto("mult","jet_cc_1",vec_jet_cc.size(),weight);
    //step 2 MET/INV MASS CUT   
        
    if(run0lep) {
       // metcut
      
      if (ev.met_pt < 170.) continue;
      METCUT++;
      mon.fillHisto("eventflow","histo",2,weight);
      //mon.fillHisto("raw_eventflow","histo",2);
      hr->Fill(3);
      mon.fillHisto("mult","fjet_2",vec_fjet.size(),weight);
      mon.fillHisto("mult","jet_cc_2",vec_jet_cc.size(),weight);
      mon.fillHisto("mult","bjet_cc_2",vec_bjets_cc.size(),weight);
    }
    else {
      // inv mass cuts
      
      float dileptonmass= (vec_leptons[0]+vec_leptons[1]).M();
      if (dileptonmass>100. || dileptonmass < 80.) continue;
      INVM++;
      mon.fillHisto("eventflow","histo",2,weight);
      //mon.fillHisto("raw_eventflow","histo",2);
      hr->Fill(3);
      mon.fillHisto("mult","fjet_2",vec_fjet.size(),weight);
      mon.fillHisto("mult","jet_cc_2",vec_jet_cc.size(),weight);
      mon.fillHisto("mult","bjet_cc_2",vec_bjets_cc.size(),weight);
    }
   
    
    //step 3 : // 
    Float_t m4b(-1) ,pt4b(-1.),m4b_SR1(-1) ,pt4b_SR1(-1.),ptf1(-1.),ptb1(-1.),ptb1_SR1(-1.),ptb2(-1),ht(-1.),dilep_pt(-1.),drll(-1.),dphiHZ(-1.),n_ad_j(-1.),met(-1.),
      dilep_pt_SR1(-1.),drll_SR1(-1.),dphiHZ_SR1(-1.),n_ad_j_SR1(-1.),ht_SR1(-1),met_SR1(-1.),btag1(-1),
      btag3(-1), btag3_SR1(-1),sd_mass1(-1.),xbb1(-2.),xbbccqq1(-2.),drjj(-1),dphi_met_j(-1),dphi_met_l(-1),drjj_SR1(-1),dphi_met_j_SR1(-1),dphi_met_l_SR1(-1);
    float mvaBDT1(-10.0),mvaBDT2(-10.0);
    bool isSR2(false);
    bool isSR3(false);
    bool isSR1(false);
    if(vec_bjets.size()==2)isSR2=true;
    if(vec_bjets.size()>=3)isSR3=true;
    if (vec_fjet.size()>=1 && vec_bjets_cc.size()>=1)isSR1=true;
    
    if(isSR2||isSR3)
      {
	

	bool matched(false);
        bool matched1(false);bool matched2(false);

        if(isSignal)
          {
            if(fjet_matched(vec_bjets[0],vec_genbb1_truth)>0)
              {
      
                matched1=true;
              }
	    else if( fjet_matched(vec_bjets[0],vec_genbb2_truth)>0)
	      {
		matched2=true;
	      }
            if(matched1||matched2)matched=true;
          }
        else { // backgrounds:                                                                                                                                                                     
          if( fjet_matched(vec_bjets[0],vec_genqg)>0){
       
            matched=true;
	  }
        }
        if(matched){
          mon.fillHisto("pt","bb_AK4",vec_bjets[0].Pt(),weight);
        }

	ptb1=vec_bjets[0].Pt();
	ptb2=vec_bjets[1].Pt();
	TLorentzVector tot_bjets;
	for(std::vector<TLorentzVector>::size_type i=0; i<vec_bjets.size(); i++){
	  if(i==4) break;
	   tot_bjets += vec_bjets[i];
	}
	if(!run0lep){
	  mon.fillHisto("inv_mass","dilept",(vec_leptons[0]+vec_leptons[1]).M(),weight);
	  mon.fillHisto("dR","dilept",getDeltaR(vec_leptons[0],vec_leptons[1]),weight);
	  drll=getDeltaR(vec_leptons[0],vec_leptons[1]);
	  dilep_pt=(vec_leptons[0]+vec_leptons[1]).Pt();
	  dphi_met_l=min(fabs(vec_leptons[0].Phi()-ev.met_phi),fabs(vec_leptons[1].Phi()-ev.met_phi));
	  if(fabs( tot_bjets.Phi()-(vec_leptons[0]+vec_leptons[1]).Phi())<TMath::Pi())
		  {
		    dphiHZ=std::fabs( tot_bjets.Phi()-(vec_leptons[0]+vec_leptons[1]).Phi());
		  }
		else{
		  dphiHZ=2*TMath::Pi()-std::fabs( tot_bjets.Phi()-(vec_leptons[0]+vec_leptons[1]).Phi());
		}
	  dphi_met_j=999;
	  for(std::vector<TLorentzVector>::size_type i=0; i<vec_jet.size(); i++)if (fabs( deltaPhi(vec_jet[i].Phi(),ev.met_phi))<dphi_met_j)dphi_met_j=fabs( deltaPhi(vec_jet[i].Phi(),ev.met_phi));
	   
	}
	drjj=getDeltaR(vec_bjets[0],vec_bjets[1]);
	met=ev.met_pt;
	n_ad_j=vec_jet.size();
	btag1=vec_jet_and_btag[0].second;
	pt4b = tot_bjets.Pt();
	m4b = tot_bjets.M();
	ht = 0;
	
	for(std::vector<TLorentzVector>::size_type i=0; i<vec_jet.size(); i++){
	  ht += vec_jet[i].Pt();
	}
	
	if(isSR2){
	  btag3=vec_jet_and_btag[1].second;
	  mon.fillHisto("eventflow","histo",3,weight);
	  hr->Fill(4);
	}
	else if(isSR3){
	  btag3=vec_jet_and_btag[2].second;
	  mon.fillHisto("eventflow","histo",4,weight);
          hr->Fill(5);
	}
	//MVAREADER
	
	if(isSR2)
	  {
	    if(!run0lep)
	      {
		mvaBDT2 = SR2Reader.GenReMVAReader
		  (
		   dilep_pt, drll,dphiHZ,dphi_met_l,dphi_met_j,
		   m4b,pt4b, met,ht,drjj,n_ad_j,
	       ptb1,ptb2 ,btag1,btag3,
		   "SR2Class"
		   );
	      }
	    else if(run0lep)
	      {
		mvaBDT2 = SR2Reader.GenReMVAReader
                  (
                   
                   m4b,pt4b, met,ht,drjj,n_ad_j,
		   ptb1,ptb2 ,btag1,btag3,
                   "SR2Class"
                   );
	      }
	    mon.fillHisto("bdt","sr2",mvaBDT2,weight);
	 }
	else if (isSR3)
	  {
	    if(!run0lep){
	      mvaBDT2 = SR3Reader.GenReMVAReader
		(
		 dilep_pt, drll,dphiHZ,dphi_met_l,dphi_met_j,
		 m4b,pt4b, met,ht,n_ad_j,
		 ptb1,ptb2 ,btag1,btag3,
		 "SR3Class"
		 );
	    }
	    else if(run0lep)
	      {
		mvaBDT2 = SR3Reader.GenReMVAReader
		  (
		   
		   m4b,pt4b, met,ht,n_ad_j,
		   ptb1,ptb2 ,btag1,btag3,
                 "SR3Class"
		   );
	      }
	    mon.fillHisto("bdt","sr3",mvaBDT2,weight);
	  }
	 
      }
    
    if(isSR1)
      {
	
     	mon.fillHisto("mult","fjet_3",vec_fjet.size(),weight);
	mon.fillHisto("mult","bjet_cc_3",vec_bjets_cc.size(),weight);
	mon.fillHisto("mult","jet_cc_3",vec_jet_cc.size(),weight);
	mon.fillHisto("Ht","total",Ht(vec_jet_cc,vec_fjet),weight);
	ht_SR1=Ht(vec_jet_cc,vec_fjet);
	//met
	met_SR1=ev.met_pt;
	//mon.fillHisto("met","step_3",ev.met_pt,weight);
	//mon.fillHisto("met_phi","step_3",ev.met_phi,weight);
	if(hasMETtrigger1||hasMETtrigger2){
	  mon.fillHisto("met","step3_trig",ev.met_pt,weight);
	  mon.fillHisto("met_phi","step3_trig",ev.met_phi,weight);
	}
	//2lept inv mass
	if(!run0lep){
	  mon.fillHisto("inv_mass","dilept",(vec_leptons[0]+vec_leptons[1]).M(),weight);
	  mon.fillHisto("dR","dilept",getDeltaR(vec_leptons[0],vec_leptons[1]),weight);
	  drll_SR1=getDeltaR(vec_leptons[0],vec_leptons[1]);
	  dilep_pt_SR1=(vec_leptons[0]+vec_leptons[1]).Pt();
	  dphi_met_l_SR1=min(fabs(vec_leptons[0].Phi()-ev.met_phi),fabs(vec_leptons[1].Phi()-ev.met_phi));       
	  
	  mon.fillHisto("pt","lept1",vec_leptons[0].Pt(),weight);
	  mon.fillHisto("eta","lept1",vec_leptons[0].Eta(),weight);
	  mon.fillHisto("phi","lept1",vec_leptons[0].Phi(),weight);
	  //lepton2 kinematics
	  
	  mon.fillHisto("pt","lep2",vec_leptons[1].Pt(),weight);
	  mon.fillHisto("eta","lep2",vec_leptons[1].Eta(),weight);
	  mon.fillHisto("phi","lep2",vec_leptons[1].Phi(),weight);
	}
	
	// Fat-jet analysis:
	//fj1 kinematics
	ptf1=vec_fjet[0].Pt();
	mon.fillHisto("pt","fjet1",vec_fjet[0].Pt(),weight);
	mon.fillHisto("eta","fjet1",vec_fjet[0].Eta(),weight);
	mon.fillHisto("phi","fjet1",vec_fjet[0].Phi(),weight);
	
        // fj1 sdmass
	sd_mass1=ev.fjet_softdropM[fjet_index[0].second];
	mon.fillHisto("fj_sd_mass","1",ev.fjet_softdropM[fjet_index[0].second],weight);
	mon.fillHisto("fj_pt_sd_mass","1",vec_fjet[0].Pt(),ev.fjet_softdropM[fjet_index[0].second],weight);
	mon.fillHisto("subcount","1",ev.fjet_subjet_count[fjet_index[0].second],weight);
	
	bool matched(false);
	bool matched1(false);bool matched2(false);
	
	if(isSignal) 
	  {
	    if(fjet_matched(vec_fjet[0],vec_genbb1_truth)>0)
	      {
		mon.fillHisto("fj_pt_bb_pt","1",vec_fjet[0].Pt(),(vec_genbb1_truth[0]+vec_genbb1_truth[1]).Pt(),weight);
		matched1=true;
	      } 
	else if( fjet_matched(vec_fjet[0],vec_genbb2_truth)>0)
	  {
	    mon.fillHisto("fj_pt_bb_pt","1",vec_fjet[0].Pt(),(vec_genbb2_truth[0]+vec_genbb2_truth[1]).Pt(),weight);
	    matched2=true;  
	  }
	    if(matched1||matched2)matched=true;
	  }
	else { // backgrounds:
	  
	  if( fjet_matched(vec_fjet[0],vec_genqg)>0){
	    mon.fillHisto("fj_pt_bb_pt","1",vec_fjet[0].Pt(),fjet_matched(vec_fjet[0],vec_genqg),weight);
	    matched=true;
      }
	}
	if(matched){
	  mon.fillHisto("pt","bb_AK8",vec_fjet[0].Pt(),weight);
	  mon.fillHisto("fj_sd_mass","1_mat",ev.fjet_softdropM[fjet_index[0].second],weight);
	  mon.fillHisto("subcount","1_mat",ev.fjet_subjet_count[fjet_index[0].second],weight);
	}
											      
	//xbb discriminants
	xbb1=ev.fjet_btag10[fjet_index[0].second]/(ev.fjet_btag10[fjet_index[0].second]+ev.fjet_btag13[fjet_index[0].second]+ev.fjet_btag14[fjet_index[0].second]+ev.fjet_btag15[fjet_index[0].second]+ev.fjet_btag16[fjet_index[0].second]+ev.fjet_btag17[fjet_index[0].second]);
	if(TMath::IsNaN(xbb1)||std::isinf(xbb1))xbb1=-1;
	mon.fillHisto("fjet_btagXbb","1",xbb1,weight);
	xbbccqq1=(ev.fjet_btag10[fjet_index[0].second]+ev.fjet_btag11[fjet_index[0].second]+ev.fjet_btag12[fjet_index[0].second])/(ev.fjet_btag10[fjet_index[0].second]+ev.fjet_btag11[fjet_index[0].second]+ev.fjet_btag12[fjet_index[0].second]+ev.fjet_btag13[fjet_index[0].second]+ev.fjet_btag14[fjet_index[0].second]+ev.fjet_btag15[fjet_index[0].second]+ev.fjet_btag16[fjet_index[0].second]+ev.fjet_btag17[fjet_index[0].second]);
	if(TMath::IsNaN(xbbccqq1)||std::isinf(xbbccqq1))xbbccqq1=-1;
	mon.fillHisto("fjet_btagXbbXccXqq","1",xbbccqq1,weight);
	mon.fillHisto("fjet_XbbXccXqq_pt","1",vec_fjet[0].Pt(),xbbccqq1);
	mon.fillHisto("fjet_Xbb_pt","1",vec_fjet[0].Pt(),xbb1,weight);
	
	if(verbose)
          {
            if(ev.fjet_softdropM[fjet_index[0].second]<1)
              {
		std::cout << " this f-jet has pt " <<vec_fjet[0].Pt()  << ", no of sujets  " << ev.fjet_subjet_count[fjet_index[0].second]<<"sd_mass" <<ev.fjet_softdropM[fjet_index[0].second]<< " , Xbbccqq: " <<xbbccqq1<<",and matched;: "<< matched<<std::endl;
	      }
	  }
       
	btag3_SR1=ev.jet_btag1[bjet_index_cc[0].second];
	ptb1_SR1=vec_bjets_cc[0].Pt();
	n_ad_j_SR1=vec_jet_cc.size();
	drjj_SR1=getDeltaR(vec_fjet[0],vec_bjets_cc[0]);
	mon.fillHisto("phi","met_sr1",ev.met_phi,weight);
	mon.fillHisto("eventflow","histo",5,weight);
	hr->Fill(6);
	mon.fillHisto("btag1","3", ev.jet_btag1[bjet_index_cc[0].second],weight);
	mon.fillHisto("pt","bjet1",vec_bjets_cc[0].Pt(),weight);
	mon.fillHisto("eta","bjet1",vec_bjets_cc[0].Eta(),weight);
	mon.fillHisto("phi","bjet1",vec_bjets_cc[0].Phi(),weight);
	if(vec_bjets_cc.size()==1){
	  m4b_SR1=(vec_bjets_cc[0]+vec_fjet[0]).M();
	  pt4b_SR1=(vec_bjets_cc[0]+vec_fjet[0]).Pt();
	  
	  if(!run0lep){
	    dphi_met_j_SR1=min(fabs(deltaPhi(vec_fjet[0].Phi(),ev.met_phi)),fabs(deltaPhi(vec_bjets_cc[0].Phi(),ev.met_phi)));
	    if(fabs((vec_fjet[0]+vec_bjets_cc[0]).Phi()-(vec_leptons[0]+vec_leptons[1]).Phi())<TMath::Pi())
	      {
		dphiHZ_SR1=std::fabs((vec_fjet[0]+vec_bjets_cc[0]).Phi()-(vec_leptons[0]+vec_leptons[1]).Phi());
	      }
	    else{
	      dphiHZ_SR1=2*TMath::Pi()-std::fabs((vec_fjet[0]+vec_bjets_cc[0]).Phi()-(vec_leptons[0]+vec_leptons[1]).Phi());
	    }
	  }
	  mon.fillHisto("m_bbb","tot",(vec_bjets_cc[0]+vec_fjet[0]).M(),weight);
	  mon.fillHisto("pt_bbb","tot",(vec_bjets_cc[0]+vec_fjet[0]).Pt(),weight);
	  mon.fillHisto("m_bbb","sr1",(vec_bjets_cc[0]+vec_fjet[0]).M(),weight);
	  mon.fillHisto("pt_bbb","sr1",(vec_bjets_cc[0]+vec_fjet[0]).Pt(),weight);
	}
	else if(vec_bjets_cc.size()>1){
	  
	  m4b_SR1=(vec_bjets_cc[0]+vec_bjets_cc[1]+vec_fjet[0]).M();
	  pt4b_SR1=(vec_bjets_cc[0]+vec_bjets_cc[1]+vec_fjet[0]).Pt();
	  
	  if(!run0lep){
	    float dphi_met_j1=min(fabs(deltaPhi(vec_fjet[0].Phi(),ev.met_phi)),fabs(deltaPhi(vec_bjets_cc[0].Phi(),ev.met_phi)));
	    float dphi_met_j2=min(fabs(deltaPhi(vec_bjets_cc[0].Phi(),ev.met_phi)),fabs(deltaPhi(vec_bjets_cc[1].Phi(),ev.met_phi)));
	    dphi_met_j_SR1=min(dphi_met_j1,dphi_met_j2);
	    if(std::fabs((vec_fjet[0]+vec_bjets_cc[0]+vec_bjets_cc[1]).Phi()-(vec_leptons[0]+vec_leptons[1]).Phi())< TMath::Pi())
	      {
		dphiHZ_SR1=std::fabs((vec_fjet[0]+vec_bjets_cc[0]+vec_bjets_cc[1]).Phi()-(vec_leptons[0]+vec_leptons[1]).Phi());
	      }
	    else {
	      dphiHZ_SR1=2*TMath::Pi()-std::fabs((vec_fjet[0]+vec_bjets_cc[0]+vec_bjets_cc[1]).Phi()-(vec_leptons[0]+vec_leptons[1]).Phi());
	    }
	  }
	  mon.fillHisto("pt","bjet2",vec_bjets_cc[1].Pt(),weight);
	  mon.fillHisto("m_bbb","tot",(vec_bjets_cc[0]+vec_bjets_cc[1]+vec_fjet[0]).M(),weight);
	  mon.fillHisto("pt_bbb","tot",(vec_bjets_cc[0]+vec_bjets_cc[1]+vec_fjet[0]).Pt(),weight);
	  mon.fillHisto("m_bbb","sr1",(vec_bjets_cc[0]+vec_bjets_cc[1]+vec_fjet[0]).M(),weight);
	  mon.fillHisto("pt_bbb","sr1",(vec_bjets_cc[0]+vec_bjets_cc[1]+vec_fjet[0]).Pt(),weight);
	  
	}
	if(TMath::IsNaN(sd_mass1)|| isinf(sd_mass1))continue;
	if(run0lep)
	  {	  
	    mvaBDT1 = SR1Reader.GenReMVAReader
	      (
	       
	    m4b_SR1,pt4b_SR1, met_SR1,ht_SR1,
	    ptf1,sd_mass1,xbb1, xbbccqq1,n_ad_j_SR1,ptb1_SR1,btag3_SR1,drjj_SR1,"S\
R1Class");
	  }
	else if(! run0lep)
	  {
	    mvaBDT1 = SR1Reader.GenReMVAReader
	      (
	       dilep_pt_SR1, drll_SR1,dphiHZ_SR1,dphi_met_l_SR1,dphi_met_j_SR1,
	       m4b_SR1,pt4b_SR1, met_SR1,ht_SR1,
	       ptf1,sd_mass1,xbb1, xbbccqq1,n_ad_j_SR1,ptb1_SR1,btag3_SR1,drjj_SR1,"SR1Class");
	    mon.fillHisto("bdt","sr1",mvaBDT1,weight);
	  }
	
      }
     //##############################################################################
  
    
    if (runMVA) {
      
       float mvaweight = 1.0; // weight; //xsecWeight; //1.0;
       
       myMVAHandler_.getEntry (isSR1,isSR2,isSR3,m4b,pt4b,m4b_SR1,pt4b_SR1,ptf1,sd_mass1,ht,met,ht_SR1,met_SR1,xbb1,xbbccqq1,ptb1,ptb1_SR1,ptb2,n_ad_j,n_ad_j_SR1,btag1,btag3,drjj,dphi_met_j,dilep_pt,drll,dphiHZ,dphi_met_l,btag3_SR1,drjj_SR1,dphi_met_j_SR1,dilep_pt_SR1,drll_SR1,dphiHZ_SR1,dphi_met_l_SR1,weight);
       
       myMVAHandler_.fillTree();
       //	  }
       //	}
     } // runMVA
     
    }// end event loop
  printf("\n");
  file->Close();
  
  //write MVA files
 
  TString mvaout = TString ( runProcess.getParameter<std::string>("outdir") ) + "/mva_" + outFileUrl + ".root";
  if (runMVA)myMVAHandler_.writeTree(mvaout);
  std::cout << "reconstruction level: "  << std::endl;
  
  std::cout<< "Cuts," << "No. of events," << "Efficiency" <<std::endl;
  std::cout << "n_events(raw) " <<  totalEntries << std::endl;
  std::cout << " lept cut," << n_event_lepton_test << ","<< (float)n_event_lepton_test/(float)totalEntries*100<<std::endl;
  std::cout << "inv mass cut," <<INVM<< ","<< (float)INVM/(float)totalEntries*100<<std::endl;
  std::cout << "MET cut," <<METCUT<< ","<< (float)METCUT/(float)totalEntries*100<<std::endl;
  std::cout << "1fjet  cut," << sr1<< ","<< (float)sr1/(float)totalEntries*100<<std::endl;
  std::cout << ">2fjet  cut," << sr2<< ","<< (float)sr2/(float)totalEntries*100<<std::endl;
  std::cout << "number of events after at least 2 fjet cuts and z->vv: " << st4_vv << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and Z->qq_light: " << st4_qq << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and Z->bb: " << st4_bb << std::endl;
  std::cout << "number of events after at least 2 fjet cuts and Z->ll: " << st4_ll << std::endl;
  std::cout << "generator level: "  << std::endl;
  std::cout << "number of events with z->qq_light: " << zqq_light << std::endl;
  std::cout << "number of events with z->bb: " << zbb << std::endl;
  std::cout << "number of events with z->vv: " << zvv << std::endl;
  std::cout << "number of events with z->ll: " << zll << std::endl;
  // std::cout << "n_expected" <<nev_exp << std::endl;
  // std::cout << "n_expected_final" <<nev_exp*n_fjets/totalEntries << std::endl;
  // std::cout << "weight " << weight << std::endl;

  //##############################################
  //########     SAVING HISTO TO FILE     ########
  //##############################################
  //save control plots to file
  outUrl += "/";
  outUrl += outFileUrl + ".root";
  //    outUrl = outFileUrl + ".root";
  printf("Results saved in %s\n", outUrl.Data());

  
  // in the end: save all to the file
  int nTrial = 0;
  TFile *ofile=TFile::Open(outUrl, "recreate");
  while( !ofile->IsOpen() || ofile->IsZombie() ){
    if(nTrial > 3){
      printf("Output file open failed!");
      if ( outTxtFile_final ) fclose(outTxtFile_final);
      return -1;
    }
    nTrial++;
    usleep(1000000*nTrial);
    ofile=TFile::Open(outUrl, "update");
  }
  mon.Write();
  hr->Write();
  drbf->Write();
  drjf->Write();
  ofile->Close();
  
  if ( outTxtFile_final ) fclose(outTxtFile_final);
  
} // end main



//functions
//deltar
double getDeltaR(TLorentzVector vec_1, TLorentzVector vec_2)
{
  double delta_phi;
  double delta_eta;

  delta_phi = vec_1.Phi() - vec_2.Phi();
  delta_eta = vec_1.Eta() - vec_2.Eta();

  return std::sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
}

bool Jet_matched(TLorentzVector jet,std::vector<TLorentzVector> vecb) 
{
  bool matched=false;
  float dR_min=999;
  float dR=0;
   for (int i=0; i<vecb.size(); i++) {
      dR=getDeltaR(jet, vecb[i]);
      if (dR<dR_min) dR_min=dR;
   }
   if (dR_min<0.4) matched=true;
   return matched; 
}

float fjet_matched(TLorentzVector fjet,std::vector<TLorentzVector> vecb)
{
  float qg_pt=-999.;

  std::vector<TLorentzVector> vecmin;
  TLorentzVector vecqgpt;
  
  bool ismatched(false);

  for (int i=0; i<vecb.size(); i++)
   {
     ismatched=false;
      float dR=getDeltaR(fjet, vecb[i]);
      if (dR<0.8) {
	ismatched=true;
	vecmin.push_back(vecb[i]);
      }
   }
  for (int i=0; i<vecmin.size(); i++)
    {
      vecqgpt+=vecmin[i];
    }
  
  if(ismatched) qg_pt=vecqgpt.Pt();
  return qg_pt;
  
}

double dR_j_fj_min(TLorentzVector jet,std::vector<TLorentzVector> vecfjet)
{
      double dR[vecfjet.size()];
      double dR_min=0;
      int  min_ind=0;
      
      for (int fj = 0; fj < vecfjet.size(); fj++)
      {
         dR[fj] = getDeltaR(jet, vecfjet[fj]);
         if(dR[fj] <dR[min_ind]) min_ind = fj; 
      }
      dR_min=dR[min_ind];
      return dR_min;
}
double dphi_met_j_min(float metphi,std::vector<TLorentzVector> jet,std::vector<TLorentzVector> fjet)
{
  
  double dphif_min=999;
  double dphij_min=999;
  double dphi_min=999;
  for (int fj = 0; fj < fjet.size(); fj++)
    {
      double dphif = fabs(fjet[fj].Phi()- metphi);
      if(dphif <dphif_min) dphif_min=dphif;
    }
  for (int j = 0; j < jet.size(); j++)
    {
      double dphij = fabs(jet[j].Phi()- metphi);
      if(dphij <dphij_min) dphij_min=dphij;
    }
  dphi_min=min(dphif_min,dphij_min);
  return dphi_min;
}

double Ht(std::vector<TLorentzVector> vecjet,std::vector<TLorentzVector> vecfjet)
{
   double Ht=0;
   if (vecjet.size()>0){
   for(int i=0;i< vecjet.size();i++) Ht+=vecjet[i].Pt();}
   if(vecfjet.size()>0){
   for(int j=0;j< vecfjet.size();j++) Ht+=vecfjet[j].Pt();}
   return Ht;
}

    
std::vector<TLorentzVector> sort_vec_pt(std::vector<TLorentzVector> vec)
{
if (vec.size()>0)
    {
      for (int i=0;i<vec.size()-1;i++)
      {
          int maxIndex = i;  
          for (int j = i + 1; j < vec.size(); j++) 
            {  
              if (vec[j].Pt() > vec[maxIndex].Pt())
          {  
            maxIndex = j;  
          }  
            }   
          TLorentzVector temp;
          temp = vec[i];  
          vec[i] = vec[maxIndex];  
          vec[maxIndex] = temp; 
      }
    }
  return vec;
}
