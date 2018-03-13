#include <iostream>
#include <map>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"
#include "UserCode/bsmhiggs_fwk/interface/TMVAReader.h"
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/bsmhiggs_fwk/interface/PDFInfo.h"
#include "UserCode/bsmhiggs_fwk/interface/rochcor2016.h" 
#include "UserCode/bsmhiggs_fwk/interface/muresolution_run2.h" 
#include "UserCode/bsmhiggs_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/bsmhiggs_fwk/interface/BTagCalibrationStandalone.h"
#include "UserCode/bsmhiggs_fwk/interface/BtagUncertaintyComputer.h"
//#include "UserCode/bsmhiggs_fwk/interface/METUtils.h"
//#include "UserCode/bsmhiggs_fwk/interface/BTagUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/statWgt.h"

#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

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


using namespace std;


namespace LHAPDF {
void initPDFSet(int nset, const std::string& filename, int member=0);
int numberPDF(int nset);
void usePDFMember(int nset, int member);
double xfx(int nset, double x, double Q, int fl);
double getXmin(int nset, int member);
double getXmax(int nset, int member);
double getQ2min(int nset, int member);
double getQ2max(int nset, int member);
void extrapolate(bool extrapolate=true);
}

struct stPDFval {
    stPDFval() {}
    stPDFval(const stPDFval& arg) :
        qscale(arg.qscale),
        x1(arg.x1),
        x2(arg.x2),
        id1(arg.id1),
        id2(arg.id2) {
    }

    double qscale;
    double x1;
    double x2;
    int id1;
    int id2;
};

struct ptsort: public std::binary_function<LorentzVector, LorentzVector, bool> 
{
  bool operator () (const LorentzVector & x, const LorentzVector & y) 
  { return  ( x.pt() > y.pt() ) ; }
};

struct ptsortinpair: public std::binary_function<std::pair<int,LorentzVector>, std::pair<int,LorentzVector>, bool>  
{
  bool operator () (const std::pair<int,LorentzVector> & x, const std::pair<int,LorentzVector> & y) 
  { return (x.second.pt() > y.second.pt() ); }
};

struct btagsort: public std::binary_function<PhysicsObject_Jet, PhysicsObject_Jet, float> 
{
  bool operator () (const PhysicsObject_Jet & x, PhysicsObject_Jet & y) 
  { return  ( x.btag0 > y.btag0 ) ; }
}; 

//bool runDBversion = false;

// Physics objects offline thresholds
//const float lep_threshold_=25.; 
const float mu_threshold_=25.; 
const float ele_threshold_=30.; 
const float jet_threshold_=20.; 

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
const float CSVLooseWP = 0.5426;  // Updated to 80X Moriond17 Loose
const float CSVMediumWP = 0.800;
const float CSVTightWP = 0.935;

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/Hbbtagging#8_0_X
// https://indico.cern.ch/event/543002/contributions/2205058/attachments/1295600/1932387/cv-doubleb-tagging-btv-approval.pdf (definition of WPs: slide 16)
const float DBLooseWP = 0.300;
const float DBMediumWP = 0.600;
const float DBTightWP = 0.900;

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors
const float DeepCSVLooseWP = 0.2219;
const float DeepCSVMediumWP = 0.6324;
const float DeepCSVTightWP = 0.8958;

int main(int argc, char* argv[])
{
    //##################################################################################
    //##########################    GLOBAL INITIALIZATION     ##########################
    //##################################################################################

    //check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }
   
    //load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    //AutoLibraryLoader::enable();
    FWLiteEnabler::enable();

    //configure the process
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  
    bool isMC = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");

    double xsec = runProcess.getParameter<double>("xsec");

    TString proc=runProcess.getParameter<std::string>("proc");
    TString dtag=runProcess.getParameter<std::string>("tag");
    TString suffix=runProcess.getParameter<std::string>("suffix");

    bool verbose = runProcess.getParameter<bool>("verbose");

    // will reweight the top pt in TT+jets sample (optional)
    bool reweightTopPt = runProcess.getParameter<bool>("reweightTopPt");

    // will produce the input root trees to BDT training (optional)
    bool runMVA = runProcess.getParameter<bool>("runMVA");

    // Will set DeepCSV as the default b-tagger automaticaly:
    bool use_DeepCSV = runProcess.getParameter<bool>("useDeepCSV");
    
    bool usemetNoHF = runProcess.getParameter<bool>("usemetNoHF");
    
    TString url = runProcess.getParameter<std::string>("input");
    TString outFileUrl( dtag ); //gSystem->BaseName(url));
    //outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir = runProcess.getParameter<std::string>("outdir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);


    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_final.Data());
    fprintf(outTxtFile_final,"run lumi event\n");

    int fType(0);
    if(url.Contains("DoubleEG")) fType=EE;
    if(url.Contains("DoubleMuon"))  fType=MUMU;
    if(url.Contains("MuonEG"))      fType=EMU;
    if(url.Contains("SingleMuon"))  fType=MU;
    if(url.Contains("SingleElectron")) fType=E;
    bool isSingleMuPD(!isMC && url.Contains("SingleMuon"));
    bool isDoubleMuPD(!isMC && url.Contains("DoubleMuon"));
    bool isSingleElePD(!isMC && url.Contains("SingleElectron"));
    bool isDoubleElePD(!isMC && url.Contains("DoubleEG"));

    bool isMC_ZZ  = isMC && ( string(url.Data()).find("TeV_ZZ_")  != string::npos );
    
    bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZamcatnloFXFX")  != string::npos
                              || string(url.Data()).find("MC13TeV_WZpowheg")  != string::npos );

    bool isMC_VVV = isMC && ( string(url.Data()).find("MC13TeV_WZZ")  != string::npos
                              || string(url.Data()).find("MC13TeV_WWZ")  != string::npos
                              || string(url.Data()).find("MC13TeV_ZZZ")  != string::npos );

    bool isMCBkg_runPDFQCDscale = (isMC_ZZ || isMC_WZ || isMC_VVV);

    bool isMC_ttbar = isMC && (string(url.Data()).find("TeV_TTJets")  != string::npos);
    bool isMC_stop  = isMC && (string(url.Data()).find("TeV_SingleT")  != string::npos);

    bool isMC_WJets = isMC && ( (string(url.Data()).find("MC13TeV_WJets")  != string::npos) || (string(url.Data()).find("MC13TeV_W1Jets")  != string::npos) || (string(url.Data()).find("MC13TeV_W2Jets")  != string::npos) || (string(url.Data()).find("MC13TeV_W3Jets")  != string::npos) || (string(url.Data()).find("MC13TeV_W4Jets")  != string::npos) );
    bool isMC_DY = isMC && ( (string(url.Data()).find("MC13TeV_DY")  != string::npos) );
    
    bool isMC_Wh = isMC && (string(url.Data()).find("Wh")  != string::npos); 
    bool isMC_Zh = isMC && (string(url.Data()).find("Zh")  != string::npos); 
    bool isMC_VBF = isMC && (string(url.Data()).find("VBF")  != string::npos); 

    bool isQCD = isMC && (string(url.Data()).find("QCD")  != string::npos);

    bool isSignal = (isMC_Wh || isMC_Zh || isMC_VBF );

    if (isSignal) printf("Signal url = %s\n",url.Data());

    //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
    //the scale factors are taken as average numbers from the pT dependent curves see:
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
    BTagSFUtil btsfutil;
    float beff(0.68), sfb(0.99), sfbunc(0.015);
    float leff(0.13), sfl(1.05), sflunc(0.12);

    //setup calibration readers 80X
    std::string b_tagging_name, csv_file_path;
    float LooseWP = -1, MediumWP = -1, TightWP = -1;
    if (!use_DeepCSV) 
    {
      b_tagging_name = "CSVv2";
       
      csv_file_path = std::string(std::getenv("CMSSW_BASE"))+
                      "/src/UserCode/bsmhiggs_fwk/data/weights/CSVv2_Moriond17_B_H.csv";
       
      LooseWP = CSVLooseWP;
      MediumWP = CSVMediumWP;
      TightWP = CSVTightWP;
    }
    if ( use_DeepCSV) {
      b_tagging_name = "DeepCSV";
       
      csv_file_path = std::string(std::getenv("CMSSW_BASE"))+
                      "/src/UserCode/bsmhiggs_fwk/data/weights/DeepCSV_Moriond17_B_H.csv";
       
      LooseWP = DeepCSVLooseWP;
      MediumWP = DeepCSVMediumWP;
      TightWP = DeepCSVTightWP;
    }
    
    BTagCalibration btagCalib(b_tagging_name, csv_file_path);
    // setup calibration readers 80X
    BTagCalibrationReader80X btagCal80X   (BTagEntry::OP_LOOSE, "central", {"up", "down"});
    btagCal80X.load(btagCalib, BTagEntry::FLAV_B, "comb");
    btagCal80X.load(btagCalib, BTagEntry::FLAV_C, "comb");
    btagCal80X.load(btagCalib, BTagEntry::FLAV_UDSG, "incl");


    //systematics
    bool runSystematics = runProcess.getParameter<bool>("runSystematics");
    std::vector<TString> varNames(1,"");
    if(runSystematics) {
        cout << "Systematics will be computed for this analysis: " << endl;
        varNames.push_back("_jerup"); //1
        varNames.push_back("_jerdown"); //2
        varNames.push_back("_jesup"); //3
        varNames.push_back("_jesdown"); //4
        varNames.push_back("_umetup"); //5
        varNames.push_back("_umetdown");//6
        varNames.push_back("_lesup"); //7
        varNames.push_back("_lesdown"); //8
        varNames.push_back("_puup"); //9
        varNames.push_back("_pudown"); //10
        varNames.push_back("_btagup"); //11
        varNames.push_back("_btagdown");//12
        if(isMCBkg_runPDFQCDscale) {
            varNames.push_back("_pdfup");
            varNames.push_back("_pdfdown");
            varNames.push_back("_qcdscaleup");
            varNames.push_back("_qcdscaledown");

        }
        if(isSignal) {
            varNames.push_back("_pdfup");
            varNames.push_back("_pdfdown");
            varNames.push_back("_qcdscaleup");
            varNames.push_back("_qcdscaledown");
        }
        if(isMC_ZZ) {
            varNames.push_back("_qqZZewkup");
            varNames.push_back("_qqZZewkdown");
        }

        for(size_t sys=1; sys<varNames.size(); sys++) {
            cout << varNames[sys] << endl;
        }
    }
    size_t nvarsToInclude=varNames.size();


    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    
    //jet energy scale uncertainties
    TString jecDir = runProcess.getParameter<std::string>("jecDir");
    //gSystem->ExpandPathName(uncFile);
    cout << "Loading jet energy scale uncertainties from: " << jecDir << endl;

    if     (dtag.Contains("2016B") || dtag.Contains("2016C") ||dtag.Contains("2016D")) jecDir+="Summer16_80X/Summer16_23Sep2016BCDV4_DATA/";
    else if(dtag.Contains("2016E") || dtag.Contains("2016F")) jecDir+="Summer16_80X/Summer16_23Sep2016EFV4_DATA/";
    else if(dtag.Contains("2016G")) jecDir+="Summer16_80X/Summer16_23Sep2016GV4_DATA/";
    else if(dtag.Contains("2016H")) jecDir+="Summer16_80X/Summer16_23Sep2016HV4_DATA/";
    if(isMC) {jecDir+="Summer16_80X/Summer16_23Sep2016V4_MC/";}

    gSystem->ExpandPathName(jecDir);
    FactorizedJetCorrector *jesCor = NULL;
    jesCor = utils::cmssw::getJetCorrector(jecDir,isMC);

    TString pf(isMC ? "MC" : "DATA");
    JetCorrectionUncertainty *totalJESUnc = NULL;
    totalJESUnc = new JetCorrectionUncertainty((jecDir+"/"+pf+"_Uncertainty_AK4PFchs.txt").Data());
    //JetCorrectionUncertainty jecUnc(uncFile.Data());

    //Lepton scale corrections
    EnergyScaleCorrection_class eScaler_("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_74x_pho");     
    eScaler_.doScale=true;
    eScaler_.doSmearings=true;

    std::string bit_string_stat = "001";
    std::string bit_string_syst = "010";
    std::string bit_string_gain = "100";
    std::bitset<6> bit_stat(bit_string_stat);
    std::bitset<6> bit_syst(bit_string_syst);
    std::bitset<6> bit_gain(bit_string_gain);

    //muon energy scale and uncertainties
    TString muscleDir = runProcess.getParameter<std::string>("muscleDir");
    gSystem->ExpandPathName(muscleDir);

    rochcor2016* muCor2016 = new rochcor2016(); //replace the MuScleFitCorrector we used at run1
    
    //pdf info
    
    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;

    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 9,0,9) );
    h->GetXaxis()->SetBinLabel(1,"Raw");
    h->GetXaxis()->SetBinLabel(2,"Trigger");
    h->GetXaxis()->SetBinLabel(3,"1 lepton");
    h->GetXaxis()->SetBinLabel(4,"E_{T}^{miss}>25");
    h->GetXaxis()->SetBinLabel(5,"50<M_{T}^{W}<250");
    h->GetXaxis()->SetBinLabel(6,">=(2-jets,2b-jets)");
    h->GetXaxis()->SetBinLabel(7,">=3b-tags");
    h->GetXaxis()->SetBinLabel(8,">=4b-tags");
 
    //generator level plots
    mon.addHistogram( new TH1F( "pileup", ";pileup;Events", 100,-0.5,99.5) );
    
    mon.addHistogram( new TH1F( "higgsMass",";m_{h} [GeV];Events",50,0.,1000.) );
    mon.addHistogram( new TH1F( "higgsPt",";p_{T}^{h} [GeV];Events",30,0.,500.));
    //mon.addHistogram( new TH1F( "higgsEta",";#eta (h);Evenets",100,-5,5) );
 
    //Top pt
    mon.addHistogram( new TH1F( "toppt",";#it{p}_{T}^{top} [GeV];Events",100,0.,2000.) );    

    //RECO level, physics objects
    mon.addHistogram( new TH1F( "dR_raw",";#Delta R(SV,b);Events",50,0.,5.));
    mon.addHistogram( new TH1F( "dRlj_raw",";#Delta R(lep,jet);Events",100,0.,5.));

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l} [GeV];Events", 30,0.,500.) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );
    
    mon.addHistogram( new TH1F( "jet_pt_raw", ";#it{p}_{T} [GeV];Events",30,0.,500.) );
    mon.addHistogram( new TH1F( "softjet_pt_raw", ";#it{p}_{T} [GeV];Events",20,0.,40.) );
    mon.addHistogram( new TH1F( "jet_eta_raw",";jet #eta;Events", 70,-3,3) );
    mon.addHistogram( new TH1F( "jet_phi_raw",";jet #phi;Events", 70,-6,6) );

    mon.addHistogram( new TH1F( "b_discrim"," ;b discriminator;",50,0,1.) );
    mon.addHistogram( new TH1F( "db_discrim"," ;double-b discriminator;",25,-1.,1.) );
    mon.addHistogram( new TH1F( "sd_mass"," ;soft-drop Mass;",60,0.,300.) );
    mon.addHistogram( new TH1F( "pruned_mass"," ;pruned Mass;",60,0.,300.) );
    mon.addHistogram( new TH1F( "softb_ntrk"," ; SV Ntrks;",21,-0.5,21.5) );
    mon.addHistogram( new TH1F( "softb_dxy"," ; SV dxy;",50,0.,20.) );
    mon.addHistogram( new TH1F( "softb_dxyz_signif"," ; SVSIP3D;",50,1.,100.) );
    mon.addHistogram( new TH1F( "softb_cos"," ; SV cos((PV,SV),p_{SV});",25,-1.,1.) );


    //for MC normalization (to 1/pb)
    TH1F* Hcutflow = (TH1F*) mon.addHistogram( new TH1F ("cutflow" , "cutflow" ,6,0,6) ) ;
    
    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nleptons", ";Lepton multiplicity;Events", 4,0,4) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h1->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h2 = (TH1F *)mon.addHistogram( new TH1F("njets_raw",  ";Jet multiplicity (#it{p}_{T}>20 GeV);Events",6,0,6) );
    for(int ibin=1; ibin<=h2->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h2->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h2->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h3 = (TH1F *)mon.addHistogram( new TH1F("nbjets_raw",    ";b-tag Jet multiplicity (#it{p}_{T}>20 GeV);Events",6,0,6) );
    for(int ibin=1; ibin<=h3->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h3->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h3->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH2F *h4 = (TH2F *)mon.addHistogram( new TH2F("nbjets_2D", ";AK4 jet multiplicity (#it{p}_{T}>20 GeV); b-Jet multiplicity",5,0,5,5,0,5) );
    for(int ibin=1; ibin<=h4->GetXaxis()->GetNbins(); ibin++) {
      TString label("");
      if(ibin==h4->GetXaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      h4->GetXaxis()->SetBinLabel(ibin,label);
    }
    for(int ibin=1; ibin<=h4->GetYaxis()->GetNbins(); ibin++) {
      TString label("");
      if(ibin==h4->GetYaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      h4->GetYaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h5 = (TH1F *)mon.addHistogram( new TH1F("nsubjets_raw",  ";Subjet multiplicity (#it{p}_{T}>20 GeV);Events",3,0,3) );
    for(int ibin=1; ibin<=h5->GetXaxis()->GetNbins(); ibin++) {
      TString label("");
      if(ibin==h5->GetXaxis()->GetNbins()) label +="#geq";
      else                                label +="=";
      label += (ibin-1);
      h5->GetXaxis()->SetBinLabel(ibin,label);
    }

    //--------------------------------------------------------------------------
    //some strings for tagging histograms:
    const char* astr[] = {"_b1","_b2","_b3","_b4"};
    std::vector<TString> htag(astr, astr+4);
    //--------------------------------------------------------------------------
    
    //btaging efficiency

    //#################################################
    //############# CONTROL PLOTS #####################
    //#################################################
    
    mon.addHistogram( new TH1F( "nvtx_raw",";Vertices;Events",100,0,100) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",";Vertices;Events",100,0,100) );
    mon.addHistogram( new TH1F( "pfmet",    ";E_{T}^{miss} [GeV];Events", 50,0.,500.) );
    mon.addHistogram( new TH1F( "ht",    ";H_{T} (p_{T}^{j}>20) [GeV];Events", 50,0.,800.) );
    mon.addHistogram( new TH1F( "ht_b30",    ";H_{T} (p_{T}^{j}>30) [GeV];Events", 50,0.,800.) );
    mon.addHistogram( new TH1F( "mtw",       ";#it{m}_{T}^{W} [GeV];Events", 50,0.,500.) );
    mon.addHistogram( new TH1F( "ptw",       ";#it{p}_{T}^{W} [GeV];Events",30,0.,500.) );
    mon.addHistogram( new TH1F( "dphiWh", ";#Delta#it{#phi}(#it{W},h);Events", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dRave",";#Delta R(b,b)_{ave};Events",50,0.,5.));
    mon.addHistogram( new TH1F( "dmmin",";#Delta m_{b,b}^{min};Events",25,0.,250.));
    mon.addHistogram( new TH1F( "dphijmet", ";|#Delta#it{#phi}(jet,E_{T}^{miss})|;#jet", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphilepmet", ";|#Delta#it{#phi}(lep,E_{T}^{miss})|;Events", 20,0,TMath::Pi()) );

    //MVA
    mon.addHistogram( new TH1F( "bdt", ";BDT;Events", 30, -0.3, 0.3) );
    
    //##################################################################################
    //########################## STUFF FOR CUTS OPTIMIZATION  ##########################
    //##################################################################################

    TH1F* Hoptim_systs     =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;;", nvarsToInclude,0,nvarsToInclude) ) ;
    for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
      Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);
    }

    std::vector<double> optim_Cuts1_bdt;
    optim_Cuts1_bdt.push_back(-0.4); //add a bin in the shapes with a BDT cut of -0.4
    for(double bdt=-0.30;bdt<0.30;bdt+=0.02) { optim_Cuts1_bdt.push_back(bdt); }

    TH2F* Hoptim_cuts =(TH2F*)mon.addHistogram(new TProfile2D("optim_cut", ";cut index;variable", optim_Cuts1_bdt.size(),0,optim_Cuts1_bdt.size(), 1, 0, 1)) ;
    Hoptim_cuts->GetYaxis()->SetBinLabel(1, "BDT>");
    for(unsigned int index=0;index<optim_Cuts1_bdt.size();index++){ Hoptim_cuts->Fill(index, 0.0, optim_Cuts1_bdt[index]); }
    for(size_t ivar=0; ivar<nvarsToInclude; ivar++){
      mon.addHistogram( new TH2F (TString("bdt_shapes")+varNames[ivar],";cut index;BDT;Events",optim_Cuts1_bdt.size(),0,optim_Cuts1_bdt.size(), 30,-0.3,0.3) );
    }



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

    //MC normalization (to 1/pb)
    double xsecWeight = 1.0;
    float cnorm=1.0;
    if (isMC) {
      //xsecWeight = 0.; //disable MC sample if not present in the map
      
      double totalNumberofEvents(0.);
      
      TH1F* nevtH = (TH1F *) file->Get("mainNtuplizer/nevents");
      totalNumberofEvents = nevtH->GetBinContent(1);
      TH1F* posH = (TH1F *) file->Get("mainNtuplizer/n_posevents");
      TH1F* negH = (TH1F *) file->Get("mainNtuplizer/n_negevents");
      if(posH && negH) cnorm = posH->GetBinContent(1) - negH->GetBinContent(1);
      if(rescaleFactor>0) cnorm /= rescaleFactor;
      printf("cnorm = %f and totalNumberOfEvents= %f\n", cnorm, totalNumberofEvents);
      //xsecWeight=xsec/totalNumberofEvents;
      xsecWeight=xsec/cnorm; // effective luminosity}
      /*
      std::map<std::string, int> xsec_map = mStat;

      //std::string myproc = proc.Data();
      //std::cout << "Runnin process " << myproc << std::endl;
      std::map<std::string, int>::iterator it;
      for ( it = xsec_map.begin(); it != xsec_map.end(); it++ ) {
        if (it->first == proc.Data()) {
          xsecWeight = (xsec/(float)it->second);
          //totalNumberofEvents = it->second;
          std::cout << "weight = " << (xsecWeight*35866.9) << std::endl;
        }
      }
      */
      //float pereventwgt=(xsecWeight*35866.9);
      //printf("\n Running process with xSec = %f , and totalNumEvents = %d  . Per event weight is (L=35.9 fb-1): %f \n\n",
      //xsec, totalNumberofEvents, pereventwgt );
    }
    //Hcutflow->SetBinContent(1,cnorm);

    //pileup weighting
    TString PU_Central = runProcess.getParameter<std::string>("pu_central");
    gSystem->ExpandPathName(PU_Central);
    cout << "Loading PU weights Central: " << PU_Central << endl;
    TFile *PU_Central_File = TFile::Open(PU_Central);
    
    TH1F* PU_intended = (TH1F *) PU_Central_File->Get("pileup"); // Data pileup distribution
    TH1F* PU_generated=NULL; // MC pileup distribution 

    TH1F* PU_weight=new TH1F("hPUweight","",100,-0.5,99.5);
    
    if (isMC) {
      //PU_generated = (TH1F*)file->Get("mainNtuplizer/pileuptrue"); // MC pileup distribution
      PU_generated = (TH1F*)file->Get("mainNtuplizer/pileuptrue");
      
      PU_intended->Scale(1./PU_intended->Integral());
      TH1F* PUnorm = PU_intended;
      PU_generated->Scale(1./PU_generated->Integral());
      TH1F *PUgen = PU_generated;
      
      TH1F* Quotient = PUnorm; //->Clone("quotient");
      Quotient->Divide(PUgen); 

      for(int ibin=0; ibin<100; ibin++){
        float x = Quotient->GetBinContent(ibin);
        PU_weight->SetBinContent(ibin,x);
        if ( verbose ) printf("pu = %3d has weight = %7.3f \n",ibin,x);
      }
    }//is MC
    
    //event categorizer
    //EventCategory eventCategoryInst(1);   

    //Lepton scale factors
    LeptonEfficiencySF lepEff;

    //####################################################################################################################
    //###########################################           MVAHandler         ###########################################
    //####################################################################################################################
    //construct MVA out put file name
    TString mvaout = TString ( runProcess.getParameter<std::string>("outdir") ) + "/mva_" + outFileUrl + ".root";
    MVAHandler myMVAHandler_;
    if (runMVA) { myMVAHandler_.initTree(mvaout); }

    //####################################################################################################################
    //###########################################           TMVAReader         ###########################################
    //####################################################################################################################

    TMVAReader myTribTMVAReader;
    myTribTMVAReader.InitTMVAReader();
    std::string TribMVA_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/Haa4bSBClassificationTribMVA_BDT.weights.xml";
    myTribTMVAReader.SetupMVAReader( "Haa4bSBClassificationTribMVA", TribMVA_xml_path );

    TMVAReader myQuabTMVAReader;
    myQuabTMVAReader.InitTMVAReader();
    std::string QuabMVA_xml_path = std::string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/mva/Haa4bSBClassificationQuabMVA_BDT.weights.xml"; 
    myQuabTMVAReader.SetupMVAReader( "Haa4bSBClassificationQuabMVA", QuabMVA_xml_path );
    
    
    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################

    //loop on all the events
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    DuplicatesChecker duplicatesChecker;
    int nDuplicates(0);
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");

    for( int iev=evStart; iev<evEnd; iev++) 
    {
        if((iev-evStart)%treeStep==0) {
          printf("."); fflush(stdout);
        }

        if ( verbose ) printf("\n\n Event info %3d: \n",iev);


        //##############################################   EVENT LOOP STARTS   ###########################################
        //load the event content from tree
        summaryHandler_.getEntry(iev);
        DataEvtSummary_t &ev=summaryHandler_.getEvent();
        if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) {
            nDuplicates++;
            cout << "nDuplicates: " << nDuplicates << endl;
            continue;
        }

        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev); 

        std::vector<TString> tags(1,"all");
        //genWeight
        float genWeight = 1.0;
        if (isMC) {
          if(ev.genWeight<0) { genWeight = -1.0; }
        }
        //systematical weight
        float weight = 1.0; //xsecWeight;
        if(isMC) 
        {
          weight *= genWeight;
          //Here is the tricky part.,... rewrite xsecWeight for WJets/WXJets and DYJets/DYXJets
          if( isMC_WJets )
          {
            xsecWeight = xsecWeightCalculator::xsecWeightCalcLHEJets(0, ev.lheNJets);
          }
          else if( isMC_DY )
          {
            if (string(url.Data()).find("10to50")  != string::npos)
            {
              xsecWeight = xsecWeightCalculator::xsecWeightCalcLHEJets(1, ev.lheNJets);
            }
            else
            {
              xsecWeight = xsecWeightCalculator::xsecWeightCalcLHEJets(2, ev.lheNJets);
            }
          }
          weight *= xsecWeight; 
        }

	// Apply Top pt-reweighting
	double top_wgt(1.0);    

	if(reweightTopPt && isMC_ttbar){
	  PhysicsObjectCollection &partons = phys.genpartons;
	  double SFtop(0.);
	  double SFantitop(0.);

	  int itop(0);
	  for (auto & top : partons) {

	    printf("Parton : ID=%6d, m=%5.1f, momID=%6d : pt=%6.1f, status=%d\n",
		   top.id,
		   top.mass(),
		   top.momid,
		   top.pt(),
		   top.status
		   );

	    if (top.id==6 && top.status==62) {
	      SFtop=exp(0.0615-0.0005*top.pt());
	      itop++;

	      mon.fillHisto("toppt","top",top.pt(),weight);
	    }
	    if (top.id==-6 && top.status==62) {
	      SFantitop=exp(0.0615-0.0005*top.pt());
	      itop++;

	      mon.fillHisto("toppt","antitop",top.pt(),weight);
	    }
	  }
  
	  if (itop<2) { 
	    printf("Did not found tt pair!!\n"); 
	  } else if (itop==2) {
	    top_wgt=sqrt(SFtop*SFantitop);
	  } else {
	    printf("More than 2 top particles found. Please check\n");
	  }
    
	  //printf("weight= %3f and top weight= %3f\n",weight,top_wgt);
	  weight *= top_wgt;
	  //printf("Final weight is : %3f\n\n",weight);
	}

        //only take up and down from pileup effect
        double TotalWeight_plus = 1.0;
        double TotalWeight_minus = 1.0;

        if(isMC) mon.fillHisto("pileup", tags, ev.ngenTruepu, 1.0);
        
        float puWeight(1.0);
        if(isMC) {
          puWeight = getSFfrom1DHist(ev.ngenTruepu, PU_weight) ;
          if ( verbose ) printf("pu = %3d has weight = %7.3f \n",ev.ngenTruepu,puWeight);
          weight *= puWeight;
          //TotalWeight_plus  *= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Up);
          //TotalWeight_minus *= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Down);
        }
        
        Hcutflow->Fill(1,genWeight);
        Hcutflow->Fill(2,xsecWeight);
        Hcutflow->Fill(3,puWeight);
        Hcutflow->Fill(4,weight);
        //Hcutflow->Fill(3,weight*TotalWeight_minus);
        //Hcutflow->Fill(4,weight*TotalWeight_plus);

        //add PhysicsEvent_t class, get all tree to physics objects
        //PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        // FIXME need to have a function: loop all leptons, find a Z candidate,
        // can have input, ev.mn, ev.en
        // assign ee,mm,emu channel
        // check if channel name is consistent with trigger
        // store dilepton candidate in lep1 lep2, and other leptons in 3rdleps


        bool hasMMtrigger = ev.triggerType & 0x1;
        bool hasMtrigger  = (ev.triggerType >> 1 ) & 0x1;
        bool hasEEtrigger = (ev.triggerType >> 2 ) & 0x1;
        // type 3 is high-pT eeTrigger (safety)
        bool hasEtrigger  = (ev.triggerType >> 4 ) & 0x1;
        bool hasEMtrigger = (ev.triggerType >> 5 ) & 0x1;


        //#########################################################################
        //#####################      Objects Selection       ######################
        //#########################################################################

        //
        // MET ANALYSIS
        //
        //apply Jet Energy Resolution corrections to jets (and compute associated variations on the MET variable)
        //std::vector<PhysicsObjectJetCollection> variedJets;
        //LorentzVectorCollection variedMET;

        //METUtils::computeVariation(phys.jets, phys.leptons, (usemetNoHF ? phys.metNoHF : phys.met), variedJets, variedMET, &jecUnc);

        LorentzVector metP4=phys.met; //variedMET[0];
        PhysicsObjectJetCollection &corrJets = phys.jets; //variedJets[0];
        PhysicsObjectFatJetCollection &fatJets = phys.fatjets;
        PhysicsObjectSVCollection &secVs = phys.svs;

        //
        // LEPTON ANALYSIS
        //
        PhysicsObjectLeptonCollection &leps = phys.leptons;

	int nGoodLeptons(0);
	std::vector<std::pair<int,LorentzVector> > goodLeptons;
	int nExtraLeptons(0);
	std::vector<LorentzVector> extraLeptons;

	std::vector<LorentzVector> vetoLeptons; // ---> Use this collection to remove jet-lepton overlaps for e/mu below 30(25) GeV. 

	float lep_threshold(25.);
	float eta_threshold=2.5;

	for (auto &ilep : leps) {
	  //if ( ilep.pt()<5. ) continue;

	  int lepid = ilep.id;
	  if (abs(lepid)==13) eta_threshold=2.4;
	  if (fabs(ilep.eta())>eta_threshold) continue;
  
	  bool hasTightIdandIso(true);
	  if (abs(lepid)==11) {
	    lep_threshold=ele_threshold_;
	    hasTightIdandIso &= (ilep.passIdEl && ilep.passIsoEl);
	  } else if (abs(lepid)==13) {
	    lep_threshold=mu_threshold_;
	    hasTightIdandIso &= (ilep.passIdMu && ilep.passIsoMu);
	  } else continue;

	  bool hasExtraLepton(false);
	  
	  if ( hasTightIdandIso && (ilep.pt()>lep_threshold) ) {

	    if(abs(lepid)==11) { // ele scale corrections
	      double et = ilep.en_cor_en / cosh(fabs(ilep.en_EtaSC));
	      if (isMC) {
		double sigma= eScaler_.getSmearingSigma(phys.run,(fabs(ilep.en_EtaSC)<=1.447),ilep.en_R9, ilep.en_EtaSC, et,ilep.en_gainSeed,0,0);
		//Now smear the MC energy
		TRandom3 *rgen_ = new TRandom3(0);
		double smearValue = rgen_->Gaus(1, sigma) ;
		//TLorentzVector p4        
		ilep.SetPxPyPzE(ilep.Px()*smearValue, ilep.Py()*smearValue, ilep.Pz()*smearValue, ilep.E()*smearValue);
	      } else {
		double scale_corr=eScaler_.ScaleCorrection(phys.run,(fabs(ilep.en_EtaSC)<=1.447),ilep.en_R9, ilep.en_EtaSC, et,ilep.en_gainSeed); 
		//TLorentzVector p4
		ilep.SetPxPyPzE(ilep.Px()*scale_corr, ilep.Py()*scale_corr, ilep.Pz()*scale_corr, ilep.E()*scale_corr); 
	      }
	    } else if (abs(lepid)==13) { // mu scale corrections    
	      if(muCor2016){
		float qter =1.0;
		int ntrk = ilep.mn_trkLayersWithMeasurement;

		TLorentzVector p4(ilep.Px(),ilep.Py(),ilep.Pz(),ilep.E());
		//printf("Muon P4 (before roch): px=%f, py=%f, pz=%f, e=%f\n",p4.Px(),p4.Py(),p4.Pz(),p4.E());
		if (isMC) { muCor2016->momcor_mc(p4, lepid<0 ? -1 :1, ntrk, qter);}
		else { muCor2016->momcor_data(p4, lepid<0 ? -1 :1, 0, qter); }

		ilep.SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
		//printf("Muon P4 (AFTER roch): px=%f, py=%f, pz=%f, e=%f\n\n",ilep.Px(),ilep.Py(),ilep.Pz(),ilep.E());
	      }
	    }

	    nGoodLeptons++;
	    std::pair <int,LorentzVector> goodlep;
	    goodlep = std::make_pair(lepid,ilep);
	    goodLeptons.push_back(goodlep);
	  } else { // extra loose leptons

	    if (abs(lepid)==11) {
	      hasExtraLepton = (ilep.passIdLooseEl && ilep.passIsoEl && ilep.pt()>10.);
	    } else if (abs(lepid)==13) {
	      hasExtraLepton = ( (ilep.passIdLooseMu && ilep.passIsoMu && ilep.pt()>10.) || (ilep.passSoftMuon) );
	    }

	  }

	  if ( hasTightIdandIso && ilep.pt()>20. && ilep.pt()<30.) { 
            vetoLeptons.push_back(ilep); 
          }   

	  if (hasExtraLepton) {
	    nExtraLeptons++;
	    extraLeptons.push_back(ilep);
	  }
	} // leptons

	sort(goodLeptons.begin(), goodLeptons.end(), ptsortinpair());
	sort(vetoLeptons.begin(), vetoLeptons.end(), ptsort());

	mon.fillHisto("nleptons","raw", goodLeptons.size(),weight);
	mon.fillHisto("nleptons","raw_extra", extraLeptons.size(),weight);
	
	std::vector<TString> tag_cat;
	//TString tag_cat;
        int evcat=-1;
	if (goodLeptons.size()==1) evcat = getLeptonId(abs(goodLeptons[0].first));
	if (goodLeptons.size()>1) evcat = getDileptonId(abs(goodLeptons[0].first),abs(goodLeptons[1].first)); 
        switch(evcat) {
        case MUMU :
	  tag_cat.push_back("mumu");
            break;
        case EE   :
	  tag_cat.push_back("ee");
            break;
	case MU  :
	  tag_cat.push_back("mu");
	  tag_cat.push_back("l");
	    break;
	case E   :
	  tag_cat.push_back("e");
	  tag_cat.push_back("l");
	    break;
        case EMU  :
	  tag_cat.push_back("emu");
            break;
	    //default :
	    //continue;
        }
        /*
        //split inclusive DY sample into DYToLL and DYToTauTau
        if(isMC && mctruthmode==1) {
            //if(phys.genleptons.size()!=2) continue;
            if(phys.genleptons.size()==2 && isDYToTauTau(phys.genleptons[0].id, phys.genleptons[1].id) ) continue;
        }

        if(isMC && mctruthmode==2) {
            if(phys.genleptons.size()!=2) continue;
            if(!isDYToTauTau(phys.genleptons[0].id, phys.genleptons[1].id) ) continue;
        }
        */

	// All: "Raw"
	mon.fillHisto("eventflow","all",0,weight);
	mon.fillHisto("eventflow","bdt",0,weight);      

        bool hasTrigger(false);

        if(!isMC) {
	    //if(evcat!=fType) continue;

            //if(evcat==EE   && !(hasEEtrigger||hasEtrigger) ) continue;
            //if(evcat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
            //if(evcat==EMU  && !hasEMtrigger ) continue;

            //this is a safety veto for the single mu PD
            if(isSingleMuPD) {
	      if(!hasMtrigger) continue;
		//                if(hasMtrigger && hasMMtrigger) continue;
            }
            // if(isDoubleMuPD) {
            //     if(!hasMMtrigger) continue;
            // }

            //this is a safety veto for the single Ele PD
            if(isSingleElePD) {
	      if(!hasEtrigger) continue;
		//    if(hasEtrigger && hasEEtrigger) continue;
            }
            // if(isDoubleElePD) {
            //     if(!hasEEtrigger) continue;
            // }

            hasTrigger=true;

        } else { // isMC
	  // if(evcat==E   && hasEtrigger ) hasTrigger=true;
	  // if(evcat==MU && hasMtrigger ) hasTrigger=true;
	  // if(evcat==EMU  && ( hasEtrigger || hasMtrigger)) hasTrigger=true;
	  hasTrigger = (hasEtrigger || hasMtrigger);
	  if(!hasTrigger) continue;
        }

	//tags.push_back(tag_cat); //add ee, mumu, emu category

        //prepare the tag's vectors for histo filling
	//for(size_t ich=0; ich<tag_cat.size(); ich++) {
	//  tags.push_back( tag_cat[ich] );
	//}

        // pielup reweightiing
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      xsecWeight*genWeight);
        mon.fillHisto("nvtxwgt_raw",tags, phys.nvtx,      weight);

    
	// Trigger
	mon.fillHisto("eventflow","all",1,weight);
	mon.fillHisto("eventflow","bdt",1,weight); 

	// -------------------------------------------------------------------------
	// Exactly 1 good lepton
	bool passOneLepton(goodLeptons.size()==1); 
	if (!passOneLepton) continue;
	// -------------------------------------------------------------------------
	//if(goodLeptons.size()!=1) continue; // at least 1 tight leptons

        // lepton ID + ISO scale factors 
        if(isMC) {
	  if (evcat==E) {
	    weight *= lepEff.getRecoEfficiency( goodLeptons[0].second.eta(), 11).first;
	    weight *= lepEff.getLeptonEfficiency( goodLeptons[0].second.pt(), goodLeptons[0].second.eta(), 11, "tight" ,patUtils::CutVersion::ICHEP16Cut ).first ; //ID
	    
	  } else if (evcat==MU) {
	    weight *= lepEff.getTrackingEfficiency( goodLeptons[0].second.eta(), 13).first; //Tracking eff
	    weight *= lepEff.getLeptonEfficiency( goodLeptons[0].second.pt(), goodLeptons[0].second.eta(), 13, "tight" ,patUtils::CutVersion::ICHEP16Cut ).first ; //ID
	    weight *= lepEff.getLeptonEfficiency( goodLeptons[0].second.pt(), goodLeptons[0].second.eta(), 13, "tightiso",patUtils::CutVersion::ICHEP16Cut ).first; //ISO w.r.t ID
	  }
        }

	mon.fillHisto("eventflow","all",2,weight);
	mon.fillHisto("eventflow","bdt",2,weight); 
	// // -------------------------------------------------------------------------
	// // 2nd lepton veto
	// // -------------------------------------------------------------------------
	// bool pass2ndlepVeto(extraLeptons.size()==0);
	// if (!pass2ndlepVeto) continue;
	// mon.fillHisto("eventflow","all",3,weight);

	// Lepton kinematics
	if (abs(goodLeptons[0].first==11)) {
	  mon.fillHisto("leadlep_pt_raw","e",goodLeptons[0].second.pt(),weight);
	  mon.fillHisto("leadlep_eta_raw","e",goodLeptons[0].second.eta(),weight);
	} else if (abs(goodLeptons[0].first==13)) {
	  mon.fillHisto("leadlep_pt_raw","mu",goodLeptons[0].second.pt(),weight);
	  mon.fillHisto("leadlep_eta_raw","mu",goodLeptons[0].second.eta(),weight);
	}

	//Dphi(lep, MET) ?
	float dphilepmet=fabs(deltaPhi(goodLeptons[0].second.phi(),metP4.phi()));
	mon.fillHisto("dphilepmet","raw",dphilepmet,weight);

        //
        //MET AND MT ANALYSIS
        //

	LorentzVector wsum=metP4+goodLeptons[0].second;
	// mtW
	double tMass = 2.*goodLeptons[0].second.pt()*metP4.pt()*(1.-TMath::Cos(deltaPhi(goodLeptons[0].second.phi(),metP4.phi())));

	mon.fillHisto("pfmet","raw",metP4.pt(),weight);
	mon.fillHisto("mtw","raw",sqrt(tMass),weight);
	mon.fillHisto("ptw","raw",wsum.pt(),weight);

	//-------------------------------------------------------------------
	//MET>25 GeV 
	bool passMet25(metP4.pt()>25);
	if (!passMet25) continue;
	mon.fillHisto("eventflow","all",3,weight); // MEt cut
	mon.fillHisto("eventflow","bdt",3,weight);
	//-------------------------------------------------------------------

	//mtW >50 GeV
	bool passMt(sqrt(tMass)>50. && sqrt(tMass)<250.);
	if (!passMt) continue;
	mon.fillHisto("eventflow","all",4,weight); // MT cut
	mon.fillHisto("eventflow","bdt",4,weight);
	
        //
        //JET AND BTAGGING ANALYSIS
        //

	//###########################################################
	// The AK8 fat jets configuration
	//###########################################################

	//AK8 + double-b tagger fat-jet collection
	PhysicsObjectFatJetCollection DBfatJets; // collection of AK8 fat jets

	int ifjet(0);
	for (auto & ijet : fatJets ) {

	  if(ijet.pt()<jet_threshold_) continue;
	  if(fabs(ijet.eta())>2.4) continue;

	  double dR = deltaR( ijet, goodLeptons[0].second );
	  mon.fillHisto("dRlj_raw","all_fjet",dR,weight);
    
	  if (dR<0.4) continue;
  
	  ifjet++;

	  int count_sbj(0);
	    // Examine soft drop subjets in AK8 jet:
	  count_sbj = ijet.subjets.size(); // count subjets above 20 GeV only
  
	  if ( verbose ) printf("\n\n Print info for subjets in AK8 %3d : ", ifjet);

	  if ( verbose ) {
	    // loop over subjets of AK8 jet
	    for (auto & it : ijet.subjets ) {
	      printf("\n subjet in Ntuple has : pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f",   
		     it.pt(),
		     it.eta(),
		     it.Phi(),
		     it.M()
		     );
	    }
	    printf("\n\n");
	  } // verbose
  
	  mon.fillHisto("db_discrim","fjet",ijet.btag0,weight);
	  if (ijet.motherid == 36) mon.fillHisto("db_discrim","fjet_true",ijet.btag0,weight);
	  mon.fillHisto("nsubjets_raw","fjet",count_sbj,weight);
	  if (ijet.motherid == 36) mon.fillHisto("nsubjets_raw","fjet_true",count_sbj,weight);
	  mon.fillHisto("sd_mass","fjet",ijet.softdropM,weight);
	  if (ijet.motherid == 36) mon.fillHisto("sd_mass","fjet_true",ijet.softdropM,weight);
	  mon.fillHisto("pruned_mass","fjet",ijet.prunedM,weight);
	  if (ijet.motherid == 36) mon.fillHisto("pruned_mass","fjet_true",ijet.prunedM,weight);

	  if (ijet.softdropM<=40.) {
	    if (ijet.softdropM>=7.)  {
	      mon.fillHisto("db_discrim","fjet_lowm",ijet.btag0,weight);
	      if (ijet.motherid == 36)  mon.fillHisto("db_discrim","fjet_true_lowm",ijet.btag0,weight);
	    }
	  } else {
	    mon.fillHisto("db_discrim","fjet_highm",ijet.btag0,weight);
	    if (ijet.motherid == 36) mon.fillHisto("db_discrim","fjet_true_highm",ijet.btag0,weight);
	  }

	  bool hasDBtag(ijet.btag0>DBMediumWP);
  
	  if (ijet.softdropM>=7. && ijet.softdropM<=40. && count_sbj>0) {
	    if (hasDBtag) DBfatJets.push_back(ijet);
	  }
	  
	} // AK8 fatJets loop

	//--------------------------------------------------------------------------
	// AK8 + double-b jets
	sort(DBfatJets.begin(), DBfatJets.end(), ptsort());

	mon.fillHisto("nbjets_raw","nfatJet", DBfatJets.size(),weight);

	int is(0);
	for (auto & jet : DBfatJets) {
	   mon.fillHisto("jet_pt_raw", "fat"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "fat"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}
	//--------------------------------------------------------------------------


        //###########################################################
        // AK4 jets ,
        // AK4 jets + CSVloose b-tagged configuration
        //###########################################################

        PhysicsObjectJetCollection GoodIdJets;
        PhysicsObjectJetCollection CSVLoosebJets; // used to define the SRs

        int nJetsGood30(0);
        int nCSVLtags(0),nCSVMtags(0),nCSVTtags(0);
        double BTagWeights(1.0);

	float mindphijmet(999.);
	for(size_t ijet=0; ijet<corrJets.size(); ijet++) {

	  if(corrJets[ijet].pt()<jet_threshold_) continue;
	  if(fabs(corrJets[ijet].eta())>2.4) continue;
  
	  //jet ID
	  if(!corrJets[ijet].isPFLoose) continue;
	  //if(corrJets[ijet].pumva<0.5) continue;
  
	  // //check overlaps with selected leptons
	  bool hasOverlap(false);
	  for(size_t ilep=0; ilep<goodLeptons.size(); ilep++) {
	    double dR = deltaR( corrJets[ijet], goodLeptons[ilep].second );
	    mon.fillHisto("dRlj_raw","all",dR,weight);
    
	    if (abs(goodLeptons[ilep].first)==11) hasOverlap = (dR<0.2); // within 0.2 for electrons
	    if (abs(goodLeptons[ilep].first)==13) hasOverlap = (dR<0.4); // within 0.4 for muons
	  }
	  if(hasOverlap) continue;
  
	  // if (vetoLeptons.size()>0) {  
	  //   double dR_thr=deltaR(corrJets[ijet],vetoLeptons[0]); 
	  //   if (dR_thr<0.4) continue; // reject jet if found close to e/mu below 30 GeV 
	  // }

	  GoodIdJets.push_back(corrJets[ijet]);
	  if(corrJets[ijet].pt()>30) nJetsGood30++;

	  // Dphi (j,met)
	  float dphijmet=fabs(deltaPhi(corrJets[ijet].phi(),metP4.phi()));
	  if (dphijmet<mindphijmet) mindphijmet=dphijmet;

  
	  // B-tagging
	  bool hasCSVtag;
          double btag_dsc = -1;
          if ( use_DeepCSV ) {btag_dsc = corrJets[ijet].btag1;} else {btag_dsc = corrJets[ijet].btag0;}
	  nCSVLtags += (btag_dsc>LooseWP);
	  nCSVMtags += (btag_dsc>MediumWP);
          nCSVTtags += (btag_dsc>TightWP);
          mon.fillHisto("b_discrim",b_tagging_name,btag_dsc,weight);
          if (corrJets[ijet].motherid == 36) mon.fillHisto("b_discrim",b_tagging_name+"_true",btag_dsc,weight);
          hasCSVtag = btag_dsc>LooseWP;

	  if (isMC && (corrJets[ijet].pt()>30.) ) {
	    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
	    btsfutil.SetSeed(ev.event*10 + ijet*10000);
    
	    if(abs(corrJets[ijet].flavid)==5) {
	      //  80X recommendation
	      btsfutil.modifyBTagsWithSF(hasCSVtag , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_B ,
										 corrJets[ijet].eta(), corrJets[ijet].pt()), beff);
	    } else if(abs(corrJets[ijet].flavid)==4) {
	      //  80X recommendation
	      btsfutil.modifyBTagsWithSF(hasCSVtag , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_C ,
										 corrJets[ijet].eta(), corrJets[ijet].pt()), beff);
	    } else {
	      //  80X recommendation
	      btsfutil.modifyBTagsWithSF(hasCSVtag , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_UDSG ,
										 corrJets[ijet].eta(), corrJets[ijet].pt()), leff);
	    }
	  } // isMC
  
	    // Fill b-jet vector:

	  if (hasCSVtag) {
	      /*
	      if (runDBversion) {
		float dRmin(999.);
		for (auto & ifb : DBfatJets) {
		  for (auto & it : ifb.subjets) { // subjets loop
		    float dR = deltaR(corrJets[ijet], it);
		    if (dR<dRmin) dRmin=dR;
		  }//subjets
		} // AK8
		if (dRmin>0.4) CSVLoosebJets.push_back(corrJets[ijet]);
	      }
              else {
	      */ 
	    CSVLoosebJets.push_back(corrJets[ijet]); 
	      //}
	  }
	  //} // b-jet loop
	} // jet loop
    

	//--------------------------------------------------------------------------
	// AK4 jets:
	sort(GoodIdJets.begin(), GoodIdJets.end(), ptsort());
	// Fill Histograms with AK4,AK4 + CVS, AK8 + db basics:
	mon.fillHisto("njets_raw","nj", GoodIdJets.size(),weight);

	// AK4 jets pt:
	is=0;
	for (auto & jet : GoodIdJets) {
	  mon.fillHisto("jet_pt_raw", "jet"+htag[is], jet.pt(),weight); 
	  mon.fillHisto("jet_eta_raw", "jet"+htag[is], jet.eta(),weight); 
	  mon.fillHisto("jet_phi_raw","jet"+htag[is], jet.phi(),weight); 
	
	  if (jet.pt()<30.) {
	    mon.fillHisto("jet_pt_raw", "pt_20to30_"+htag[is], jet.pt(),weight); 
	    mon.fillHisto("jet_eta_raw", "pt_20to30_"+htag[is], jet.eta(),weight); 
	    mon.fillHisto("jet_phi_raw", "pt_20to30_"+htag[is], jet.phi(),weight);
	  }
	  is++; 
	  if (is>3) break; // plot only up to 4 b-jets ?                                                                                                                                                                                   
        }

	//--------------------------------------------------------------------------
	// AK4 + CSV jets:
	sort(CSVLoosebJets.begin(), CSVLoosebJets.end(), ptsort());
	mon.fillHisto("nbjets_raw","nb", CSVLoosebJets.size(),weight);


        //-------------------------------------------------------------------
        // AK4 + CSV jets 
        is=0; 
        for (auto & jet : CSVLoosebJets) {
          mon.fillHisto("jet_pt_raw", b_tagging_name+htag[is], jet.pt(),weight); 
          mon.fillHisto("jet_eta_raw", b_tagging_name+htag[is], jet.eta(),weight); 
	  mon.fillHisto("jet_phi_raw", b_tagging_name+htag[is], jet.phi(),weight);

          if (jet.pt()<30.) { 
            mon.fillHisto("jet_pt_raw", "pt_20to30_"+b_tagging_name+htag[is], jet.pt(),weight); 
            mon.fillHisto("jet_eta_raw", "pt_20to30_"+b_tagging_name+htag[is], jet.eta(),weight); 
	    mon.fillHisto("jet_phi_raw", "pt_20to30_"+b_tagging_name+htag[is], jet.phi(),weight); 
          } 
	  is++;
          if (is>3) break; // plot only up to 4 b-jets ?
	}
	
	//--------------------------------------------------------------------------
	// dphi(jet,MET)
	mon.fillHisto("dphijmet","raw",mindphijmet,weight);

	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	// minDR between a b-jet and AK8 jet

	int ibs(0);
	for (auto & ib : CSVLoosebJets) {
  
	  float dRmin(999.);
	  float dRmin_sub(999.);
  
	  for (auto & it : DBfatJets) {
	    float dR = deltaR(ib, it);
	    if (dR<dRmin) dRmin=dR;

	    // loop in subjets
	    for (auto & isub : it.subjets){
	      if (isub.pt()<20.) continue;
	      float dRsub = deltaR(ib, isub);
	      if (dRsub<dRmin_sub) dRmin_sub=dRsub;
	    }

	  }
	  mon.fillHisto("dR_raw","drmin"+htag[ibs],dRmin, weight);
	  mon.fillHisto("dR_raw","drmin_sub"+htag[ibs],dRmin_sub, weight);

	  ibs++;
	  if (ibs>3) break; // plot only up to 4 b-jets ?
	} 


	//###########################################################
	// Soft b-jets from SVs configuration
	//###########################################################

	// SVs collection
	PhysicsObjectSVCollection SVs;
	PhysicsObjectSVCollection SVs_raw; // non-cross-cleaned secondary vertices
	
	
	for (auto & isv : secVs) {

	  if (isv.pt()>=jet_threshold_) continue; // SV pT>20 GeV
  
	  mon.fillHisto("softb_ntrk","raw",isv.ntrk,weight);
	  if (isv.ntrk<3) continue; // nTrks associated to SV >= 3

	  if ( verbose ) {

	    printf("\n SV has : pt=%6.1f, ntrk=%3d, dxy=%7.3f, dxyz_signif=%7.3f, isv.cos_dxyz_p=%7.3f",
		   isv.pt(),   
		   isv.ntrk,
		   isv.dxy,
		   isv.dxyz_signif,
		   isv.cos_dxyz_p
		   );
	    
	  } // verbose 
  
	  // check overlap with any other jet
	  bool hasOverlap(false);
 
	  mon.fillHisto("softb_dxy","raw",isv.dxy,weight);
	  if (isv.sv_mc_mcbh_ind>0)mon.fillHisto("softb_dxy","raw_true",isv.dxy,weight);
	  mon.fillHisto("softb_dxyz_signif","raw",isv.dxyz_signif,weight);
	  if (isv.sv_mc_mcbh_ind>0) mon.fillHisto("softb_dxyz_signif","raw_true",isv.dxyz_signif,weight);
	  mon.fillHisto("softb_cos","raw",isv.cos_dxyz_p,weight);
	  if (isv.sv_mc_mcbh_ind>0) mon.fillHisto("softb_cos","raw_true",isv.cos_dxyz_p,weight);
  
	  if (isv.dxy>3.) continue;
	  if (isv.dxyz_signif<4.) continue;
	  if (isv.cos_dxyz_p<0.98) continue;

	  SVs_raw.push_back(isv);
	  
	  // plot minDR(SV,b)
	  float dRmin_csv(999.);
	  for (auto & it : CSVLoosebJets) {
	    double dR=deltaR(it, isv);
	    if (dR<dRmin_csv) dRmin_csv=dR;
	  }
	  mon.fillHisto("dR_raw","sv_b",dRmin_csv,weight);

	  //if (!runDBversion) { // use soft-b tags only if AK8 jets are not used
	  hasOverlap=(dRmin_csv<0.4);
	  if (!hasOverlap) {// continue;
	    // Fill final soft-bs from SVs
	    SVs.push_back(isv);
	  }
	  //}

	}

	//--------------------------------------------------------------------------
	// Soft-bs properties
	//--------------------------------------------------------------------------

	sort(SVs.begin(), SVs.end(), ptsort());
	sort(SVs_raw.begin(), SVs_raw.end(), ptsort());

	mon.fillHisto("nbjets_raw","nb_soft",SVs.size(),weight);

	is=0;
	for (auto & isv : SVs) {
	   mon.fillHisto("softjet_pt_raw", "softb"+htag[is], isv.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "softb"+htag[is], isv.eta(),weight);
	   is++;
	   if (is>3) break;
	}

	// DR between 2 SVs
	if (SVs.size()>1) {
	  double dR=deltaR(SVs[0], SVs[1]);
	  mon.fillHisto("dR_raw","svs",dR,weight);
	}


	//-------------------------------------------------------------------
	//-------------------------------------------------------------------
	
        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

	//At least 2 jets
	if (GoodIdJets.size()<2) continue;
	sort(GoodIdJets.begin(), GoodIdJets.end(), btagsort());

	// Use highest b-tagged jet 0.55<b-tag_high<0.8
	double btag_high(-1.);
	if (use_DeepCSV) {
	  btag_high=GoodIdJets[0].btag1;
	} else {
	  btag_high=GoodIdJets[0].btag0;
	}

	bool btag_sideband((btag_high<=0.7 && btag_high>=0.57) || (btag_high<=0.5 && btag_high>=0.3));
	
	is=0;
	for (auto & jet : GoodIdJets) {
	  if (use_DeepCSV) {
	    mon.fillHisto("b_discrim",b_tagging_name+htag[is],jet.btag1,weight);
	  } else {
	    mon.fillHisto("b_discrim",b_tagging_name+htag[is],jet.btag0,weight);
	  }
	  is++;
	  if (is>3) break;
	}
	//-------------------------------------------------------------------
	//-------------------------------------------------------------------
	// First, set all b-jets (x-cleaned) in one vector<LorentzVector>
	vector<LorentzVector> GoodIdbJets;
	
	// SRs: (nj,3b), (nj, 4b), ... or
	// Top CRs: (5j, 2b), (4j, 2b), (3j, 2b), ...
	if (nCSVMtags>=1) {
	  if (CSVLoosebJets.size()>=2 ) {
	    
	    if (CSVLoosebJets.size()==2 && SVs.size()==0) { // Top Control Regions
	      
	      for (auto & i : GoodIdJets) {
		GoodIdbJets.push_back(i);
	      }
	      
	    } else { // Signal Region GoodIdbJets
	      for (auto & i : CSVLoosebJets) 
		{
		  GoodIdbJets.push_back(i);
		} // AK4 + CSV
	      for (auto & i : SVs) {
		GoodIdbJets.push_back(i);
	      } // soft-b from SV
	    }
	  
	    // At least 2 jets and 2 b-jets
	    mon.fillHisto("eventflow","all",5,weight); 
	    mon.fillHisto("eventflow","bdt",5,weight);

	  } else { continue; } // At least 2 CSVv2 b-jets

	} else if ( btag_sideband) { // non-TT backgrouns (W, DY, QCD) CRs:
	  
	  for (auto & i : GoodIdJets) {
	    GoodIdbJets.push_back(i);
	  }

	} else {
	  //printf("\n Unknown category, please check \n");
	  continue;
	}



	//-------------------------------------------------------------------
	//At least 2 jets and 2 b-jets
	
	//if (GoodIdJets.size()<2 || CSVLoosebJets.size()<2) continue;
	//if (nCSVMtags<1) continue; // At least 2 CSVv2 b-jets with LooseWP(0.54) and at least 1 satisfying the MediumWP(0.80)
	// mon.fillHisto("eventflow","all",5,weight); 
	// mon.fillHisto("eventflow","bdt",5,weight);
	
	//-------------------------------------------------------------------
	/*
	if (runDBversion) {
	  for (auto & i : DBfatJets) {
	    for (auto & it : i.subjets ) {
	      GoodIdbJets.push_back(it);
	    } // subjets
	  }// AK8 jet
	}
	*/

	// // 2D plots
	// //-------------------------------------------------------------------
	// mon.fillHisto("nbjets_2D","cat_raw",GoodIdJets.size(),GoodIdbJets.size(),weight);
	// mon.fillHisto("nbjets_raw","merged",GoodIdbJets.size(),weight);
	// //-------------------------------------------------------------------
	// mon.fillHisto("nbjets_2D","cat2_raw",GoodIdbJets.size(),DBfatJets.size(),weight);
	// mon.fillHisto("nbjets_2D","cats_raw",CSVLoosebJets.size(),SVs.size(),weight);
	
	// //-------------------------------------------------------------------
	// //-------------------------------------------------------------------
	// // Merged CSV jets + SV
	// is=0;
	// for (auto & jet : GoodIdbJets) {
	//    mon.fillHisto("jet_pt_raw", "merged"+htag[is], jet.pt(),weight);
	//    mon.fillHisto("jet_eta_raw", "merged"+htag[is], jet.eta(),weight);
	//    is++;
	//    if (is>3) break; // plot only up to 4 b-jets ?
	// }

	//-------------------------------------------------------------------                                                                                                  

        //##############################################
        //########  Main Event Selection        ########
        //##############################################

        //-------------------------------------------------------------------
        // At least 3 b-tags
	if (GoodIdbJets.size()<3) continue;
	
	TString ch;
	if (evcat==E) { ch="E_";}
	else if (evcat==MU) { ch="MU_"; }
	else { printf("UNKNOWN lepton category - please check\n"); }

	bool isSignalRegion(true);

	if (nCSVMtags>=1) {
	  if (CSVLoosebJets.size()>2 || SVs.size()>0) {
	    // SR categories
	    mon.fillHisto("eventflow","all",6,weight); 
	    
	    // Cats: 3b
	    if (GoodIdbJets.size()==3) { 
	      tags.push_back("SR_3b");
	      tags.push_back(ch+"SR_3b"); 
	    }
	    else {
	      tags.push_back("SR_geq4b"); 
	      tags.push_back(ch+"SR_geq4b"); 

	      mon.fillHisto("eventflow","all",7,weight);
	      
	      if (GoodIdbJets.size()==4) { 
		tags.push_back("SR_4b"); 
		tags.push_back(ch+"SR_4b"); 
	      }
	      else {
		if (GoodIdbJets.size()==5) { 
		  tags.push_back("SR_5b"); 
		  tags.push_back(ch+"SR_5b");
		}
		tags.push_back("SR_geq5b");
		tags.push_back(ch+"SR_geq5b");  
	      }
	    }
	  } else { 
	    // thats the Top Control Regions
	    isSignalRegion=false;
	    
	    // Top Control Region categories
	    if (GoodIdbJets.size()==3) { 
	      tags.push_back("CR_3b");
	      tags.push_back(ch+"CR_3b");  
	    }
	    else {
	      tags.push_back("CR_geq4b");
	      tags.push_back(ch+"CR_geq4b");
	      if (GoodIdbJets.size()==4) { 
		tags.push_back("CR_4b"); 
		tags.push_back(ch+"CR_4b"); 
	      }
	      else {
		if (GoodIdbJets.size()==5) { 
		  tags.push_back("CR_5b"); 
		  tags.push_back(ch+"CR_5b");
		}
		tags.push_back("CR_geq5b");
		tags.push_back(ch+"CR_geq5b"); 
	      }
	    }
	  }
	} else if (btag_sideband) {  // thats the non-TT (W,DY,QCD) Control Regions
	  
	  isSignalRegion=false;

	  // Non-TT Control Region categories
	  if (GoodIdbJets.size()==3) { 
	    tags.push_back("CR_nonTT_3b");
	    tags.push_back(ch+"CR_nonTT_3b");
	  }
	  else {
	    tags.push_back("CR_nonTT_geq4b");
	    tags.push_back(ch+"CR_nonTT_geq4b");
	    if (GoodIdbJets.size()==4) { 
	      tags.push_back("CR_nonTT_4b"); 
	      tags.push_back(ch+"CR_nonTT_4b"); 
	    }
	    else {
	      if (GoodIdbJets.size()==5) { 
		tags.push_back("CR_nonTT_5b"); 
		tags.push_back(ch+"CR_nonTT_5b");
	      }
	      tags.push_back("CR_nonTT_geq5b");
	      tags.push_back(ch+"CR_nonTT_geq5b");
	    }
	  }
	} else {
	  tags.push_back("UNKNOWN");
	  printf("\n Unknown category, please check \n");
        }


	// Here define all variables 
        LorentzVector allHadronic;
        //std::pair <int,LorentzVector> pairHadronic;

	// HT from all CSV + soft b's
	float ht(0.); float ht_b30(0.);
        
	// Hadronic vector sum:
	int countb(0);
	for (auto & thisb : GoodIdbJets) {
	  allHadronic+=thisb;
	  countb++; if (countb>3) break;
	}
	// Hadronic scalar sum (HT):
	for (auto & thisb : GoodIdbJets) {
	  ht+=thisb.pt();
	  if (thisb.pt()>30.) ht_b30+=thisb.pt();
	}
 
        //-----------------------------------------------------------
        // Control plots
        //----------------------------------------------------------

        // 3,4 b's pT
        mon.fillHisto("nbjets_raw",tags,GoodIdbJets.size(),weight);
        is=0;
        for (auto & jet : GoodIdbJets) 
        {
            mon.fillHisto("jet_pt_raw", "merged_final"+htag[is], jet.pt(),weight);
            mon.fillHisto("jet_eta_raw", "merged_final"+htag[is], jet.eta(),weight);
            is++;
            if (is>3) break; // plot only up to 4 b-jets ?
        }
 
        // higgs mass
        mon.fillHisto("higgsMass",tags,allHadronic.mass(),weight);
        // higgs pT
        mon.fillHisto("higgsPt",tags,allHadronic.pt(),weight);
        // HT from all CSV + soft b's
        mon.fillHisto("ht",tags,ht,weight);
	mon.fillHisto("ht_b30",tags,ht_b30,weight);
        // MET
        mon.fillHisto("pfmet",tags,metP4.pt(),weight);
        // dphi(jet,MET)
        mon.fillHisto("dphijmet",tags,mindphijmet,weight);
        // pTW
        //LorentzVector wsum=metP4+goodLeptons[0].second;
        mon.fillHisto("ptw",tags,wsum.pt(),weight);
        // mtW 
        mon.fillHisto("mtw",tags,sqrt(tMass),weight);
        // Dphi(W,h) instead of DRmin(l,b)
        double dphi_Wh=fabs(deltaPhi(allHadronic.phi(),wsum.phi()));
        mon.fillHisto("dphiWh",tags,dphi_Wh,weight);

        // DR(bb)_average
        vector<float> dRs;
        float dm(0.);

	if (GoodIdbJets.size()==3) {
	  dRs.push_back(deltaR(GoodIdbJets[0],GoodIdbJets[1]));
	  dRs.push_back(deltaR(GoodIdbJets[0],GoodIdbJets[2]));
	  dRs.push_back(deltaR(GoodIdbJets[1],GoodIdbJets[2]));
	}
        else if (GoodIdbJets.size()>=4) 
        {
            dRs.push_back(deltaR(GoodIdbJets[0],GoodIdbJets[3]));
            dRs.push_back(deltaR(GoodIdbJets[1],GoodIdbJets[3]));
            dRs.push_back(deltaR(GoodIdbJets[2],GoodIdbJets[3]));

            float dm1 = fabs( (GoodIdbJets[0]+GoodIdbJets[1]).mass() - (GoodIdbJets[2]+GoodIdbJets[3]).mass() );
            float dm2 = fabs( (GoodIdbJets[0]+GoodIdbJets[2]).mass() - (GoodIdbJets[1]+GoodIdbJets[3]).mass() );

            dm1 = min(dm1, dm2);
            dm2 = fabs( (GoodIdbJets[0]+GoodIdbJets[3]).mass() - (GoodIdbJets[1]+GoodIdbJets[2]).mass() );
            dm = min(dm1, dm2);
        }

        float dRave_(0.);
        for (auto & it : dRs)
        {
            dRave_+=it;
        }
        dRave_/=dRs.size();
        mon.fillHisto("dRave",tags,dRave_,weight);
        mon.fillHisto("dmmin",tags,dm, weight);

	//##############################################################################
        //############ MVA Reader #####################################################
	//##############################################################################

        float mvaBDT(-10.0);
        if (GoodIdbJets.size() == 3)
        {
            mvaBDT = myTribTMVAReader.GenReMVAReader
                     (
                      wsum.pt(),
                      allHadronic.mass(), allHadronic.pt(), dRave_, dm, ht,
                      dphi_Wh,
                      "Haa4bSBClassificationTribMVA"
                     );
        }
        else if (GoodIdbJets.size() >= 4)
        {
            mvaBDT = myQuabTMVAReader.GenReMVAReader
                     (
                      wsum.pt(),
                      allHadronic.mass(), allHadronic.pt(), dRave_, dm, ht,
                      dphi_Wh,
                      "Haa4bSBClassificationQuabMVA"
                     );
        }
        //else continue;
        mon.fillHisto("bdt", tags, mvaBDT, weight);

	if (GoodIdbJets.size() == 3) 
	{
	  if (mvaBDT>0.19) mon.fillHisto("eventflow","bdt",6,weight);  
	}
	else if (GoodIdbJets.size() >= 4) 
	{
	  if (mvaBDT>0.14) mon.fillHisto("eventflow","bdt",7,weight);  
	}

	//##############################################################################
        //############ MVA Handler ####################################################
	//##############################################################################

	if (runMVA)
        {
	  float mvaweight = 1.0;
	  genWeight > 0 ? mvaweight = weight/xsecWeight : mvaweight = -weight / xsecWeight; // Include all weights except for the xsecWeight
	  if ( isSignalRegion && GoodIdbJets.size() >= 3 )
	  {
	    myMVAHandler_.getEntry
	    (
		GoodIdbJets.size() == 3, GoodIdbJets.size() >= 4, // 3b cat, 4b cat
		wsum.pt(), //W only, w pt
		allHadronic.mass(), allHadronic.pt(), dRave_, dm, ht, //Higgs only, higgs mass, higgs pt, bbdr average, bb dm min, sum pt from all bs
		dphi_Wh, //W and H, dr 
                mvaweight, //note, since weight is not the weight we want, we store all others except xSec weigh
                ev.lheNJets //AUX variable for weight calculation
	    );
	    myMVAHandler_.fillTree();
	  }
	}

        //##############################################################################
        //### HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //##############################################################################

	// LOOP ON SYSTEMATIC VARIATION FOR THE STATISTICAL ANALYSIS
	for(size_t ivar=0; ivar<nvarsToInclude; ivar++)
        {
	  if(!isMC && ivar>0 ) continue; //loop on variation only for MC samples

	  //scan the BDT cut and fill the shapes
	  for(unsigned int index=0;index<optim_Cuts1_bdt.size();index++){
	    if(mvaBDT>optim_Cuts1_bdt[index]){
	      mon.fillHisto(TString("bdt_shapes")+varNames[ivar],tags,index, mvaBDT,weight);
	    }
	  }
	}

        //##############################################
        // recompute MET/MT if JES/JER was varied
        //##############################################
        //LorentzVector vMET = variedMET[ivar>8 ? 0 : ivar];
        //PhysicsObjectJetCollection &vJets = ( ivar<=4 ? variedJets[ivar] : variedJets[0] );
    } // loop on all events END

    printf("\n");
    file->Close();

    //write MVA files
    if (runMVA) { myMVAHandler_.writeTree(); }

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################
    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results saved in %s\n", outUrl.Data());

    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();
    ofile->Close();

    if ( outTxtFile_final ) fclose(outTxtFile_final);
}
