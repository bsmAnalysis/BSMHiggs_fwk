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
//#include "UserCode/bsmhiggs_fwk/interface/EventCategory.h"
#include "UserCode/bsmhiggs_fwk/interface/statWgt.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
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

struct btagsort: public std::binary_function<PhysicsObject_Jet, PhysicsObject_Jet, float> 
{
  bool operator () (const PhysicsObject_Jet & x, PhysicsObject_Jet & y) 
  { return  ( x.btag0 > y.btag0 ) ; }
}; 

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

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }
   
    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    //AutoLibraryLoader::enable();
    FWLiteEnabler::enable();

    // configure the process
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  
    bool isMC = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");

    double xsec = runProcess.getParameter<double>("xsec");

    TString proc=runProcess.getParameter<std::string>("proc");
    TString dtag=runProcess.getParameter<std::string>("tag");
    TString suffix=runProcess.getParameter<std::string>("suffix");

    bool verbose = runProcess.getParameter<bool>("verbose");

    bool runCR = runProcess.getParameter<bool>("runControl");
    
    bool use_DeepCSV = runProcess.getParameter<bool>("useDeepCSV"); // Will set DeepCSV as the default b-tagger automaticaly
    bool usemetNoHF = runProcess.getParameter<bool>("usemetNoHF");
    
    TString url = runProcess.getParameter<std::string>("input");
    TString outFileUrl( dtag ); //gSystem->BaseName(url));
    //    outFileUrl.ReplaceAll(".root","");
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

    bool isMC_WJets = isMC && (string(url.Data()).find("MC13TeV_WJets")  != string::npos);
    bool isMC_DY = isMC && (string(url.Data()).find("MC13TeV_DY")  != string::npos);
    
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

    // setup calibration readers 80X
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
    //    gSystem->ExpandPathName(uncFile);
    cout << "Loading jet energy scale uncertainties from: " << jecDir << endl;

    if(dtag.Contains("2016B") || dtag.Contains("2016C") ||dtag.Contains("2016D")) jecDir+="Summer16_80X/Summer16_23Sep2016BCDV4_DATA/";
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
 
     // generator level plots
    mon.addHistogram( new TH1F( "pileup", ";pileup;Events", 100,-0.5,99.5) );
    
    mon.addHistogram( new TH1F( "higgsMass",";m_{h} [GeV];Events",30,0.,600.) );
    mon.addHistogram( new TH1F( "higgsPt",";p_{T}^{h} [GeV];Events",30,0.,600.));
    // mon.addHistogram( new TH1F( "higgsEta",";#eta (h);Evenets",100,-5,5) );
    // if (dogen) {
    //   mon.addHistogram( new TH1F( "a1mass",";m_{a1} [GeV];Events",400,0.,200.) );
    //   mon.addHistogram( new TH1F( "a2mass",";m_{a2} [GeV];Events",400,0.,200.) );
    //   mon.addHistogram( new TH1F( "a1pt",";p_{T}^{a1};Events",100,0.,600.));
    //   mon.addHistogram( new TH1F( "a2pt",";p_{T}^{a2};Events",100,0.,600.));
    //   mon.addHistogram( new TH1F( "aabalance",";p_{T}^{a1}/p_{T}^{a2};Events",200,0.,10.));
    //   mon.addHistogram( new TH1F( "a1DR",";#Delta R1(b,#bar{b});Events",100,0.,5.));
    //   mon.addHistogram( new TH1F( "a2DR",";#Delta R2(b,#bar{b});Events",100,0.,5.)); 
    //   mon.addHistogram( new TH1F( "aaDR",";#Delta R(a_{1},a_{2});Events",100,0.,5.)); 
      
    //   mon.addHistogram( new TH1F( "b1pt",";p_{T}^{b1};Events",60,0.,600.));   
    //   mon.addHistogram( new TH1F( "b2pt",";p_{T}^{b2};Events",60,0.,600.)); 
    //   mon.addHistogram( new TH1F( "b3pt",";p_{T}^{b3};Events",60,0.,600.)); 
    //   mon.addHistogram( new TH1F( "b4pt",";p_{T}^{b4};Events",60,0.,600.)); 
    // }
    /*

    mon.addHistogram( new TH1F( "zpt_raw",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "mt_raw",       ";#it{m}_{T} [GeV];Events", 100,0,2000) );
    double MTBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000,2000};
    const int nBinsMT = sizeof(MTBins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "mt2_raw",       ";#it{m}_{T} [GeV];Events", nBinsMT,MTBins) );
    mon.addHistogram( new TH1F( "zmass_raw",    ";#it{m}_{ll} [GeV];Events", 100,40,250) );
    */
    
     // RECO level
    mon.addHistogram( new TH1F( "dR_raw",";#Delta R(SV,b);Events",50,0.,5.));
    mon.addHistogram( new TH1F( "dRlj_raw",";#Delta R(lep,jet);Events",100,0.,5.));

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l} [GeV];Events", 30,0.,600.) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );
    
    mon.addHistogram( new TH1F( "jet_pt_raw", ";#it{p}_{T} [GeV];Events",30,0.,600.) );
    mon.addHistogram( new TH1F( "softjet_pt_raw", ";#it{p}_{T} [GeV];Events",20,0.,40.) );
    mon.addHistogram( new TH1F( "jet_eta_raw",";#eta;Events", 70,-3,3) );

    mon.addHistogram( new TH1F( "b_discrim"," ;b discriminator;",50,0,1.) );
    mon.addHistogram( new TH1F( "db_discrim"," ;double-b discriminator;",25,-1.,1.) );
    mon.addHistogram( new TH1F( "sd_mass"," ;soft-drop Mass;",60,0.,300.) );
    mon.addHistogram( new TH1F( "pruned_mass"," ;pruned Mass;",60,0.,300.) );
    mon.addHistogram( new TH1F( "softb_ntrk"," ; SV Ntrks;",21,-0.5,21.5) );
    mon.addHistogram( new TH1F( "softb_dxy"," ; SV dxy;",50,0.,20.) );
    mon.addHistogram( new TH1F( "softb_dxyz_signif"," ; SVSIP3D;",50,1.,100.) );
    mon.addHistogram( new TH1F( "softb_cos"," ; SV cos((PV,SV),p_{SV});",25,-1.,1.) );

    mon.addHistogram( new TH1F( "nvtx_raw",";Vertices;Events",100,0,100) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",";Vertices;Events",100,0,100) );
    mon.addHistogram( new TH1F( "pfmet",    ";E_{T}^{miss} [GeV];Events", 60,0.,600.) );
    mon.addHistogram( new TH1F( "ht",    ";H_{T} [GeV];Events", 30,0.,600.) );
    mon.addHistogram( new TH1F( "mtw",       ";#it{m}_{T}^{W} [GeV];Events", 60,0.,600.) );
    mon.addHistogram( new TH1F( "ptw",       ";#it{p}_{T}^{W} [GeV];Events",30,0.,600.) );
    mon.addHistogram( new TH1F( "dphiWh", ";#Delta#it{#phi}(#it{W},h);Events", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dRave",";#Delta R(b,b)_{ave};Events",50,0.,5.));
    mon.addHistogram( new TH1F( "dmmin",";#Delta m_{b,b}^{min};Events",25,0.,250.));
    mon.addHistogram( new TH1F( "dphijmet", ";|#Delta#it{#phi}(jet,E_{T}^{miss})|;#jet", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphilepmet", ";|#Delta#it{#phi}(lep,E_{T}^{miss})|;Events", 20,0,TMath::Pi()) );

    //MVA
    mon.addHistogram( new TH1F( "MVABDT", "BDT", 100, -0.5, 0.5) );

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
    
    /*
    // preselection plots
    double METBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
    const int nBinsMET = sizeof(METBins)/sizeof(double) - 1;

    double METBins2[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000};
    const int nBinsMET2 = sizeof(METBins2)/sizeof(double) - 1;

    double MET2Bins[]= {0,80,160,240,320,400,480,560,640,800,1200};
    const int xnBinsMET2 = sizeof(MET2Bins)/sizeof(double) - 1;

    double MT2Bins[]= {0,100,200,300,400,500,600,700,800,1000,1200};
    const int xnBinsMT2 = sizeof(MT2Bins)/sizeof(double) - 1;
    */

    //--------------------------------------------------------------------------
    // some strings for tagging histograms:
    const char* astr[] = {"_b1","_b2","_b3","_b4"};
    std::vector<TString> htag(astr, astr+4);
    //--------------------------------------------------------------------------
    
    // btaging efficiency

    //#################################################
    //############# CONTROL PLOTS #####################
    //#################################################


    //##################################################################################
    //########################## STUFF FOR CUTS OPTIMIZATION  ##########################
    //##################################################################################


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
      xsecWeight = 0.; // disable MC sample if not present in the map
      /*
      double totalNumberofEvents(0.);
      
      TH1F* nevtH = (TH1F *) file->Get("mainNtuplizer/nevents");
      totalNumberofEvents = nevtH->GetBinContent(1);
      TH1F* posH = (TH1F *) file->Get("mainNtuplizer/n_posevents");
      TH1F* negH = (TH1F *) file->Get("mainNtuplizer/n_negevents");
      if(posH && negH) cnorm = posH->GetBinContent(1) - negH->GetBinContent(1);
      if(rescaleFactor>0) cnorm /= rescaleFactor;
      printf("cnorm = %f and totalNumberOfEvents= %f\n",cnorm, totalNumberofEvents);
      
      //xsecWeight=xsec/totalNumberofEvents;
      xsecWeight=xsec/cnorm; // effective luminosity
      */
      
      std::map<std::string, int> xsec_map = mStat;

      // std::string myproc = proc.Data();
      //   std::cout << "Runnin process " << myproc << std::endl;
      std::map<std::string, int>::iterator it;
      for ( it = xsec_map.begin(); it != xsec_map.end(); it++ ) {
	if (it->first == proc.Data()) {
	  xsecWeight = (xsec/(float)it->second);
	  //	  totalNumberofEvents = it->second;
	  std::cout << "weight = " << (xsecWeight*35866.9) << std::endl;
	}
      }
      
      //      float pereventwgt=(xsecWeight*35866.9);
      // printf("\n Running process with xSec = %f , and totalNumEvents = %d  . Per event weight is (L=35.9 fb-1): %f \n\n",
      //	     xsec, totalNumberofEvents, pereventwgt );
    }
    //  Hcutflow->SetBinContent(1,cnorm);

    //pileup weighting
    TString PU_Central = runProcess.getParameter<std::string>("pu_central");
    gSystem->ExpandPathName(PU_Central);
    cout << "Loading PU weights Central: " << PU_Central << endl;
    TFile *PU_Central_File = TFile::Open(PU_Central);
    
    TH1F* PU_intended = (TH1F *) PU_Central_File->Get("pileup"); // Data pileup distribution
    TH1F* PU_generated=NULL; // MC pileup distribution 

    TH1F* PU_weight=new TH1F("hPUweight","",100,-0.5,99.5);
    
    if (isMC) {
      // PU_generated = (TH1F*)file->Get("mainNtuplizer/pileuptrue"); // MC pileup distribution
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
    } // is MC
    
    // event categorizer
    //    EventCategory eventCategoryInst(1);   

    // Lepton scale factors
    LeptonEfficiencySF lepEff;

    //####################################################################################################################
    //###########################################           MVAHandler         ###########################################
    //####################################################################################################################
    //construct MVA out put file name
    TString mvaout = TString ( runProcess.getParameter<std::string>("outdir") ) + "/mva_" + outFileUrl + ".root";
    MVAHandler myMVAHandler_;
    myMVAHandler_.initTree(mvaout);

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

    // loop on all the events
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

        std::vector<TString> tags(1,"all");
        //genWeight
        float genWeight = 1.0;
        if (isMC) {
            if(ev.genWeight<0) { genWeight = -1.0; }
        }
        //systematical weight
        float weight = 1.0; //xsecWeight;
        if(isMC) {
	  weight *= genWeight;
	  weight *= xsecWeight; 
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

        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

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

	float lep_threshold(25.);
	float eta_threshold=2.5;
	
	for (auto &ilep : leps) {
	  if ( ilep.pt()<3. ) continue;

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

	  if (hasExtraLepton) {
	    nExtraLeptons++;
	    extraLeptons.push_back(ilep);
	  }
	} // leptons

	// sort(goodLeptons.begin(), goodLeptons.end(), ptsort());
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
	    //    default   :
	    // continue;
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
	  //   if(evcat!=fType) continue;

            // if(evcat==EE   && !(hasEEtrigger||hasEtrigger) ) continue;
            // if(evcat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
            // if(evcat==EMU  && !hasEMtrigger ) continue;

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
	// for(size_t ich=0; ich<tag_cat.size(); ich++){
	//   tags.push_back( tag_cat[ich] );
	// }

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

	// Dphi(lep, MET) ?
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
	// MET>25 GeV 
	bool passMet25(metP4.pt()>25);
	if (!passMet25) continue;
	mon.fillHisto("eventflow","all",3,weight); // MEt cut
	mon.fillHisto("eventflow","bdt",3,weight);
	//-------------------------------------------------------------------

	// mtW >50 GeV
	bool passMt(sqrt(tMass)>50. && sqrt(tMass)<250.);
	if (!passMt) continue;
	mon.fillHisto("eventflow","all",4,weight); // MT cut
	mon.fillHisto("eventflow","bdt",4,weight); 
        //
        //JET AND BTAGGING ANALYSIS
        //

        //###########################################################
        //  AK4 jets ,
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
	  if(fabs(corrJets[ijet].eta())>2.5) continue;
  
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

	  if (isMC) {
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
	  if (runCR) {
	    // To start with, simply use all AK4 jets in the pseudo b-jet collection:
	    CSVLoosebJets.push_back(corrJets[ijet]);
	  } else {
	    if (hasCSVtag) {
	      CSVLoosebJets.push_back(corrJets[ijet]);
	    }
	  }

	  //} // b-jet loop
	} // jet loop
    

	//--------------------------------------------------------------------------
	// AK4 jets:
	sort(GoodIdJets.begin(), GoodIdJets.end(), ptsort());
	// Fill Histograms with AK4,AK4 + CVS, AK8 + db basics:
	mon.fillHisto("njets_raw","nj", GoodIdJets.size(),weight);


	//--------------------------------------------------------------------------
	// AK4 + CSV jets:
	sort(CSVLoosebJets.begin(), CSVLoosebJets.end(), ptsort());
	mon.fillHisto("nbjets_raw","nb", CSVLoosebJets.size(),weight);

	if (runCR) {
	  // CR: Fake b-jets are sorted in b-tag discriminator rather than pt
	  sort(CSVLoosebJets.begin(), CSVLoosebJets.end(), btagsort());
	}

	//--------------------------------------------------------------------------
	// dphi(jet,MET)
	mon.fillHisto("dphijmet","raw",mindphijmet,weight);


	//###########################################################
	// The AK8 fat jets configuration
	//###########################################################

	// AK8 + double-b tagger fat-jet collection
	PhysicsObjectFatJetCollection DBfatJets; // collection of AK8 fat jets

	int ifjet(0);
	for (auto & ijet : fatJets ) {

	  if(ijet.pt()<jet_threshold_) continue;
	  if(fabs(ijet.eta())>2.5) continue;

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
  
	  if (ijet.softdropM>=7. && ijet.softdropM<=40.) {
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
	//  Configure cleaned AK4 and AK4 + CSVloose Jets : 
	//       Require DR separation with AK8 subjets
	//###########################################################

	PhysicsObjectJetCollection cleanedCSVLoosebJets;

	// AK4 + CSV jets
	for (auto & ib : CSVLoosebJets) {

	  bool hasOverlap(false);

	  float dRmin(999.);
	  // Loop over AK8 jets and find subjets
	  for (auto & ifb : DBfatJets) {
	    for (auto & it : ifb.subjets) { // subjets loop
	      //   hasOverlap = (deltaR(ib, it)<0.4);
	      float dR=deltaR(ib, it);
	      if (dR<dRmin) dRmin=dR;

	      hasOverlap=(dRmin<0.4);
	      
	      if ( verbose ) {
		if (hasOverlap) {

		  printf("\n Found overlap of b-jet : pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f \n with subjet in AK8: pt=%6.1f, eta=%7.3f, phi=%7.3f, mass=%7.3f \n",
			 ib.pt(),
			 ib.eta(),
			 ib.Phi(),
			 ib.M(),

			 it.pt(),
			 it.eta(),
			 it.Phi(),
			 it.M()
			 );
		  
		} // hasOverlap
	      } //verbose
      
	    } // subjets loop
	  } // AK8 loop

	  if (dRmin>0.4) {
	    cleanedCSVLoosebJets.push_back(ib);
	  }
	} // CSV b-jet loop


	//--------------------------------------------------------------------------
	// Cross-cleaned AK4 CSV b-jets:
	//sort(cleanedGoodIdJets.begin(), cleanedGoodIdJets.end(), ptsort());
	sort(cleanedCSVLoosebJets.begin(), cleanedCSVLoosebJets.end(), ptsort());

	//mon.fillHisto("njets_raw","cleaned", cleanedGoodIdJets.size(),weight);
	mon.fillHisto("nbjets_raw","cleaned", cleanedCSVLoosebJets.size(),weight);


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
  
	  hasOverlap=(dRmin_csv<0.4);
	  if (!hasOverlap) {// continue;
	  // Fill final soft-bs from SVs
	    SVs.push_back(isv);
	  }

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

 
	//--------------------------------------------------------------------------
	// First , set all b-jets (x-cleaned) in one vector<LorentzVector>
	vector<LorentzVector> GoodIdbJets;

	for (auto & i : CSVLoosebJets) 
	  {
	    GoodIdbJets.push_back(i);
	  } // AK4 + CSV
	if (!runCR) 
	  { // disable soft b's in the CR for the moment
	    for (auto & i : SVs) {
	      GoodIdbJets.push_back(i);
	    } // soft-b from SV
	  }
	mon.fillHisto("nbjets_2D","cat_raw",GoodIdJets.size(),GoodIdbJets.size(),weight);
	//mon.fillHisto("nbjets_2D","cat_cleaned_raw",cleanedGoodIdJets.size(),GoodIdbJets.size(),weight);
	mon.fillHisto("nbjets_raw","merged",GoodIdbJets.size(),weight);

	
	//--------------------------------------------------------------------------
	vector<LorentzVector> cleanedGoodIdbJets;

	for (auto & ib : GoodIdbJets) {
	  float dRmin(999.);
	  // Loop over AK8 jets and find subjets
	  for (auto & ifb : DBfatJets) {
	    for (auto & it : ifb.subjets) { // subjets loop
	      //   hasOverlap = (deltaR(ib, it)<0.4);
	      float dR=deltaR(ib, it);
	      if (dR<dRmin) dRmin=dR;
      
	    } // subjets loop
	  } // AK8 loop
  
	  if (dRmin>0.4) {
	    cleanedGoodIdbJets.push_back(ib);
	  }
	}

        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

	// LorentzVector wsum=metP4+goodLeptons[0].second;
	// // mtW
	// double tMass = 2.*goodLeptons[0].second.pt()*metP4.pt()*(1.-TMath::Cos(deltaPhi(goodLeptons[0].second.phi(),metP4.phi())));

	// mon.fillHisto("pfmet","raw",metP4.pt(),weight);
	// mon.fillHisto("mtw","raw",sqrt(tMass),weight);
	// mon.fillHisto("ptw","raw",wsum.pt(),weight);

	// //-------------------------------------------------------------------
	//  // MET>25 GeV 
	// bool passMet25(metP4.pt()>25);
	// if (!passMet25) continue;
	// mon.fillHisto("eventflow","all",3,weight); // MEt cut
	// //-------------------------------------------------------------------

	// // mtW >50 GeV
	// bool passMt(sqrt(tMass)>50. && sqrt(tMass)<200.);
	// if (!passMt) continue;
	// mon.fillHisto("eventflow","all",4,weight); // MT cut

	//-------------------------------------------------------------------
	mon.fillHisto("nbjets_2D","cat2_raw",GoodIdbJets.size(),DBfatJets.size(),weight);
	mon.fillHisto("nbjets_2D","cat3_raw",cleanedGoodIdbJets.size(),DBfatJets.size(),weight);
	mon.fillHisto("nbjets_2D","cats_raw",CSVLoosebJets.size(),SVs.size(),weight);

	//-------------------------------------------------------------------
	// At least 2 jets and 2 b-jets
	if (GoodIdJets.size()<2 || CSVLoosebJets.size()<2) continue;
	mon.fillHisto("eventflow","all",5,weight); 
	mon.fillHisto("eventflow","bdt",5,weight);
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// AK4 jets pt
	is=0;
	for (auto & jet : GoodIdJets) {
	   mon.fillHisto("jet_pt_raw", "jet"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "jet"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}
	//-------------------------------------------------------------------
	// AK4 + CSV jets
	is=0;
	for (auto & jet : CSVLoosebJets) {
	     mon.fillHisto("jet_pt_raw", b_tagging_name+htag[is], jet.pt(),weight);
	     mon.fillHisto("jet_eta_raw", b_tagging_name+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}
	//-------------------------------------------------------------------
	// x-cleaned AK4 + CSV jets
	is=0;
	for (auto & jet : cleanedCSVLoosebJets) {
	   mon.fillHisto("jet_pt_raw", "cleaned"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "cleaned"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}
	//-------------------------------------------------------------------
	// Merged CSV jets + SV
	is=0;
	for (auto & jet : GoodIdbJets) {
	   mon.fillHisto("jet_pt_raw", "merged"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "merged"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}


        //##############################################
        //########  Main Event Selection        ########
        //##############################################

        //-------------------------------------------------------------------
        // At least 3 b-tags
	if (runCR) 
	  {
	    if (CSVLoosebJets[1].btag0<CSVLooseWP) continue;
	    if (CSVLoosebJets.size()>=3) 
	      {
		if (CSVLoosebJets[2].btag0>CSVLooseWP) continue;    
	      }
	  }
	if (GoodIdbJets.size()<3) continue;
	mon.fillHisto("eventflow","all",6,weight); 
	
        //-------------------------------------------------------------------
        // if (GoodIdbJets.size()==1 && DBfatJets.size()==0) continue; // only allow =1b cat. if a fat-jet is present (in 3b cat)
        // if (GoodIdbJets.size()==2 && DBfatJets.size()==0) continue; // only allow =2b cat. if a fat-jet is present (in 4b cat)

        //----------------------------------------------------------------------------------------------------------//
        // Event categories according to (n-j, m-b, k-fat) jet multiplicities [nj>=2, (nb==1 + kf=1), nb>=2, kf>=0 ]
        //----------------------------------------------------------------------------------------------------------//

	// Here define all variables 
        LorentzVector allHadronic;
        //std::pair <int,LorentzVector> pairHadronic;

	// HT from all CSV + soft b's
	float ht(0.);

	if (GoodIdbJets.size()==3) // 3b category
        {// 3b cat.
            tags.push_back("SR_3b");
	    // Hadronic vector sum:
	    for (auto & thisb : GoodIdbJets) 
	      {
		allHadronic+=thisb;
	      }
	    // Hadronic scalar sum (HT):
	    for (auto & thisb : GoodIdbJets) 
	      {
		ht+=thisb.pt();
	      }
	    //mon.fillHisto("eventflow","all",6,weight); 
        } 
        else if (GoodIdbJets.size()>=4) // 4b category
        {// 4b cat.
	  tags.push_back("SR_geq4b"); 

	  if (GoodIdbJets.size()==4) 
	    { 
	      tags.push_back("SR_4b"); 
	    }
	  else
	    {
	      tags.push_back("SR_geq5b");
	    }
	  // Hadronic vector sum:
	  int countb(0);
	  for (auto & thisb : GoodIdbJets) 
	    {
	      allHadronic+=thisb;
	      countb++; if (countb>3) break;
	    }
	  // Hadronic scalar sum (HT):
	  for (auto & thisb : GoodIdbJets) 
	    {
	      ht+=thisb.pt();
	    }
	  
	  mon.fillHisto("eventflow","all",7,weight);
        } 
        else 
        {
            tags.push_back("UNKNOWN");
            printf("\n Unknown category, please check \n");
        }
 
        //-----------------------------------------------------------
        // Control plots
        // ----------------------------------------------------------

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
        // // HT from all CSV + soft b's
        mon.fillHisto("ht",tags,ht,weight);
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
        mon.fillHisto("MVABDT", tags, mvaBDT, weight);

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
	
        float mvaweight = 1.0;
        genWeight > 0 ? mvaweight = puWeight : mvaweight = -puWeight; // absorb the negative sign 
        if ( GoodIdbJets.size() >= 3 )
        {
          myMVAHandler_.getEntry
          (
	   GoodIdbJets.size() == 3, GoodIdbJets.size() >= 4, // 3b cat, 4b cat
	   wsum.pt(), //W only, w pt
            allHadronic.mass(), allHadronic.pt(), dRave_, dm, ht, //Higgs only, higgs mass, higgs pt, bbdr average, bb dm min, sum pt from all bs
	   dphi_Wh, //W and H, dr 
	   mvaweight //note, since weight is not the weight we want, we store all others except xSec weight
	   );
          myMVAHandler_.fillTree();
        }
	
        //##############################################################################
        //### HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //##############################################################################

        //##############################################
        // recompute MET/MT if JES/JER was varied
        //##############################################
        //LorentzVector vMET = variedMET[ivar>8 ? 0 : ivar];
        //PhysicsObjectJetCollection &vJets = ( ivar<=4 ? variedJets[ivar] : variedJets[0] );
	
    } // loop on all events END

    printf("\n");
    file->Close();

    //write MVA files
    myMVAHandler_.writeTree();

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
