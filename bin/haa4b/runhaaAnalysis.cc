#include <iostream>

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/bsmhiggs_fwk/interface/PatUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"
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

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
const float CSVLooseWP = 0.5426;  // Updated to 80X Moriond17 Loose
const float CSVMediumWP = 0.800;
const float CSVTightWP = 0.935;

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/Hbbtagging#8_0_X
// https://indico.cern.ch/event/543002/contributions/2205058/attachments/1295600/1932387/cv-doubleb-tagging-btv-approval.pdf (definition of WPs: slide 16)
const float DBLooseWP = 0.300;
const float DBMediumWP = 0.600;
const float DBTightWP = 0.900;

// Physics objects offline thresholds
const float lep_threshold_=25.;
const float mu_threshold_=25.;
const float ele_threshold_=30.;
const float jet_threshold_=20.;


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
   
    TString dtag=runProcess.getParameter<std::string>("tag");
    TString suffix=runProcess.getParameter<std::string>("suffix");

    bool verbose = runProcess.getParameter<bool>("verbose");

    bool usemetNoHF = runProcess.getParameter<bool>("usemetNoHF");
    
    TString url = runProcess.getParameter<std::string>("input");
    TString outFileUrl( dtag ); //gSystem->BaseName(url));
    //    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");
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
    if(url.Contains("SingleMuon"))  fType=MUMU;
    if(url.Contains("SingleElectron")) fType=EE;
    bool isSingleMuPD(!isMC && url.Contains("SingleMuon"));
    bool isDoubleMuPD(!isMC && url.Contains("DoubleMuon"));
    bool isSingleElePD(!isMC && url.Contains("SingleElectron"));
    bool isDoubleElePD(!isMC && url.Contains("DoubleEG"));

    bool isMC_ZZ2L2Nu  = isMC && ( string(url.Data()).find("TeV_ZZTo2L2Nu")  != string::npos);
    bool isMC_ZZTo4L   = isMC && ( string(url.Data()).find("TeV_ZZTo4L")  != string::npos);
    bool isMC_ZZTo2L2Q = isMC && ( string(url.Data()).find("TeV_ZZTo2L2Q")  != string::npos);

    bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZamcatnloFXFX")  != string::npos
                              || string(url.Data()).find("MC13TeV_WZpowheg")  != string::npos );

    bool isMC_VVV = isMC && ( string(url.Data()).find("MC13TeV_WZZ")  != string::npos
                              || string(url.Data()).find("MC13TeV_WWZ")  != string::npos
                              || string(url.Data()).find("MC13TeV_ZZZ")  != string::npos );

    bool isMCBkg_runPDFQCDscale = (isMC_ZZ2L2Nu || isMC_ZZTo4L || isMC_ZZTo2L2Q ||
                                   isMC_WZ || isMC_VVV);

    bool isMC_ttbar = isMC && (string(url.Data()).find("TeV_TT")  != string::npos);
    bool isMC_stop  = isMC && (string(url.Data()).find("TeV_SingleT")  != string::npos);

    bool isMC_Wh = isMC && (string(url.Data()).find("Wh")  != string::npos); 
    bool isMC_Zh = isMC && (string(url.Data()).find("Zh")  != string::npos); 
    bool isMC_VBF = isMC && (string(url.Data()).find("VBF")  != string::npos); 

    bool isSignal = (isMC_Wh || isMC_Zh || isMC_VBF );
    if (isSignal) printf("Signal url = %s\n",url.Data());

    //b-tagging: beff and leff must be derived from the MC sample using the discriminator vs flavor
    //the scale factors are taken as average numbers from the pT dependent curves see:
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
    BTagSFUtil btsfutil;
    float beff(0.68), sfb(0.99), sfbunc(0.015);
    float leff(0.13), sfl(1.05), sflunc(0.12);

    // setup calibration readers 80X
    BTagCalibration btagCalib("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/weights/CSVv2_Moriond17_B_H.csv");
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
        varNames.push_back("_jerup"); 	//1
        varNames.push_back("_jerdown"); //2
        varNames.push_back("_jesup"); 	//3
        varNames.push_back("_jesdown"); //4
        varNames.push_back("_umetup"); 	//5
        varNames.push_back("_umetdown");//6
        varNames.push_back("_lesup"); 	//7
        varNames.push_back("_lesdown"); //8
        varNames.push_back("_puup"); 	//9
        varNames.push_back("_pudown"); 	//10
        varNames.push_back("_btagup"); 	//11
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
//            varNames.push_back("_pdfacceptup");
//            varNames.push_back("_pdfacceptdown");
//            varNames.push_back("_qcdscaleacceptup");
//            varNames.push_back("_qcdscaleacceptdown");
        }
        if(isMC_ZZ2L2Nu) {
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
    //    JetCorrectionUncertainty jecUnc(uncFile.Data());

    //pdf info	
    
    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;

    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 8,0,8) );
    h->GetXaxis()->SetBinLabel(1,"Raw");
    h->GetXaxis()->SetBinLabel(2,"Trigger");
    h->GetXaxis()->SetBinLabel(3,"1 lepton");
    // h->GetXaxis()->SetBinLabel(3,"#Delta#it{#phi}(jet,E_{T}^{miss})>0.5");
    h->GetXaxis()->SetBinLabel(4,"E_{T}^{miss}>25");
    h->GetXaxis()->SetBinLabel(5,"M_{T}^{W}>0");
    h->GetXaxis()->SetBinLabel(6,">=(2-jets,2b-jets)");
    h->GetXaxis()->SetBinLabel(7,">=3b-tags");
    h->GetXaxis()->SetBinLabel(8,">=4b-tags");
 
     
     // GEN level kinematics
    mon.addHistogram( new TH1F( "higgsMass",";m_{h} [GeV];Events",200,0.,600.) );
    mon.addHistogram( new TH1F( "higgsPt",";p_{T}^{h} [GeV];Events",300,0,1500));
    mon.addHistogram( new TH1F( "higgsEta",";#eta (h);Evenets",100,-5,5) );
    mon.addHistogram( new TH1F( "a1mass",";m_{a1} [GeV];Events",800,0.,200.) );
    mon.addHistogram( new TH1F( "a2mass",";m_{a2} [GeV];Events",800,0.,200.) );
    mon.addHistogram( new TH1F( "a1pt",";p_{T}^{a1};Events",500,0,1000));
    mon.addHistogram( new TH1F( "a2pt",";p_{T}^{a2};Events",500,0,1000));
    mon.addHistogram( new TH1F( "aabalance",";p_{T}^{a1}/p_{T}^{a2};Events",200,0.,10.));
    mon.addHistogram( new TH1F( "a1DR",";#Delta R1(b,#bar{b});Events",100,0.,5.));
    mon.addHistogram( new TH1F( "a2DR",";#Delta R2(b,#bar{b});Events",100,0.,5.)); 
    //   mon.addHistogram( new TH1F( "dphiaa",";#Delta #phi(a1,a2);Events",40, 0, 4) ); 
    mon.addHistogram( new TH1F( "aaDR",";#Delta R(a_{1},a_{2});Events",100,0.,5.)); 

    mon.addHistogram( new TH1F( "b1pt",";p_{T}^{b1};Events",500,0,1000));   
    mon.addHistogram( new TH1F( "b2pt",";p_{T}^{b2};Events",500,0,1000)); 
    mon.addHistogram( new TH1F( "b3pt",";p_{T}^{b3};Events",500,0,1000)); 
    mon.addHistogram( new TH1F( "b4pt",";p_{T}^{b4};Events",500,0,1000)); 
    
    /*

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    mon.addHistogram( new TH1F( "nvtx_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "zpt_raw",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "mt_raw",       ";#it{m}_{T} [GeV];Events", 100,0,2000) );
    double MTBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000,2000};
    const int nBinsMT = sizeof(MTBins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "mt2_raw",       ";#it{m}_{T} [GeV];Events", nBinsMT,MTBins) );
    mon.addHistogram( new TH1F( "zmass_raw",    ";#it{m}_{ll} [GeV];Events", 100,40,250) );

    mon.addHistogram( new TH2F( "ptlep1vs2_raw",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );
    mon.addHistogram( new TH1F( "trailep_pt_raw", ";Trailing lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "trailep_eta_raw",";Trailing lepton #eta^{l};Events", 52,-2.6,2.6) );
    */

     // RECO level
    mon.addHistogram( new TH1F( "dR_raw",";#Delta R(SV,b);Events",100,0.,5.));

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l} [GeV];Events", 200,0,500) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );
    mon.addHistogram( new TH1F( "jet_pt_raw", ";#it{p}_{T}^{b} [GeV];Events", 500,0,500) );
    mon.addHistogram( new TH1F( "jet_eta_raw",";#eta;Events", 100,-5,5) );

    mon.addHistogram( new TH1F( "b_discrim"," ;b discriminator;",200,0,1.) );
    mon.addHistogram( new TH1F( "db_discrim"," ;double-b discriminator;",100,-1.,1.) );
    mon.addHistogram( new TH1F( "sd_mass"," ;soft-drop Mass;",100,-5.,295.) );
    mon.addHistogram( new TH1F( "pruned_mass"," ;pruned Mass;",100,-5.,295.) );
    mon.addHistogram( new TH1F( "softb_ntrk"," ; SV Ntrks;",11,-0.5,11.5) );
    mon.addHistogram( new TH1F( "softb_dxy"," ; SV dxy;",1000,0.,10.) );
    mon.addHistogram( new TH1F( "softb_dxyz_signif"," ; SVSIP3D;",100,1.,20.) );
    mon.addHistogram( new TH1F( "softb_cos"," ; SV cos((PV,SV),p_{SV});",250,-1.,1.) );

    mon.addHistogram( new TH1F( "nvtx_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "pfmet",    ";E_{T}^{miss} [GeV];Events", 100,0,500) );
    mon.addHistogram( new TH1F( "ht",    ";H_{T} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "mtw",       ";#it{m}_{T}^{W} [GeV];Events", 200,0,1000) );
    mon.addHistogram( new TH1F( "ptw",       ";#it{p}_{T}^{W} [GeV];Events", 200,0,1000) );
    mon.addHistogram( new TH1F( "dphiWh", ";#Delta#it{#phi}(#it{W},h);Events", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dRave",";#Delta R(b,b)_{ave};Events",100,0.,5.));
    mon.addHistogram( new TH1F( "dphijmet", ";#Delta#it{#phi}(jet,E_{T}^{miss});Events", 20,0,TMath::Pi()) );
    
    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nleptons_raw", ";Lepton multiplicity;Events", 4,0,4) );
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
    if(file==0) {
      return -1;
      printf("file is 0");
    }
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
    /*
    float cnorm=1.0;
    if(isMC) {
        //TH1F* cutflowH = (TH1F *) file->Get("mainAnalyzer/llvv/nevents");
        //if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        TH1F* posH = (TH1F *) file->Get("mainAnalyzer/llvv/n_posevents");
        TH1F* negH = (TH1F *) file->Get("mainAnalyzer/llvv/n_negevents");
        if(posH && negH) cnorm = posH->GetBinContent(1) - negH->GetBinContent(1);
        if(rescaleFactor>0) cnorm /= rescaleFactor;
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);
    */
    // muon trigger efficiency SF
    //Electron ID RECO SF

    // event categorizer
    //    EventCategory eventCategoryInst(1);   //jet(0,1,>=2) binning

    // Lepton scale factors
    LeptonEfficiencySF lepEff;

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

    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
	  printf("."); fflush(stdout);
        }

	if ( verbose ) printf("\n\n Event info %3d: \n",iev);

	
        //##############################################   EVENT LOOP STARTS   ##############################################
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
        if(isMC && ev.genWeight<0) genWeight = -1.0;

        //systematical weight
        float weight = 1.0;
        if(isMC) weight *= genWeight;

	//pileup re-weighting
	//	if(isMC) weight *= ev.puWeight;
	
        //only take up and down from pileup effect
        double TotalWeight_plus = 1.0;
        double TotalWeight_minus = 1.0;

        if(isMC) mon.fillHisto("pileup", tags, ev.ngenTruepu, 1.0);
	/*
        if(isMC) {
            weight 		*= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Central);
            TotalWeight_plus 	*= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Up);
            TotalWeight_minus 	*= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Down);
        }

        Hcutflow->Fill(1,genWeight);
        Hcutflow->Fill(2,weight);
        Hcutflow->Fill(3,weight*TotalWeight_minus);
        Hcutflow->Fill(4,weight*TotalWeight_plus);
	*/

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
        //#####################      Gen Selection       ##########################
        //#########################################################################
	
	bool iswithinAcceptance(true);
	PhysicsObjectCollection genbs;
	
	// Fill a PhysicsObjectCollection with GEN particles (from hard process)
	if (isSignal) {


	  if ( verbose ) {
	    PhysicsObjectCollection &particles = phys.genparticles;
	    
	    // dump the MC particle list
	    printf("\n\n Gen particles:\n" );
	    int ipar(0);
	    for (auto & it : particles) {

	        printf("  %3d : ID=%6d, m=%5.1f, momID=%6d , momIndx=%6d : pt=%6.1f, eta=%7.3f, phi=%7.3f\n",
		       ipar,
		       it.id,
		       it.mass(),
		       it.momid,
		       it.momidx,
		       it.pt(),
		       it.eta(),
		       it.phi()
		       ) ;
		ipar++;
	    }
	   
	  } // verbose

	  //#### Find a1 and a2 positions in mcparticle list
	  //#### for Wh, Zh : 6 and 7 ####################### 
	  //#### for VBF: 5 and 6 ###########################
	  
	  PhysicsObjectCollection &higgses = phys.genHiggs;
	  PhysicsObjectCollection &partons = phys.genpartons;
	  PhysicsObjectCollection &leptons = phys.genleptons;
	  
	  const int pos1(6);
	  const int pos2(7);

	  LorentzVector higgs, a1, a2;
       
	  int as(0);
	  for (auto & par : higgses ) {
	    // 125 GeV Higgs:
	    if (par.id==25) {
	      mon.fillHisto("higgsMass","raw",par.mass(),weight);
	      mon.fillHisto("higgsPt","raw",par.pt(),weight);
	      mon.fillHisto("higgsEta","raw",par.eta(), weight);
	    }

	    // a's from h decay
	    if (par.id==36) {
	      if (as==0) { a1 += par; as++; 
	      } else {	a2 += par; }
	    }
	    
	  } // higgses

	  mon.fillHisto("a1mass","raw",a1.mass(),weight);
	  mon.fillHisto("a2mass","raw",a2.mass(),weight);
	  double raw_aaDR=deltaR(a1,a2);
	  mon.fillHisto("aaDR","raw",raw_aaDR,weight);  
	  
	  // h->aa->4b (all)
	
	  PhysicsObjectCollection genbFromA1;
	  PhysicsObjectCollection genbFromA2;
	  
	  for (auto & par : partons) {
	    if (par.momid!=36) continue; // look only b's from a->bb decay

	    // Apply acceptance eta cuts:
	    if (fabs(par.eta())>2.5) continue;
	    
	    if ( abs(par.id)==5 ) { 
	      genbs.push_back(par);
	      higgs += par; 
	    } 

	    if (abs(par.id)==5) {
	      if (par.momidx==pos1) {
		genbFromA1.push_back(par);
	      } else if (par.momidx==pos2) {
		genbFromA2.push_back(par);
	      }
	      
	    } // b-partons
	  } // partons

	  //sort gen b's in pt
	  sort(genbs.begin(), genbs.end(), ptsort());

	  // // Filtering on leptons from the W decay
	  // bool hasWlep(false);
	  // for (auto & it : leptons) {
	  //   hasWlep = (abs(it.momid)==24);
	  // }
	  // if (!hasWlep) continue;
	  
	  if (genbFromA1.size()==2 && genbFromA2.size()==2) {
	     // Higgs pT
	    mon.fillHisto("higgsMass","gen",higgs.mass(),weight);
	    mon.fillHisto("higgsPt","gen",higgs.pt(),weight);

	    mon.fillHisto("a1mass","gen",(genbFromA1[0]+genbFromA1[1]).mass(),weight);
	    mon.fillHisto("a2mass","gen",(genbFromA2[0]+genbFromA2[1]).mass(),weight);

	    double aaDR=deltaR( (genbFromA1[0]+genbFromA1[1]),(genbFromA2[0]+genbFromA2[1]) );
	    mon.fillHisto("aaDR","gen",aaDR,weight);

	    double balance = (genbFromA1[0]+genbFromA1[1]).pt()/(genbFromA2[0]+genbFromA2[1]).pt();
	    mon.fillHisto("aabalance","gen",balance,weight);
	    
	    double pt1,pt2;
	    double dR1,dR2;

	    pt1=(genbFromA1[0]+genbFromA1[1]).pt();
	    pt2=(genbFromA2[0]+genbFromA2[1]).pt();

	    dR1=deltaR(genbFromA1[0],genbFromA1[1]);
	    dR2=deltaR(genbFromA2[0],genbFromA2[1]);

	    mon.fillHisto("a1pt","gen",pt1,weight);
	    mon.fillHisto("a2pt","gen",pt2,weight);
	    mon.fillHisto("a1DR","gen",dR1,weight);
	    mon.fillHisto("a2DR","gen",dR2,weight);

	     // pT of the 4bs
	    if (genbs.size()>=4) {
	      mon.fillHisto("b1pt","gen",genbs[0].pt(),weight);
	      mon.fillHisto("b2pt","gen",genbs[1].pt(),weight);
	      mon.fillHisto("b3pt","gen",genbs[2].pt(),weight);
	      mon.fillHisto("b4pt","gen",genbs[3].pt(),weight);
	    }
	    
	  } else {
	    if ( verbose ) {
	      printf("Not all 4-bs found in h -> aa decays\n");
	    }
	    iswithinAcceptance=false;
	    //	    continue; // Look at signal events in detector acceptance
	  }

	  // Test Jet collection / jet multiplicities
	  PhysicsObjectJetCollection &recoJets = phys.jets; //
	  PhysicsObjectJetCollection rawJets, rawJets_eta2p5, rawJets_eta2p5_noLep;
	  
	  for(size_t ijet=0; ijet<recoJets.size(); ijet++) {
	    
	    if(recoJets[ijet].pt()<jet_threshold_) continue;
	    rawJets.push_back(recoJets[ijet]);
	    if (fabs(recoJets[ijet].eta())>2.5) continue;
	    rawJets_eta2p5.push_back(recoJets[ijet]);
	    
	    
	    // //check overlaps with gen leptons
	    bool hasOverlap(false);
	    for (auto & ilep : leptons) { 
	      double dR = deltaR( recoJets[ijet], ilep);
	      if (abs(ilep.id)==11) hasOverlap = (dR<0.2); // within 0.2 for electrons
	      if (abs(ilep.id)==13) hasOverlap = (dR<0.4); // within 0.4 for muons
	      if (hasOverlap) break;
	    }
	    if(hasOverlap) continue;
	    
	    rawJets_eta2p5_noLep.push_back(recoJets[ijet]);
	    
	  }
	  mon.fillHisto("njets_raw","nj_all",rawJets.size(),weight);
	  mon.fillHisto("njets_raw","nj_eta2p5_all",rawJets_eta2p5.size(),weight);
	  mon.fillHisto("njets_raw","nj_eta2p5_noLep_all",rawJets_eta2p5_noLep.size(),weight);

	  
	} // isSignal
	  
	
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
        // looping leptons (electrons + muons)
        int nGoodLeptons(0);
        std::vector<std::pair<int,LorentzVector> > goodLeptons;
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
	    //          if(lep.pt()<lep_threshold_) continue;
	    if(abs(lepid)==11 && lep.pt()<ele_threshold_) continue;
	    if(abs(lepid)==13 && lep.pt()<mu_threshold_) continue;
	    
            if(abs(lepid)==13 && fabs(lep.eta())> 2.4) continue;
            if(abs(lepid)==11 && fabs(lep.eta())> 2.5) continue;
            if(abs(lepid)==11 && fabs(lep.eta()) > 1.442 && fabs(lep.eta()) < 1.556) continue;

            bool hasTightIdandIso(true);
            if(abs(lepid)==13) { //muon
	      hasTightIdandIso &= (phys.leptons[ilep].passIdMu && phys.leptons[ilep].passIsoMu);
            } else if(abs(lepid)==11) { //electron
	      hasTightIdandIso &= (phys.leptons[ilep].passIdEl && phys.leptons[ilep].passIsoEl);
            } else continue;


            if(!hasTightIdandIso) continue;
            nGoodLeptons++;
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(lepid,lep);
            goodLeptons.push_back(goodlep);

        }

	// sort goodLeptons in pT
	//	sort(goodLeptons.begin(), goodLeptons.end(), ptsort());
	mon.fillHisto("nleptons_raw","all", goodLeptons.size(),weight);
	
	/*
        // ID + ISO scale factors 
        if(isMC) {
        }
	*/
	std::vector<TString> tag_cat;
	//	TString tag_cat;
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
	  if(evcat==E   && hasEtrigger ) hasTrigger=true;
	  if(evcat==MU && hasMtrigger ) hasTrigger=true;
	  if(evcat==EMU  && ( hasEtrigger || hasMtrigger)) hasTrigger=true;
	  if(!hasTrigger) continue;
        }
	
	//	tags.push_back(tag_cat); //add ee, mumu, emu category
	
        //prepare the tag's vectors for histo filling
	// for(size_t ich=0; ich<tag_cat.size(); ich++){
	//   tags.push_back( tag_cat[ich] );
	// }
	
        // pielup reweightiing
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_raw",tags, phys.nvtx,      weight);

        //
        //apply muon trigger efficiency scale factors
        //
	
	// Trigger
	mon.fillHisto("eventflow","all",1,weight);
	
	// Exactly 1 good lepton
	if(goodLeptons.size()!=1) continue; // at least 1 tight leptons
	mon.fillHisto("eventflow","all",2,weight);
	
	if (abs(goodLeptons[0].first==11)) {
	  mon.fillHisto("leadlep_pt_raw","e",goodLeptons[0].second.pt(),weight);
	  mon.fillHisto("leadlep_eta_raw","e",goodLeptons[0].second.eta(),weight);
	} else if (abs(goodLeptons[0].first==13)) {
	  mon.fillHisto("leadlep_pt_raw","mu",goodLeptons[0].second.pt(),weight);
	  mon.fillHisto("leadlep_eta_raw","mu",goodLeptons[0].second.eta(),weight);
	}
	
        //
        //JET AND BTAGGING ANALYSIS
        //

	//###########################################################
	//  AK4 jets ,
	// AK4 jets + CSVloose b-tagged configuration
	//###########################################################

        PhysicsObjectJetCollection GoodIdJets;
	PhysicsObjectJetCollection CSVLoosebJets;

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
	  
	  nCSVLtags += (corrJets[ijet].btag0>CSVLooseWP);
	  nCSVMtags += (corrJets[ijet].btag0>CSVMediumWP);
	  nCSVTtags += (corrJets[ijet].btag0>CSVTightWP);
	  
	  mon.fillHisto("b_discrim","csv",corrJets[ijet].btag0,weight);
	  if (corrJets[ijet].motherid == 36) mon.fillHisto("b_discrim","csv_true",corrJets[ijet].btag0,weight);
	  
	  bool hasCSVtag(corrJets[ijet].btag0>CSVLooseWP);
	  if (isMC) {
	    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
	    btsfutil.SetSeed(ev.event*10 + ijet*10000);
	    
	    if(abs(corrJets[ijet].flavid)==5) {
	      //  80X recommendation
	      btsfutil.modifyBTagsWithSF(hasCSVtag , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_B , corrJets[ijet].eta(), corrJets[ijet].pt()), beff);
	    } else if(abs(corrJets[ijet].flavid)==4) {
	      //  80X recommendation
	      btsfutil.modifyBTagsWithSF(hasCSVtag , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_C , corrJets[ijet].eta(), corrJets[ijet].pt()), beff);
	    } else {
	      //  80X recommendation
	      btsfutil.modifyBTagsWithSF(hasCSVtag , btagCal80X.eval_auto_bounds("central", BTagEntry::FLAV_UDSG , corrJets[ijet].eta(), corrJets[ijet].pt()), leff);
	    }
	  } // isMC
	  
	    // Fill b-jet vector:
	  if (hasCSVtag) {
	    CSVLoosebJets.push_back(corrJets[ijet]);
	  }
	  
	  //	} // b-jet loop
	} // jet loop
    

	//###########################################################
	// The AK8 fat jets configuration
	//###########################################################
	
	// AK8 + double-b tagger fat-jet collection
	PhysicsObjectFatJetCollection DBfatJets; // collection of AK8 fat jets
	
	int ifjet(0);
	for (auto & ijet : fatJets ) {
	
	  if(ijet.pt()<jet_threshold_) continue;
	  if(fabs(ijet.eta())>2.5) continue;

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
	  
	  //	  bool hasDBtag(ijet.btag0>DBMediumWP);
	  
	  mon.fillHisto("db_discrim","fjet",ijet.btag0,weight);
	  if (ijet.motherid == 36) mon.fillHisto("db_discrim","fjet_true",ijet.btag0,weight);
	  mon.fillHisto("nsubjets_raw","fjet",count_sbj,weight);
	  if (ijet.motherid == 36) mon.fillHisto("nsubjets_raw","fjet_true",count_sbj,weight);
	  mon.fillHisto("sd_mass","fjet",ijet.softdropM,weight);
	  if (ijet.motherid == 36) mon.fillHisto("sd_mass","fjet_true",ijet.softdropM,weight);
	  mon.fillHisto("pruned_mass","fjet",ijet.prunedM,weight);
	  if (ijet.motherid == 36) mon.fillHisto("pruned_mass","fjet_true",ijet.prunedM,weight);
	  
	  // double-b tagger + at least 1 subjet in AK8
	  //   if (hasDBtag && count_sbj>0) {
	  //   if (ijet.prunedM>10 && ijet.prunedM<60) {
	  DBfatJets.push_back(ijet);
	  //  }
	  //}
	  
	} // AK8 fatJets loop
	
		
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


	//###########################################################
	// Soft b-jets from SVs configuration
	//###########################################################

	// SVs collection
	PhysicsObjectSVCollection SVs;
	
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
	  
	  float dRmin(999.);
	  for (auto & it : GoodIdJets) {
	    double dR=deltaR(it, isv);
	    if (dR<dRmin) dRmin=dR;
	  }
	  mon.fillHisto("dR_raw","sv_jet",dRmin,weight);
	  
	  // plot minDR(SV,b)
	  float dRmin_csv(999.);
	  for (auto & it : CSVLoosebJets) {
	    double dR=deltaR(it, isv);
	    if (dR<dRmin_csv) dRmin_csv=dR;
	  }
	  mon.fillHisto("dR_raw","sv_b",dRmin_csv,weight);
	  
	  hasOverlap=(dRmin_csv<0.4);
	  if (hasOverlap) continue;
	  
	  // Fill final soft-bs from SVs
	  SVs.push_back(isv);
	  
	}

	
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	// some strings for tagging histograms:
	const char* astr[] = {"_b1","_b2","_b3","_b4"};
        std::vector<TString> htag(astr, astr+4);

	//--------------------------------------------------------------------------
	// AK4 jets:
	sort(GoodIdJets.begin(), GoodIdJets.end(), ptsort());
	
	// Fill Histograms with AK4,AK4 + CVS, AK8 + db basics:
	mon.fillHisto("njets_raw","nj", GoodIdJets.size(),weight);
	
	int is(0);
	for (auto & jet : GoodIdJets) {
	   mon.fillHisto("jet_pt_raw", "jet"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "jet"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}
	
	//--------------------------------------------------------------------------
	// AK4 + CSV jets:
	sort(CSVLoosebJets.begin(), CSVLoosebJets.end(), ptsort());
	
	mon.fillHisto("nbjets_raw","nb", CSVLoosebJets.size(),weight);

	is=0;
	for (auto & jet : CSVLoosebJets) {
	   mon.fillHisto("jet_pt_raw", "csv"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "csv"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}

	//--------------------------------------------------------------------------
	// AK8 + double-b jets
	sort(DBfatJets.begin(), DBfatJets.end(), ptsort());

	mon.fillHisto("nbjets_raw","nfatJet", DBfatJets.size(),weight);

	is=0;
	for (auto & jet : DBfatJets) {
	   mon.fillHisto("jet_pt_raw", "fat"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "fat"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}

	//--------------------------------------------------------------------------
	// Cross-cleaned AK4 CSV b-jets:
	//sort(cleanedGoodIdJets.begin(), cleanedGoodIdJets.end(), ptsort());
	sort(cleanedCSVLoosebJets.begin(), cleanedCSVLoosebJets.end(), ptsort());

	//	mon.fillHisto("njets_raw","cleaned", cleanedGoodIdJets.size(),weight);
	mon.fillHisto("nbjets_raw","cleaned", cleanedCSVLoosebJets.size(),weight);

	is=0;
	for (auto & jet : cleanedCSVLoosebJets) {
	   mon.fillHisto("jet_pt_raw", "cleaned"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "cleaned"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}
	
	//--------------------------------------------------------------------------
	// Soft-bs properties
	//--------------------------------------------------------------------------

	sort(SVs.begin(), SVs.end(), ptsort());

	mon.fillHisto("nbjets_raw","nb_soft",SVs.size(),weight);

	is=0;
	for (auto & isv : SVs) {
	   mon.fillHisto("jet_pt_raw", "softb"+htag[is], isv.pt(),weight);
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

	//--------------------------------------------------------------------------
	// First , set all b-jets (x-cleaned) in one vector<LorentzVector>
	vector<LorentzVector> GoodIdbJets;

	for (auto & i : CSVLoosebJets) {
	  GoodIdbJets.push_back(i);
	} // AK4 + CSV
	for (auto & i : SVs) {
	  GoodIdbJets.push_back(i);
	} // soft-b from SV

	mon.fillHisto("nbjets_2D","cat_raw",GoodIdJets.size(),GoodIdbJets.size(),weight);
	//	mon.fillHisto("nbjets_2D","cat_cleaned_raw",cleanedGoodIdJets.size(),GoodIdbJets.size(),weight);
	mon.fillHisto("nbjets_raw","merged",GoodIdbJets.size(),weight);
	
	is=0;
	for (auto & jet : GoodIdbJets) {
	   mon.fillHisto("jet_pt_raw", "merged"+htag[is], jet.pt(),weight);
	   mon.fillHisto("jet_eta_raw", "merged"+htag[is], jet.eta(),weight);
	   is++;
	   if (is>3) break; // plot only up to 4 b-jets ?
	}

	// Fill true b-jet multiplicity for signal eff
	mon.fillHisto("nbjets_raw","true",genbs.size(),weight);
	
        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

	LorentzVector wsum=metP4+goodLeptons[0].second;
	// mtW
	double tMass = 2.*goodLeptons[0].second.pt()*metP4.pt()*(1.-TMath::Cos(deltaPhi(goodLeptons[0].second.phi(),metP4.phi())));

	mon.fillHisto("pfmet","raw",metP4.pt(),weight);
	mon.fillHisto("mtw","raw",sqrt(tMass),weight);
	mon.fillHisto("ptw","raw",wsum.pt(),weight);


	 // MET>25 GeV 
	bool passMet25(metP4.pt()>25);
	if (!passMet25) continue;
	mon.fillHisto("eventflow","all",3,weight); // MEt cut
		
	// mtW >50 GeV
	bool passMt(sqrt(tMass)>50);
	//if (!passMt) continue;
	mon.fillHisto("eventflow","all",4,weight); // MT cut
	
	// At least 2 jets and 2 b-jets
	if (GoodIdJets.size()<2 || CSVLoosebJets.size()<2) continue;
	mon.fillHisto("eventflow","all",5,weight); 
	
	// At least 3 b-tags
	if (GoodIdbJets.size()<3) continue;
	mon.fillHisto("eventflow","all",6,weight); 

	// if (GoodIdbJets.size()==1 && DBfatJets.size()==0) continue; // only allow =1b cat. if a fat-jet is present (in 3b cat)
	// if (GoodIdbJets.size()==2 && DBfatJets.size()==0) continue; // only allow =2b cat. if a fat-jet is present (in 4b cat)
	
	//----------------------------------------------------------------------------------------------------------//
	// Event categories according to (n-j, m-b, k-fat) jet multiplicities [nj>=2, (nb==1 + kf=1), nb>=2, kf>=0 ]
	//----------------------------------------------------------------------------------------------------------//

	 LorentzVector allHadronic;
	 //	 std::pair <int,LorentzVector> pairHadronic;
 
	 if (GoodIdbJets.size()==3) {// 3b cat.
	    tags.push_back("3b");
	   for (auto & thisb : GoodIdbJets) {
	     allHadronic+=thisb;
	   }
	   //  mon.fillHisto("eventflow",tags,4,weight);
	 } else if (GoodIdbJets.size()>=4) {// 4b cat.
	   tags.push_back("4b");
	   int countb(0);
	   for (auto & thisb : GoodIdbJets) {
	     allHadronic+=thisb;
	     countb++;
	     if (countb>3) break;
	   }
	   mon.fillHisto("eventflow","all",7,weight);
	 } else {
	   tags.push_back("UNKNOWN");
	   printf("\n Unknown category, please check \n");
	 }
	 
	 //-----------------------------------------------------------
	 // Control plots
	 // ----------------------------------------------------------

	 // 3,4 b's pT
	 mon.fillHisto("nbjets_raw",tags,GoodIdbJets.size(),weight);
	 is=0;
	 for (auto & jet : GoodIdbJets) {
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
	 float ht(0.);
	 for (auto & thisb : GoodIdbJets) {
	   ht+=thisb.pt();
	 }
	 mon.fillHisto("ht",tags,ht,weight);
	 // MET
	 mon.fillHisto("pfmet",tags,metP4.pt(),weight);
	 // dphi(jet,MET)
	 mon.fillHisto("dphijmet",tags,mindphijmet,weight);
	 // pTW
	 //	 LorentzVector wsum=metP4+goodLeptons[0].second;
	 mon.fillHisto("ptw",tags,wsum.pt(),weight);
	 // // mtW 
	 // double tMass = pow(sqrt(pow(goodLeptons[0].second.pt(),2))+sqrt(pow(metP4.pt(),2)+pow(80.,2)),2);
	 // tMass-=pow(wsum.pt(),2);
	 mon.fillHisto("mtw",tags,sqrt(tMass),weight);
	 // Dphi(W,h) instead of DRmin(l,b)
	 double dphi_Wh=fabs(deltaPhi(allHadronic.phi(),wsum.phi()));
	 mon.fillHisto("dphiWh",tags,dphi_Wh,weight);
	 // DR(bb)_average
	 vector<float> dRs;
	 dRs.push_back(deltaR(GoodIdbJets[0],GoodIdbJets[1]));
	 dRs.push_back(deltaR(GoodIdbJets[0],GoodIdbJets[2]));
	 dRs.push_back(deltaR(GoodIdbJets[1],GoodIdbJets[2]));
	 if (GoodIdbJets.size()>=4) {
	   dRs.push_back(deltaR(GoodIdbJets[0],GoodIdbJets[3]));
	   dRs.push_back(deltaR(GoodIdbJets[1],GoodIdbJets[3]));
	   dRs.push_back(deltaR(GoodIdbJets[2],GoodIdbJets[3]));
	 }

	 float dRave_(0.);
	 for (auto & it : dRs) {
	   dRave_+=it;
	 }
	 dRave_/=dRs.size();
	 mon.fillHisto("dRave",tags,dRave_,weight);
	 
        //##############################################
        //########  Main Event Selection        ########
        //##############################################


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

    if(outTxtFile_final)fclose(outTxtFile_final);
}

