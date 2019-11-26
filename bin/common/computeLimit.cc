#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include "UserCode/bsmhiggs_fwk/interface/tdrstyle.h"
#include "UserCode/bsmhiggs_fwk/interface/JSONWrapper.h"
#include "UserCode/bsmhiggs_fwk/interface/RootUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/MacroUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/HxswgUtils.h"
#include "UserCode/bsmhiggs_fwk/interface/th1fmorph.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>
#include <regex>


using namespace std;

TString signalSufix="";

TString histo(""), histoVBF("");

int rebinVal = 1;
double MCRescale = 1.0;
double SignalRescale = 1.0;

double datadriven_qcd_Syst = 0.50;    

bool postfit=false;

double norm_top=0.89;
double enorm_3b_w=1.37;
double enorm_4b_w=1.29;
double munorm_3b_w=1.37;
double munorm_4b_w=1.44;

int mass;
bool shape = true;
TString postfix="";
TString systpostfix="";
bool runSystematics = true; 

bool runZh = false;
bool modeDD = false;
bool simfit = false;

TString vh_tag;

std::vector<TString> Channels;
std::vector<string> AnalysisBins;

double DDRescale = 1.0;
TString DYFile ="";
TString FREFile="";
string signalTag="";

bool BackExtrapol  = false;
bool subNRB        = false;
bool MCclosureTest = false;
bool scaleVBF      = false;

bool mergeWWandZZ = false;
bool skipWW = true;
bool skipGGH = false;
bool skipQQH = false;
bool subDY = false;
bool subWZ = false;
bool subFake = false;
bool blindData = false;
bool blindWithSignal = false; 

TString inFileUrl(""),jsonFile("");

double shapeMin =-9999;
double shapeMax = 9999;
double shapeMinVBF =-9999;
double shapeMaxVBF = 9999;
bool doInterf = false;
double minSignalYield = 0;
float statBinByBin = -1;
bool useLogy = true;
bool blindSR = false;
double lumi = -1;
int signalScale = 1;

bool dirtyFix1 = false;
bool dirtyFix2 = false;

std::vector<int> shapeBinToConsider;

std::vector<int> indexcutV;
std::vector<int> indexcutVL;
std::vector<int> indexcutVR;

std::map<string, int> indexcutM;
std::map<string, int> indexcutML;
std::map<string, int> indexcutMR;

std::vector<string> keywords;


int indexvbf = -1;
int massL=-1, massR=-1;

double dropBckgBelow=0.01; 

/*
 *Case Sensitive Implementation of startsWith()
 *It checks if the string 'mainStr' starts with given string 'toMatch'
 */
bool startsWith(std::string mainStr, std::string toMatch){
  // std::string::find returns 0 if toMatch is found at starting
  if(mainStr.find(toMatch) == 0)
    return true;
  else
    return false;
}

bool matchKeyword(JSONWrapper::Object& process, std::vector<string>& keywords){
  if(keywords.size()<=0)return true;
  if(process.isTag("keys")){
    std::vector<JSONWrapper::Object> dsetkeywords = process["keys"].daughters();
    for(size_t ikey=0; ikey<dsetkeywords.size(); ikey++){
      for(unsigned int i=0;i<keywords.size();i++){if(std::regex_match(dsetkeywords[ikey].toString(),std::regex(keywords[i])))return true;}
    }
  }else{
    return true;
  }
  return false;
}



void filterBinContent(TH1* histo){
  if(shapeBinToConsider.size()<=0)return;
  for(int i=0;i<=histo->GetNbinsX()+1;i++){
    bool toBeConsidered=false;  for(unsigned int j=0;j<shapeBinToConsider.size();j++){if(shapeBinToConsider[j]==i){toBeConsidered=true;break;}}
    if(!toBeConsidered){histo->SetBinContent(i,0); histo->SetBinError(i,0);}
  }
}


//wrapper for a projected shape for a given proc
class ShapeData_t
{
  public:
  	std::map<string, double> uncScale;
  	std::map<string, TH1*  > uncShape;
  	TH1* fit;

  	ShapeData_t(){
	  	fit=NULL;
  	}
  	~ShapeData_t(){}

  	TH1* histo(){
     	if(uncShape.find("")==uncShape.end())return NULL;
     	return uncShape[""];
  	}

  	void clearSyst(){
     	TH1* nominal = histo();
     	uncScale.clear();
     	uncShape.clear();
     	uncShape[""] = nominal;
  	}


  	void removeStatUnc(){
     	for(auto unc = uncShape.begin(); unc!= uncShape.end(); unc++){
        TString name = unc->first.c_str();
        if(name.Contains("stat") && (name.Contains("Up") || name.Contains("Down"))){
          uncShape.erase(unc);
          unc--;
        }
     	}
  	}

  	void makeStatUnc(string prefix="", string suffix="", string suffix2="", bool noBinByBin=false){
	  if(!histo() || histo()->Integral()<=0)return;
	  string delimiter = "_";
	  unsigned firstDelimiter = suffix.find(delimiter);
	  unsigned lastDelimiter = suffix.find_last_of(delimiter);
	  unsigned endPosOfFirstDelimiter = firstDelimiter + delimiter.length();
	  string channel_and_bin = suffix.substr(endPosOfFirstDelimiter, lastDelimiter-endPosOfFirstDelimiter);
	  
	  if(suffix.find("instrmet") != std::string::npos){
	    TH1* h = (TH1*) histo()->Clone("TMPFORSTAT");
	    int BIN=0;
	    std::vector<unsigned int > v_lowStatBin;
	    v_lowStatBin.clear();
	    TString InstrMET_gammaStats_Url(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/bsmhiggs_fwk/data/InstrMET_systematics/InstrMET_systematics_GAMMASTATS.root");
	    TFile* f_InstrMET_gammaStats = TFile::Open(InstrMET_gammaStats_Url);
	    TH1* h_InstrMET_Up_gammaStats = (TH1*)utils::root::GetObjectFromPath(f_InstrMET_gammaStats, (channel_and_bin+"_mt_InstrMET_absolute_shape_up").c_str() );
	    TH1* h_InstrMET_Down_gammaStats = (TH1*)utils::root::GetObjectFromPath(f_InstrMET_gammaStats, (channel_and_bin+"_mt_InstrMET_absolute_shape_down").c_str() );   			
	    for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){           
	      if( true /*h->GetBinContent(ibin)/h->Integral()>0.01*/){ //This condition is removed for the moment, we may put it back in the future. 
		
		char ibintxt[255]; sprintf(ibintxt, "_b%i", BIN);BIN++;
		TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU"+ibintxt);//  statU->Reset();
		TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD"+ibintxt);//  statD->Reset();           
		
		statU->SetBinContent(ibin, std::max(0.0, h_InstrMET_Up_gammaStats->GetBinContent(ibin))>0 ? h_InstrMET_Up_gammaStats->GetBinContent(ibin) : 0.115);
		statD->SetBinContent(ibin, std::max(0.0, h_InstrMET_Down_gammaStats->GetBinContent(ibin))); 
		uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Up"  ] = statU;
		uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Down"] = statD;
		
	      }
	      else{
		v_lowStatBin.push_back(ibin);
	      }
	    }
	    
	    TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
	    TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");
	    
	    if(v_lowStatBin.size()>0){
	      for(unsigned int j=0; j < v_lowStatBin.size(); j++){
		statU->SetBinContent(v_lowStatBin[j], std::max(0.0, h_InstrMET_Up_gammaStats->GetBinContent(v_lowStatBin[j])));   
		statD->SetBinContent(v_lowStatBin[j], std::max(0.0, h_InstrMET_Down_gammaStats->GetBinContent(v_lowStatBin[j])));   
	      }
	      uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
	      uncShape[prefix+"stat"+suffix+"Down"] = statD;
	    }	
	    
	    //f_InstrMET_gammaStats->Close();
	    delete h; //all done with this copy
	    
	  }
	  else{
	    
	    TH1* h = (TH1*) histo()->Clone("TMPFORSTAT");
	    
	    //bin by bin stat uncertainty
	    if(statBinByBin>0 && shape==true && !noBinByBin){
	      int BIN=0;
	      for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){           
		//		if(h->Integral()<=0.01) continue;
		//		if(h->GetBinContent(ibin)/h->Integral()<0.01 ) continue; // printf("Found bin with cont/Integral < 0.01\n");
		//		if(h->GetBinContent(ibin)<=0.) printf("Found bin with 0 contnent\n");
		if( !(h->GetBinContent(ibin)<=0 && h->GetBinError(ibin)>0) &&  (h->GetBinContent(ibin)<=0 || h->GetBinContent(ibin)/h->Integral()<0.01 || h->GetBinError(ibin)/h->GetBinContent(ibin)<statBinByBin))continue;
		//		if(h->GetBinContent(ibin)<=0.)continue;
		char ibintxt[255]; sprintf(ibintxt, "_b%i", BIN);BIN++;
		TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU"+ibintxt);//  statU->Reset();
		TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD"+ibintxt);//  statD->Reset();           
		if(h->GetBinContent(ibin)>0){
		  statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), h->GetBinContent(ibin) + h->GetBinError(ibin))));   statU->SetBinError(ibin, 0.0);
		  statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), h->GetBinContent(ibin) - h->GetBinError(ibin))));   statD->SetBinError(ibin, 0.0);
		  // statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.0, h->GetBinContent(ibin) + h->GetBinError(ibin))));   statU->SetBinError(ibin, 0);
		  // statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.0, h->GetBinContent(ibin) - h->GetBinError(ibin))));   statD->SetBinError(ibin, 0);
		}else{
		  statU->SetBinContent(ibin,              statU->GetBinContent(ibin) + statU->GetBinError(ibin));
		  statD->SetBinContent(ibin,std::max(0.0, statD->GetBinContent(ibin) - statD->GetBinError(ibin)));
		}
		uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Up"  ] = statU;
		uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Down"] = statD;
		/*h->SetBinContent(ibin, 0);*/  h->SetBinError(ibin, 0);  //remove this bin from shape variation for the other ones
		//printf("%s --> %f - %f - %f\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Up").c_str(), statD->Integral(), h->GetBinContent(ibin), statU->Integral() );
	      }
	    }
	    
	    //after this line, all bins with large stat uncertainty have been considered separately
	    //so now it remains to consider all the other bins for which we assume a total correlation bin by bin
	    if(h->Integral()<=0)return; //all non empty bins have already bin variated
	    TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
	    TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");
	    for(int ibin=1; ibin<=statU->GetXaxis()->GetNbins(); ibin++){
	      if(h->GetBinContent(ibin)>0){
		statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statU->GetBinContent(ibin) + statU->GetBinError(ibin))));
		statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statD->GetBinContent(ibin) - statD->GetBinError(ibin))));
	      }else{
		statU->SetBinContent(ibin,              statU->GetBinContent(ibin) + statU->GetBinError(ibin));
		statD->SetBinContent(ibin,std::min(0.0, statD->GetBinContent(ibin) - statD->GetBinError(ibin)));
	      }
	    }
	    uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
	    uncShape[prefix+"stat"+suffix+"Down"] = statD;
	    
	    delete h; //all done with this copy
	  }
	}
  
  	double getScaleUncertainty(){
     	double Total=0;
     	for(std::map<string, double>::iterator unc=uncScale.begin();unc!=uncScale.end();unc++){
        if(unc->second<0)continue;
        Total+=pow(unc->second,2);
     	}     
     	return Total>0?sqrt(Total):-1;
  	}

        double getIntegratedShapeUncertainty(string name, string upORdown){
     	double Total=0;
     	//this = ch->second.shapes[histoName.Data()]
     	for(std::map<string, TH1*>::iterator var = uncShape.begin(); var!=uncShape.end(); var++){
       	TString systName = var->first.c_str();
       	if(var->first=="")continue; //Skip Nominal shape
       	if(!systName.Contains(upORdown))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"

       	TH1* hvar = (TH1*)(var->second->Clone((name+var->first).c_str()));

	double varYield = hvar->Integral();
	TH1* h = NULL;
	if (this->histo()!=NULL) h = (TH1*)(this->histo()->Clone((name+"Nominal").c_str()));
	double yield = 0.; 
	if (h!=NULL) yield = h->Integral();
       	Total+=pow(varYield-yield,2); //the total shape unc is the sqrt of the quadratical sum of the difference between the nominal and the variated yields.
     	}     
     	return Total>0?sqrt(Total):-1;
  	}

        double getBinShapeUncertainty(string name, int bin, string upORdown){
     	double Total=0;
     	//this = ch->second.shapes[histoName.Data()]
     	for(std::map<string, TH1*>::iterator var = uncShape.begin(); var!=uncShape.end(); var++){
       	TString systName = var->first.c_str();

       	if(var->first=="")continue; //Skip Nominal shape
       	if(!systName.Contains(upORdown))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"
       	TH1* hvar = (TH1*)(var->second->Clone((name+var->first).c_str()));

	double varYield = hvar->GetBinContent(bin);
	TH1* h = (TH1*)(this->histo()->Clone((name+"Nominal").c_str()));
	double yield = h->GetBinContent(bin);
       	Total+=pow(varYield-yield,2); //the total shape unc is the sqrt of the quadratical sum of the difference between the nominal and the variated yields.
     	}     
     	return Total>0?sqrt(Total):-1;
  	}


  	void rescaleScaleUncertainties(double StartIntegral, double EndIntegral){
     	for(std::map<string, double>::iterator unc=uncScale.begin();unc!=uncScale.end();unc++){
        printf("%E/%E = %E = %E/%E\n", unc->second, StartIntegral, unc->second/StartIntegral, EndIntegral * (unc->second/StartIntegral), EndIntegral); 
        if(StartIntegral!=0){unc->second = EndIntegral * (unc->second/StartIntegral);}else{unc->second = unc->second * EndIntegral;}
     	}     
  	}


};

class ChannelInfo_t
{
  public:
    string bin;
    string channel;

  std::map<string, ShapeData_t> shapes;

  ChannelInfo_t(){}
  ~ChannelInfo_t(){}
};

class ProcessInfo_t
{
  public:
  	bool isData;
  	bool isBckg;
  	bool isSign;
  	double xsec;
  	double br;
  	double mass;
  	string shortName;
  	std::map<string, ChannelInfo_t> channels;
  	JSONWrapper::Object jsonObj;

  	ProcessInfo_t(){xsec=0;}
  	~ProcessInfo_t(){}
};

class AllInfo_t
{
  public:
  
    std::map<string, ProcessInfo_t> procs;
    std::vector<string> sorted_procs;

		AllInfo_t(){};
		~AllInfo_t(){};

    // reorder the procs to get the backgrounds; total bckg, signal, data 
    void sortProc();

    // Sum up all background processes and add this as a total process
    void addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst = false);

    // Sum up all background processes and add this as a total process
    void addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst = false);

    // Sum up all background processes and add this as a total process
    void computeTotalBackground();

    // Replace the Data process by TotalBackground
    void blind();

    // Print the Yield table
    void getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc=NULL);

    // Dump efficiencies
    void getEffFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName);

    // drop background process that have a negligible yield
    void dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold);

    // drop control channels
    void dropCtrlChannels(std::vector<TString>& selCh);

   // Subtract nonQCD MC processes from A,C,D regions in data
    void doBackgroundSubtraction(FILE* pFile, std::vector<TString>& selCh,TString mainHisto);

    // Make a summary plot
    void showShape(std::vector<TString>& selCh , TString histoName, TString SaveName);

    // Make a summary plot of the uncertainties
    void showUncertainty(std::vector<TString>& selCh , TString histoName, TString SaveName);

    // Turn to cut&count (rebin all histo to 1 bin only)
    void turnToCC(string histoName);

    // Make a summary plot
    void saveHistoForLimit(string histoName, TFile* fout);

    // Add hardcoded uncertainties 
    void addHardCodedUncertainties(string histoName);

    // produce the datacards 
    void buildDataCards(string histoName, TString url);

    // Load histograms from root file and json to memory
    void getShapeFromFile(TFile* inF, std::vector<string> channelsAndShapes, int cutBin, JSONWrapper::Object &Root,  double minCut=0, double maxCut=9999, bool onlyData=false);

    // Rebin histograms to make sure that high mt/met region have no empty bins
    void rebinMainHisto(string histoName);

    //Merge bins together
    void mergeBins(std::vector<string>& binsToMerge, string NewName);

    // Handle empty bins
    void HandleEmptyBins(string histoName);

};


void printHelp();
void printHelp()
{
  printf("Options\n");
  printf("--in        --> input file with from plotter\n");
  printf("--json      --> json file with the sample descriptor\n");
  printf("--histoVBF  --> name of histogram to be used for VBF\n");
  printf("--histo     --> name of histogram to be used\n");
  printf("--shapeMin  --> left cut to apply on the shape histogram\n");
  printf("--shapeMax  --> right cut to apply on the shape histogram\n");
  printf("--shapeMinVBF  --> left cut to apply on the shape histogram for Vbf bin\n");
  printf("--shapeMaxVBF  --> right cut to apply on the shape histogram for Vbf bin\n");
  printf("--indexvbf  --> index of selection to be used for the vbf bin (if unspecified same as --index)\n");
  printf("--index     --> index of selection to be used (Xbin in histogram to be used); different comma separated values can be given for each analysis bin\n");
  printf("--indexL    --> index of selection to be used (Xbin in histogram to be used) used for interpolation;  different comma separated values can be given for each analysis bin\n");
  printf("--indexR    --> index of selection to be used (Xbin in histogram to be used) used for interpolation;  different comma separated values can be given for each analysis bin\n");
  printf("--m         --> higgs mass to be considered\n");
  printf("--mL        --> higgs mass on the left  of the mass to be considered (used for interpollation\n");
  printf("--mR        --> higgs mass on the right of the mass to be considered (used for interpollation\n");
  printf("--syst      --> use this flag if you want to run systematics, default is no systematics\n");
  printf("--shape     --> use this flag if you want to run shapeBased analysis, default is cut&count\n");
  printf("--subNRB    --> use this flag if you want to subtract non-resonant-backgounds similarly to what was done in 2011 (will also remove H->WW)\n");
  printf("--subNRB12  --> use this flag if you want to subtract non-resonant-backgounds using a new technique that keep H->WW\n");
  printf("--subDY     --> histogram that contains the Z+Jets background estimated from Gamma+Jets)\n");
  printf("--subWZ     --> use this flag if you want to subtract WZ background by the 3rd lepton SB)\n");
  printf("--DDRescale --> factor to be used in order to multiply/rescale datadriven estimations\n");
  printf("--closure   --> use this flag if you want to perform a MC closure test (use only MC simulation)\n");
  printf("--bins      --> list of bins to be used (they must be comma separated without space)\n");
  printf("--HWW       --> use this flag to consider HWW signal)\n");
  printf("--skipGGH   --> use this flag to skip GGH signal)\n");
  printf("--skipQQH   --> use this flag to skip GGH signal)\n");
  printf("--blind     --> use this flag to replace observed data by total predicted background)\n");
  printf("--blindWithSignal --> use this flag to replace observed data by total predicted background+signal)\n");
  printf("--postfix    --> use this to specify a postfix that will be added to the process names)\n");
  printf("--systpostfix    --> use this to specify a syst postfix that will be added to the process names)\n");
  printf("--MCRescale    --> use this to rescale the cross-section of all MC processes by a given factor)\n");
  printf("--postfit  ---> use this to apply postfit Normalization values for W and Top processes in the Signal + Control regions \n");
  printf("--signalRescale    --> use this to rescale signal cross-section by a given factor)\n");
  printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
  printf("--minSignalYield   --> use this to specify the minimum Signal yield you want in each channel)\n");
  printf("--signalSufix --> use this flag to specify a suffix string that should be added to the signal 'histo' histogram\n");
  printf("--signalTag   --> use this flag to specify a tag that should be present in signal sample name\n");
  printf("--signalScale   --> use this flag to specify a Scale applied on signal\n");
  printf("--rebin         --> rebin the histogram\n");
  printf("--statBinByBin --> make bin by bin statistical uncertainty\n");
  printf("--inclusive  --> merge bins to make the analysis inclusive\n");
  printf("--dropBckgBelow --> drop all background processes that contributes for less than a threshold to the total background yields\n");
  printf("--scaleVBF    --> scale VBF signal by ggH/VBF\n");
  printf("--key        --> provide a key for sample filtering in the json\n");  
  printf("--noLogy        --> use this flag to make y-axis linear scale\n");  
}

//
int main(int argc, char* argv[])
{
  setTDRStyle();
  gStyle->SetPadTopMargin   (0.06);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadRightMargin (0.16);
  gStyle->SetPadLeftMargin  (0.14);
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.45);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
  gStyle->SetOptStat(0);  
  gStyle->SetOptFit(0);

  //get input arguments
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")          !=string::npos) { printHelp(); return -1;} 
    else if(arg.find("--minSignalYield") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&minSignalYield ); i++; printf("minSignalYield = %f\n", minSignalYield);}
    else if(arg.find("--scaleVBF") !=string::npos) { scaleVBF=true; printf("scaleVBF = True\n");}
    else if(arg.find("--subNRB")   !=string::npos) { subNRB=true; skipWW=true; printf("subNRB = True\n");}
    else if(arg.find("--subDY")    !=string::npos) { subDY=true; DYFile=argv[i+1];  i++; printf("Z+Jets will be replaced by %s\n",DYFile.Data());}
    else if(arg.find("--subFake")  !=string::npos) { subFake=true; printf("Fake lepton QCD procs will be replaced by DD\n");}
    else if(arg.find("--subWZ")    !=string::npos) { subWZ=true; printf("WZ will be estimated from 3rd lepton SB\n");}
    else if(arg.find("--DDRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&DDRescale); i++;}
    else if(arg.find("--MCRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&MCRescale); i++;}
    else if(arg.find("--signalRescale")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&SignalRescale); i++;}
    else if(arg.find("--HWW")      !=string::npos) { skipWW=false; printf("HWW = True\n");}
    else if(arg.find("--skipGGH")  !=string::npos) { skipGGH=true; printf("skipGGH = True\n");}
    else if(arg.find("--skipQQH")  !=string::npos) { skipQQH=true; printf("skipQQH = True\n");}
    else if(arg.find("--blindWithSignal")    !=string::npos) { blindData=true; blindWithSignal=true; printf("blindData = True; blindWithSignal = True\n");}
    else if(arg.find("--blind")    !=string::npos) { blindData=true; printf("blindData = True\n");}
    else if(arg.find("--closure")  !=string::npos) { MCclosureTest=true; printf("MCclosureTest = True\n");}
    else if(arg.find("--shapeBinToConsider")    !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); shapeBinToConsider.push_back(C);  pch = strtok(NULL,",");} i++; printf("Only the following histo bins will be considered: "); for(unsigned int i=0;i<shapeBinToConsider.size();i++)printf(" %i ", shapeBinToConsider[i]);printf("\n");}
    else if(arg.find("--shapeMinVBF") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMinVBF); i++; printf("Min cut on shape for VBF = %f\n", shapeMinVBF);}
    else if(arg.find("--shapeMaxVBF") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMaxVBF); i++; printf("Max cut on shape for VBF = %f\n", shapeMaxVBF);}
    else if(arg.find("--shapeMin") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMin); i++; printf("Min cut on shape = %f\n", shapeMin);}
    else if(arg.find("--shapeMax") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&shapeMax); i++; printf("Max cut on shape = %f\n", shapeMax);}
    else if(arg.find("--interf")    !=string::npos) { doInterf=true; printf("doInterf = True\n");}
    else if(arg.find("--indexvbf") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&indexvbf); i++; printf("indexVBF = %i\n", indexvbf);}
    else if(arg.find("--index" )   !=string::npos && i+1<argc)   { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); indexcutV .push_back(C);  pch = strtok(NULL,",");} i++; printf("index  = "); for(unsigned int i=0;i<indexcutV .size();i++)printf(" %i ", indexcutV [i]);printf("\n");}
    else if(arg.find("--indexL")    !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); indexcutVL.push_back(C);  pch = strtok(NULL,",");} i++; printf("indexL = "); for(unsigned int i=0;i<indexcutVL.size();i++)printf(" %i ", indexcutVL[i]);printf("\n");}
    else if(arg.find("--indexR")    !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");while (pch!=NULL){int C;  sscanf(pch,"%i",&C); indexcutVR.push_back(C);  pch = strtok(NULL,",");} i++; printf("indexR = "); for(unsigned int i=0;i<indexcutVR.size();i++)printf(" %i ", indexcutVR[i]);printf("\n");}
    else if(arg.find("--in")       !=string::npos && i+1<argc)  { inFileUrl = argv[i+1];  i++;  printf("in = %s\n", inFileUrl.Data());  }
    else if(arg.find("--json")     !=string::npos && i+1<argc)  { jsonFile  = argv[i+1];  i++;  printf("json = %s\n", jsonFile.Data()); }
    else if(arg.find("--histoVBF") !=string::npos && i+1<argc)  { histoVBF  = argv[i+1];  i++;  printf("histoVBF = %s\n", histoVBF.Data()); }
    else if(arg.find("--histo")    !=string::npos && i+1<argc)  { histo     = argv[i+1];  i++;  printf("histo = %s\n", histo.Data()); }
    else if(arg.find("--mL")       !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&massL ); i++; printf("massL = %i\n", massL);}
    else if(arg.find("--mR")       !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&massR ); i++; printf("massR = %i\n", massR);}
    else if(arg.find("--m")        !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&mass ); i++; printf("mass = %i\n", mass);}
    else if(arg.find("--bins")     !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("bins are : ");while (pch!=NULL){printf(" %s ",pch); AnalysisBins.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--channels") !=string::npos && i+1<argc)  { char* pch = strtok(argv[i+1],",");printf("channels are : ");while (pch!=NULL){printf(" %s ",pch); Channels.push_back(pch);  pch = strtok(NULL,",");}printf("\n"); i++; }
    else if(arg.find("--postfix")   !=string::npos && i+1<argc)  { postfix = argv[i+1]; systpostfix = argv[i+1]; i++;  printf("postfix '%s' will be used\n", postfix.Data());  }
    else if(arg.find("--systpostfix")   !=string::npos && i+1<argc)  { systpostfix = argv[i+1];  i++;  printf("systpostfix '%s' will be used\n", systpostfix.Data());  }
    else if(arg.find("--shape")  !=string::npos) { shape=true; printf("shapeBased = True\n");}   
    else if(arg.find("--syst")  !=string::npos) { runSystematics=true; printf("syst = True\n");}     
    else if(arg.find("--simfit")  !=string::npos) { simfit=true; printf("simfit = True\n");}    
    //    else if(arg.find("--postfit")  !=string::npos) { postfit=true; printf("postfit = True\n");}    
    else if(arg.find("--dirtyFix2")    !=string::npos) { dirtyFix2=true; printf("dirtyFix2 = True\n");}
    else if(arg.find("--dirtyFix1")    !=string::npos) { dirtyFix1=true; printf("dirtyFix1 = True\n");}
    else if(arg.find("--signalSufix") !=string::npos) { signalSufix = argv[i+1]; i++; printf("signalSufix '%s' will be used\n", signalSufix.Data()); }
    else if(arg.find("--signalTag") !=string::npos) { signalTag = argv[i+1]; i++; printf("signalTag '%s' will be used\n", signalTag.c_str()); }
    else if(arg.find("--signalScale") !=string::npos) { sscanf(argv[i+1],"%d",&signalScale); i++; printf("signalScale = %d\n", signalScale);}
    else if(arg.find("--rebin")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&rebinVal); i++; printf("rebin = %i\n", rebinVal);}
    else if(arg.find("--BackExtrapol")    !=string::npos) { BackExtrapol=true; printf("BackExtrapol = True\n");}
    else if(arg.find("--statBinByBin")    !=string::npos) { sscanf(argv[i+1],"%f",&statBinByBin); i++; printf("statBinByBin = %f\n", statBinByBin);}
    else if(arg.find("--dropBckgBelow")   !=string::npos) { sscanf(argv[i+1],"%lf",&dropBckgBelow); i++; printf("dropBckgBelow = %f\n", dropBckgBelow);}
    else if(arg.find("--key"          )   !=string::npos && i+1<argc){ keywords.push_back(argv[i+1]); printf("Only samples matching this (regex) expression '%s' are processed\n", argv[i+1]); i++;  }
    else if(arg.find("--noLogy")    !=string::npos) { useLogy=false; printf("useLogy = False\n");}
    else if(arg.find("--SRblind")    !=string::npos) { blindSR=true; printf("blindSR = True\n");}
    else if(arg.find("--lumi") !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&lumi); i++; printf("Lumi = %lf\n", lumi);}
    if(arg.find("--runZh") !=string::npos) { runZh=true; printf("runZh = True\n");}
    if(arg.find("--modeDD") !=string::npos) { modeDD=true; printf("modeDD = True\n");} 
    if(arg.find("--postfit")  !=string::npos) { postfit=true; printf("postfit = True\n");}   
  }
  if(jsonFile.IsNull()) { printf("No Json file provided\nrun with '--help' for more details\n"); return -1; }
  if(inFileUrl.IsNull()){ printf("No Inputfile provided\nrun with '--help' for more details\n"); return -1; }
  if(histo.IsNull())    { printf("No Histogram provided\nrun with '--help' for more details\n"); return -1; }
  if(mass==-1)          { printf("No massPoint provided\nrun with '--help' for more details\n"); return -1; }
  if(indexcutV.size()<=0){printf("INDEX CUT SIZE IS NULL\n"); printHelp(); return -1; }
  if(AnalysisBins.size()==0)AnalysisBins.push_back("all");
  if(Channels.size()==0){ 
    if (modeDD) {
      Channels.push_back("e_A_SR"); Channels.push_back("mu_A_SR");
      Channels.push_back("e_B_SR"); Channels.push_back("mu_B_SR");  
      Channels.push_back("e_C_SR"); Channels.push_back("mu_C_SR");  
      Channels.push_back("e_D_SR"); Channels.push_back("mu_D_SR");  
    } else {
      if(runZh){
	Channels.push_back("ee_A_SR"); Channels.push_back("mumu_A_SR");  
      } else {
	Channels.push_back("e_A_SR"); Channels.push_back("mu_A_SR");
      }
    }
    if(simfit){
      if (modeDD) {
	Channels.push_back("e_A_CR");Channels.push_back("mu_A_CR"); // Top CR
	Channels.push_back("e_B_CR");Channels.push_back("mu_B_CR"); // Top CR 
	Channels.push_back("e_C_CR");Channels.push_back("mu_C_CR"); // Top CR 
	Channels.push_back("e_D_CR");Channels.push_back("mu_D_CR"); // Top CR

	Channels.push_back("e_A_CR5j");Channels.push_back("mu_A_CR5j"); // tt+bb CR
	Channels.push_back("e_B_CR5j");Channels.push_back("mu_B_CR5j"); // tt+bb CR 
	Channels.push_back("e_C_CR5j");Channels.push_back("mu_C_CR5j"); // tt+bb CR 
	Channels.push_back("e_D_CR5j");Channels.push_back("mu_D_CR5j"); // tt+bb CR 
      } else {
	if(runZh){ // Zh
	  Channels.push_back("ee_A_CR");Channels.push_back("mumu_A_CR"); // DY CR
	  Channels.push_back("emu_A_SR");Channels.push_back("emu_A_CR"); // Top CR     
	}else{ // Wh
	  Channels.push_back("e_A_CR");Channels.push_back("mu_A_CR"); // Top/W CR
	  Channels.push_back("e_A_CR5j");Channels.push_back("mu_A_CR5j"); // tt+bb CR   
	}
      }
    }
  }

  vh_tag = runZh ? "_zh" : "_wh";

  //make sure that the index vector are well filled
  if(indexcutVL.size()==0) indexcutVL.push_back(indexcutV [0]);
  if(indexcutVR.size()==0) indexcutVR.push_back(indexcutV [0]);
  while(indexcutV .size()<AnalysisBins.size()){indexcutV .push_back(indexcutV [0]);}
  while(indexcutVL.size()<AnalysisBins.size()){indexcutVL.push_back(indexcutVL[0]);}
  while(indexcutVR.size()<AnalysisBins.size()){indexcutVR.push_back(indexcutVR[0]);}
  if(indexvbf>=0){for(unsigned int i=0;i<AnalysisBins.size();i++){if(AnalysisBins[i].find("vbf")!=string::npos){indexcutV[i]=indexvbf; indexcutVL[i]=indexvbf; indexcutVR[i]=indexvbf;} }}



  //handle merged bins
  std::vector<std::vector<string> > binsToMerge;
  for(unsigned int b=0;b<AnalysisBins.size();b++){
    if(AnalysisBins[b].find('+')!=std::string::npos){
      std::vector<string> subBins;
      char* pch = strtok(&AnalysisBins[b][0],"+"); 
      while (pch!=NULL){
        indexcutV.push_back(indexcutV[b]);
        indexcutVL.push_back(indexcutVL[b]);
        indexcutVR.push_back(indexcutVR[b]);
        AnalysisBins.push_back(pch);
        subBins.push_back(pch);
        pch = strtok(NULL,"+");
      }
      binsToMerge.push_back(subBins);
      AnalysisBins.erase(AnalysisBins.begin()+b);
      indexcutV .erase(indexcutV .begin()+b);
      indexcutVL.erase(indexcutVL.begin()+b);
      indexcutVR.erase(indexcutVR.begin()+b);
      b--;
    }
  }


  //fill the index map
  for(unsigned int i=0;i<AnalysisBins.size();i++){indexcutM[AnalysisBins[i]] = indexcutV[i]; indexcutML[AnalysisBins[i]] = indexcutVL[i]; indexcutMR[AnalysisBins[i]] = indexcutVR[i];}


  ///////////////////////////////////////////////


  //init the json wrapper
  JSONWrapper::Object Root(jsonFile.Data(), true);


  //init globalVariables
  TString massStr(""); if(mass>0)massStr += mass;
  std::vector<TString> allCh,allProcs;

  std::vector<TString> ch;
  if (modeDD) {
    ch.push_back("e_A_SR"); ch.push_back("mu_A_SR");  
    ch.push_back("e_B_SR"); ch.push_back("mu_B_SR");   
    ch.push_back("e_C_SR"); ch.push_back("mu_C_SR");   
    ch.push_back("e_D_SR"); ch.push_back("mu_D_SR");   
    if (simfit) {
      ch.push_back("e_A_CR"); ch.push_back("mu_A_CR");  
      ch.push_back("e_B_CR"); ch.push_back("mu_B_CR"); 
      ch.push_back("e_C_CR"); ch.push_back("mu_C_CR"); 
      ch.push_back("e_D_CR"); ch.push_back("mu_D_CR");

      ch.push_back("e_A_CR5j"); ch.push_back("mu_A_CR5j");  
      ch.push_back("e_B_CR5j"); ch.push_back("mu_B_CR5j"); 
      ch.push_back("e_C_CR5j"); ch.push_back("mu_C_CR5j"); 
      ch.push_back("e_D_CR5j"); ch.push_back("mu_D_CR5j"); 
    }
  } else {
    if(runZh){ // Zh
      ch.push_back("ee_A_SR"); ch.push_back("mumu_A_SR");     
      if (simfit) {
	ch.push_back("ee_A_CR"); ch.push_back("mumu_A_CR");   
	ch.push_back("emu_A_SR"); ch.push_back("emu_A_CR");       
      }

    } else { // Wh
      ch.push_back("e_A_SR"); ch.push_back("mu_A_SR");
      if (simfit) { 
	ch.push_back("e_A_CR"); ch.push_back("mu_A_CR");
	ch.push_back("e_A_CR5j"); ch.push_back("mu_A_CR5j"); 
      }
    }

  }
  //TString ch[]={"SR"}; //"mumu","ee","emu"};
  const size_t nch=ch.size(); //sizeof(ch)/sizeof(TString);
  std::vector<TString> sh;
  sh.push_back(histo);
  if(subNRB)sh.push_back(histo+"_NRBctrl");
  if(subWZ)sh.push_back(histo+"_3rdLepton");

  AllInfo_t allInfo;

  //open input file
  TFile* inF = TFile::Open(inFileUrl);
  if( !inF || inF->IsZombie() ){ printf("Invalid file name : %s\n", inFileUrl.Data());}
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  //LOAD shapes
  const size_t nsh=sh.size();
  for(size_t b=0; b<AnalysisBins.size(); b++){
    std::vector<string> channelsAndShapes;
    for(size_t i=0; i<nch; i++){
      for(size_t j=0; j<nsh; j++){
	channelsAndShapes.push_back((ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
	printf("Adding shape %s\n",(ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
      }
    }
    double cutMin=shapeMin; double cutMax=shapeMax;
    allInfo.getShapeFromFile(inF, channelsAndShapes, indexcutM[AnalysisBins[b]], Root, cutMin, cutMax   );     
  }


  inF->Close();
  printf("Loading all shapes... Done\n");

  allInfo.computeTotalBackground();
  if(MCclosureTest)allInfo.blind();

  FILE* pFile;

  //define vector for search
  std::vector<TString>& selCh = Channels;

  if(modeDD) {
    pFile = fopen("datadriven_qcd.tex","w");
    if(subFake)allInfo.doBackgroundSubtraction(pFile,selCh,histo);
  }

  //replace data by total MC background
  if(blindData)allInfo.blind();

  //extrapolate backgrounds toward higher BDT region to make sure that there is no empty bins
  if(shape && BackExtrapol)allInfo.rebinMainHisto(histo.Data());

  //drop backgrounds with rate<1%
  allInfo.dropSmallBckgProc(selCh, histo.Data(), dropBckgBelow);

  //drop control channels
  allInfo.dropCtrlChannels(selCh);

  //merge bins  
  for(unsigned int B=0;B<binsToMerge.size();B++){
    std::string NewBinName = string("["); binsToMerge[B][0];  for(unsigned int b=1;b<binsToMerge[B].size();b++){NewBinName += "+"+binsToMerge[B][b];} NewBinName+="]";
    allInfo.mergeBins(binsToMerge[B],NewBinName);
  }


  //turn to CC analysis eventually
  if(!shape)allInfo.turnToCC(histo.Data());

  allInfo.HandleEmptyBins(histo.Data()); //needed for negative bin content --> May happens due to NLO interference for instance

  // Blind data in Signal Regions only
  //  if(blindData)allInfo.blind();

  //print event yields from the histo shapes
  pFile = fopen(runZh?"Yields_zh.tex":"Yields_wh.tex","w");  FILE* pFileInc = fopen(runZh?"YieldsInc_zh.tex":"YieldsInc_wh.tex","w");
  allInfo.getYieldsFromShape(pFile, selCh, histo.Data(), pFileInc);
  fclose(pFile); fclose(pFileInc);

  //print signal efficiency
  pFile = fopen(runZh?"Efficiency_zh.tex":"Efficiency_wh.tex","w");
  allInfo.getEffFromShape(pFile, selCh, histo.Data());
  fclose(pFile);

  //add by hand the hard coded uncertainties
  allInfo.addHardCodedUncertainties(histo.Data());

  //produce a plot
  allInfo.showShape(selCh,histo,"plot"); //this produce the final global shape

  //produce a plot
  if(runSystematics) allInfo.showUncertainty(selCh,histo,"plot"); //this produces all the plots with the syst

  //prepare the output
  string limitFile=("haa4b_"+massStr+systpostfix+vh_tag+".root").Data();
  TFile *fout=TFile::Open(limitFile.c_str(),"recreate");

  allInfo.saveHistoForLimit(histo.Data(), fout);

  allInfo.buildDataCards(histo.Data(), limitFile);

  //all done
  fout->Close();
}

//
// reorder the procs to get the backgrounds; total bckg, signal, data 
//
void AllInfo_t::sortProc(){
  std::vector<string>bckg_procs;
  std::vector<string>sign_procs;
  bool isTotal=false, isData=false;
  for(unsigned int p=0;p<sorted_procs.size();p++){
    string procName = sorted_procs[p];
    std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
    if(it==procs.end())continue;
    if(it->first=="total"){isTotal=true; continue;}
    if(it->first=="data"){isData=true; continue;}
    if(it->second.isSign)sign_procs.push_back(procName);
    if(it->second.isBckg)bckg_procs.push_back(procName);
  }
  sorted_procs.clear();
  sorted_procs.insert(sorted_procs.end(), bckg_procs.begin(), bckg_procs.end());
  if(isTotal)sorted_procs.push_back("total");
  if(isData)sorted_procs.push_back("data");
  sorted_procs.insert(sorted_procs.end(), sign_procs.begin(), sign_procs.end());
}

//
// Sum up all shapes from one src channel to a total shapes in the dest channel
//
void AllInfo_t::addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst){
  std::map<string, ShapeData_t>& shapesInfoDest = dest.shapes;
  std::map<string, ShapeData_t>& shapesInfoSrc  = src.shapes;

  if(!computeSyst){
    for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
      if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();
      
      //Loop on all shape systematics (including also the central value shape)
      for(std::map<string, TH1*>::iterator uncS = sh->second.uncShape.begin();uncS!= sh->second.uncShape.end();uncS++){
	if(uncS->first!="") continue; //We only take nominal shapes
	if(shapesInfoDest[sh->first].uncShape.find(uncS->first)==shapesInfoDest[sh->first].uncShape.end()){
	  shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) uncS->second->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
	}else{
	  shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second);
	}
      }
      
      //take care of the scale uncertainty 
      for(std::map<string, double>::iterator unc = sh->second.uncScale.begin();unc!= sh->second.uncScale.end();unc++){
	if(shapesInfoDest[sh->first].uncScale.find(unc->first)==shapesInfoDest[sh->first].uncScale.end()){
	  shapesInfoDest[sh->first].uncScale[unc->first] = unc->second;
	}else{
	  shapesInfoDest[sh->first].uncScale[unc->first] = sqrt( pow(shapesInfoDest[sh->first].uncScale[unc->first],2) + pow(unc->second,2) );
	}
      }
    }  
  }
  else {

    for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
      if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();
 		//Loop on all shape systematics (including also the central value shape)
      for(std::map<string, TH1*>::iterator uncS = sh->second.uncShape.begin();uncS!= sh->second.uncShape.end();uncS++){
	if(uncS->first=="") continue; //We only take systematic (i.e non-nominal) shapes
	if(shapesInfoSrc[sh->first].uncShape.find("")==shapesInfoSrc[sh->first].uncShape.end()) continue;
	//1. Copy the nominal shape
	shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) shapesInfoDest[sh->first].uncShape[""]->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
	//2. we remove the nominal value of the process we are running on
	shapesInfoDest[sh->first].uncShape[uncS->first]->Add(shapesInfoSrc[sh->first].uncShape[""], -1);
	//3. and add the variation up/down
	shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second); 
      }
    }
  }

}

//
// Sum up all background processes and add this as a total process
//
void AllInfo_t::addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst){
  dest.xsec = src.xsec*src.br;
  for(std::map<string, ChannelInfo_t>::iterator ch = src.channels.begin(); ch!=src.channels.end(); ch++){
    if(dest.channels.find(ch->first)==dest.channels.end()){   //this channel does not exist, create it
      dest.channels[ch->first]         = ChannelInfo_t();
      dest.channels[ch->first].bin     = ch->second.bin;
      dest.channels[ch->first].channel = ch->second.channel;
    }

    addChannel(dest.channels[ch->first], ch->second, computeSyst);
  }
}

//
// Subtract nonQCD MC from A,C,D regions in QCD Analysis
//
void AllInfo_t::doBackgroundSubtraction(FILE* pFile,std::vector<TString>& selCh,TString mainHisto) {

  // Closure test for DD QCD predictions
  char Lcol     [1024] = "";
  char Lchan    [1024] = "";
  char Lalph1   [1024] = "";
  char Lalph2   [1024] = "";
  char Lyield [1024] = "";
  char LyieldMC [1024] = "";
  // MC closure test for QCD predictions
  char LalphMC  [1024] = "";
  char LyieldPred[1024] = "";
  char RatioMC  [1024] = "";

  //check that the data proc exist
  std::map<string, ProcessInfo_t>::iterator dataProcIt=procs.find("data");             
  if(dataProcIt==procs.end()){printf("The process 'data' was not found... can not do QCD background prediction\n"); return;}

  // create 3 new processes for A,C,D regions in data
  TString NRBProcName = "NonQCD";
  for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==NRBProcName.Data()){sorted_procs.erase(p);break;}} 
  sorted_procs.push_back(NRBProcName.Data());
  procs[NRBProcName.Data()] = ProcessInfo_t(); //reset
  ProcessInfo_t& procInfo_NRB = procs[NRBProcName.Data()];
  procInfo_NRB.shortName = "nonqcd";
  procInfo_NRB.isData = true;
  procInfo_NRB.isSign = false;
  procInfo_NRB.isBckg = true;
  procInfo_NRB.xsec = 0.0;
  procInfo_NRB.br = 1.0;

  //create an histogram containing all the nonQCD MC backgrounds
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(it->second.isData)continue;
    TString procName = it->first.c_str();
    //if(!(procName.Contains("Single Top") || procName.Contains("t#bar{t}+#gammaZW") || procName.Contains("Z#rightarrow") || procName.Contains("VV") || procName.Contains("Vh") || procName.Contains("t#bar{t}") || procName.Contains("W#rightarrow")))continue;
    if(!(procName.Contains("Other Bkgds") || procName.Contains("Z#rightarrow") || procName.Contains("t#bar{t}") || procName.Contains("W#rightarrow")))continue;
        printf("Subtracting nonQCD process from data: %s \n",procName.Data()); 
    addProc(procInfo_NRB, it->second, false);
  }

  for(std::map<string, ChannelInfo_t>::iterator chData = dataProcIt->second.channels.begin(); chData!=dataProcIt->second.channels.end(); chData++){
    if(std::find(selCh.begin(), selCh.end(), chData->second.channel)==selCh.end())continue;

    // if(chData->first.find("CR5j"))continue;
    if(chData->first.find("_A_")!=string::npos)continue; // do not subtract nonQCD MC in regions A...

    // now look at the channels in the nonQCD MC process
    std::map<string, ChannelInfo_t>::iterator chNRB = procInfo_NRB.channels.find(chData->first); 
    if(chNRB==procInfo_NRB.channels.end()){  //this channel does not exist, create it
      procInfo_NRB.channels[chData->first] = ChannelInfo_t();     
      chNRB                = procInfo_NRB.channels.find(chData->first);
      chNRB->second.bin     = chData->second.bin;
      chNRB->second.channel = chData->second.channel;
    }

    // Subtract NonQCD MC from Data in regions B,C,D: 
    std::map<string, ShapeData_t>& shapesInfoDest = chData->second.shapes;
    std::map<string, ShapeData_t>& shapesInfoSrc = chNRB->second.shapes;

    for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
      if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();
    
      //Loop on all shape systematics (including also the central value shape)
      for(std::map<string, TH1*>::iterator uncS = sh->second.uncShape.begin();uncS!= sh->second.uncShape.end();uncS++){
	if(uncS->first!="") continue; //We only take nominal shapes
	//	if(shapesInfoDest[sh->first].uncShape.find(uncS->first)==shapesInfoDest[sh->first].uncShape.end()){
	//	  shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) uncS->second->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
	//	}else{
	//	printf("Here start subtracting NonQCD %s \n",uncS->second->GetName());
	shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second,-1);
	  //	}
      }
    }
  
    // Now here all the data channels B,C,D have subtracted NonQCD components
  } // end data channels

  // Replace shape qcdD*(qcdA/qcdC) with qcd shape B     

  //create a new proc for DD QCD backgrounds
  TString DDProcName = "ddqcd";
  for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==DDProcName.Data()){sorted_procs.erase(p);break;}}           
  sorted_procs.push_back(DDProcName.Data());
  procs[DDProcName.Data()] = ProcessInfo_t(); //reset
  ProcessInfo_t& procInfo_DD = procs[DDProcName.Data()];
  procInfo_DD.shortName = "ddqcd";
  procInfo_DD.isData = true;
  procInfo_DD.isSign = false;
  procInfo_DD.isBckg = true;
  procInfo_DD.xsec   = 0.0;
  procInfo_DD.br     = 1.0;

  //create an histogram containing all the QCD MC backgrounds
  std::vector<string> toBeDelete;
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(!it->second.isBckg || it->second.isData)continue;
    TString procName = it->first.c_str();

    if(!( procName.Contains("QCD"))) continue;
    addProc(procInfo_DD, it->second);

    for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==it->first){sorted_procs.erase(p);break;}}
    toBeDelete.push_back(it->first);
  }
  for(std::vector<string>::iterator p=toBeDelete.begin();p!=toBeDelete.end();p++){procs.erase(procs.find((*p)));}

  for(std::map<string, ChannelInfo_t>::iterator chData = dataProcIt->second.channels.begin(); chData!=dataProcIt->second.channels.end(); chData++){
    if(std::find(selCh.begin(), selCh.end(), chData->second.channel)==selCh.end())continue;

    if(!(chData->first.find("_A_")!=string::npos))continue; // replace only in regions A...

    //   printf("Channel 2:%s\n",chData->first.c_str());
    //  if (chData->second.channel.find("CR5j")) continue;
    
    std::map<string, ChannelInfo_t>::iterator chDD = procInfo_DD.channels.find(chData->first); 
    if(chDD==procInfo_DD.channels.end()){  //this channel does not exist, create it
      procInfo_DD.channels[chData->first] = ChannelInfo_t();     
      chDD                = procInfo_DD.channels.find(chData->first);
      chDD->second.bin     = chData->second.bin;
      chDD->second.channel = chData->second.channel;
    }

    TString binName;
    binName=chData->second.channel.c_str(); // e.g. e_A_SR
    //load data histograms in the QCD control regions
    binName.ReplaceAll("A_","D_");         
    //    printf("binName= %s\n",binName.Data());  
    TH1* hCtrl_SB = dataProcIt->second.channels[(binName+"_"+chData->second.bin.c_str()).Data()].shapes[mainHisto.Data()].histo(); // Region D  
    binName.ReplaceAll("D_","B_");    
    TH1* hCtrl_SI = dataProcIt->second.channels[(binName+"_"+chData->second.bin.c_str()).Data()].shapes[mainHisto.Data()].histo();  // Region B 
    binName.ReplaceAll("B_","C_");       
    TH1* hChan_SB = dataProcIt->second.channels[(binName+"_"+chData->second.bin.c_str()).Data()].shapes[mainHisto.Data()].histo(); // Region C

    TH1* hDD     =  chDD->second.shapes[mainHisto.Data()].histo(); // Region B
    //    if(hDD==NULL){std::cout << "hDD does not exist:" << chDD->second.bin << "_" << chDD->second.channel << mainHisto.Data() << std::endl;}
    
    // load MC histograms in the QCD control regions
    //    binName.ReplaceAll("C_","B_");
    TH1* hDD_C = procInfo_DD.channels.find((binName+"_"+chData->second.bin.c_str()).Data())->second.shapes[mainHisto.Data()].histo(); 
    binName.ReplaceAll("C_","D_");  
    TH1* hDD_D = procInfo_DD.channels.find((binName+"_"+chData->second.bin.c_str()).Data())->second.shapes[mainHisto.Data()].histo();
    binName.ReplaceAll("D_","B_");    
    TH1* hDD_B = procInfo_DD.channels.find((binName+"_"+chData->second.bin.c_str()).Data())->second.shapes[mainHisto.Data()].histo();  
    
    
    binName.ReplaceAll("B_","");   

    //compute alpha
    double alpha=0 ,alpha_err=0;
    double alphaMC=0, alphaMC_err=0;

    double errC,errD;
    double errMC_C, errMC_D;
    double valMC=0, valMC_err=0;
    double valDD_MC=0, valDD_MC_err=0;
    double ratioMC=0, ratioMC_err=0;
 
    if(hCtrl_SB->Integral()>0){
      alpha     = hChan_SB->IntegralAndError(1,hChan_SB->GetXaxis()->GetNbins(),errC) / hCtrl_SB->IntegralAndError(1,hCtrl_SB->GetXaxis()->GetNbins(),errD);
      alpha_err = ( fabs( hChan_SB->Integral() * errD ) + fabs(errC * errD )  ) / pow(hCtrl_SB->Integral(), 2);        
    }
    
    if(hDD && hDD_B && hDD_C && hDD_D){
      // alpha in MC
      if(hDD_D!=NULL && hDD_D->Integral()>0){
        alphaMC = hDD_C->IntegralAndError(1,hDD_C->GetXaxis()->GetNbins(),errMC_C) / hDD_D->IntegralAndError(1,hDD_D->GetXaxis()->GetNbins(),errMC_D);
        alphaMC_err =  ( fabs( hDD_C->Integral() * errMC_D ) + fabs(errMC_C * errMC_D )  ) / pow(hDD_D->Integral(), 2);  
      }
    
      if(hDD!=NULL) valMC = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins(),valMC_err); if(valMC<1E-6){valMC=0.0; valMC_err=0.0;}   

      TH1 *hDD_MC = (TH1*)hDD->Clone("mcobs");
      hDD_MC->Scale(alphaMC);

      valDD_MC = hDD_MC->IntegralAndError(1,hDD_MC->GetXaxis()->GetNbins()+1,valDD_MC_err); if(valDD_MC<1E-6){valDD_MC=0.0; valDD_MC_err=0.0;}

      // Compute ratio of predicted (valDD_MC) vs observed (valDD) in MC
      ratioMC=valDD_MC/valMC;
      ratioMC_err=(fabs(valDD_MC * valMC_err) + fabs(valDD_MC_err * valMC_err)) / pow(valMC, 2);
    }
    if(!hDD) hDD = (TH1*) hCtrl_SI->Clone();

    std::cout << "hDD_B hDD_C hDD_D hCtrl_SI hChan_SB hCtrl_SB" << std::endl;
    std::cout << (hDD_B!=NULL) << "     " << (hDD_C!=NULL) << "     " << (hDD_D!=NULL) << "     " << (hCtrl_SI!=NULL) << "        " << (hChan_SB!=NULL) << "        " << (hCtrl_SB!=NULL) << std::endl;
    hDD->Reset();
    hDD->Add(hCtrl_SI , 1.0);

    if(hCtrl_SB->Integral()<=0 || hCtrl_SI->Integral()<0 || hChan_SB->Integral()<0) alpha = 0;
    hDD->Scale(alpha);
    hDD->SetTitle(DDProcName.Data());

    //save values for printout
    double valDD, valDD_err;
    //valDD = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-6){valDD=0.0; valDD_err=0.0;}
    valDD = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-3){valDD=0.0; valDD_err=0.0;}
    
    //remove all syst uncertainty
    chDD->second.shapes[mainHisto.Data()].clearSyst();
    //add syst uncertainty
    chDD->second.shapes[mainHisto.Data()].uncScale[string("CMS_haa4b_sys_ddqcd_") + binName.Data() +"_"+chData->second.bin.c_str() + systpostfix.Data()] = valDD*datadriven_qcd_Syst; //:1.8*valDD;
    //    chDD->second.shapes[mainHisto.Data()].uncScale[string("CMS_haa4b_sys_ddqcd_") + binName.Data() + systpostfix.Data()] = ratioMC<0.5?valDD*datadriven_qcd_Syst:fabs(1.-ratioMC)*valDD;    

    //printout
    sprintf(Lcol    , "%s%s"  ,Lcol,    "|c");
    sprintf(Lchan   , "%s%25s",Lchan,   (string(" &") + chData->second.channel+string(" - ")+chData->second.bin).c_str());
    sprintf(Lalph1  , "%s%25s",Lalph1,  (string(" &") + utils::toLatexRounded(alpha,alpha_err)).c_str());
    sprintf(Lyield  , "%s%25s",Lyield,  (string(" &") + utils::toLatexRounded(valDD,valDD_err,valDD*datadriven_qcd_Syst)).c_str());
    sprintf(LyieldMC, "%s%25s",LyieldMC,(string(" &") + utils::toLatexRounded(valMC,valMC_err)).c_str());
    sprintf(LalphMC  , "%s%25s",LalphMC,  (string(" &") + utils::toLatexRounded(alphaMC,alphaMC_err)).c_str()); 
    sprintf(LyieldPred  , "%s%25s",LyieldPred,  (string(" &") + utils::toLatexRounded(valDD_MC,valDD_MC_err)).c_str());     
    sprintf(RatioMC, "%s%25s",RatioMC, (string(" &") + utils::toLatexRounded(ratioMC,ratioMC_err)).c_str());
  } // end data channels

  procs["NonQCD"] = ProcessInfo_t(); //reset  

  //recompute total background  
  computeTotalBackground(); 

  if(pFile){
    if (postfit){
      fprintf(pFile,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\tiny\n\\begin{center}\n\\caption{Data-driven QCD background estimation (using initial W and Top scale factors).}\n\\label{tab:table}\n");} else {
      fprintf(pFile,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\tiny\n\\begin{center}\n\\caption{Data-driven QCD background estimation.}\n\\label{tab:table}\n");
    }
    fprintf(pFile,"\\begin{tabular}{%s|}\\hline\n", Lcol);
    fprintf(pFile,"channel               %s\\\\\\hline\n", Lchan);
    fprintf(pFile,"$\\text{SF}_{qcd}$ measured    %s\\\\\n", Lalph1);
    fprintf(pFile,"QCD yield predicted (data)            %s\\\\\n", Lyield);
    fprintf(pFile,"QCD yield observed (MC)              %s\\\\\n", LyieldMC);
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"$\\text{SF}_{qcd}$ (MC)    %s\\\\\n", LalphMC);
    fprintf(pFile,"QCD yield predicted (MC)   %s\\\\\n", LyieldPred);
    fprintf(pFile,"QCD yield observed (MC)  %s\\\\\n", LyieldMC); 
    fprintf(pFile,"ratio MC (syst)        %s\\\\\n", RatioMC); 
    fprintf(pFile,"\\hline\n");   
    fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n\\end{document}\n");
  }

}


//
// Sum up all background processes and add this as a total process
//
void AllInfo_t::computeTotalBackground(){
  for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)=="total"){sorted_procs.erase(p);break;}}           
  sorted_procs.push_back("total");
  procs["total"] = ProcessInfo_t(); //reset
  ProcessInfo_t& procInfo_Bckgs = procs["total"];
  procInfo_Bckgs.shortName = "total";
  procInfo_Bckgs.isData = false;
  procInfo_Bckgs.isSign = false;
  procInfo_Bckgs.isBckg = true;
  procInfo_Bckgs.xsec   = 0.0;
  procInfo_Bckgs.br     = 1.0;
  //Compute total background nominal
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(it->first=="total" || it->second.isBckg!=true)continue;
    addProc(procInfo_Bckgs, it->second, false);
  }
  //Compute total background systematics
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(it->first=="total" || it->second.isBckg!=true)continue;
    addProc(procInfo_Bckgs, it->second, true);
  }

}


//
// Replace the Data process by TotalBackground
//
void AllInfo_t::blind() {
  if(procs.find("total")==procs.end())computeTotalBackground();

  if(true){ //always replace data
    //if(procs.find("data")==procs.end()){ //true only if there is no "data" samples in the json file
    sorted_procs.push_back("data");           
    procs["data"] = ProcessInfo_t(); //reset
    ProcessInfo_t& procInfo_Data = procs["data"];
    procInfo_Data.shortName = "data";
    procInfo_Data.isData = true;
    procInfo_Data.isSign = false;
    procInfo_Data.isBckg = false;
    procInfo_Data.xsec   = 0.0;
    procInfo_Data.br     = 1.0;
    for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
      if(it->first!="total")continue;
      /*
      for (std::map<string, ChannelInfo_t>::iterator ch=it->second.channels.begin(); ch!=it->second.channels.end();ch++){ 
	printf(" ---> Blind in channel %s :\n",ch->second.channel.c_str()); //find("CR")==string::npos));

      }
      */
      addProc(procInfo_Data, it->second);
    }
  }
}

void AllInfo_t::getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc){
  if(!pFileInc)pFileInc=pFile;

  std::vector<string> VectorProc;
  std::map<string, bool> MapChannel;
  std::map<string, std::map<string, string> > MapProcChYields;         
  std::map<string, bool> MapChannelBin;
  std::map<string, std::map<string, string> > MapProcChYieldsBin;         

  std::map<string, string> rows;
  std::map<string, string> rowsBin;
  string rows_header = "\\begin{tabular}{|c|";
  string rows_title  = "channel";

  //order the proc first
  sortProc();

  for(unsigned int p=0;p<sorted_procs.size();p++){
    string procName = sorted_procs[p];
    std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
    if(it==procs.end())continue;
    rows_header += "c|";
    rows_title  += "& " + it->second.shortName;
    std::map<string, double> bin_valerr;
    std::map<string, double> bin_val;
    std::map<string, double> bin_systUp;
    std::map<string, double> bin_systDown;

    VectorProc.push_back(it->first);
    for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
      if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
      if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;

      // Only get yields from shapes for regions A 
      if(modeDD && ((ch->first.find("_B_")!=string::npos) || (ch->first.find("_C_")!=string::npos) ||(ch->first.find("_D_")!=string::npos)))continue;   

      printf("Get yields from shapes:\n");
      printf("Process: %s , channel: %s \n",procName.c_str(),ch->first.c_str());

      TH1* h = ch->second.shapes[histoName].histo();
      double valerr = 0.;
      double val  = 0.;
      if (h!=NULL) val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);

      double syst_scale = std::max(0.0, ch->second.shapes[histoName].getScaleUncertainty());
      double syst_shapeUp = std::max(0.0, ch->second.shapes[histoName].getIntegratedShapeUncertainty((it->first+ch->first).c_str(), "Up"));
      double syst_shapeDown = std::max(0.0, ch->second.shapes[histoName].getIntegratedShapeUncertainty((it->first+ch->first).c_str(), "Down"));
      double systUp = sqrt(pow(syst_scale,2)+pow(syst_shapeUp,2));
      double systDown = sqrt(pow(syst_scale,2)+pow(syst_shapeDown,2));
      systUp= (systUp >0)?systUp: -1; //Set to -1 if no syst, to be coherent with other convention in this file
      systDown= (systDown >0)?systDown: -1; //Set to -1 if no syst, to be coherent with other convention in this file
      if(val<1E-5 && valerr>=10*val && procName.find("ww")!=std::string::npos){val=0.0;}
      else if(val<1E-5 && valerr>=10*val){val=0.0; systUp=-1;systDown=-1;}
      else if(val<1E-6){val=0.0; valerr=0.0; systUp=-1;systDown=-1;}
      if(it->first=="data"){valerr=-1.0; systUp=-1;systDown=-1;}
      string YieldText = "";

      if(it->first=="data" || it->first=="total")YieldText += "\\boldmath ";
      if(it->first=="data"){char tmp[256];sprintf(tmp, "$%.0f$", val); YieldText += tmp;
      }else{                YieldText += utils::toLatexRounded(val,valerr, systUp, true, systDown);     }


      printf("%f %f %f %f --> %s\n", val, valerr, systUp, systDown, utils::toLatexRounded(val,valerr, systUp, true, systDown).c_str());

      if(rows.find(ch->first)==rows.end())rows[ch->first] = string("$ ")+ch->first+" $";
      rows[ch->first] += string("&") + YieldText;

      TString LabelText = TString("$") + ch->second.channel+ " " +ch->second.bin + TString("$");
      LabelText.ReplaceAll("eq"," ="); LabelText.ReplaceAll("g =","\\geq"); LabelText.ReplaceAll("l =","\\leq"); 
      //      LabelText.ReplaceAll("_OS","OS "); LabelText.ReplaceAll("el","e"); LabelText.ReplaceAll("mu","\\mu");  LabelText.ReplaceAll("ha","\\tau_{had}");

      TString BinText = TString("$") + ch->second.bin + TString("$");
      BinText.ReplaceAll("eq"," ="); BinText.ReplaceAll("g =","\\geq"); BinText.ReplaceAll("l =","\\leq");
      //      BinText.ReplaceAll("_OS","OS "); BinText.ReplaceAll("el","e"); BinText.ReplaceAll("mu","\\mu");  BinText.ReplaceAll("ha","\\tau_{had}");


      bin_val   [BinText.Data()] = val;
      bin_valerr[BinText.Data()] = pow(valerr,2);
      bin_systUp  [BinText.Data()] = systUp>=0?pow(systUp,2):-1;
      bin_systDown  [BinText.Data()] = systDown>=0?pow(systDown,2):-1;
      if(systUp<0)bin_systUp  [BinText.Data()]=-1;
      if(systDown<0)bin_systDown  [BinText.Data()]=-1;

      bin_val   [" Inc."] += val;
      bin_valerr[" Inc."] += pow(valerr,2);
      bin_systUp  [" Inc."] += systUp>=0?pow(systUp,2):0; //We are doing a quadratic sum here, have to add 0 if we have negative value
      bin_systDown  [" Inc."] += systDown>=0?pow(systDown,2):0; //We are doing a quadratic sum here, have to add 0 if we have negative value

      MapChannel[LabelText.Data()] = true;
      MapProcChYields[it->first][LabelText.Data()] = YieldText;
    }
    if(bin_systUp  [" Inc."] <= 0) bin_systUp  ["Inc"]=-1; //If negative value, or 0, set it to -1
    if(bin_systDown  [" Inc."] <= 0) bin_systDown  ["Inc"]=-1; //If negative value, or 0, set it to -1

    for(std::map<string, double>::iterator bin=bin_val.begin(); bin!=bin_val.end(); bin++){
      string YieldText = "";                 
      if(it->first=="data" || it->first=="total" || bin->first==" Inc.")YieldText += "\\boldmath ";
      if(it->first=="data"){char tmp[256];sprintf(tmp, "%.0f", bin_val[bin->first]); YieldText += tmp;  //unblinded
				//                 if(it->first=="data"){char tmp[256];sprintf(tmp, "-"); rowsBin[bin->first] += tmp;  //blinded
      }else{                YieldText += utils::toLatexRounded(bin_val[bin->first],sqrt(bin_valerr[bin->first]), bin_systUp[bin->first]<0?-1:sqrt(bin_systUp[bin->first]), true, bin_systDown[bin->first]<0?-1:sqrt(bin_systDown[bin->first]));   }

      if(rowsBin.find(bin->first)==rowsBin.end())rowsBin[bin->first] = string("$ ")+bin->first+" $";
      rowsBin[bin->first] += string("&") + YieldText;

      MapChannelBin[bin->first] = true;
      MapProcChYieldsBin[it->first][bin->first] = YieldText;                

      if(bin->first==" Inc."){
        MapChannel[bin->first] = true;
        MapProcChYields[it->first][bin->first] = YieldText;                
      }
    }
  }


		//           //All Channels
		//           fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
		//           fprintf(pFile, "%s}\\\\\n", rows_header.c_str());
		//           fprintf(pFile, "%s\\\\\n", rows_title .c_str());
		//           for(std::map<string, string>::iterator row = rows.begin(); row!= rows.end(); row++){
		//              fprintf(pFile, "%s\\\\\n", row->second.c_str());
		//           }
		//           fprintf(pFile,"\\hline\n");
		//           fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");

		//           //All Bins
		//           fprintf(pFileInc,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
		//           fprintf(pFileInc, "%s}\\\\\n", rows_header.c_str());
		//           fprintf(pFileInc, "%s\\\\\n", rows_title .c_str());
		//           for(std::map<string, string>::iterator row = rowsBin.begin(); row!= rowsBin.end(); row++){
		//              fprintf(pFileInc, "%s\\\\\n", row->second.c_str());
		//           }
		//           fprintf(pFileInc,"\\hline\n");
		//           fprintf(pFileInc,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");



    //All Channels
  fprintf(pFile,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
  fprintf(pFile, "\\begin{tabular}{|c|"); for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ fprintf(pFile, "c|"); } fprintf(pFile, "}\\\\\n");
  fprintf(pFile, "channel");   for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ fprintf(pFile, " & %s", ch->first.c_str()); } fprintf(pFile, "\\\\\\hline\n");
  for(auto proc = VectorProc.begin();proc!=VectorProc.end(); proc++){
    if(*proc=="total")fprintf(pFile, "\\hline\n");
    auto ChannelYields = MapProcChYields.find(*proc);
    if(ChannelYields == MapProcChYields.end())continue;
    fprintf(pFile, "%s ", proc->c_str()); 
    for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ 
      fprintf(pFile, " & ");
      if(ChannelYields->second.find(ch->first)!=ChannelYields->second.end()){
	fprintf(pFile, " %s", (ChannelYields->second)[ch->first].c_str());
      }
    }
    fprintf(pFile, "\\\\\n");
    if(*proc=="data")fprintf(pFile, "\\hline\n");             
  }
  fprintf(pFile,"\\hline\n");
  fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n\\end{document}\n");
  
    //All Bins
  fprintf(pFileInc,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
  fprintf(pFileInc, "\\begin{tabular}{|c|"); for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ fprintf(pFileInc, "c|"); } fprintf(pFileInc, "}\\\\\\hline\n");
  fprintf(pFileInc, "channel");   for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ fprintf(pFileInc, " & %s", ch->first.c_str()); } fprintf(pFileInc, "\\\\\\hline\n");
  for(auto proc = VectorProc.begin();proc!=VectorProc.end(); proc++){
    if(*proc=="total")fprintf(pFileInc, "\\hline\n");
    auto ChannelYields = MapProcChYieldsBin.find(*proc);
    if(ChannelYields == MapProcChYieldsBin.end())continue;
    fprintf(pFileInc, "%s ", proc->c_str()); 
    for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ 
      fprintf(pFileInc, " & ");
      if(ChannelYields->second.find(ch->first)!=ChannelYields->second.end()){
	fprintf(pFileInc, " %s", (ChannelYields->second)[ch->first].c_str());
      }
    }
    fprintf(pFileInc, "\\\\\n");
    if(*proc=="data")fprintf(pFileInc, "\\hline\n");             
  }
  fprintf(pFileInc,"\\hline\n");
  fprintf(pFileInc,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n\\end{document}\n");
  
  
  
}

  //
  // Dump efficiencies
  //
  void AllInfo_t::getEffFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName)
  {
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      if(!it->second.isSign)continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        TH1* h = ch->second.shapes[histoName].histo();
        double valerr = 0.;
        double val = 0.;
	if (h!=NULL) val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
        fprintf(pFile,"%30s %30s %4.0f %6.2E %6.2E %6.2E %6.2E\n",ch->first.c_str(), it->first.c_str(), it->second.mass, it->second.xsec, it->second.br, val/(it->second.xsec*it->second.br), valerr/(it->second.xsec*it->second.br));
      }
    }
  }




  //
  // drop control channels
  //
  void AllInfo_t::dropCtrlChannels(std::vector<TString>& selCh)
  {
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end()){it->second.channels.erase(ch); ch=it->second.channels.begin();}
      }
    }
  }



  //
  // drop background process that have a negligible yield
  //
  void AllInfo_t::dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold)
  {
    std::map<string, double> map_yields;          
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      if(!it->second.isBckg)continue;
      map_yields[it->first] = 0;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
	//	map_yields[it->first] += ch->second.shapes[histoName].histo()->Integral();
	TH1 *h=ch->second.shapes[histoName].histo();
        if (h!=NULL) { map_yields[it->first] += h->Integral();}
	//	else {map_yields[it->first] += 0.;}
      }
    }

    double total = map_yields["total"];
    for(std::map<string, double>::iterator Y=map_yields.begin();Y!=map_yields.end();Y++){
      //      if(Y->first.find("ddqcd")<std::string::npos)continue;//never drop this background
      //      if(Y->first.find("VV")<std::string::npos)continue;//never drop this background
      // if(Y->first.find("Vh")<std::string::npos)continue;//never drop this background 
      //      if(Y->first.find("t#bar{t}+#gammaZW")<std::string::npos)continue;//never drop this background    
      if(Y->second/total<threshold){
        printf("Drop %s from the list of backgrounds because of negligible rate (%f%% of total bckq)\n", Y->first.c_str(), Y->second/total);
        for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==Y->first ){sorted_procs.erase(p);break;}}
        procs.erase(procs.find(Y->first ));
      }
    }
  }

  //
  // Make a summary plot
  //
  void AllInfo_t::showShape(std::vector<TString>& selCh , TString histoName, TString SaveName)
  {
    int NLegEntry = 0;
    
    std::map<string, THStack*          > map_stack;
    std::map<string, TH1*              > map_mc;
    std::map<string, TGraphAsymmErrors*     > map_unc;
    std::map<string, TH1*              > map_uncH;
    std::map<string, TH1*              > map_data;
    std::map<string, TGraphAsymmErrors*> map_dataE;
    std::map<string, std::vector<TH1*> > map_signals;
    std::map<string, int               > map_legend;
		//           TLegend* legA  = new TLegend(0.6,0.5,0.99,0.85, "NDC");
		//           TLegend* legA  = new TLegend(0.03,0.00,0.97,0.70, "NDC");
		//           TLegend* legA  = new TLegend(0.03,0.99,0.97,0.89, "NDC");
    //  TLegend* legA  = new TLegend(0.08,0.89,0.97,0.95, "");
    TLegend *legA = new TLegend(0.30,0.74,0.93,0.96, "NDC");
    legA->SetHeader("");
    legA->SetNColumns(3);   
    legA->SetBorderSize(0);
    legA->SetTextFont(42);   legA->SetTextSize(0.03);
    legA->SetLineColor(0);   legA->SetLineStyle(1);   legA->SetLineWidth(1);
    legA->SetFillColor(0); legA->SetFillStyle(0);//blind>-1E99?1001:0);
    std::vector<TLegendEntry*> legEntries;  //needed to have the entry in reverse order

    //order the proc first
    sortProc();

    //loop on sorted proc
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      TString process(procName.c_str());
      //      if( process.Contains("BOnly_B") || process.Contains("SandBandInterf_SBI") ) continue;

      if(it==procs.end())continue;
      //loop on channels for each process
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
	if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;

        if(ch->second.shapes.find(histoName.Data())==(ch->second.shapes).end())continue;

        TH1* h = ch->second.shapes[histoName.Data()].histo();
	if (!h) continue;
        //if(process.Contains("SOnly_S") ) h->Scale(10);  
        if(it->first=="total"){
          //double Uncertainty = std::max(0.0, ch->second.shapes[histoName.Data()].getScaleUncertainty() / h->Integral() );;
          double syst_scale = std::max(0.0, ch->second.shapes[histoName.Data()].getScaleUncertainty());
          //double syst_shape = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str()));
          //double syst = sqrt(pow(syst_scale,2)+pow(syst_shape,2)); //On perd de l'info ici car on considere l'ecart constant au lieu d'y aller bin par bin
          //double Uncertainty = syst / h->Integral();
          double Uncertainty_scale=syst_scale / h->Integral();

          double Maximum = 0;
          TGraphAsymmErrors* errors = new TGraphAsymmErrors(h->GetXaxis()->GetNbins());
					//                    errors->SetFillStyle(3427);
					//                    errors->SetFillColor(kGray+1);
          errors->SetFillStyle(3005);
          errors->SetFillColor(kGray+3);                    
          errors->SetLineStyle(1);
          errors->SetLineColor(1);
          int icutg=0;
          for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){
            if(h->GetBinContent(ibin)>0)
              errors->SetPoint(icutg,h->GetXaxis()->GetBinCenter(ibin), h->GetBinContent(ibin));
            //This is the part where we define which errors will be shown on the shape plot
            double syst_shape_binUp = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str(), ibin, "Up"));
            double syst_shape_binDown = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str(), ibin, "Down"));
	    double syst_binUp = sqrt(pow(Uncertainty_scale*h->GetBinContent(ibin), 2) + pow(syst_shape_binUp,2));
	    double syst_binDown = sqrt(pow(Uncertainty_scale*h->GetBinContent(ibin), 2) + pow(syst_shape_binDown,2));
	    double Uncertainty_binUp = syst_binUp / h->GetBinContent(ibin);
	    double Uncertainty_binDown = syst_binDown / h->GetBinContent(ibin);

            //errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0, sqrt(pow(h->GetBinContent(ibin)*Uncertainty_bin,2) + pow(h->GetBinError(ibin),2) ) );
	    errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0,h->GetXaxis()->GetBinWidth(ibin)/2.0, sqrt(pow(h->GetBinContent(ibin)*Uncertainty_binDown,2) + pow(h->GetBinError(ibin),2) ), sqrt(pow(h->GetBinContent(ibin)*Uncertainty_binUp,2) + pow(h->GetBinError(ibin),2) ) );
	    
	    //                        printf("Unc=%6.2f  X=%6.2f Y=%6.2f+-%6.2f+-%6.2f=%6.2f\n", Uncertainty, h->GetXaxis()->GetBinCenter(ibin), h->GetBinContent(ibin), h->GetBinContent(ibin)*Uncertainty, h->GetBinError(ibin), sqrt(pow(h->GetBinContent(ibin)*Uncertainty,2) + pow(h->GetBinError(ibin),2) ) );
	    //                        errors->SetPointError(icutg,h->GetXaxis()->GetBinWidth(ibin)/2.0, 0 );
            Maximum =  std::max(Maximum , h->GetBinContent(ibin) + errors->GetErrorYhigh(icutg));
            icutg++;
          }errors->Set(icutg);
          errors->SetMaximum(Maximum);
          map_unc[ch->first] = errors;

	  map_uncH[ch->first] = (TH1D*)h->Clone((ch->first+"histPlusSyst").c_str()); //utils::root::checkSumw2((ch->first+"histPlusSyst").c_str());
	  // loop over hist to set the stat+syst errors on map_uncH
	  icutg=0;
	  for(int ibin=1; ibin<=map_uncH[ch->first]->GetXaxis()->GetNbins(); ibin++){
	    double ierr=errors->GetErrorY(icutg);
	    map_uncH[ch->first]->SetBinError(ibin,ierr);
	    icutg++;
	  }

	  continue;//otherwise it will fill the legend
        }else if(it->second.isBckg){                 
          if(map_stack.find(ch->first)==map_stack.end()){
            map_stack[ch->first] = new THStack((ch->first+"stack").c_str(),(ch->first+"stack").c_str());
	    map_mc   [ch->first] = (TH1D*)h->Clone((ch->first+"mc").c_str()); //utils::root::checkSumw2((ch->first+"mc").c_str());
          }
          map_stack   [ch->first]->Add(h,"HIST");
	  map_mc [ch->first]->Add(h); 
	  //if(h!=NULL && h->Integral()>0){map_mc   [ch->first] = (TH1D*)h->Clone((ch->first+"mc").c_str());utils::root::checkSumw2((ch->first+"mc").c_str());}else{map_mc   [ch->first]->Add(h);}

        }else if(it->second.isSign){                    
          map_signals [ch->first].push_back(h);

        }else if(it->first=="data"){
          h->SetFillStyle(0);
          h->SetFillColor(0);
          h->SetMarkerSize(0.7);
          h->SetMarkerStyle(20);
          h->SetMarkerColor(1);
          h->SetBinErrorOption(TH1::kPoisson);
          map_data[ch->first] = h;

          //poisson error bars
          const double alpha = 1 - 0.6827;
          TGraphAsymmErrors * g = new TGraphAsymmErrors(h);
          g->SetMarkerSize(0.7);
          g->SetMarkerStyle (20);

          for (int i = 0; i < g->GetN(); ++i) {
            int N = g->GetY()[i];
            double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
            double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
            g->SetPointEYlow(i, N-L);
            g->SetPointEYhigh(i, U-N);
          }
          map_dataE[ch->first] = g;
        }

        if(map_legend.find(it->first)==map_legend.end()){
          map_legend[it->first]=1;
          if(it->first=="data"){
            legA->AddEntry(h,it->first.c_str(),"PE0");
          }else if(it->second.isSign){
            legEntries.insert(legEntries.begin(), new TLegendEntry(h, it->first.c_str(), "L") );
          }else{
            legEntries.push_back(new TLegendEntry(h, it->first.c_str(), "F") );
          }
          NLegEntry++;
        }
      }
    }
    //fill the legend in reverse order
    while(!legEntries.empty()){
      legA->AddEntry(legEntries.back()->GetObject(), legEntries.back()->GetLabel(), legEntries.back()->GetOption());
      legEntries.pop_back();      
    }
    if(map_unc.begin()!=map_unc.end())legA->AddEntry(map_unc.begin()->second, "Syst. + Stat.", "F");

   
    TCanvas* c[50];
    int I=1;
    for(std::map<string, THStack*>::iterator p = map_stack.begin(); p!=map_stack.end(); p++){
      //init tab
      string ires;
      ostringstream convert;
      convert << I;   
      
      int NBins = map_data.size()/selCh.size();
      c[I] = new TCanvas("c_"+(char)I,"c_",800,800); //selCh.size());
      //   c[I]->SetTopMargin(0.00); c[I]->SetRightMargin(0.00); c[I]->SetBottomMargin(0.00);  c[I]->SetLeftMargin(0.00);
      //   TPad* t1 = new TPad("t1","t1", 0.03, 0.03, 1.00, 0.90, 4, 1);  t1->Draw();  t1->cd();
      //  t1->SetTopMargin(0.07); t1->SetRightMargin(0.07); t1->SetBottomMargin(0.07);  t1->SetLeftMargin(0.12);
      //t1->SetFillColor(0);
      TPad* t1 = new TPad("t1","t1", 0.0, 0.2, 1.0, 1.0);
      t1->SetFillColor(0);
      t1->SetBorderMode(0);
      t1->SetBorderSize(2);
      t1->SetTickx(1);
      t1->SetTicky(1);
      t1->SetLeftMargin(0.10);
      t1->SetRightMargin(0.05);
      t1->SetTopMargin(0.05);
      t1->SetBottomMargin(0.10);
      t1->SetFrameFillStyle(0);
      t1->SetFrameBorderMode(0);
      t1->SetFrameFillStyle(0);
      t1->SetFrameBorderMode(0);
      
      t1->Draw();
      t1->cd();

      t1->SetLogy(useLogy); 
 
      //print histograms
      TH1* axis = (TH1*)map_data[p->first]->Clone("axis");
      axis->Reset();      
      //      axis->GetXaxis()->SetRangeUser(axis->GetXaxis()->FindBin(0.), axis->GetXaxis()->GetXmax());
      axis->GetXaxis()->SetRangeUser(axis->GetXaxis()->GetXmin(),axis->GetXaxis()->GetXmax());
      //double signalHeight=0;
      //for(unsigned int s=0;s<map_signals[p->first].size();s++){signalHeight = std::max(signalHeight, map_signals[p->first][s]->GetMaximum());}
      //axis->SetMaximum(1.5*std::max(signalHeight , std::max( map_unc[p->first]->GetMaximum(), map_data[p->first]->GetMaximum())));
      axis->SetMaximum(1.5*std::max(map_unc[p->first]->GetMaximum(), map_data[p->first]->GetMaximum()));       
      //axis->SetMaximum(5000.5*std::max(map_unc[p->first]->GetMaximum(), map_data[p->first]->GetMaximum()));
      
      //hard code range
      if(useLogy){
        if(procs["data"].channels[p->first].bin.find("vbf")!=string::npos){
          axis->SetMinimum(1E-1);
          axis->SetMaximum(std::max(axis->GetMaximum(), 5E1));
        }else{
          axis->SetMinimum(1E-2);
          axis->SetMaximum(std::max(axis->GetMaximum(), 1E9));
        }
      }

      axis->GetXaxis()->SetLabelOffset(0.007);
      axis->GetXaxis()->SetLabelSize(0.04);
      axis->GetXaxis()->SetTitleOffset(1.2);
      axis->GetXaxis()->SetTitleFont(42);
      axis->GetXaxis()->SetTitleSize(0.04);
      axis->GetYaxis()->SetLabelFont(42);
      axis->GetYaxis()->SetLabelOffset(0.007);
      axis->GetYaxis()->SetLabelSize(0.04);
      axis->GetYaxis()->SetTitleOffset(1.35);
      axis->GetYaxis()->SetTitleFont(42);
      axis->GetYaxis()->SetTitleSize(0.04);

      if ( startsWith(p->first,"mu_A_SR_3b") || startsWith(p->first,"mu_A_SR_4b") ||
        startsWith(p->first,"e_A_SR_3b") || startsWith(p->first,"e_A_SR_4b") ||
	startsWith(p->first,"mumu_A_SR_3b") || startsWith(p->first,"mumu_A_SR_4b") || 
	startsWith(p->first,"ee_A_SR_3b") || startsWith(p->first,"ee_A_SR_4b") ) {
	
	int bbin=axis->FindBin(0.1); //map_data[p->first]->FindBin(0.1);
       	for(unsigned int i=bbin;i<axis->GetNbinsX()+1; i++){   
       	  axis->SetBinContent(i, 0); axis->SetBinError(i, 0);
	}
      }
      //if((I-1)%NBins!=0)
      axis->GetYaxis()->SetTitle("Events");
      //if(I<=NBins)
      //      axis->GetXaxis()->SetTitle("BDT");
      axis->Draw();

      t1->Update();

      p->second->Draw("same"); // MC stack histogram
      map_unc [p->first]->Draw("2 same");
      for(unsigned int i=0;i<map_signals[p->first].size();i++){
	TH1* hs= map_signals[p->first][i];
        if (hs) {
	  if ( startsWith(p->first,"mu_A_SR_3b") || startsWith(p->first,"mu_A_SR_4b") ||
		  startsWith(p->first,"e_A_SR_3b") || startsWith(p->first,"e_A_SR_4b") ||
		  startsWith(p->first,"mumu_A_SR_3b") || startsWith(p->first,"mumu_A_SR_4b") ||
		  startsWith(p->first,"ee_A_SR_3b") || startsWith(p->first,"ee_A_SR_4b") ) 
	    hs->Scale(signalScale);
	  hs->Draw("HIST same");
	}
      }
      
      if(blindSR){
	 if ( startsWith(p->first,"mu_A_SR_3b") || startsWith(p->first,"mu_A_SR_4b") ||
       	   startsWith(p->first,"e_A_SR_3b") || startsWith(p->first,"e_A_SR_4b") ||
	   startsWith(p->first,"mumu_A_SR_3b") || startsWith(p->first,"mumu_A_SR_4b") || 
	   startsWith(p->first,"ee_A_SR_3b") || startsWith(p->first,"ee_A_SR_4b") ) {

	 std::cout << "Channel name: " << p->first << std::endl;
	 int bbin=axis->FindBin(0.1);
	 int totNbins = map_dataE[p->first]->GetN();
       	 for (int i=bbin-1; i<totNbins; i++){
       	 //for (int i=0; i<map_dataE[p->first]->GetN()+1; i++){
       	   map_dataE[p->first]->RemovePoint(bbin-1);
	 }
	
       	 //TH1 *hist=(TH1*)p->second->GetHistogram();
       	 TPave* blinding_box = new TPave(axis->GetBinLowEdge(axis->FindBin(0.1)), axis->GetMinimum(),axis->GetXaxis()->GetXmax(), axis->GetMaximum(), 0, "NB" );  
       	 blinding_box->SetFillColor(15); blinding_box->SetFillStyle(3013); blinding_box->Draw("same F");
        }
      }

      if(!blindData) map_dataE[p->first]->Draw("P0 same");
      
      bool printBinContent = false;
      if(printBinContent){
        TLatex* tex = new TLatex();
        tex->SetTextSize(0.04); tex->SetTextFont(42);
        tex->SetTextAngle(60);
        TH1* histdata = map_data[p->first];
        double Xrange = axis->GetXaxis()->GetXmax()-axis->GetXaxis()->GetXmin();
        for(int xi=1;xi<=histdata->GetNbinsX();++xi){
          double x=histdata->GetBinCenter(xi);
          double y=histdata->GetBinContent(xi);
          double yData=histdata->GetBinContent(xi);

          int graphBin=-1;  for(int k=0;k<map_unc[p->first]->GetN();k++){if(fabs(map_unc[p->first]->GetX()[k] - x)<histdata->GetBinWidth(xi)){graphBin=k;}} 
          if(graphBin<0){
            printf("MC bin not found for X=%f\n", x);
          }else{
            double yMC =  map_unc[p->first]->GetY()[graphBin];
            double yMCerr = map_unc[p->first]->GetErrorY(graphBin);
            y = std::max(y, yMC+yMCerr);
            if(yMC>=1){tex->DrawLatex(x-0.02*Xrange,y*1.15,Form("#color[4]{B=%.1f#pm%.1f}",yMC, yMCerr));
            }else{     tex->DrawLatex(x-0.02*Xrange,y*1.15,Form("#color[4]{B=%.2f#pm%.2f}",yMC, yMCerr));
            }
          }
          tex->DrawLatex(x+0.02*Xrange,y*1.15,Form("D=%.0f",yData));                 
        }
      }

      //print tab channel header
      TPaveText* Label = new TPaveText(0.2,0.81,0.84,0.89, "NDC");
      Label->SetFillColor(0);  Label->SetFillStyle(0);  Label->SetLineColor(0); Label->SetBorderSize(0);  Label->SetTextAlign(31);
      TString LabelText = procs["data"].channels[p->first].channel+"  "+procs["data"].channels[p->first].bin;
      LabelText.ReplaceAll("e_","e "); 
      //LabelText.ReplaceAll("l=","#leq");LabelText.ReplaceAll("g=","#geq"); 
      //      LabelText.ReplaceAll("_OS","OS "); LabelText.ReplaceAll("el","e"); LabelText.ReplaceAll("mu","#mu");  LabelText.ReplaceAll("ha","#tau_{had}");
      //  Label->AddText(LabelText);  Label->Draw();

      gPad->RedrawAxis();

    //print legend
    // c1->cd(0);
      //      legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);  legA->SetBorderSize(0); legA->SetHeader(LabelText);
      //      legA->SetNColumns((NLegEntry/2) + 1);
      legA->Draw("same");    legA->SetTextFont(42);
      
    //print canvas header
		/*
       t2->cd(0);
		//           TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
    TPaveText* T = new TPaveText(0.1,0.7,0.9,1.0, "NDC");
    T->SetFillColor(0);  T->SetFillStyle(0);  T->SetLineColor(0); T->SetBorderSize(0);  T->SetTextAlign(22);
    if(systpostfix.Contains('3'))      { T->AddText("CMS preliminary, #sqrt{s}=13.0 TeV");
    }else if(systpostfix.Contains('8')){ T->AddText("CMS preliminary, #sqrt{s}=8.0 TeV");
    }else{                               T->AddText("CMS preliminary, #sqrt{s}=7.0 TeV");
    }T->Draw();

*/

    // c1->cd(0);
      /*   
      double L=0.03, R=0.03, T=0.02, B=0.0;
      char LumiText[1024];
      if(systpostfix.Contains('3'))      { double iLumi= 35914;sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", 13.0);
      }else if(systpostfix.Contains('8')){ double iLumi=20000;sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", 8.0);
      }else{                               double iLumi= 5000;sprintf(LumiText, "%.1f %s^{-1} (%.0f TeV)", iLumi>100?iLumi/1000:iLumi, iLumi>100?"fb":"pb", 7.0); 
      }
      */
      double iLumi=35.9;
      double iEcm=13;
      if(lumi > 0) iLumi = lumi;
      /*
      TPaveText* T1 = new TPaveText(1.0-R-0.50, 1.0-T-0.05, 1.02-R, 1.0-T-0.005, "NDC");
      T1->SetTextFont(43); T1->SetTextSize(22);   T1->SetTextAlign(31);
      T1->SetFillColor(0); T1->SetFillStyle(0);   T1->SetBorderSize(0);
      T1->AddText(LumiText);  T1->Draw();
      
      //TOP LEFT IN-FRAME
      TPaveText* T2 = new TPaveText(L+0.005, 1.0-T-0.05, L+0.20, 1.0-T-0.005, "NDC");
      T2->SetTextFont(63); T2->SetTextSize(22);   T2->SetTextAlign(11);
      T2->SetFillColor(0); T2->SetFillStyle(0);   T2->SetBorderSize(0);
      T2->AddText("CMS"); T2->Draw();
      
      if(true){ //Right to CMS
	TPaveText* T3 = new TPaveText(L+0.095, 1.0-T-0.05, L+0.50, 1.0-T-0.005, "NDC");
	T3->SetTextFont(53); T3->SetTextSize(22);   T3->SetTextAlign(11);
	T3->SetFillColor(0); T3->SetFillStyle(0);   T3->SetBorderSize(0);
	T3->AddText("Preliminary"); T3->Draw();
      }
      */
   
      c[I]->cd();
      
      TPad *t2 = new TPad("t2", "t2",0.0,0.0, 1.0,0.2);
      t2->SetFillColor(0);
      t2->SetBorderMode(0);
      t2->SetBorderSize(2);
      t2->SetGridy();
      t2->SetTickx(1);
      t2->SetTicky(1);
      t2->SetLeftMargin(0.10);
      t2->SetRightMargin(0.05);
      t2->SetTopMargin(0.0);
      t2->SetBottomMargin(0.20);
      t2->SetFrameFillStyle(0);
      t2->SetFrameBorderMode(0);
      t2->SetFrameFillStyle(0);
      t2->SetFrameBorderMode(0);
      t2->Draw();
      t2->cd();
      t2->SetGridy(true);
      t2->SetPad(0,0.0,1.0,0.2);


      TH1D* denSystUncH = (TH1D*)map_uncH[p->first];
      utils::root::checkSumw2(denSystUncH);
      
      if(blindSR){
	 if ( startsWith(p->first,"mu_A_SR_3b") || startsWith(p->first,"mu_A_SR_4b") ||
       	   startsWith(p->first,"e_A_SR_3b") || startsWith(p->first,"e_A_SR_4b") ||
	   startsWith(p->first,"mumu_A_SR_3b") || startsWith(p->first,"mumu_A_SR_4b") || 
	   startsWith(p->first,"ee_A_SR_3b") || startsWith(p->first,"ee_A_SR_4b") ) {
	
	  int bbin=denSystUncH->FindBin(0.1); 
       	  for(unsigned int i=bbin;i<=denSystUncH->GetNbinsX()+1; i++){   
       	    denSystUncH->SetBinContent(i, 0); denSystUncH->SetBinError(i, 0);
	  }
	}
      }

      int GPoint=0;
      TGraphErrors *denSystUnc=new TGraphErrors(denSystUncH->GetXaxis()->GetNbins()); 
      for(int xbin=1; xbin<=denSystUncH->GetXaxis()->GetNbins(); xbin++){
	denSystUnc->SetPoint(GPoint, denSystUncH->GetBinCenter(xbin), 1.0);
	denSystUnc->SetPointError(GPoint,  denSystUncH->GetBinWidth(xbin)/2, denSystUncH->GetBinContent(xbin)!=0?denSystUncH->GetBinError(xbin)/denSystUncH->GetBinContent(xbin):0);
	GPoint++;
	
      	//if(denSystUncH->GetBinContent(xbin)==0) {std::cout << "!!!!!!! bin i: " << xbin << " is 0" << std::endl;continue;}
      	if(denSystUncH->GetBinContent(xbin)==0) {continue;}
      	Double_t err=denSystUncH->GetBinError(xbin)/denSystUncH->GetBinContent(xbin);
      	denSystUncH->SetBinContent(xbin,1);
      	denSystUncH->SetBinError(xbin,err);
      }denSystUnc->Set(GPoint);
      
      denSystUnc->SetLineColor(1);
      denSystUnc->SetFillStyle(3004);
      denSystUnc->SetFillColor(kGray+2);
      denSystUnc->SetMarkerColor(1);
      denSystUnc->SetMarkerStyle(1);
      denSystUncH->Reset("ICE");       
      denSystUncH->SetTitle("");
      denSystUncH->SetStats(kFALSE);

      denSystUncH->Draw();
      denSystUnc->Draw("2 0 SAME");
      float yscale = (1.0-0.2)/(0.2);       
      denSystUncH->GetYaxis()->SetTitle("Data/#Sigma Bkg.");
      denSystUncH->GetXaxis()->SetTitle(""); //drop the tile to gain space
      //denSystUncH->GetYaxis()->CenterTitle(true);
      denSystUncH->SetMinimum(0.4);
      denSystUncH->SetMaximum(1.6);
      
      denSystUncH->GetXaxis()->SetLabelFont(42);
      denSystUncH->GetXaxis()->SetLabelOffset(0.007);
      denSystUncH->GetXaxis()->SetLabelSize(0.04 * yscale);
      denSystUncH->GetXaxis()->SetTitleFont(42);
      denSystUncH->GetXaxis()->SetTitleSize(0.035 * yscale);
      denSystUncH->GetXaxis()->SetTitleOffset(0.8);
      denSystUncH->GetYaxis()->SetLabelFont(42);
      denSystUncH->GetYaxis()->SetLabelOffset(0.007);
      denSystUncH->GetYaxis()->SetLabelSize(0.03 * yscale);
      denSystUncH->GetYaxis()->SetTitleFont(42);
      denSystUncH->GetYaxis()->SetTitleSize(0.035 * yscale);
      denSystUncH->GetYaxis()->SetTitleOffset(0.3);
       


         //add comparisons
      //  for(size_t icd=0; icd<compDists.size(); icd++){
      TString name("CompHistogram"); //name+=icd;
      TH1 *dataToObsH = (TH1D*)map_data[p->first]->Clone(name);
      utils::root::checkSumw2(dataToObsH);

      TH1 *mc = (TH1D*)map_mc[p->first]->Clone("mc");
      utils::root::checkSumw2(mc);
      dataToObsH->Divide(mc);
      TGraphErrors* dataToObs = new TGraphErrors(dataToObsH);
      dataToObs->SetMarkerColor(1);
      dataToObs->SetMarkerStyle(20);
      dataToObs->SetMarkerSize(0.7);
      
      if(blindSR){
	if ( startsWith(p->first,"mu_A_SR_3b") || startsWith(p->first,"mu_A_SR_4b") ||
       	startsWith(p->first,"e_A_SR_3b") || startsWith(p->first,"e_A_SR_4b") ||
	startsWith(p->first,"mumu_A_SR_3b") || startsWith(p->first,"mumu_A_SR_4b") || 
	startsWith(p->first,"ee_A_SR_3b") || startsWith(p->first,"ee_A_SR_4b") ) {

	 int bbin=axis->FindBin(0.1);
	 int totNbins = dataToObs->GetN();
       	 for (int i=bbin-1; i<totNbins; i++){
       	   dataToObs->RemovePoint(bbin-1);
	 }
	
        }
      }
      
      dataToObs->Draw("P 0 SAME");

      
      if(blindSR){
        if ( startsWith(p->first,"mu_A_SR_3b") || startsWith(p->first,"mu_A_SR_4b") ||
       	startsWith(p->first,"e_A_SR_3b") || startsWith(p->first,"e_A_SR_4b") ||
	startsWith(p->first,"mumu_A_SR_3b") || startsWith(p->first,"mumu_A_SR_4b") || 
	startsWith(p->first,"ee_A_SR_3b") || startsWith(p->first,"ee_A_SR_4b") ) {
	
       	  TPave* blinding_box = new TPave(axis->GetBinLowEdge(axis->FindBin(0.1)), 0.4,axis->GetXaxis()->GetXmax(), 1.6, 0, "NB" );  
       	  blinding_box->SetFillColor(15); blinding_box->SetFillStyle(3013); blinding_box->Draw("same F");
        }
      }
      
      TLegend *legR = new TLegend(0.56,0.78,0.93,0.96, "NDC");
      legR->SetHeader("");
      legR->SetNColumns(2);
      legR->SetBorderSize(1);
      legR->SetTextFont(42);   legR->SetTextSize(0.03 * yscale);
      legR->SetLineColor(1);   legR->SetLineStyle(1);   legR->SetLineWidth(1);
      legR->SetFillColor(0);   legR->SetFillStyle(1001);//blind>-1E99?1001:0);
      //legR->AddEntry(denRelUnc, "Stat. Unc.", "F");
      legR->AddEntry(denSystUnc, "Syst. + Stat.", "F");
      //legR->Draw("same");
      //  gPad->RedrawAxis();
      t1->cd();
      c[I]->cd();

      utils::root::DrawPreliminary(iLumi, iEcm, t1);
      c[I]->Modified();  
      c[I]->Update();

      //save canvas
      LabelText.ReplaceAll(" ","_"); 
      c[I]->SaveAs(LabelText+"_Shape.png");
      c[I]->SaveAs(LabelText+"_Shape.pdf");
      c[I]->SaveAs(LabelText+"_Shape.C");
      delete c[I];
      
      I++;
    }

  }


  //
  // Make a summary plot
  //
  void AllInfo_t::showUncertainty(std::vector<TString>& selCh , TString histoName, TString SaveName)
  {
    string UncertaintyOnYield="";  char txtBuffer[4096];

    //loop on sorted proc
    for(unsigned int p=0;p<sorted_procs.size();p++){
      int NLegEntry = 0;
      std::map<string, int               > map_legend;
      std::vector<TH1*>                    toDelete;             
      TLegend* legA  = new TLegend(0.03,0.89,0.97,0.95, "");

      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      if(it->first=="total" || it->first=="data")continue;  //only do samples which have systematics
      //      if(it->first=="data")continue;  //only do samples which have systematics  

      std::map<string, bool> mapUncType;
      std::map<string, std::map< string, double> > mapYieldPerBin;
      std::map<string, std::pair< double, double> > mapYieldInc;


      int NBins = it->second.channels.size()/selCh.size();
      TCanvas* c1 = new TCanvas("c1","c1",300*NBins,300*selCh.size());
      c1->SetTopMargin(0.00); c1->SetRightMargin(0.00); c1->SetBottomMargin(0.00);  c1->SetLeftMargin(0.00);
      TPad* t2 = new TPad("t2","t2", 0.03, 0.90, 1.00, 1.00, -1, 1);  t2->Draw();  c1->cd();
      t2->SetTopMargin(0.00); t2->SetRightMargin(0.00); t2->SetBottomMargin(0.00);  t2->SetLeftMargin(0.00);
      TPad* t1 = new TPad("t1","t1", 0.03, 0.03, 1.00, 0.90, 4, 1);  t1->Draw();  t1->cd();
      t1->SetTopMargin(0.00); t1->SetRightMargin(0.00); t1->SetBottomMargin(0.00);  t1->SetLeftMargin(0.00);
      t1->Divide(NBins, selCh.size(), 0, 0);

      int I=1;
      mapYieldInc[""].first = 0;  mapYieldInc[""].second = 0;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++, I++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
        if(ch->second.shapes.find(histoName.Data())==(ch->second.shapes).end())continue;

        //add the stat uncertainty is there;
        //ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );//add stat uncertainty to the uncertainty map;
        if((it->second.shortName).find("wh")!=std::string::npos)ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_wh")).Data(),systpostfix.Data(), false );// attention
	//	else if((it->second.shortName).find("qqH")!=std::string::npos)ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_qqH")).Data(),systpostfix.Data(), false );
        else ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );
        TVirtualPad* pad = t1->cd(I); 
        pad->SetTopMargin(0.06); pad->SetRightMargin(0.03); pad->SetBottomMargin(0.07);  pad->SetLeftMargin(0.06);
	//pad->SetLogy(true); 
	//TH1* h = (TH1*)(ch->second.shapes[histoName.Data()].histo()->Clone((it->first+ch->first+"Nominal").c_str())); 
	TH1* hh = ch->second.shapes[histoName.Data()].histo();
	if (hh==NULL) continue;
	
        TH1* h = (TH1*)(hh->Clone((it->first+ch->first+"Nominal").c_str())); 
	
        double yield = h->Integral();
        toDelete.push_back(h);
        mapYieldPerBin[""][ch->first] = yield;
        mapYieldInc[""].first  += yield;
        mapYieldInc[""].second = 1;

        //print histograms
        TH1* axis = (TH1*)h->Clone("axis");
        axis->Reset();      
        axis->GetXaxis()->SetRangeUser(axis->GetXaxis()->GetXmin(), axis->GetXaxis()->GetXmax());
        axis->GetYaxis()->SetRangeUser(0.5, 1.5); //100% uncertainty
        if((I-1)%NBins!=0)axis->GetYaxis()->SetTitle("");
        axis->Draw();
        toDelete.push_back(axis);


        //print tab channel header
        TPaveText* Label = new TPaveText(0.1,0.81,0.94,0.89, "NDC");
        Label->SetFillColor(0);  Label->SetFillStyle(0);  Label->SetLineColor(0); Label->SetBorderSize(0);  Label->SetTextAlign(31);
        TString LabelText = ch->second.channel+"  -  "+ch->second.bin;
        LabelText.ReplaceAll("eq","="); LabelText.ReplaceAll("l=","#leq");LabelText.ReplaceAll("g=","#geq"); 
        LabelText.ReplaceAll("_OS","OS "); LabelText.ReplaceAll("el","e"); LabelText.ReplaceAll("mu","#mu");  LabelText.ReplaceAll("ha","#tau_{had}");
        Label->AddText(LabelText);  Label->Draw();

        TLine* line = new TLine(axis->GetXaxis()->GetXmin(), 1.0, axis->GetXaxis()->GetXmax(), 1.0);
        toDelete.push_back((TH1*)line);
        line->SetLineWidth(2);  line->SetLineColor(1); line->Draw("same");
        if(I==1){legA->AddEntry(line,"Nominal","L");  NLegEntry++;}

        int ColorIndex=3;

        //draw scale uncertainties
        for(std::map<string, double>::iterator var = ch->second.shapes[histoName.Data()].uncScale.begin(); var!=ch->second.shapes[histoName.Data()].uncScale.end(); var++){
          if(h->Integral()<=0)continue;
          double ScaleChange   = var->second/h->Integral();
          double ScaleUp   = 1 + ScaleChange;
          double ScaleDn   = 1 - ScaleChange;

          TString systName = var->first.c_str();
          systName.ToLower();
          systName.ReplaceAll("cms","");
          systName.ReplaceAll("haa4b","");
          systName.ReplaceAll("sys","");
          systName.ReplaceAll("13tev","");
          systName.ReplaceAll("_","");
          systName.ReplaceAll("up","");
          systName.ReplaceAll("down","");

          TLine* lineUp = new TLine(axis->GetXaxis()->GetXmin(), ScaleUp, axis->GetXaxis()->GetXmax(), ScaleUp);
          TLine* lineDn = new TLine(axis->GetXaxis()->GetXmin(), ScaleDn, axis->GetXaxis()->GetXmax(), ScaleDn);
          toDelete.push_back((TH1*)lineUp);  toDelete.push_back((TH1*)lineDn);


          int color = ColorIndex;
          if(map_legend.find(systName.Data())==map_legend.end()){
            map_legend[systName.Data()]=color;
            legA->AddEntry(lineUp,systName.Data(),"L");
            NLegEntry++;
            ColorIndex++;
          }else{
            color = map_legend[systName.Data()];
          }

          if(mapYieldPerBin[systName.Data()].find(ch->first)==mapYieldPerBin[systName.Data()].end()){
            mapYieldPerBin[systName.Data()][ch->first] = fabs( ScaleChange );                       
            mapUncType[systName.Data()] = false;
          }else{
            mapYieldPerBin[systName.Data()][ch->first] = std::max( ScaleChange , mapYieldPerBin[systName.Data()][ch->first]);
          }

          if(mapYieldInc.find(systName.Data())==mapYieldInc.end()){ mapYieldInc[systName.Data()].first = 0; mapYieldInc[systName.Data()].second = 0;  }
          mapYieldInc[systName.Data()].first += ScaleChange*h->Integral();
          mapYieldInc[systName.Data()].second += h->Integral();

          lineUp->SetLineWidth(2);  lineUp->SetLineColor(color); lineUp->Draw("same");
          lineDn->SetLineWidth(2);  lineDn->SetLineColor(color); lineDn->Draw("same");
        }

        //draw shape uncertainties
        //double syst_shape = std::max(0.0, ch->second.shapes[histoName.Data()].getShapeUncertainty((it->first+ch->first).c_str()));
        for(std::map<string, TH1*>::iterator var = ch->second.shapes[histoName.Data()].uncShape.begin(); var!=ch->second.shapes[histoName.Data()].uncShape.end(); var++){
          if(var->first=="")continue;

          TH1* hvar = (TH1*)(var->second->Clone((it->first+ch->first+var->first).c_str())); 
          double varYield = hvar->Integral();
          hvar->Divide(h);
          toDelete.push_back(hvar);

          TString systName = var->first.c_str();
          systName.ToLower();
          systName.ReplaceAll("cms","");
          systName.ReplaceAll("haa4b","");
          systName.ReplaceAll("sys","");
          systName.ReplaceAll("13tev","");
          systName.ReplaceAll("_","");
          systName.ReplaceAll("up","");
          systName.ReplaceAll("down","");

          int color = ColorIndex;
          if(systName.Contains("stat")){systName = "stat"; color=2;}
          if(map_legend.find(systName.Data())==map_legend.end()){
            map_legend[systName.Data()]=color;
            legA->AddEntry(hvar,systName.Data(),"L");
            NLegEntry++;
            ColorIndex++;
          }else{
            color = map_legend[systName.Data()];
          }

          if(yield>0){
            if(mapYieldPerBin[systName.Data()].find(ch->first)==mapYieldPerBin[systName.Data()].end()){
              mapYieldPerBin[systName.Data()][ch->first] = fabs( 1 - (varYield/yield));
              mapUncType[systName.Data()] = true;                        
            }else{
              mapYieldPerBin[systName.Data()][ch->first] = std::max(fabs( 1 - (varYield/yield) ), mapYieldPerBin[systName.Data()][ch->first]);
            }

            if(mapYieldInc.find(systName.Data())==mapYieldInc.end()){ mapYieldInc[systName.Data()].first = 0; mapYieldInc[systName.Data()].second = 0;  }
            mapYieldInc[systName.Data()].first +=  fabs( 1 - (varYield/yield))*yield;
            mapYieldInc[systName.Data()].second += yield;
          }

          hvar->SetFillColor(0);                  
          hvar->SetLineStyle(1);
          hvar->SetLineColor(color);
          hvar->SetLineWidth(2);
          hvar->Draw("HIST same");                   
        }
        //remove the stat uncertainty
	ch->second.shapes[histoName.Data()].removeStatUnc(); 
      }
      //print legend
      c1->cd(0);
      legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);  legA->SetBorderSize(0); legA->SetHeader("");
      legA->SetNColumns((NLegEntry/2) + 1);
      legA->Draw("same");    legA->SetTextFont(42);

      //print canvas header
      t2->cd(0);
      TPaveText* T = new TPaveText(0.1,0.7,0.9,1.0, "NDC");
      T->SetFillColor(0);  T->SetFillStyle(0);  T->SetLineColor(0); T->SetBorderSize(0);  T->SetTextAlign(22);
      if(systpostfix.Contains('3'))      { T->AddText((string("CMS preliminary, #sqrt{s}=13.0 TeV,   ")+it->first).c_str());
      }else if(systpostfix.Contains('8')){ T->AddText((string("CMS preliminary, #sqrt{s}=8.0 TeV,   ")+it->first).c_str());
      }else{                               T->AddText((string("CMS preliminary, #sqrt{s}=7.0 TeV,   ")+it->first).c_str());
      }T->Draw();

      //save canvas
      c1->SaveAs(SaveName+vh_tag+"_Uncertainty_"+it->second.shortName+".png");
      c1->SaveAs(SaveName+vh_tag+"_Uncertainty_"+it->second.shortName+".pdf");
      c1->SaveAs(SaveName+vh_tag+"_Uncertainty_"+it->second.shortName+".C");
      delete c1;             

      for(unsigned int i=0;i<toDelete.size();i++){delete toDelete[i];} //clear the objects


      //add inclusive uncertainty as a channel
      //
      for(auto systIt=mapYieldInc.begin(); systIt!=mapYieldInc.end(); systIt++){ mapYieldPerBin[systIt->first][" Inc"] = systIt->second.first/systIt->second.second;  }
      //print uncertainty on yield            
      //
      sprintf(txtBuffer, "\\multicolumn{%i}{'c'}{\\bf{%s}}\\\\ \n", I+1, it->first.c_str());  UncertaintyOnYield+= txtBuffer;
      sprintf(txtBuffer, "%10s & %25s", "Type", "Uncertainty");
      for(auto chIt=mapYieldPerBin[""].begin();chIt!=mapYieldPerBin[""].end();chIt++){ sprintf(txtBuffer, "%s & %12s ", txtBuffer, chIt->first.c_str()); } sprintf(txtBuffer, "%s\\\\ \\hline\n", txtBuffer);  UncertaintyOnYield += txtBuffer;
      sprintf(txtBuffer, "%10s & %25s", "", "Nominal yields ");
      for(auto chIt=mapYieldPerBin[""].begin();chIt!=mapYieldPerBin[""].end();chIt++){ sprintf(txtBuffer, "%s & %12.4E ", txtBuffer, chIt->second); } sprintf(txtBuffer, "%s\\\\ \n", txtBuffer);  UncertaintyOnYield += txtBuffer;
      for(auto varIt=mapYieldPerBin.begin();varIt!=mapYieldPerBin.end();varIt++){
        if(varIt->first == "")continue;
        sprintf(txtBuffer, "%10s & %25s", mapUncType[varIt->first.c_str()]?"shape":"scale", varIt->first.c_str());
        for(auto chIt=mapYieldPerBin[""].begin();chIt!=mapYieldPerBin[""].end();chIt++){
          if(varIt->second.find(chIt->first)==varIt->second.end()){ sprintf(txtBuffer, "%s & %10s   "    , txtBuffer, "-" );                        
          }else{                                                    sprintf(txtBuffer, "%s & %+10.1f\\%% ", txtBuffer,100.0 * varIt->second[chIt->first] ); 
          }
        }sprintf(txtBuffer, "%s\\\\ \n", txtBuffer);   UncertaintyOnYield += txtBuffer;
      }sprintf(txtBuffer, "\\hline \n"); UncertaintyOnYield += txtBuffer;
    }

    FILE* pFile = fopen(SaveName+vh_tag+"_Uncertainty.txt", "w");
    if(pFile){ fprintf(pFile, "%s\n", UncertaintyOnYield.c_str()); fclose(pFile);}
  }




  //
  // Turn to cut&count (rebin all histo to 1 bin only)
  //
  void AllInfo_t::turnToCC(string histoName){
    //order the proc first
    sortProc();

    //Loop on processes and channels
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
	//   ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
	// TH1* h = shapeInfo.histo();
	ShapeData_t& shapeInfo = ch->second.shapes[histoName];
	TH1 * hh = shapeInfo.histo();
        double integral = 0.; if (hh!=NULL) integral = hh->Integral();
	
        TString proc = it->second.shortName.c_str();
        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          TString syst   = unc->first.c_str();
          TH1*    hshape = unc->second;
          hshape->SetDirectory(0);

          hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins()); 
          //make sure to also count the underflow and overflow
          double bin  = hshape->GetBinContent(0) + hshape->GetBinContent(1) + hshape->GetBinContent(2);
          double bine = sqrt(hshape->GetBinError(0)*hshape->GetBinError(0) + hshape->GetBinError(1)*hshape->GetBinError(1) + hshape->GetBinError(2)*hshape->GetBinError(2));
          hshape->SetBinContent(0,0);              hshape->SetBinError  (0,0);
          hshape->SetBinContent(1,bin);            hshape->SetBinError  (1,bine);
          hshape->SetBinContent(2,0);              hshape->SetBinError  (2,0);
        }
      }
    }
  }

  //
  // Make a summary plot
  //
  void AllInfo_t::saveHistoForLimit(string histoName, TFile* fout){
    //order the proc first
    sortProc();

    //Loop on processes and channels
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;
        if(!fout->GetDirectory(chbin)){fout->mkdir(chbin);}fout->cd(chbin);

        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        TH1* h = shapeInfo.histo();
	if(h==NULL)continue;
	//                 shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), it->second.isSign );//add stat uncertainty to the uncertainty map;
        //shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );//add stat uncertainty to the uncertainty map;

        //Li Fix
        if((it->second.shortName).find("ggH")!=std::string::npos)shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_ggH")).Data(),systpostfix.Data(), false );// attention
	else if((it->second.shortName).find("qqH")!=std::string::npos)shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_qqH")).Data(),systpostfix.Data(), false );
        else shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );
	fout->cd(chbin);

        TString proc = it->second.shortName.c_str();
	//	printf("\n\n");
	//	printf("Process: %s and channel : %s , with N shapes= %d\n",proc.Data(),ch->first.c_str(),(int)shapeInfo.uncShape.size());

        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          TString syst   = unc->first.c_str();
          TH1*    hshape = unc->second;
          hshape->SetDirectory(0);

	  //	  printf("Shape= %s\n",unc->first.c_str());

          if(syst==""){
            //central shape (for data call it data_obs)
            hshape->SetName(proc); 
            if(it->first=="data"){
              hshape->Write("data_obs");
            }else{
              hshape->Write(proc+postfix);
            }
          }else if(runSystematics && proc!="data" && (syst.Contains("Up") || syst.Contains("Down"))){
            //if empty histogram --> no variation is applied except for stat
	    
            if(!syst.Contains("stat") && (hshape->Integral()<h->Integral()*0.01 || isnan((float)hshape->Integral()))){hshape->Reset(); hshape->Add(h,1); }

//	    std::cout << "*****************" << proc << postfix << syst << "******************" << std::endl;
	    if(hshape->Integral()<=0){
	      for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) //hmirrorshape->SetBinContent(ibin, 1E-10);	    
	        hshape->SetBinContent(ibin, 1E-10);
	    }
            //write variation to file
	    hshape->SetName(proc+syst);
	    hshape->Write(proc+postfix+syst);
          }else if(runSystematics){
            //for one sided systematics the down variation mirrors the difference bin by bin
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst+"Up");
            TH1 *hmirrorshape=(TH1 *)hshape->Clone(proc+syst+"Down");
            for(int ibin=1; ibin<=hmirrorshape->GetXaxis()->GetNbins(); ibin++){
              double bin = 2*h->GetBinContent(ibin)-hmirrorshape->GetBinContent(ibin);
              if(bin<0)bin=0;
              hmirrorshape->SetBinContent(ibin,bin);
            }
            if(hmirrorshape->Integral()<=0) hmirrorshape->SetBinContent(1, 1E-10);
            hmirrorshape->Write(proc+postfix+syst+"Down");
          }

          if(runSystematics && syst!=""){
            TString systName(syst); 
            systName.ReplaceAll("Up",""); systName.ReplaceAll("Down","");//  systName.ReplaceAll("_","");
	    systName.ReplaceAll("Pile","PileUp");
            if(systName.First("_")==0)systName.Remove(0,1);

            TH1 *temp=(TH1*) hshape->Clone();
            temp->Add(h,-1);
            if(temp->Integral()!=0){
              if(shape){
                shapeInfo.uncScale[systName.Data()]=-1.0;
              }else{
                double Unc = fabs(temp->Integral());
                if(shapeInfo.uncScale.find(systName.Data())==shapeInfo.uncScale.end()){
                  shapeInfo.uncScale[systName.Data()]=Unc;
                }else{
                  shapeInfo.uncScale[systName.Data()]=(shapeInfo.uncScale[systName.Data()] + Unc)/2.0;
                }
              }
            }
            delete temp;
          }else if(syst==""){
            shapeInfo.uncScale[syst.Data()]=hshape->Integral();
          }
        }
      }
      fout->cd("..");
    }
  }


  //
  // add hardcoded uncertainties
  //
  void AllInfo_t::addHardCodedUncertainties(string histoName){
    

    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end() || it->first=="total")continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
	ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        double integral = 0.;

	TH1* h=shapeInfo.histo(); if (h) integral=h->Integral();

        //lumi
	if( !(it->second.shortName.find("ttbarbba")!=string::npos) && //!(it->second.shortName.find("ttbarcba")!=string::npos) &&
	    !(it->second.shortName.find("wlnu")!=string::npos) ) {
	  if(!it->second.isData && systpostfix.Contains('3'))shapeInfo.uncScale["lumi_13TeV"] = integral*0.025;
	  if(!it->second.isData && systpostfix.Contains('8'))shapeInfo.uncScale["lumi_8TeV" ] = integral*0.026;
	  if(!it->second.isData && systpostfix.Contains('7'))shapeInfo.uncScale["lumi_7TeV" ] = integral*0.022;
	}
        //Id+Trigger efficiencies combined
	
        if(!it->second.isData){
          if(chbin.Contains("e"  ))  shapeInfo.uncScale["CMS_eff_e"] = integral*0.02; //0.072124;
          if(chbin.Contains("mu"))  shapeInfo.uncScale["CMS_eff_m"] = integral*0.02; //0.061788;
        }
	
	//Normalization uncertainties 
	//	if(it->second.shortName.find("wlnu")!=string::npos){shapeInfo.uncScale["norm_wlnu"] = integral*0.10;}  
	//      if(it->second.shortName.find("ttbarjet")!=string::npos){shapeInfo.uncScale["norm_tt"] = integral*0.10;}  

	//if(it->second.shortName.find("zvv")!=string::npos){shapeInfo.uncScale["norm_zvv"] = integral*0.50;}    
	//if(it->second.shortName.find("vhbb")!=string::npos){shapeInfo.uncScale["norm_vhbb"] = integral*0.10;}     
	//if(it->second.shortName.find("singleto")!=string::npos){shapeInfo.uncScale["norm_singletop"] = integral*0.05;}
	//if(it->second.shortName.find("ttbargam")!=string::npos){shapeInfo.uncScale["norm_topgzw"] = integral*0.15;} 
	if(it->second.shortName.find("otherbkg")!=string::npos){shapeInfo.uncScale["norm_otherbkgds"] = integral*0.53;} 
	if(runZh) {
	  if(it->second.shortName.find("wlnu")!=string::npos){shapeInfo.uncScale["norm_wjet"] = integral*0.02;}     
	} else {
	  if(it->second.shortName.find("zll")!=string::npos){shapeInfo.uncScale["norm_zll"] = integral*0.02;}
	}

	if(it->second.shortName.find("ttbarlig")!=string::npos){shapeInfo.uncScale["norm_toplight"] = integral*0.06;} 
	//	if(it->second.shortName.find("ttbarcba")!=string::npos){shapeInfo.uncScale["norm_topcc"] = integral*0.50;} 
	if (!subFake){
	  if(it->second.shortName.find("qcd")!=string::npos){shapeInfo.uncScale["norm_qcd"] = integral*0.50;} 
	}
        //uncertainties to be applied only in higgs analyses
        if(mass>0){
          //bin migration at theory level
	  // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#WH_Process
	  if (runZh) {
	    if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["QCDscale_wh"]  = integral*0.038;} //QCD scale
	    if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["PDFscale_wh"]  = integral*0.024;} //PDF+as scale 
	  } else {
	    if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["QCDscale_wh"]  = integral*0.006;} //QCD scale
	    if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["PDFscale_wh"]  = integral*0.019;} //PDF+as scale
	  }
        }//end of uncertainties to be applied only in higgs analyses

      }
    }
  }


  //
  // produce the datacards 
  //
  void AllInfo_t::buildDataCards(string histoName, TString url)
  {
    std::vector<string>clean_procs;
    std::vector<string>sign_procs;
    //make a map of all systematics considered and say if it's shape-based or not.
    std::map<string, bool> allChannels;
    std::map<string, bool> allSysts;
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end() || it->first=="total")continue;
      if(it->second.isSign) {
	sign_procs.push_back(procName);
      }
      if(it->second.isBckg) {clean_procs.push_back(procName); }

      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        allChannels[ch->first] = true;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        for(std::map<string, double>::iterator unc=shapeInfo.uncScale.begin();unc!=shapeInfo.uncScale.end();unc++){
          if(unc->first=="")continue;
          allSysts[unc->first] = unc->second==-1?true:false;
        }
      }
    }
    clean_procs.insert(clean_procs.begin(), sign_procs.begin(), sign_procs.end());
    int nsign = sign_procs.size();
    printf("\n\n Now building the datacards . We have %d signal points available \n",nsign);

    TString eecard = "";
    TString mumucard = "";
    TString combinedcard = "";

    TString ecrcard = "";
    TString mucrcard = "";

    for(std::map<string, bool>::iterator C=allChannels.begin(); C!=allChannels.end();C++){
      // REmove B,C,D control regions from the datacard
      if(modeDD && ((C->first.find("B_")!=string::npos) || (C->first.find("C_")!=string::npos) ||(C->first.find("D_")!=string::npos)))continue;  
      // Remove CR5j_3b channel
      //      if (C->first.find("CR5j_3b")!=string::npos) continue;

      TString dcName=url;              
      dcName.ReplaceAll(".root","_"+TString(C->first.c_str())+".dat");
      //      dcName.ReplaceAll("_qcdB","_qcdB");  

      combinedcard += (C->first+"=").c_str()+dcName+" ";
      if(runZh) {
	if(C->first.find("ee"  )!=string::npos)eecard   += (C->first+"=").c_str()+dcName+" ";
	if(C->first.find("mumu")!=string::npos)mumucard += (C->first+"=").c_str()+dcName+" ";

	if(C->first.find("ee_A_CR"  )!=string::npos)ecrcard   += (C->first+"=").c_str()+dcName+" ";      
	if(C->first.find("mumu_A_CR")!=string::npos)mucrcard += (C->first+"=").c_str()+dcName+" ";  
	
	if(C->first.find("emu_A"  )!=string::npos){
	  eecard   += (C->first+"=").c_str()+dcName+" ";mumucard += (C->first+"=").c_str()+dcName+" ";
	}        

      } else {
	if(C->first.find("e"  )!=string::npos)eecard   += (C->first+"=").c_str()+dcName+" ";    
	if(C->first.find("mu")!=string::npos)mumucard += (C->first+"=").c_str()+dcName+" ";  

	if(C->first.find("e_A_CR"  )!=string::npos)ecrcard   += (C->first+"=").c_str()+dcName+" ";  
	if(C->first.find("mu_A_CR")!=string::npos)mucrcard += (C->first+"=").c_str()+dcName+" "; 
      }

      bool TTcontrolregion(false);  
      bool nonTTcontrolregion(false); 

      if(C->first.find("CR5j"  )!=string::npos){ // this is the non-Top Control Region 
	nonTTcontrolregion=true;
      } else if(C->first.find("CR"  )!=string::npos){ // this is the Top Control Region
	TTcontrolregion=true;
      }

      dcName.ReplaceAll("[", "");
      dcName.ReplaceAll("]", "");
      dcName.ReplaceAll("+", "");

      FILE* pFile = fopen(dcName.Data(),"w");
      //header
      fprintf(pFile, "imax 1\n");
      fprintf(pFile, "jmax *\n");
      fprintf(pFile, "kmax *\n");
      fprintf(pFile, "-------------------------------\n");
      if(shape){
	fprintf(pFile, "shapes * * %s %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n",url.Data(), C->first.c_str(), C->first.c_str() );
	fprintf(pFile, "-------------------------------\n");
      }
      //observations
      fprintf(pFile, "bin bin1\n");
      fprintf(pFile, "Observation %f\n", procs["data"].channels[C->first].shapes[histoName].histo()->Integral());
      fprintf(pFile, "-------------------------------\n");

      //yields
      fprintf(pFile,"%55s  ", "bin");     for(unsigned int j=0; j<clean_procs.size(); j++){ 
	fprintf(pFile,"%8s ", "bin1")                     ;}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<clean_procs.size(); j++){ 
	fprintf(pFile,"%8s ", procs[clean_procs[j]].shortName.c_str());
      }  
      fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<clean_procs.size(); j++){ 
	fprintf(pFile,"%8i ", ((int)j)-(nsign-1)    );}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "rate");    for(unsigned int j=0; j<clean_procs.size(); j++){ 
	double fval=0;      
	if (procs[clean_procs[j]].channels[C->first].shapes[histoName].histo()!=NULL) { 
	  fval=procs[clean_procs[j]].channels[C->first].shapes[histoName].histo()->Integral();}    
	fprintf(pFile,"%8f ", fval);     
	//	fprintf(pFile,"%8f ", procs[clean_procs[j]].channels[C->first].shapes[histoName].histo()->Integral() );
      }
      fprintf(pFile,"\n");
      fprintf(pFile, "-------------------------------\n");

      for(std::map<string, bool>::iterator U=allSysts.begin(); U!=allSysts.end();U++){
	//	if(TTcontrolregion || nonTTcontrolregion){
	  //	  if(U->first=="lumi_13TeV") continue; //skip lumi unc in the Control Regions  
	//	}
        char line[2048];
        sprintf(line,"%-45s %-10s ", U->first.c_str(), U->second?"shapeN2":"lnN");
        bool isNonNull = false;
        for(unsigned int j=0; j<clean_procs.size(); j++){
          ShapeData_t& shapeInfo = procs[clean_procs[j]].channels[C->first].shapes[histoName];
          double integral = 0.;
	  if (shapeInfo.histo()!=NULL) integral=shapeInfo.histo()->Integral();
          if(shapeInfo.uncScale.find(U->first)!=shapeInfo.uncScale.end()){   isNonNull = true;   
	    if(U->second) sprintf(line,"%s%8s ",line,"       1");
            else if(integral>0)sprintf(line,"%s%8f ",line,1+(shapeInfo.uncScale[U->first]/integral));
            else sprintf(line,"%s%8f ",line,1+(shapeInfo.uncScale[U->first]));
          }else{ sprintf(line,"%s%8s ",line,"       -");        }
        }
        if(isNonNull)fprintf(pFile, "%s\n", line);
      }

      // Add lines for simultaneous fit in the CRs:  
      fprintf(pFile, "-------------------------------\n");  
      fprintf(pFile,"\n");

      if(runZh) {
	if(C->first.find("ee" )!=string::npos) {         
	  fprintf(pFile,"tt_norm_e rateParam bin1 ttbarbba 1\n"); 
	  fprintf(pFile,"tt_norm_e rateParam bin1 ttbarcba 1\n");    
	  if (C->first.find("3b")!=string::npos) fprintf(pFile,"z_norm_3b_e rateParam bin1 zll 1 \n");    
	  if (C->first.find("4b")!=string::npos) fprintf(pFile,"z_norm_4b_e rateParam bin1 zll 1 \n");
	} else if (C->first.find("mumu" )!=string::npos) {  
	    fprintf(pFile,"tt_norm_mu rateParam bin1 ttbarbba 1\n"); 
	    fprintf(pFile,"tt_norm_mu rateParam bin1 ttbarcba 1\n");    
	    if (C->first.find("3b")!=string::npos) fprintf(pFile,"z_norm_3b_mu rateParam bin1 zll 1 \n");  
	    if (C->first.find("4b")!=string::npos) fprintf(pFile,"z_norm_4b_mu rateParam bin1 zll 1 \n");   
	} else if  (C->first.find("emu" )!=string::npos) {     
	  fprintf(pFile,"tt_norm_e rateParam bin1 ttbarbba 1\n");   
	  fprintf(pFile,"tt_norm_e rateParam bin1 ttbarcba 1\n");  
	  fprintf(pFile,"tt_norm_mu rateParam bin1 ttbarbba 1\n");   
	  fprintf(pFile,"tt_norm_mu rateParam bin1 ttbarcba 1\n");  
	} 
      } else {
	if(C->first.find("e" )!=string::npos) {             
	  fprintf(pFile,"tt_norm_e rateParam bin1 ttbarbba 1\n");      
	  fprintf(pFile,"tt_norm_e rateParam bin1 ttbarcba 1\n");   
	  // if (C->first.find("3b")!=string::npos) 
	  fprintf(pFile,"w_norm_e rateParam bin1 wlnu 1 \n");     
	    // if (C->first.find("4b")!=string::npos) fprintf(pFile,"v_norm_4b_e rateParam bin1 wlnu 1 \n");         
	} else if (C->first.find("mu" )!=string::npos) {    
	  fprintf(pFile,"tt_norm_mu rateParam bin1 ttbarbba 1\n");            
	  fprintf(pFile,"tt_norm_mu rateParam bin1 ttbarcba 1\n"); 
	  //	  if (C->first.find("3b")!=string::npos) 
	  fprintf(pFile,"w_norm_mu rateParam bin1 wlnu 1 \n"); 
	  //	  if (C->first.find("4b")!=string::npos) fprintf(pFile,"v_norm_4b_mu rateParam bin1 wlnu 1 \n"); 
	}                  
      }
      
      fprintf(pFile, "-------------------------------\n");  
      fclose(pFile);

      

    }

    FILE* pFile = fopen("combineCards"+vh_tag+".sh","w");   
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcard + " > " + "card_combined"+vh_tag+".dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + eecard       + " > " + "card_e"+vh_tag+".dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + mumucard     + " > " + "card_mu"+vh_tag+".dat").Data());
    
    fprintf(pFile,"%s;\n",(TString("sed -i '/tt_norm_mu/d' card_e"+vh_tag+".dat")).Data());   
    fprintf(pFile,"%s;\n",(TString("sed -i '/tt_norm_e/d' card_mu"+vh_tag+".dat")).Data());         

    fclose(pFile);         

  }


  //
  // Load histograms from root file and json to memory
  //
  void AllInfo_t::getShapeFromFile(TFile* inF, std::vector<string> channelsAndShapes, int cutBin, JSONWrapper::Object &Root,  double minCut, double maxCut, bool onlyData){
    std::vector<TString> BackgroundsInSignal;

    //iterate over the processes required
    std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
    for(unsigned int i=0;i<Process.size();i++){
      string matchingKeyword="";
      if(!utils::root::getMatchingKeyword(Process[i], keywords, matchingKeyword))continue; //only consider samples passing key filtering


      TString procCtr(""); procCtr+=i;
      TString proc=Process[i].getString("tag", "noTagFound");

      printf("Process (get shape from file) = %s\n",proc.Data());
      printf("1:channelsANdShapes size= %d\n",(int)channelsAndShapes.size());    

      string dirName = proc.Data(); 
	    //std::<TString> keys = Process[i].getString("keys", "noKeysFound");
      if(Process[i].isTagFromKeyword(matchingKeyword, "mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i].getIntFromKeyword(matchingKeyword, "mctruthmode", 0)); dirName += buf; }
      string procSuffix = Process[i].getStringFromKeyword(matchingKeyword, "suffix", "");
      if(procSuffix!=""){dirName += "_" + procSuffix;}
      while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");         

      TDirectory *pdir = (TDirectory *)inF->Get(dirName.c_str());         
      if(!pdir){printf("Directory (%s) for proc=%s is not in the file!\n", dirName.c_str(), proc.Data()); continue;}

      bool isData = Process[i].getBool("isdata", false);
      if(onlyData && !isData) { continue; }//just here to speedup the NRB prediction     
      //      if(proc.Contains(")cp0"))continue; // skip those samples

      bool isSignal = Process[i].getBool("issignal", false);
      if(Process[i].getBool("spimpose", false) && (proc.Contains("ggH") || proc.Contains("qqH") || proc.Contains("Wh") ))isSignal=true;
      //LQ bool isInSignal = Process[i].getBool("isinsignal", false);
      int color = Process[i].getInt("color", 1);
      int lcolor = Process[i].getInt("lcolor", 1);
      int mcolor = Process[i].getInt("mcolor", color);
      int lwidth = Process[i].getInt("lwidth", 1);
      int lstyle = Process[i].getInt("lstyle", 1);
      int fill   = Process[i].getInt("fill"  , 1001);
      int marker = Process[i].getInt("marker", 20);
      /*
      if(isSignal && signalTag!=""){
        if(!proc.Contains(signalTag.c_str()) )continue;
      }
      */

      double procMass=0;  char procMassStr[128] = "";
      if(isSignal &&  mass>0 && (proc.Contains("Wh (")) ){
        if(proc.Contains("H(") && proc.Contains("A(")){sscanf(proc.Data()+proc.First("A")+2,"%lf",&procMass);
        }else if(proc.Contains("H(")){sscanf(proc.Data()+proc.First("H")+2,"%lf",&procMass);
        }else if(proc.Contains("A(")){sscanf(proc.Data()+proc.First("A")+2,"%lf",&procMass);
        }else if(proc.Contains("Wh (")){sscanf(proc.Data()+proc.First("(")+1,"%lf",&procMass);}
	  
	printf("%s --> %f\n",  proc.Data(), procMass);
	
	//skip signal sample not needed
	if(massL!=-1 && massR!=-1){
	  if(procMass!=massL && procMass!=massR)continue; 
	}else{
	  if(procMass!=mass)continue;
	}
	sprintf(procMassStr,"%i",(int)procMass);
	//printf("found signal to be %s\n",  proc.Data());
      }
      
      /*
      if(!isSignal &&  mass>0 && proc.Contains("XH(") && proc.Contains(")#rightarrow WW")){
        sscanf(proc.Data()+proc.First("H(")+2,"%lf",&procMass);
        if(!(procMass==mass || procMass==massL || procMass==massR))continue; //skip XH-->WW background sample not concerned
      }
	*/

      TString procSave = proc;
      //      if(isSignal && mass>0 && proc.Contains("ggH") && proc.Contains("ZZ"))proc = TString("ggH")  +procMassStr;
      if(isSignal && mass>0 && proc.Contains("Wh")) proc = TString("Wh")  +procMassStr;
      //      else if(isSignal && mass>0 && proc.Contains("qqH") && proc.Contains("ZZ"))proc = TString("qqH")  +procMassStr;

      //      if(skipGGH && isSignal && mass>0 && proc.Contains("ggH") )continue;
      ///      if(skipQQH && isSignal && mass>0 && (proc.Contains("qqH") || proc.Contains("VBF")) )continue;

      TString shortName = proc;
      shortName.ToLower();
      shortName.ReplaceAll(procMassStr,"");
      shortName.ReplaceAll("#bar{t}","tbar");
      shortName.ReplaceAll("Z#rightarrow ll","dy");
      shortName.ReplaceAll("#rightarrow","");
      shortName.ReplaceAll("(",""); shortName.ReplaceAll(")","");    shortName.ReplaceAll("+","");    shortName.ReplaceAll(" ","");   shortName.ReplaceAll("/","");  shortName.ReplaceAll("#",""); 
      shortName.ReplaceAll("=",""); shortName.ReplaceAll(".","");    shortName.ReplaceAll("^","");    shortName.ReplaceAll("}","");   shortName.ReplaceAll("{","");  shortName.ReplaceAll(",","");
      shortName.ReplaceAll("ggh", "ggH");
      shortName.ReplaceAll("qqh", "qqH");
      if(shortName.Length()>8)shortName.Resize(8);

      if(procs.find(proc.Data())==procs.end()){sorted_procs.push_back(proc.Data());}
      ProcessInfo_t& procInfo = procs[proc.Data()];
      procInfo.jsonObj = Process[i]; 
      procInfo.isData = isData;
      procInfo.isSign = isSignal;
      procInfo.isBckg = !procInfo.isData && !procInfo.isSign;
      procInfo.mass   = procMass;
      procInfo.shortName = shortName.Data();

      if(procInfo.isSign){
        procInfo.xsec = procInfo.jsonObj["data"].daughters()[0].getDouble("xsec", 1);
        if(procInfo.jsonObj["data"].daughters()[0].isTag("br")){
          std::vector<JSONWrapper::Object> BRs = procInfo.jsonObj["data"].daughters()[0]["br"].daughters();
          double totalBR=1.0; for(size_t ipbr=0; ipbr<BRs.size(); ipbr++){totalBR*=BRs[ipbr].toDouble();}   
          procInfo.br = totalBR;
        }
      }

      //Loop on all channels, bins and shape to load and store them in memory structure
      TH1* syst = (TH1*)pdir->Get("all_optim_systs");
      if(syst==NULL){
	printf("Please check all_optim_systs histo is missing!\n\n");
	syst=new TH1F("all_optim_systs","all_optim_systs",1,0,1);syst->GetXaxis()->SetBinLabel(1,"");
      }

      for(unsigned int c=0;c<channelsAndShapes.size();c++){
	TString chName    = (channelsAndShapes[c].substr(0,channelsAndShapes[c].find(";"))).c_str();
        TString binName   = (channelsAndShapes[c].substr(channelsAndShapes[c].find(";")+1, channelsAndShapes[c].rfind(";")-channelsAndShapes[c].find(";")-1)).c_str();
        TString shapeName = (channelsAndShapes[c].substr(channelsAndShapes[c].rfind(";")+1)).c_str();
        TString ch        = chName+TString("_")+binName;

	//	printf("channel= %s, bin= %s, shape name= %s, ch name= %s\n",chName.Data(), binName.Data(), shapeName.Data(), ch.Data());

        ChannelInfo_t& channelInfo = procInfo.channels[ch.Data()];
        channelInfo.bin        = binName.Data();
        channelInfo.channel    = chName.Data();
        ShapeData_t& shapeInfo = channelInfo.shapes[shapeName.Data()];

	//   printf("%s SYST SIZE=%i\n", (ch+"_"+shapeName).Data(), syst->GetNbinsX() );
        for(int ivar = 1; ivar<=syst->GetNbinsX();ivar++){                
          TH1D* hshape   = NULL;
          TString varName   = syst->GetXaxis()->GetBinLabel(ivar);
          TString histoName = ch+"_"+shapeName+(isSignal?signalSufix:"")+varName ;
	  //  if(shapeName==histo && histoVBF!="" && ch.Contains("vbf"))histoName = ch+"_"+histoVBF+(isSignal?signalSufix:"")+varName ;
          //if(isSignal && ivar==1)printf("Syst %i = %s\n", ivar, varName.Data()); 

          TH2* hshape2D = (TH2*)pdir->Get(histoName ); 
	  //	  else {hshape = (TH1D*)pdir->Get(histoName ); }
          if(!hshape2D){
	    //	    printf("Histo %s is NULL for syst:%s\n", histoName.Data(), varName.Data());   
            if(shapeName==histo && histoVBF!="" && ch.Contains("vbf")){   hshape2D = (TH2*)pdir->Get(TString("all_")+histoVBF+(isSignal?signalSufix:"")+varName);
            }else{                                                        hshape2D = (TH2*)pdir->Get(TString("all_")+shapeName+varName);
            }
	  
            if(hshape2D){
              hshape2D->Reset();
            }else{  //if still no histo, skip this proc...
	      //              printf("Histo %s does not exist for syst:%s\n", histoName.Data(), varName.Data());
              continue;
            }
	  }

          //special treatment for side mass points
          int cutBinUsed = cutBin;
          if(shapeName == histo && !ch.Contains("vbf") && procMass==massL)cutBinUsed = indexcutML[channelInfo.bin];
          if(shapeName == histo && !ch.Contains("vbf") && procMass==massR)cutBinUsed = indexcutMR[channelInfo.bin];

          histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
          hshape   = hshape2D->ProjectionY(histoName,cutBinUsed,cutBinUsed);
	  //  else {hshape = hshape2D; }
          filterBinContent(hshape);

          if(isnan((float)hshape->Integral())){hshape->Reset();}
          hshape->SetDirectory(0);
          hshape->SetTitle(proc);
	  utils::root::fixExtremities(hshape,false,true);
          hshape->SetFillColor(color); hshape->SetLineColor(lcolor); hshape->SetMarkerColor(mcolor);
          hshape->SetFillStyle(fill);  hshape->SetLineWidth(lwidth); hshape->SetMarkerStyle(marker); hshape->SetLineStyle(lstyle);

          //if current shape is the one to cut on, then apply the cuts
          if(shapeName == histo){
            //if(ivar==1 && isSignal)printf("A %s %s Integral = %f\n", ch.Data(), shortName.Data(), hshape->Integral() );
            for(int x=0;x<=hshape->GetXaxis()->GetNbins()+1;x++){
              if(hshape->GetXaxis()->GetBinCenter(x)<=minCut || hshape->GetXaxis()->GetBinCenter(x)>=maxCut){ hshape->SetBinContent(x,0); hshape->SetBinError(x,0); }
            }

            if(rebinVal>1){ hshape->Rebin(rebinVal); }
            hshape->GetYaxis()->SetTitle("Entries");// (/25GeV)");
          }
          hshape->Scale(MCRescale);
	  if (postfit) {
	    if(!(ch.Contains("_A_"))) {
	      printf("W/Top NORMALIZATIONs: Process = %s and channel = %s\n\n",proc.Data(),ch.Data());

	      // Top normalization
	      if ( (proc.Contains("t#bar{t} + b#bar{b}")!=std::string::npos) ||
		   (proc.Contains("t#bar{t} + c#bar{c}")!=std::string::npos) ){   
		if (ch.Contains("e")) {  
		  if(ch.Contains("3b")) hshape->Scale(norm_top);
		  if(ch.Contains("4b")) hshape->Scale(norm_top);   
		}
		if (ch.Contains("mu")) { 
		 if(ch.Contains("3b")) hshape->Scale(norm_top);
		 if(ch.Contains("4b")) hshape->Scale(norm_top);   
		}    
	      }
	    
	    // W normalization
	      if (proc.Contains("W#rightarrow l#nu")!=std::string::npos) {
		if (ch.Contains("e")) {  
		  if(ch.Contains("3b")) hshape->Scale(enorm_3b_w);    
		  if(ch.Contains("4b")) hshape->Scale(enorm_4b_w);     
		}
		if (ch.Contains("mu")) {          
		  if(ch.Contains("3b")) hshape->Scale(munorm_3b_w);      
		  if(ch.Contains("4b")) hshape->Scale(munorm_4b_w);          
		}
	      }
	    }
	  }
          if(isSignal)hshape->Scale(SignalRescale);

          //Do Renaming and cleaning
          varName.ReplaceAll("down","Down");
          varName.ReplaceAll("up","Up");

          if(varName==""){//does nothing
	    //	  }else if(varName.EndsWith("_jes")){varName.ReplaceAll("_jes","_CMS_scale_j");
	    //	  }else if(varName.BeginsWith("_umet")) { continue; //skip this one for now
          }else if(varName.BeginsWith("_jer")){varName.ReplaceAll("_jer","_CMS_res_j"); // continue;//skip res for now
	  }else if(varName.BeginsWith("_les")){
	    continue; // skip this one for now
	    //	    if(ch.Contains("e"  ))varName.ReplaceAll("_les","_CMS_scale_e");
	    //	    if(ch.Contains("mu"))varName.ReplaceAll("_les","_CMS_scale_m");
          }else if(varName.BeginsWith("_btag"  )){varName.ReplaceAll("_btag","_CMS_eff_b");
	  }else if(varName.BeginsWith("_ctag"  )){varName.ReplaceAll("_ctag","_CMS_eff_c");
	  }else if(varName.BeginsWith("_ltag"  )){varName.ReplaceAll("_ltag","_CMS_eff_mistag");
	    //          }else if(varName.BeginsWith("_pu"    )){varName.ReplaceAll("_pu", "_CMS_haa4b_pu");
	    //	  }else if(varName.BeginsWith("_pdf" )){
	    //	    if (proc.Contains("wh")!=std::string::npos) {continue; }
          }else if(varName.BeginsWith("_bnorm"  )){continue; //skip this one
	  }else{ varName="_CMS_haa4b"+varName;}

	  hshape->SetTitle(proc+varName);
	  if(shapeInfo.uncShape.find(varName.Data())==shapeInfo.uncShape.end()){
	    shapeInfo.uncShape[varName.Data()] = hshape;
	  }else{
	    shapeInfo.uncShape[varName.Data()]->Add(hshape);
	  }
        }
      }
    }
  }

  //
  // Rebin histograms to make sure that high BDT region have no empty bins
  //
  void AllInfo_t::rebinMainHisto(string histoName)
  {
    //Loop on processes and channels
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          TH1* histo = unc->second;
          if(!histo)continue;
	    TString jetBin = ch->second.bin.c_str();

	    //	    printf("Now going to rebinMainHisto of: %s\n",jetBin.Data());

	    if(jetBin.Contains("3b")){
	      //	      double xbins[] = {-0.30, -0.18, -0.06, 0.06, 0.18, 0.30};  
	      double xbins[] = {-0.30, -0.18, 0.0, 0.16, 0.2, 0.30};     
	      int nbins=sizeof(xbins)/sizeof(double);    
	      unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);  
	      utils::root::fixExtremities(unc->second, false, true); 
	    }else if(jetBin.Contains("4b")){ 
	      double xbins[] = {-0.30, -0.18, 0.0, 0.10, 0.14, 0.30}; 
	      int nbins=sizeof(xbins)/sizeof(double);
	      unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins); 
	      utils::root::fixExtremities(unc->second, false, true); 
	    }

	      /*
          if(jetBin.Contains("vbf")){
            double xbins[] = {150, 225, 300, 375, 450, 600, 750, 1100, 3000};
            int nbins=sizeof(xbins)/sizeof(double);
            unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);
            utils::root::fixExtremities(unc->second, false, true);
          }else if( jetBin.Contains("eq0jets") || jetBin.Contains("geq1jets") ){
            double xbins[] = {150, 225, 300, 375, 450, 525, 600, 725, 850, 975, 1100, 1350, 1600, 2100, 3000};
            int nbins=sizeof(xbins)/sizeof(double);
            unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);
            utils::root::fixExtremities(unc->second, false, true);
          }
	      */

        }
      }
    }
  }

    //
    // merge histograms from different bins together... but keep the channel separated 
    //
    void AllInfo_t::mergeBins(std::vector<string>& binsToMerge, string NewName){
      printf("Merge the following bins of the same channel together: "); for(unsigned int i=0;i<binsToMerge.size();i++){printf("%s ", binsToMerge[i].c_str());}
      printf("The resulting bin will be called %s\n", NewName.c_str());
      for(unsigned int p=0;p<sorted_procs.size();p++){
        string procName = sorted_procs[p];
        std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
        if(it==procs.end())continue;

        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          if(find(binsToMerge.begin(), binsToMerge.end(), ch ->second.bin)==binsToMerge.end())continue;  //make sure this bin should be merged
          for(std::map<string, ChannelInfo_t>::iterator ch2 = ch; ch2!=it->second.channels.end(); ch2++){
            if(ch->second.channel != ch2->second.channel)continue; //make sure we merge bin in the same channel
            if(ch->second.bin     == ch2->second.bin    )continue; //make sure we do not merge with itself
            if(find(binsToMerge.begin(), binsToMerge.end(), ch2->second.bin)==binsToMerge.end())continue;  //make sure this bin should be merged
            addChannel(ch->second, ch2->second); //FIXME this only adds the nominal shapes, not also the syst
            it->second.channels.erase(ch2);  
            ch2=ch;
          }
          ch->second.bin = NewName;
        }

        //also update the map keys
        std::map<string, ChannelInfo_t> newMap;
        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          newMap[ch->second.channel+ch->second.bin] = ch->second;
        }
        it->second.channels = newMap;
      }
    }


    void AllInfo_t::HandleEmptyBins(string histoName){
      for(unsigned int p=0;p<sorted_procs.size();p++){
        string procName = sorted_procs[p];
        //if(procName!="FakeLep")continue; //only do this for the FakeLepbackground right now
        std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
        if(it==procs.end())continue;
        if(it->second.isData && procName!="Instr. MET")continue; //only do this for MC or InstrMET (which also contains MC and negative bins)
        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          ShapeData_t& shapeInfo = ch->second.shapes[histoName];
          TH1* histo = (TH1*)shapeInfo.histo();
          if(!histo){printf("Histo does not exist... skip it \n"); fflush(stdout); continue;}

          double StartIntegral = histo->Integral();
          for(int binx=1;binx<=histo->GetNbinsX();binx++){
            if(histo->GetBinContent(binx)<=0){histo->SetBinContent(binx, 1E-6); histo->SetBinError(binx, 1E-6);  }
          }
          double EndIntegral = histo->Integral();                 
          shapeInfo.rescaleScaleUncertainties(StartIntegral, EndIntegral);


          for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
            for(int binx=1;binx<=unc->second->GetNbinsX();binx++){
              if(unc->second->GetBinContent(binx)<=0){unc->second->SetBinContent(binx, 1E-6); }; //histo->SetBinError(binx, 1.8);  }
	    }
	  }
	}
      }

    //recompute total background
    computeTotalBackground();
  }


