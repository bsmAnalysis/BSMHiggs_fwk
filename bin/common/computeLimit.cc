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

#include "RooFitResult.h"
#include "RooRealVar.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>
#include <regex>

using namespace std;

bool verbose = false ;
float minErrOverSqrtNBGForBinByBin = 0.35 ;

bool autoMCStats = false ;

TString signalSufix="";

TString histo(""), histoVBF("");

int rebinVal = 1;
double MCRescale = 1.0;
double SignalRescale = 1.0;

double datadriven_qcd_Syst = 0.50;    

// Add externally produced systematic shape variations: 
//bool addsyst=false;
// Use postfit normalizations in W/DY/Top bkg components:
bool postfit=false;

double rfr_tt_norm_e ;
double rfr_tt_norm_mu ;
double rfr_w_norm_e ;
double rfr_w_norm_mu ;
double rfr_z_norm_3b_e ;
double rfr_z_norm_3b_mu ;
double rfr_z_norm_4b_e ;
double rfr_z_norm_4b_mu ;

int mass;
bool shape = true;
TString postfix="";
TString systpostfix="";
bool runSystematics = false; 

bool runZh = false;
bool modeDD = false;
bool simfit = false;

TString vh_tag("");
TString year("");

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
bool replaceHighSensitivityBinsWithBG = false ;
bool noCorrelatedStatUnc = false ;
bool plotsOnly = false ;

TString inFileUrl(""),jsonFile("");

TString sumFileUrl("") ;
TString fdInputFile("") ;
TString rfrInputFile("") ;

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
double sstyCut = -1;
bool docut = false;

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



//--------------

void resetNegativeBinsAndErrors( TH1* hp, double val_for_reset = 1., double min_error = 1. ) {
   
   if ( hp == 0x0 ) { printf("\n\n *** resetNegativeBinsAndErrors: null pointer.\n\n") ; gSystem -> Exit(-1) ; }

   for ( int hbi=1; hbi<= hp->GetNbinsX(); hbi++ ) {
      double val, err ;
      val = hp -> GetBinContent( hbi ) ;
      err = hp -> GetBinError( hbi ) ;
      if ( err < min_error ) {
         if ( verbose ) { printf("  resetNegativeBinsAndErrors : hist %s, bin %d, err = %.1f, reset err to %.1f\n", hp->GetName(), hbi, err, min_error ) ; }
         hp -> SetBinError( hbi, min_error ) ;
         err = min_error ;
      }
      if ( val <= 0 ) {
         if ( verbose ) {
            printf("  resetNegativeBinsAndErrors : hist %s, bin %d, val = %.1f, reset val to %.1f and err to %.1f.\n",
             hp->GetName(), hbi, val, val_for_reset, sqrt( pow( err, 2. ) + pow( val, 2. ) ) ) ;
         }
         hp -> SetBinContent( hbi, val_for_reset ) ;
         hp -> SetBinError( hbi, sqrt( pow( err, 2. ) + pow( val, 2. ) ) ) ;
      }
   } // hbi

} // resetNegativeBinsAndErrors

//--------------





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
	//	if(name.Contains("stat") && (name.Contains("Up") || name.Contains("Down"))){  
	if(name.Contains("stat") && (name.EndsWith("Up") || name.EndsWith("Down"))){
          //-- owen: do I need to delete the histogram before calling erase to avoid a memory leak???
          //delete unc->second ;
          uncShape.erase(unc);
          unc--;
        }
        }
        }





      //----------------------------------------------------------------------------------
       void makeStatUnc(string prefix="", string suffix="", string suffix2="", bool noBinByBin=false, TH1* total_hist = 0x0 ){

          if ( verbose ) { printf("\n --- verbose : makeStatUnc : begin.  prefix = %s, suffix = %s, suffix2 = %s\n", prefix.c_str(), suffix.c_str(), suffix2.c_str() ) ; fflush(stdout) ; }

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

                if ( verbose ) { 
                   TString hname( TString(h->GetName())+"StatU"+ibintxt ) ;
                   printf(" --- verbose : makeStatUnc : making clones 1.  name = %s\n", hname.Data() ) ; fflush(stdout) ;
                }


                TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU"+ibintxt);//  statU->Reset();
                TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD"+ibintxt);//  statD->Reset();           
                
                statU->SetBinContent(ibin, std::max(0.0, h_InstrMET_Up_gammaStats->GetBinContent(ibin))>0 ? h_InstrMET_Up_gammaStats->GetBinContent(ibin) : 0.115);
                statD->SetBinContent(ibin, std::max(0.0, h_InstrMET_Down_gammaStats->GetBinContent(ibin))); 

                if ( verbose ) { printf(" --- verbose : makeStatUnc : 1 adding to uncShape with key %s\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Up").c_str() ) ; fflush(stdout) ; }

                uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Up"  ] = statU;
                uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Down"] = statD;
                
              }
              else{
                v_lowStatBin.push_back(ibin);
              }
            }
            
            if ( verbose ) {
               TString hname( TString(h->GetName())+"StatU" ) ;
               printf(" --- verbose : makeStatUnc : making clones 2.  name = %s\n", hname.Data() ) ; fflush(stdout) ;
            }


            TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
            TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");
            
            if(v_lowStatBin.size()>0){
              for(unsigned int j=0; j < v_lowStatBin.size(); j++){
                statU->SetBinContent(v_lowStatBin[j], std::max(0.0, h_InstrMET_Up_gammaStats->GetBinContent(v_lowStatBin[j])));   
                statD->SetBinContent(v_lowStatBin[j], std::max(0.0, h_InstrMET_Down_gammaStats->GetBinContent(v_lowStatBin[j])));   
              }

              if ( verbose ) { printf(" --- verbose : makeStatUnc : 2 adding to uncShape with key %s\n", (prefix+"stat"+suffix+"Up").c_str() ) ; fflush(stdout) ; }

              uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
              uncShape[prefix+"stat"+suffix+"Down"] = statD;
            }   

            delete h; //all done with this copy
            
          }
          else{

         // if ( verbose ) {
         //    printf(" --- verbose: makeStatUnc : %s   %s   %s  histo name: %s\n", prefix.c_str(), suffix.c_str(), suffix2.c_str(), histo()->GetName() ) ;
         //    if ( total_hist != 0x0 ) {
         //       printf("     --- verbose: makeStatUnc : have total_hist histogram,  name = %s,  n bins = %d\n", total_hist -> GetName(), total_hist -> GetNbinsX() ) ;
         //    } else {
         //       printf("     --- verbose: makeStatUnc : total_hist pointer is null.\n") ;
         //    }
         //    fflush(stdout) ;
         // }

            
            TH1* h = (TH1*) histo()->Clone("TMPFORSTAT");
            
            //bin by bin stat uncertainty
            if(statBinByBin>0 && shape==true && !noBinByBin){

              int BIN=0;
              for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++){           

                /////////////if( !(h->GetBinContent(ibin)<=0 && h->GetBinError(ibin)>0) &&  (h->GetBinContent(ibin)<=0 || h->GetBinContent(ibin)/h->Integral()<0.01 || h->GetBinError(ibin)/h->GetBinContent(ibin)<statBinByBin))continue;

                if ( h->GetBinContent(ibin) <= 0 ) continue ;
                if ( total_hist == 0x0 ) {
                   if ( verbose ) { printf("  *** verbose: makeStatUnc :  statBinByBin requested bu no total_hist.  Skipping statBinByBin.\n") ; fflush(stdout) ; }
                   continue ;
                }
                if ( total_hist -> GetNbinsX() != h -> GetNbinsX() ) {
                   printf("\n\n *** makeStatUnc : inconsistent histogram binnings:  %d for this hist, %d for total_hist.\n\n", h -> GetNbinsX(), total_hist -> GetNbinsX() ) ;
                   gSystem -> Exit(-1) ;
                }
                double Nbg = total_hist -> GetBinContent( ibin ) ;
                double err = h->GetBinError( ibin ) ;
                if ( verbose ) { printf("  --- verbose: makeStatUnc :  bin %d,  Nbg = %9.1f, err = %.1f.  ", ibin, Nbg, err ) ; fflush(stdout) ; }
                if ( Nbg <= 0 ) {
                   if ( verbose ) { printf("\n") ; fflush(stdout) ; }
                   continue ;
                }
                if ( verbose ) { printf("  err / sqrt(Nbg) = %7.3f , threshold = %7.3f\n", (err / sqrt(Nbg) ), minErrOverSqrtNBGForBinByBin ) ; fflush(stdout) ; }
                if ( (err / sqrt(Nbg) ) < minErrOverSqrtNBGForBinByBin ) continue ;



                char ibintxt[255]; sprintf(ibintxt, "_b%i", BIN);BIN++;

                if ( verbose ) {
                   TString hname( TString(h->GetName())+"StatU"+ibintxt ) ;
                   printf(" --- verbose : makeStatUnc : making clones 3.  name = %s\n", hname.Data() ) ; fflush(stdout) ;
                }

                TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU"+ibintxt);//  statU->Reset();
                TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD"+ibintxt);//  statD->Reset();           
                if(h->GetBinContent(ibin)>0){
                  statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), h->GetBinContent(ibin) + h->GetBinError(ibin))));   statU->SetBinError(ibin, 0.0);
                  statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), h->GetBinContent(ibin) - h->GetBinError(ibin))));   statD->SetBinError(ibin, 0.0);
                }else{
                  statU->SetBinContent(ibin,std::max(0.0, statU->GetBinContent(ibin) + statU->GetBinError(ibin)));
                  statD->SetBinContent(ibin,std::max(0.0, statD->GetBinContent(ibin) - statD->GetBinError(ibin)));
                }

                if ( verbose ) { printf(" --- verbose: makeStatUnc : setting uncShape[%s] to hist with name %s\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Up").c_str(), statU -> GetName() ) ; fflush(stdout) ; }
                if ( verbose ) { printf(" --- verbose: makeStatUnc : setting uncShape[%s] to hist with name %s\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Down").c_str(), statD -> GetName() ) ; fflush(stdout) ; }

                if ( verbose ) { printf(" --- verbose : makeStatUnc : 3 adding to uncShape with key %s\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Up").c_str() ) ; fflush(stdout) ; }

                uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Up"  ] = statU;
                uncShape[prefix+"stat"+suffix+ibintxt+suffix2+"Down"] = statD;
                /*h->SetBinContent(ibin, 0);*/  h->SetBinError(ibin, 0);  //remove this bin from shape variation for the other ones
                //printf("%s --> %f - %f - %f\n", (prefix+"stat"+suffix+ibintxt+suffix2+"Up").c_str(), statD->Integral(), h->GetBinContent(ibin), statU->Integral() );
              }
            }
            
            //after this line, all bins with large stat uncertainty have been considered separately
            //so now it remains to consider all the other bins for which we assume a total correlation bin by bin
            if(h->Integral()<=0)return; //all non empty bins have already bin variated


            if ( verbose ) {
               TString hname( TString(h->GetName())+"StatU" ) ;
               printf(" --- verbose : makeStatUnc : making clones 4.  name = %s\n", hname.Data() ) ; fflush(stdout) ;
            }

            TH1* statU=(TH1 *)h->Clone(TString(h->GetName())+"StatU");
            TH1* statD=(TH1 *)h->Clone(TString(h->GetName())+"StatD");
            for(int ibin=1; ibin<=statU->GetXaxis()->GetNbins(); ibin++){
              if(h->GetBinContent(ibin)>0){
                statU->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statU->GetBinContent(ibin) + statU->GetBinError(ibin))));
                statD->SetBinContent(ibin,std::min(2*h->GetBinContent(ibin), std::max(0.01*h->GetBinContent(ibin), statD->GetBinContent(ibin) - statD->GetBinError(ibin))));
              }else{
                //statU->SetBinContent(ibin,              statU->GetBinContent(ibin) + statU->GetBinError(ibin));
                statU->SetBinContent(ibin,std::min(0.0, statU->GetBinContent(ibin) + statU->GetBinError(ibin)));
                statD->SetBinContent(ibin,std::min(0.0, statD->GetBinContent(ibin) - statD->GetBinError(ibin)));
              }
            }

            if ( verbose ) { printf(" --- verbose : makeStatUnc : 4 adding to uncShape with key %s\n", (prefix+"stat"+suffix+"Up").c_str() ) ; fflush(stdout) ; }

           //-- owen: first check if this is already set.  If so, delete the existing one to avoid memory leak.
	    
            if ( uncShape.find( prefix+"stat"+suffix+"Up" ) != uncShape.end() ) {
               if ( uncShape[prefix+"stat"+suffix+"Up"  ] != 0x0 ) {
                  if ( verbose ) { printf(" --- verbose : makeStatUnc : 4 deleting existing hist with name %s before assignment.\n", uncShape[prefix+"stat"+suffix+"Up"  ] -> GetName() ) ; fflush(stdout) ; }
                  delete uncShape[prefix+"stat"+suffix+"Up"  ] ;
               }
            }
            if ( uncShape.find( prefix+"stat"+suffix+"Down" ) != uncShape.end() ) {
               if ( uncShape[prefix+"stat"+suffix+"Down"  ] != 0x0 ) {
                  if ( verbose ) { printf(" --- verbose : makeStatUnc : 4 deleting existing hist with name %s before assignment.\n", uncShape[prefix+"stat"+suffix+"Down"  ] -> GetName() ) ; fflush(stdout) ; }
                  delete uncShape[prefix+"stat"+suffix+"Down"  ] ;
               }
            }
	    
            uncShape[prefix+"stat"+suffix+"Up"  ] = statU;
            uncShape[prefix+"stat"+suffix+"Down"] = statD;
            
            delete h; //all done with this copy
          }
        } // makeStatUnc
  


      //----------------------------------------------------------------------------------





        double getScaleUncertainty(){
        double Total=0;
        TH1* h = NULL;
        double integral = 0;
        int start_bin = 1;
        if (this->histo()!=NULL){
          h = (TH1*)(this->histo()->Clone("nominal"));
          integral = h->Integral();
          if(docut) start_bin = h->FindBin(sstyCut);
        }
        for(std::map<string, double>::iterator unc=uncScale.begin();unc!=uncScale.end();unc++){
        if(unc->second<0)continue;
        double unc_val = unc->second;
        if(docut) {
          if(h!=NULL && integral>0) unc_val = unc_val/integral*h->Integral(start_bin, h->GetXaxis()->GetNbins());
          else unc_val = 0;
        }
        Total+=pow(unc_val,2);
//      std::cout << "scale Unc: " << unc->first << ", value: " << unc->second << std::endl;
        }
        return Total>0?sqrt(Total):-1;
        }

        double getIntegratedShapeUncertainty(string name, string upORdown){
        double Total=0;
        //this = ch->second.shapes[histoName.Data()]
        for(std::map<string, TH1*>::iterator var = uncShape.begin(); var!=uncShape.end(); var++){
        TString systName = var->first.c_str();
        if(var->first=="")continue; //Skip Nominal shape
        //if(!systName.Contains(upORdown))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down", buggy code
        if(!systName.EndsWith(upORdown.c_str()))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"

        TH1* hvar = (TH1*)(var->second->Clone((name+var->first).c_str()));

        double varYield = 0.;
        int start_bin = 1;
        if(hvar) {
          if(docut) start_bin = hvar->FindBin(sstyCut);
          varYield = hvar->Integral(start_bin, hvar->GetXaxis()->GetNbins());
        }
        TH1* h = NULL;
        if (this->histo()!=NULL) h = (TH1*)(this->histo()->Clone((name+"Nominal").c_str()));
        double yield = 0.; 
        if (h!=NULL) {
          if(docut) start_bin = h->FindBin(sstyCut);
          yield = h->Integral(start_bin, h->GetXaxis()->GetNbins());
        }
        Total+=pow(varYield-yield,2); //the total shape unc is the sqrt of the quadratical sum of the difference between the nominal and the variated yields.
        }     
        return Total>0?sqrt(Total):-1;
        }

      //------------------------------------

        double getBinShapeUncertainty(string name, int bin, string upORdown){
	  double Total=0;
	  //this = ch->second.shapes[histoName.Data()]
	  for(std::map<string, TH1*>::iterator var = uncShape.begin(); var!=uncShape.end(); var++){
	    TString systName = var->first.c_str();
	    
	    if(var->first=="") continue; //Skip Nominal shape
	    //if(!systName.Contains(upORdown))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"
	    if(!systName.EndsWith(upORdown.c_str()))continue; //only look for syst up or down at a time (upORdown should be either "Up" or "Down"
	    
	    //--owen: aug 15, 2020:  Why are the histograms cloned?  This is super slow.  Is it necessary???  Are they ever used later???
	    //                       Looks safe to not clone, so removing this.  Speeds it up a lot.
	    
            //--- with cloning
	      //	      TH1* hvar = (TH1*)(var->second->Clone((name+var->first).c_str()));
	      //	      double varYield = hvar->GetBinContent(bin);
	      //	      TH1* h = (TH1*)(this->histo()->Clone((name+"Nominal").c_str()));
	      //	      double yield = h->GetBinContent(bin);

            //--- no cloning
	    double varYield = var->second->GetBinContent(bin);
	    double yield = this->histo()->GetBinContent(bin);
	    
	    //	      if(systName.Contains("dydR") && name.find("ee_A_CR_3b")) {
	    if ( verbose ) {
	      printf("--- verbose : getBinShapeUncertainty : name = %s , bin = %d , upORdown = %s , hvar clone name = %s , h clone name = %s, varYield = %.1f, yield = %.1f, diff = %.1f \n",
		     name.c_str(), bin, upORdown.c_str(),
		     (name+var->first).c_str(),
		     (name+"Nominal").c_str(),
		     varYield, yield, (varYield-yield)
		     );
	      fflush(stdout) ;
	    }
	    Total+=pow(varYield-yield,2); 
	    //the total shape unc is the sqrt of the quadratical sum of the difference between the nominal and the variated yields.
	  } // var loop     
	  return Total>0?sqrt(Total):-1;
        } // getBinShapeUncertainty
  
      //------------------------------------

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

        ProcessInfo_t(){xsec=0;isSign=false;}
        ~ProcessInfo_t(){}

        void printProcess() {

           for ( std::map<string, ChannelInfo_t>::iterator ic = channels.begin(); ic!= channels.end(); ic++ ) {

              string chan_key = ic -> first ;
              ChannelInfo_t chan = ic -> second ;

              // printf("       proc key = %s , chan key = %s ,  bin = %s , channel = %s\n", proc_key.c_str(), chan_key.c_str(), chan.bin.c_str(), chan.channel.c_str() ) ;

              for ( std::map<string, ShapeData_t>::iterator is = chan.shapes.begin(); is!= chan.shapes.end(); is++ ) {

                 string shape_key = is -> first ;
                 ShapeData_t shape = is -> second ;

                 //printf("            proc key = %s , chan key = %s , shape key = %s , hist pointer = %p, uncScale has %lu entries, uncShape has %lu entries\n",
                 //  proc_key.c_str(), chan_key.c_str(), shape_key.c_str(), shape.histo(), shape.uncScale.size(), shape.uncShape.size() ) ;

                 int nshapesyst = shape.uncShape.size()-1 ;
                 if ( nshapesyst < 0 ) nshapesyst = 0 ;
                 printf("    printProcess:  %15s :  %18s  : N syst, scale = %2lu, shape = %2d : hist ",  shortName.c_str(), chan_key.c_str(), shape.uncScale.size(), nshapesyst ) ;
                 TH1* hp = shape.histo() ;
                 if ( hp != 0x0 ) {
                    //printf("                  hist name = %s, n bins = %d\n", hp -> GetName(), hp -> GetNbinsX() ) ;
                    printf(" %2d bins | ", hp -> GetNbinsX() ) ;
                    bool is_data_sr = false ;
                    if ( shortName == "data" && chan_key.find("SR")!=string::npos && chan_key.find("emu")==string::npos ) {
                       is_data_sr = true ;
                    }
                    if ( !is_data_sr ) {
                       printf(" entries = %9.1f integral = %9.1f | ", hp -> GetEntries(), hp -> Integral() ) ;
                    } else {
                       printf(" entries = *******.* integral = *******.* | " ) ;
                    }
                    if ( hp -> GetNbinsX() < 10 ) {
                       for ( int bi=1; bi<= hp -> GetNbinsX(); bi++ ) {
                          if ( !(is_data_sr && bi >= 4) ) {
                             printf(" b%d %9.1f |", bi, hp -> GetBinContent( bi ) ) ;
                          } else {
                             printf(" b%d *******.* |", bi ) ;
                          }
                       } // bi
                    }
                 } else {
                    printf(" *** no histogram ***" ) ;
                 }
                 printf("\n") ;
              } // is
           } // ic

        } // printProcess
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
    //////void addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst = false, bool addDiffProcs = true);
    void addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst = false, bool addDiffProcs = true, double scale_factor = 1.);

    // Sum up all background processes and add this as a total process
    ///////void addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst = false);
    void addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst = false, double scale_factor_e = 1., double scale_factor_mu = 1. );

    // Sum up all background processes and add this as a total process
    void computeTotalBackground();

    // Replace the Data process by TotalBackground
    void blind();

    // Replace high sensitivity SR bins for data with total BG.
    void replaceHighSensitivityBinsWithBG();

    // Print the Yield table
    void getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc=NULL);

    // Dump efficiencies
    void getEffFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName);

    // drop background process that have a negligible yield
    void dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold);

    // drop control channels
    void dropCtrlChannels(std::vector<TString>& selCh);

   // Subtract nonQCD MC processes from A,C,D regions in data
    void doBackgroundSubtraction(FILE* pFile, std::vector<TString>& selCh,TString mainHisto, AllInfo_t* sumAllInfo=0x0 );

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

    // Dump to understand data organization
    void printInventory();

};


void printHelp();
void printHelp()
{
  printf("Options\n");
  printf("--verbose   --> turn on a lot of extra printing\n");
  printf("--autoMCStats   --> use Combine implementation of bin-by-bin stat errors on background histograms.  Will turn of statBinByBin.\n");
  printf("--replaceHighSensitivityBinsWithBG  --> replace high-sensitivity histogram bins in SR with total BG.\n") ;
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
  //  printf("--addsyst ---> add more systematics than what was entered in syst with runhaaAnalysis, produced externally \n");
  printf("--signalRescale    --> use this to rescale signal cross-section by a given factor)\n");
  printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
  printf("--minSignalYield   --> use this to specify the minimum Signal yield you want in each channel)\n");
  printf("--signalSufix --> use this flag to specify a suffix string that should be added to the signal 'histo' histogram\n");
  printf("--signalTag   --> use this flag to specify a tag that should be present in signal sample name\n");
  printf("--signalScale   --> use this flag to specify a Scale applied on signal\n");
  printf("--rebin         --> rebin the histogram\n");
  printf("--sstyCut         --> show event yields with bdt above sstyCut\n");
  printf("--statBinByBin --> make bin by bin statistical uncertainty\n");
  printf("--inclusive  --> merge bins to make the analysis inclusive\n");
  printf("--dropBckgBelow --> drop all background processes that contributes for less than a threshold to the total background yields\n");
  printf("--scaleVBF    --> scale VBF signal by ggH/VBF\n");
  printf("--key        --> provide a key for sample filtering in the json\n");  
  printf("--noLogy        --> use this flag to make y-axis linear scale\n");  
  printf("--year        --> use this flag to indicate which year, useful when computing combined limits\n");  
  printf("--minErrOverSqrtNBGForBinByBin  --> Set minimum err / sqrt(NBG) for including a bin-by-bin stat error\n") ;
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

  if ( verbose ) { printf("  --- verbose : main :  processing %d arguments.\n", argc ) ; fflush(stdout) ; }

  //get input arguments
  for(int i=1;i<argc;i++){
    string arg(argv[i]);
    if(arg.find("--help")          !=string::npos) { printHelp(); return -1;} 
    else if(arg.find("--fitDiagnosticsInputFile") !=string::npos && i+1<argc) { fdInputFile = argv[i+1]; i++;  printf("fdInputFile = %s\n", fdInputFile.Data()); }
    else if(arg.find("--allRooFitResultsInputFile") !=string::npos && i+1<argc) { rfrInputFile = argv[i+1]; i++;  printf("rfrInputFile = %s\n", rfrInputFile.Data()); }
    else if(arg.find("--sumInputFile")       !=string::npos && i+1<argc)  { sumFileUrl = argv[i+1];  i++;  printf("sumFileUrl = %s\n", sumFileUrl.Data());  }
    else if(arg.find("--minErrOverSqrtNBGForBinByBin") !=string::npos) { sscanf(argv[i+1],"%f",&minErrOverSqrtNBGForBinByBin); printf("minErrOverSqrtNBGForBinByBin = %.3f\n", minErrOverSqrtNBGForBinByBin);}
    else if(arg.find("--replaceHighSensitivityBinsWithBG") !=string::npos) { replaceHighSensitivityBinsWithBG = true; printf("replaceHighSensitivityBinsWithBG = True\n");}
    else if(arg.find("--noCorrelatedStatUnc") !=string::npos) { noCorrelatedStatUnc = true; printf("noCorrelatedStatUnc = True\n");}
    else if(arg.find("--plotsOnly") !=string::npos) { plotsOnly = true ; printf("plotsOnly = True\n") ; }
    else if(arg.find("--autoMCStats")  !=string::npos) { autoMCStats=true; printf("autoMCStats = True\n");}
    else if(arg.find("--verbose")  !=string::npos) { verbose=true; printf("verbose = True\n");}
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
    else if(arg.find("--year")     !=string::npos && i+1<argc)  { year      = argv[i+1];  i++;  printf("year postfix = %s\n", year.Data()); }
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
    else if(arg.find("--dirtyFix2")    !=string::npos) { dirtyFix2=true; printf("dirtyFix2 = True\n");}
    else if(arg.find("--dirtyFix1")    !=string::npos) { dirtyFix1=true; printf("dirtyFix1 = True\n");}
    else if(arg.find("--signalSufix") !=string::npos) { signalSufix = argv[i+1]; i++; printf("signalSufix '%s' will be used\n", signalSufix.Data()); }
    else if(arg.find("--signalTag") !=string::npos) { signalTag = argv[i+1]; i++; printf("signalTag '%s' will be used\n", signalTag.c_str()); }
    else if(arg.find("--signalScale") !=string::npos) { sscanf(argv[i+1],"%d",&signalScale); i++; printf("signalScale = %d\n", signalScale);}
    else if(arg.find("--rebin")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%i",&rebinVal); i++; printf("rebin = %i\n", rebinVal);}
    else if(arg.find("--sstyCut")    !=string::npos && i+1<argc)  { sscanf(argv[i+1],"%lf",&sstyCut); i++; docut=true; printf("sstyCut = %f\n", sstyCut);}
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
    //    if(arg.find("--addsyst")  !=string::npos) { addsyst=true; printf("addsyst = True\n");}      
  }
  if ( autoMCStats ) {
     if ( statBinByBin > 0 ) {
        printf("\n\n *** WARNING: both autoMCStats and statBinByBin turned on.  Will use autoMCStats and turn off statBinByBin.\n\n") ;
        statBinByBin = -1 ;
     }
  }
  if ( postfit && rfrInputFile.Length() == 0 ) {
     printf("\n\n *** postfit set but no file given with --allRooFitResultsInputFile option.  Rerun with that set.\n\n") ;
     return -1 ;
  }
  if ( postfit && year.Length() == 0 ) {
     printf("\n\n *** postfit set but year not set.  Rerun with --year set.\n\n") ;
     return -1 ;
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
	/*
        Channels.push_back("e_A_CR5j");Channels.push_back("mu_A_CR5j"); // tt+bb CR
        Channels.push_back("e_B_CR5j");Channels.push_back("mu_B_CR5j"); // tt+bb CR 
        Channels.push_back("e_C_CR5j");Channels.push_back("mu_C_CR5j"); // tt+bb CR 
        Channels.push_back("e_D_CR5j");Channels.push_back("mu_D_CR5j"); // tt+bb CR 
	*/
      } else {
        if(runZh){ // Zh
          Channels.push_back("ee_A_CR");Channels.push_back("mumu_A_CR"); // DY CR
          Channels.push_back("emu_A_SR");//Channels.push_back("emu_A_CR"); // Top CR     
        }else{ // Wh
          Channels.push_back("e_A_CR");Channels.push_back("mu_A_CR"); // Top/W CR
          //Channels.push_back("e_A_CR5j");Channels.push_back("mu_A_CR5j"); // tt+bb CR   
        }
      }
    }
  }

  vh_tag = runZh ? "_zh" : "_wh";
  vh_tag = (year == "") ? vh_tag : TString("_") + year + vh_tag;
  year = (year == "") ? "" : TString("_") + year;

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
      //std::cout << "Find the string: " << AnalysisBins[b] << std::endl;
      std::vector<string> subBins;
      std::istringstream iss(AnalysisBins[b]);
      std::string token;
      while (std::getline(iss, token, '+')){
      //char* pch = strtok(&AnalysisBins[b][0],"+"); 
      //while (pch!=NULL){
        //std::cout << "subBin pushed: " << token << std::endl;
        indexcutV.push_back(indexcutV[b]);
        indexcutVL.push_back(indexcutVL[b]);
        indexcutVR.push_back(indexcutVR[b]);
        AnalysisBins.push_back(token);
        subBins.push_back(token);
        //AnalysisBins.push_back(pch);
        //subBins.push_back(pch);
      //  pch = strtok(NULL,"+");
      //}
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
      /*
      ch.push_back("e_A_CR5j"); ch.push_back("mu_A_CR5j");  
      ch.push_back("e_B_CR5j"); ch.push_back("mu_B_CR5j"); 
      ch.push_back("e_C_CR5j"); ch.push_back("mu_C_CR5j"); 
      ch.push_back("e_D_CR5j"); ch.push_back("mu_D_CR5j"); 
      */
    }
  } else {
    if(runZh){ // Zh
      ch.push_back("ee_A_SR"); ch.push_back("mumu_A_SR");     
      if (simfit) {
        ch.push_back("ee_A_CR"); ch.push_back("mumu_A_CR");   
        ch.push_back("emu_A_SR"); //ch.push_back("emu_A_CR");       
      }

    } else { // Wh
      ch.push_back("e_A_SR"); ch.push_back("mu_A_SR");
      if (simfit) { 
        ch.push_back("e_A_CR"); ch.push_back("mu_A_CR");
	//        ch.push_back("e_A_CR5j"); ch.push_back("mu_A_CR5j"); 
      }
    }

  }
  //TString ch[]={"SR"}; //"mumu","ee","emu"};
  const size_t nch=ch.size(); //sizeof(ch)/sizeof(TString);
  std::vector<TString> sh;
  sh.push_back(histo);
  if(subNRB)sh.push_back(histo+"_NRBctrl");
  if(subWZ)sh.push_back(histo+"_3rdLepton");

  if ( verbose ) {
     printf("  --- verbose : main :  contents of ch vector:\n") ;
     for ( int i=0; i<ch.size();           i++ ) { printf("     --- verbose:  ch entry %2d : %s\n", i, ch[i].Data() ) ; }
     printf("  --- verbose : main :  contents of sh vector:\n") ;
     for ( int i=0; i<sh.size();           i++ ) { printf("     --- verbose:  sh entry %2d : %s\n", i, sh[i].Data() ) ; }
     printf("  --- verbose : main :  contents of AnalysisBins vector:\n") ;
     for ( int i=0; i<AnalysisBins.size(); i++ ) { printf("     --- verbose:  AnalysisBins entry %2d : %s\n", i, AnalysisBins[i].c_str() ) ; }
     printf("  --- verbose : main :  contents of Channels vector:\n") ;
     for ( int i=0; i<Channels.size();     i++ ) { printf("     --- verbose:  Channels entry %2d : %s\n", i, Channels[i].Data() ) ; }
     fflush(stdout) ;
  }

  AllInfo_t allInfo;

  AllInfo_t* allInfoSum(0x0) ;
  TFile* inF_sum(0x0) ;
  if ( !sumFileUrl.IsNull() ) {
     printf("\n\n  sumInputFile is set to %s.  Will read it in to a separate allInfo.\n\n", sumFileUrl.Data() ) ;
     allInfoSum = new AllInfo_t() ;
     inF_sum = TFile::Open(sumFileUrl);
     if( !inF_sum || inF_sum->IsZombie() ){ printf("Invalid file name : %s\n", sumFileUrl.Data()); gSystem -> Exit(-1); }
     gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  }


  rfr_tt_norm_e = 1. ;
  rfr_tt_norm_mu = 1. ;
  rfr_w_norm_e = 1. ;
  rfr_w_norm_mu = 1. ;
  rfr_z_norm_3b_e = 1. ;
  rfr_z_norm_3b_mu = 1. ;
  rfr_z_norm_4b_e = 1. ;
  rfr_z_norm_4b_mu = 1. ;

  if ( postfit ) {

     printf("\n\n  postfit it set.  Reading in fit normalizations from %s\n\n", rfrInputFile.Data() ) ;
     TFile tf_fd( rfrInputFile, "read" ) ;
     if ( !(tf_fd.IsOpen()) ) { printf("\n\n *** bad --allRooFitResultsInputFile file %s\n\n", rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
     if ( runZh ) {

        char frname[100] ;
        RooFitResult* rfr(0x0) ;

        RooRealVar* rrv(0x0) ;
        char parname[100] ;

        sprintf( frname, "fit_b_zh_e%s", year.Data() ) ;
        rfr = (RooFitResult*) tf_fd.Get( frname ) ;
        if ( rfr == 0x0 ) {
           printf("\n\n *** postfit set but can't find RooFitResult %s \n\n", frname ) ;
           gSystem -> Exit(-1) ;
        }

        sprintf( parname, "tt_norm_e" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_tt_norm_e = rrv->getVal() ;

        sprintf( parname, "z_norm_3b_e" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_z_norm_3b_e = rrv->getVal() ;

        sprintf( parname, "z_norm_4b_e" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_z_norm_4b_e = rrv->getVal() ;

        rfr->Delete() ;


        sprintf( frname, "fit_b_zh_mu%s", year.Data() ) ;
        rfr = (RooFitResult*) tf_fd.Get( frname ) ;
        if ( rfr == 0x0 ) {
           printf("\n\n *** postfit set but can't find RooFitResult %s \n\n", frname ) ;
           gSystem -> Exit(-1) ;
        }

        sprintf( parname, "tt_norm_mu" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_tt_norm_mu = rrv->getVal() ;

        sprintf( parname, "z_norm_3b_mu" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_z_norm_3b_mu = rrv->getVal() ;

        sprintf( parname, "z_norm_4b_mu" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_z_norm_4b_mu = rrv->getVal() ;

        rfr->Delete() ;


        printf("   postfit normalizations, Zh, e :  tt_norm_e  = %6.3f , z_norm_3b_e  = %6.3f , z_norm_4b_e  = %6.3f\n", rfr_tt_norm_e , rfr_z_norm_3b_e , rfr_z_norm_4b_e  ) ;
        printf("   postfit normalizations, Zh, mu:  tt_norm_mu = %6.3f , z_norm_3b_mu = %6.3f , z_norm_4b_mu = %6.3f\n", rfr_tt_norm_mu, rfr_z_norm_3b_mu, rfr_z_norm_4b_mu ) ;
        fflush(stdout) ;

     } else {

        char frname[100] ;
        sprintf( frname, "fit_b_wh%s", year.Data() ) ;
        RooFitResult* rfr = (RooFitResult*) tf_fd.Get( frname ) ;
        if ( rfr == 0x0 ) {
           printf("\n\n *** postfit set but can't find RooFitResult %s \n\n", frname ) ;
           gSystem -> Exit(-1) ;
        }

        RooRealVar* rrv(0x0) ;
        char parname[100] ;

        sprintf( parname, "tt_norm_e" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_tt_norm_e = rrv->getVal() ;

        sprintf( parname, "tt_norm_mu" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_tt_norm_mu = rrv->getVal() ;
        
        sprintf( parname, "w_norm_e" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_w_norm_e = rrv->getVal() ;

        sprintf( parname, "w_norm_mu" ) ;
        rrv = (RooRealVar*)( rfr -> floatParsFinal()).find( parname ) ;
        if ( rrv == 0x0 ) { printf("\n\n *** postfit set but can't find %s in %s in %s\n\n", parname, frname, rfrInputFile.Data() ) ; gSystem -> Exit(-1) ; }
        rfr_w_norm_mu = rrv->getVal() ;

        rfr -> Delete() ;

        printf("   postfit normalizations, Wh:  tt_norm_e = %6.3f , tt_norm_mu = %6.3f , w_norm_e = %6.3f , w_norm_mu = %6.3f\n",
          rfr_tt_norm_e, rfr_tt_norm_mu, rfr_w_norm_e, rfr_w_norm_mu ) ;
        fflush(stdout) ;

     }

     tf_fd.Close() ;

  }



  if ( verbose ) { printf("  --- verbose : main :  Opening input root file with name inFileUrl = %s\n", inFileUrl.Data() ) ; fflush(stdout) ; }


  //open input file
  TFile* inF = TFile::Open(inFileUrl);
  if( !inF || inF->IsZombie() ){ printf("Invalid file name : %s\n", inFileUrl.Data());}
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE


  //LOAD shapes
  const size_t nsh=sh.size();
  for(size_t b=0; b<AnalysisBins.size(); b++){
    std::vector<string> channelsAndShapes;
    std::vector<string> channelsAndShapesSum;
    for(size_t i=0; i<nch; i++){
      for(size_t j=0; j<nsh; j++){
        channelsAndShapes.push_back((ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
        printf("allInfo   : Adding shape %s\n",(ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
        if ( !sumFileUrl.IsNull() && inF_sum!=0x0 ) {
          //--- only need SR 4b for sum (used in DD QCD).
           if ( ch[i].Contains("SR") && strcmp( AnalysisBins[b].c_str(), "4b" ) == 0 ) {
              channelsAndShapesSum.push_back((ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
              printf("allInfoSum: Adding shape %s\n",(ch[i]+TString(";")+AnalysisBins[b]+TString(";")+sh[j]).Data());
           }
        }
      }
    }
    double cutMin=shapeMin; double cutMax=shapeMax;
    allInfo.getShapeFromFile(inF, channelsAndShapes, indexcutM[AnalysisBins[b]], Root, cutMin, cutMax   );     
    if ( !sumFileUrl.IsNull() && inF_sum!=0x0 ) allInfoSum -> getShapeFromFile(inF_sum, channelsAndShapesSum, indexcutM[AnalysisBins[b]], Root, cutMin, cutMax   );     
  }




  inF->Close();
  printf("Loading all shapes... Done\n");
  fflush(stdout) ;

  for(unsigned int B=0;B<binsToMerge.size();B++){
    std::string NewBinName = binsToMerge[B][0]; std::cout << "binsToMerge[B][0]: " << binsToMerge[B][0]; for(unsigned int b=1;b<binsToMerge[B].size();b++){NewBinName += "_"+binsToMerge[B][b];std::cout << "binsToMerge[B][b]: " << binsToMerge[B][b] << std::endl;;}
//    std::string NewBinName = string("["); binsToMerge[B][0];  for(unsigned int b=1;b<binsToMerge[B].size();b++){NewBinName += "+"+binsToMerge[B][b];} NewBinName+="]";
    allInfo.mergeBins(binsToMerge[B],NewBinName);
    if ( !sumFileUrl.IsNull() && allInfoSum != 0x0 ) allInfoSum -> mergeBins(binsToMerge[B],NewBinName);
  }


  if ( verbose ) {
     printf("\n\n --- verbose : main :  calling allInfo.printInventory for main allInfo\n\n") ;
     allInfo.printInventory() ;
     if ( !sumFileUrl.IsNull() && allInfoSum != 0x0 ) {
        printf("\n\n --- verbose : main :  calling allInfo.printInventory for sum allInfo\n\n") ;
        allInfoSum -> printInventory() ;
     }
     fflush(stdout) ;
  }

  
  if ( verbose ) { printf("\n --- verbose : main :  calling allInfo.computeTotalBackground() for first time.\n") ; fflush(stdout) ; }
  allInfo.computeTotalBackground();
  if(MCclosureTest)allInfo.blind();


  if ( verbose ) allInfo.printInventory() ;




  if ( verbose ) { if (shape && BackExtrapol ) printf("\n  --- verbose : main :  calling allInfo.rebinMainHisto(histo.Data()) where histo = %s\n", histo.Data() ) ; fflush(stdout) ; }


  //extrapolate backgrounds toward higher BDT region to make sure that there is no empty bins
  //
  if(shape && BackExtrapol)allInfo.rebinMainHisto(histo.Data());

  if ( shape && BackExtrapol && allInfoSum != 0x0 ) {
     if ( verbose ) { printf("\n --- verbose : main :  calling rebinMainHisto(histo.Data()) for allInfoSum.\n") ; fflush(stdout) ; }
     allInfoSum -> rebinMainHisto(histo.Data());
  }


  if ( verbose ) allInfo.printInventory() ;



  FILE* pFile;

  //define vector for search
  std::vector<TString>& selCh = Channels;

  if(modeDD) {
    if ( allInfoSum != 0x0 ) {
       if ( verbose ) { printf("\n --- verbose : main :  calling doBackgroundSubtraction for sum allInfo first.\n") ; }
       pFile = fopen("datadriven_qcd"+year+"-sum.tex","w");
       if(subFake) allInfoSum -> doBackgroundSubtraction(pFile,selCh,histo);
       fclose(pFile);
    }
    if ( verbose ) { printf("\n --- verbose : main :  calling doBackgroundSubtraction for main allInfo.\n") ; }
    pFile = fopen("datadriven_qcd"+year+".tex","w");
    if(subFake)allInfo.doBackgroundSubtraction(pFile,selCh,histo, allInfoSum);
    fclose(pFile);
  }

  //replace data by total MC background
  //if(blindData)allInfo.blind();

  fflush(stdout) ;




  if ( verbose ) { printf("\n  --- verbose : main :  calling allInfo.dropSmallBckgProc(selCh, histo.Data(), dropBckgBelow) , \n") ; fflush(stdout) ; }

  //drop backgrounds with rate<1%
  allInfo.dropSmallBckgProc(selCh, histo.Data(), dropBckgBelow);




  if ( verbose ) { printf("\n  --- verbose : main :   calling allInfo.dropCtrlChannels(selCh);\n") ; fflush(stdout) ; }

  //drop control channels
  allInfo.dropCtrlChannels(selCh);



  //merge bins  
//  for(unsigned int B=0;B<binsToMerge.size();B++){
//    std::string NewBinName = string("["); binsToMerge[B][0];  for(unsigned int b=1;b<binsToMerge[B].size();b++){NewBinName += "+"+binsToMerge[B][b];} NewBinName+="]";
//    allInfo.mergeBins(binsToMerge[B],NewBinName);
//  }



  if ( verbose && !shape ) { printf("\n  --- verbose : main :  calling allInfo.turnToCC(histo.Data()); \n") ; fflush(stdout) ; }

  //turn to CC analysis eventually
  if(!shape)allInfo.turnToCC(histo.Data());





  if ( verbose ) { printf("\n  --- verbose : main :  calling allInfo.HandleEmptyBins(histo.Data()); \n") ; fflush(stdout) ; }

  allInfo.HandleEmptyBins(histo.Data()); //needed for negative bin content --> May happens due to NLO interference for instance





  if ( verbose && blindData ) { printf("\n  --- verbose : main :  calling allInfo.blind(); \n") ; fflush(stdout) ; }

  // Blind data in Signal Regions only
  if(blindData) allInfo.blind();

  if (replaceHighSensitivityBinsWithBG) allInfo.replaceHighSensitivityBinsWithBG();




  if ( verbose ) { printf("\n  --- verbose : main :   calling       allInfo.getEffFromShape(pFile, selCh, histo.Data()); \n") ; fflush(stdout) ; }

  //print signal efficiency
  pFile = fopen(vh_tag+"Efficiency.tex","w");
  allInfo.getEffFromShape(pFile, selCh, histo.Data());
  fclose(pFile);



  if ( verbose ) { printf("\n  --- verbose : main :    calling   allInfo.addHardCodedUncertainties(histo.Data()); \n") ; fflush(stdout) ; }

  //add by hand the hard coded uncertainties
  allInfo.addHardCodedUncertainties(histo.Data());



  if ( verbose ) { printf("\n  --- verbose : main :    calling allInfo.getYieldsFromShape(pFile, selCh, histo.Data(), pFileInc); \n") ; fflush(stdout) ; }

  //print event yields from the histo shapes
  pFile = fopen(vh_tag+"Yields.tex","w");  FILE* pFileInc = fopen(vh_tag+"YieldsInc.tex","w");
  allInfo.getYieldsFromShape(pFile, selCh, histo.Data(), pFileInc);
  fclose(pFile); fclose(pFileInc);
  


  if ( verbose ) { printf("\n  --- verbose : main :    calling   allInfo.showShape(selCh,histo,\"plot\"); \n") ; fflush(stdout) ; }

  //produce a plot
  allInfo.showShape(selCh,histo,"plot"); //this produce the final global shape
  
  if ( plotsOnly ) {
     printf("\n\n plotsOnly is set to true.  Bailing out now.\n\n") ; fflush(stdout) ;
     return 0 ;
  }


  if ( verbose && runSystematics ) { printf("\n  --- verbose : main :    calling allInfo.showUncertainty(selCh,histo,\"plot\"); \n") ; fflush(stdout) ; }

  //produce a plot
  //  if(runSystematics) allInfo.showUncertainty(selCh,histo,"plot"); //this produces all the plots with the syst  
  if(runSystematics && !(simfit)) allInfo.showUncertainty(selCh,histo,"plot"); //this produces all the plots with the syst  
  // georgia : now run reporting systematics only if simfit=false 
  // owen: temporarily turn this off.  Slows it down.
  
  if ( verbose ) allInfo.printInventory() ;



  //prepare the output
  string limitFile=("haa4b_"+massStr+systpostfix+vh_tag+".root").Data();
  TFile *fout=TFile::Open(limitFile.c_str(),"recreate");



  if ( verbose ) { printf("\n  --- verbose : main :   calling   allInfo.saveHistoForLimit(histo.Data(), fout);  \n") ; fflush(stdout) ; }

  allInfo.saveHistoForLimit(histo.Data(), fout);



  if ( verbose ) { printf("\n  --- verbose : main :    calling   allInfo.buildDataCards(histo.Data(), limitFile); \n") ; fflush(stdout) ; }

  allInfo.buildDataCards(histo.Data(), limitFile);

  //all done
  printf("\n\n calling fout->Close();\n\n") ; fflush(stdout) ;
  fout->Close();
  printf("\n\n At the end of main.\n\n") ; fflush(stdout) ;
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
void AllInfo_t::addChannel(ChannelInfo_t& dest, ChannelInfo_t& src, bool computeSyst, bool addDiffProcs, double scale_factor ){
  std::map<string, ShapeData_t>& shapesInfoDest = dest.shapes;
  std::map<string, ShapeData_t>& shapesInfoSrc  = src.shapes;

  if(!computeSyst){
    for(std::map<string, ShapeData_t>::iterator sh = shapesInfoSrc.begin(); sh!=shapesInfoSrc.end(); sh++){
      if(shapesInfoDest.find(sh->first)==shapesInfoDest.end())shapesInfoDest[sh->first] = ShapeData_t();
      
      //Loop on all shape systematics (including also the central value shape)
      for(std::map<string, TH1*>::iterator uncS = sh->second.uncShape.begin();uncS!= sh->second.uncShape.end();uncS++){
        if(uncS->first!="" || uncS->second == NULL) continue; //We only take nominal shapes
        if(shapesInfoDest[sh->first].uncShape.find(uncS->first)==shapesInfoDest[sh->first].uncShape.end()){
          shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) uncS->second->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
        }else{
          //////////shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second);
          shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second, scale_factor );
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
        if(uncS->first=="" || uncS->second == NULL) continue; //We only take systematic (i.e non-nominal) shapes
        if(shapesInfoSrc[sh->first].uncShape.find("")==shapesInfoSrc[sh->first].uncShape.end()) continue;
        if(addDiffProcs){ // add different procs
	    //1. Copy the nominal shape    
	  shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) shapesInfoDest[sh->first].uncShape[""]->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
          //2. we remove the nominal value of the process we are running on
	  shapesInfoDest[sh->first].uncShape[uncS->first]->Add(shapesInfoSrc[sh->first].uncShape[""], -1);
          //3. and add the variation up/down
          shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second, scale_factor ); 
        }else{ // add same proc in different channels
          if(shapesInfoDest[sh->first].uncShape.find(uncS->first)==shapesInfoDest[sh->first].uncShape.end()){
            shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) uncS->second->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
          }else{
            ////////shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second); 
            shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second, scale_factor ); 
          }
        }
      }
    }
  } // if (computeSyst)

}

//
// Sum up all background processes and add this as a total process
//
void AllInfo_t::addProc(ProcessInfo_t& dest, ProcessInfo_t& src, bool computeSyst, double scale_factor_e, double scale_factor_mu ){
  dest.xsec = src.xsec*src.br;
  for(std::map<string, ChannelInfo_t>::iterator ch = src.channels.begin(); ch!=src.channels.end(); ch++){
    if(dest.channels.find(ch->first)==dest.channels.end()){   //this channel does not exist, create it
      dest.channels[ch->first]         = ChannelInfo_t();
      dest.channels[ch->first].bin     = ch->second.bin;
      dest.channels[ch->first].channel = ch->second.channel;
    }

    //addChannel(dest.channels[ch->first], ch->second, computeSyst);
    TString ts_channel_name( ch->first ) ;
    if ( ts_channel_name.Contains( "mu_" ) ) {
       addChannel(dest.channels[ch->first], ch->second, computeSyst, true, scale_factor_mu );
    } else if ( ts_channel_name.Contains( "e_" ) ) {
       addChannel(dest.channels[ch->first], ch->second, computeSyst, true, scale_factor_e );
    } else {
       addChannel(dest.channels[ch->first], ch->second, computeSyst, true, 1.0 );
    }

  }
}

//
// Subtract nonQCD MC from A,C,D regions in QCD Analysis
//
void AllInfo_t::doBackgroundSubtraction(FILE* pFile,std::vector<TString>& selCh,TString mainHisto, AllInfo_t* sumAllInfo) {

  if ( verbose ) { printf("\n\n --- verbose :  AllInfo_t::doBackgroundSubtraction : begin\n\n") ; }

  if ( verbose ) {
     if ( sumAllInfo != 0x0 ) {
        printf("\n\n --- verbose :  AllInfo_t::doBackgroundSubtraction :  sumAllInfo is set.  Will use sum for 4b SR shape for B and C/D ratio.\n\n") ;
     } else {
        printf("\n\n --- verbose :  AllInfo_t::doBackgroundSubtraction :  sumAllInfo is NOT set.\n\n") ;
     }
     if ( fdInputFile.Length() > 0 ) {
        printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  fdInputFile is set to %s.  Will apply scale factors from fit_b from that file.\n", fdInputFile.Data() ) ;
     }
     printf("  --- verbose :  AllInfo_t::doBackgroundSubtraction :  contents of selCh vector:\n") ;
     for ( int i = 0; i<selCh.size(); i++ ) {
        printf("   %3d : %s\n", i, selCh[i].Data() ) ;
     } // i
  }

  // Closure test for DD QCD predictions
  char Lcol     [1024] = "|c";
  char Lchan    [1024] = "";
  char Lalph1   [1024] = "";
  char Lalph2   [1024] = "";
  char Lyield [1024] = "";
  char LyieldMC [1024] = "";
  // MC closure test for QCD predictions
  char LalphMC  [1024] = "";
  char LyieldPred[1024] = "";
  char RatioMC  [1024] = "";

  double w_norm_e_val(1.0) ;
  double w_norm_mu_val(1.0) ;
  double tt_norm_e_val(1.0) ;
  double tt_norm_mu_val(1.0) ;

  double w_norm_e_err(0.0) ;
  double w_norm_mu_err(0.0) ;
  double tt_norm_e_err(0.0) ;
  double tt_norm_mu_err(0.0) ;

  if ( fdInputFile.Length() > 0 ) {

     TFile tf_fd( fdInputFile, "READ" ) ;
     if ( !(tf_fd.IsOpen()) ) {
        printf("\n\n *** AllInfo_t::doBackgroundSubtraction : fdInputFile set to %s.  problem opening this file.  I quit.\n\n", fdInputFile.Data() ) ;
        gSystem -> Exit(-1) ;
     }

     RooFitResult* rfr = (RooFitResult*) tf_fd.Get( "fit_b" ) ;
     if ( rfr == 0x0 ) { printf("\n\n *** did not find fit_b RooFitResult in %s.  I quit.\n\n", fdInputFile.Data() ) ; gSystem -> Exit(-1) ; }

     RooRealVar* rrv_w_e = (RooRealVar*)(rfr -> floatParsFinal()).find( "w_norm_e" ) ;
     if ( rrv_w_e == 0x0 ) { printf("\n\n *** did not find w_norm_e in fit_b in %s.  I quit.\n\n", fdInputFile.Data() ) ; gSystem -> Exit(-1) ; }

     RooRealVar* rrv_w_mu = (RooRealVar*)(rfr -> floatParsFinal()).find( "w_norm_mu" ) ;
     if ( rrv_w_mu == 0x0 ) { printf("\n\n *** did not find w_norm_mu in fit_b in %s.  I quit.\n\n", fdInputFile.Data() ) ; gSystem -> Exit(-1) ; }

     RooRealVar* rrv_tt_e = (RooRealVar*)(rfr -> floatParsFinal()).find( "tt_norm_e" ) ;
     if ( rrv_tt_e == 0x0 ) { printf("\n\n *** did not find tt_norm_e in fit_b in %s.  I quit.\n\n", fdInputFile.Data() ) ; gSystem -> Exit(-1) ; }

     RooRealVar* rrv_tt_mu = (RooRealVar*)(rfr -> floatParsFinal()).find( "tt_norm_mu" ) ;
     if ( rrv_tt_mu == 0x0 ) { printf("\n\n *** did not find tt_norm_mu in fit_b in %s.  I quit.\n\n", fdInputFile.Data() ) ; gSystem -> Exit(-1) ; }

     w_norm_e_val = rrv_w_e -> getVal() ;
     w_norm_e_err = rrv_w_e -> getError() ;

     w_norm_mu_val = rrv_w_mu -> getVal() ;
     w_norm_mu_err = rrv_w_mu -> getError() ;

     tt_norm_e_val = rrv_tt_e -> getVal() ;
     tt_norm_e_err = rrv_tt_e -> getError() ;

     tt_norm_mu_val = rrv_tt_mu -> getVal() ;
     tt_norm_mu_err = rrv_tt_mu -> getError() ;

     if ( verbose ) {
         printf( "\n\n --- verbose :  AllInfo_t::doBackgroundSubtraction :  post-fit scale factors from %s\n", fdInputFile.Data() ) ;
         printf( "     w_norm_e = %6.3f +/- %6.3f ,    w_norm_mu = %6.3f +/- %6.3f\n", w_norm_e_val, w_norm_e_err, w_norm_mu_val, w_norm_mu_err ) ;
         printf( "    tt_norm_e = %6.3f +/- %6.3f ,   tt_norm_mu = %6.3f +/- %6.3f\n", tt_norm_e_val, tt_norm_e_err, tt_norm_mu_val, tt_norm_mu_err ) ;
         printf("\n") ;
     }

     rfr -> Delete() ;

     tf_fd.Close() ;

  }

  //check that the data proc exist
  std::map<string, ProcessInfo_t>::iterator dataProcIt=procs.find("data");             
  if(dataProcIt==procs.end()){printf("The process 'data' was not found... can not do QCD background prediction\n"); return;}

  // create 3 new processes for A,C,D regions in data
  TString NRBProcName = "NonQCD";
  for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==NRBProcName.Data()){sorted_procs.erase(p);}} 
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
        printf("Subtracting nonQCD process from data: %s, long name %s \n", it->second.shortName.c_str(), procName.Data() ); 
    if ( fdInputFile.Length() > 0 ) {
       if ( strcmp( it->second.shortName.c_str(), "wlnu" ) == 0 ) {
          if ( verbose ) { printf("  applying w_norm_e and w_norm_mu\n") ; }
          addProc(procInfo_NRB, it->second, false, w_norm_e_val, w_norm_mu_val );
       } else if ( strcmp( it->second.shortName.c_str(), "ttbarbba" ) == 0 ) {
          if ( verbose ) { printf("  applying tt_norm_e and tt_norm_mu\n") ; }
          addProc(procInfo_NRB, it->second, false, tt_norm_e_val, tt_norm_mu_val );
       } else if ( strcmp( it->second.shortName.c_str(), "ttbarcba" ) == 0 ) {
          if ( verbose ) { printf("  applying tt_norm_e and tt_norm_mu\n") ; }
          addProc(procInfo_NRB, it->second, false, tt_norm_e_val, tt_norm_mu_val );
       } else {
          addProc(procInfo_NRB, it->second, false );
       }
    } else {
       addProc(procInfo_NRB, it->second, false);
    }
  }

  if ( verbose ) { printf("  --- verbose :  AllInfo_t::doBackgroundSubtraction :  proc for all non-QCD backgrounds:\n") ; procInfo_NRB.printProcess() ; }



  for(std::map<string, ChannelInfo_t>::iterator chData = dataProcIt->second.channels.begin(); chData!=dataProcIt->second.channels.end(); chData++){
    if(std::find(selCh.begin(), selCh.end(), chData->second.channel)==selCh.end())continue;


    // if(chData->first.find("CR5j"))continue;
    if(chData->first.find("_A_")!=string::npos)continue; // do not subtract nonQCD MC in regions A...

    if ( verbose ) { printf("  --- verbose :  AllInfo_t::doBackgroundSubtraction :  loop over data channels :  %s\n", chData->first.c_str() ) ; }

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
        //      if(shapesInfoDest[sh->first].uncShape.find(uncS->first)==shapesInfoDest[sh->first].uncShape.end()){
        //        shapesInfoDest[sh->first].uncShape[uncS->first] = (TH1*) uncS->second->Clone(TString(uncS->second->GetName() + dest.channel + dest.bin ) );
        //      }else{
        //      printf("Here start subtracting NonQCD %s \n",uncS->second->GetName());
        
        if ( verbose ) { printf("  --- verbose :  AllInfo_t::doBackgroundSubtraction :  subtracting %s from %s\n", uncS->second->GetName(), shapesInfoDest[sh->first].uncShape[uncS->first]->GetName() ) ; }
        shapesInfoDest[sh->first].uncShape[uncS->first]->Add(uncS->second,-1);
          //    }
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

    if ( verbose ) { printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  ChannelInfo_t loop, chData->first = %s, chData->second.channel = %s\n", chData->first.c_str(), chData->second.channel.c_str() ) ; }

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
    TH1* hCtrl_SB = dataProcIt->second.channels[(binName+"_"+chData->second.bin.c_str()+year).Data()].shapes[mainHisto.Data()].histo(); // Region D  
    binName.ReplaceAll("D_","B_");    
    TH1* hCtrl_SI = dataProcIt->second.channels[(binName+"_"+chData->second.bin.c_str()+year).Data()].shapes[mainHisto.Data()].histo();  // Region B 
    binName.ReplaceAll("B_","C_");       
    TH1* hChan_SB = dataProcIt->second.channels[(binName+"_"+chData->second.bin.c_str()+year).Data()].shapes[mainHisto.Data()].histo(); // Region C

    TH1* hDD     =  chDD->second.shapes[mainHisto.Data()].histo(); // Region B
    if(hDD==NULL){std::cout << "hDD does not exist:" << chDD->second.bin << "_" << chDD->second.channel << mainHisto.Data() << std::endl;}
    
    // load MC histograms in the QCD control regions
    //    binName.ReplaceAll("C_","B_");
    TH1* hDD_C = procInfo_DD.channels.find((binName+"_"+chData->second.bin.c_str()+year).Data())->second.shapes[mainHisto.Data()].histo(); 
    binName.ReplaceAll("C_","D_");  
    TH1* hDD_D = procInfo_DD.channels.find((binName+"_"+chData->second.bin.c_str()+year).Data())->second.shapes[mainHisto.Data()].histo();
    binName.ReplaceAll("D_","B_");    
    TH1* hDD_B = procInfo_DD.channels.find((binName+"_"+chData->second.bin.c_str()+year).Data())->second.shapes[mainHisto.Data()].histo();  

    if ( verbose ) {

       double val, err ;

       val = hChan_SB -> IntegralAndError( 1,hChan_SB->GetXaxis()->GetNbins(), err ) ;
       printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hChan_SB (C) :  integral = %9.1f +/- %6.1f | \n", val, err ) ;

       val = hCtrl_SB -> IntegralAndError( 1,hCtrl_SB->GetXaxis()->GetNbins(), err ) ;
       printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hCtrl_SB (D) :  integral = %9.1f +/- %6.1f | \n", val, err ) ;

       val = hCtrl_SI -> IntegralAndError( 1,hCtrl_SI->GetXaxis()->GetNbins(), err ) ;
       printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hCtrl_SI (B) :  integral = %9.1f +/- %6.1f | ", val, err ) ;
       for ( int i=1; i<=hCtrl_SI->GetNbinsX(); i++ ) { printf(" bin %2d = %7.1f +/- %7.1f | ", i, hCtrl_SI->GetBinContent( i ), hCtrl_SI->GetBinError( i ) ) ; }
       printf("\n") ;

    }
    
    
    binName.ReplaceAll("B_","");   

    //compute alpha
    double alpha=0 ,alpha_err=0;
    double alphaMC=0, alphaMC_err=0;

    double errC,errD;
    double errMC_C, errMC_D;
    double valMC=0, valMC_err=0;
    double valDD_MC=0, valDD_MC_err=0;
    double ratioMC=0, ratioMC_err=0;


   //-- For every bin in B, C, and D, check if it's negative.  If it is, set error to sqrt( err^2 + val^2 ) and value to 1 event.
    resetNegativeBinsAndErrors( hChan_SB, 1. ) ; // C
    resetNegativeBinsAndErrors( hCtrl_SB, 1. ) ; // D
    resetNegativeBinsAndErrors( hCtrl_SI, 0. ) ; // B

 
    if(hCtrl_SB->Integral()>0){
      alpha     = hChan_SB->IntegralAndError(1,hChan_SB->GetXaxis()->GetNbins(),errC) / hCtrl_SB->IntegralAndError(1,hCtrl_SB->GetXaxis()->GetNbins(),errD);
      alpha_err = ( fabs( hChan_SB->Integral() * errD ) + fabs(errC * errD )  ) / pow(hCtrl_SB->Integral(), 2);        
    }
    


    //-- make it compile clean.
    ////if(hDD!=NULL) valMC = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins(),valMC_err); if(valMC<1E-6){valMC=0.0; valMC_err=0.0;}   
    if(hDD!=NULL) valMC = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins(),valMC_err);
    if(valMC<1E-6){valMC=0.0; valMC_err=0.0;}   


    if(hDD && hDD_B && hDD_C && hDD_D){
      // alpha in MC
      if(hDD_D!=NULL && hDD_D->Integral()>0){
        alphaMC = hDD_C->IntegralAndError(1,hDD_C->GetXaxis()->GetNbins(),errMC_C) / hDD_D->IntegralAndError(1,hDD_D->GetXaxis()->GetNbins(),errMC_D);
        alphaMC_err =  ( fabs( hDD_C->Integral() * errMC_D ) + fabs(errMC_C * errMC_D )  ) / pow(hDD_D->Integral(), 2);  
      }
    

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
    if ( sumAllInfo != 0x0 && ( startsWith( chData->first, "e_A_SR_4b") || startsWith( chData->first, "mu_A_SR_4b")  ) ) {

       printf("\n sumAllInfo set and this is a SR 4b channel: %s.  Will use sum for C/D ratio and B shape.\n", chData->first.c_str() ) ;

       std::map<string, ProcessInfo_t>::iterator iSumDataProc = (*sumAllInfo).procs.find("data") ;
       if ( iSumDataProc==(*sumAllInfo).procs.end() ) { printf("\n\n *** can't find data proc in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

       ProcessInfo_t sumDataProc = iSumDataProc -> second ;


      //--- B
       
       TH1* h_B(0x0) ;
       {
          string ch_name ;
          ch_name = "e_B_SR_4b" ;
          if ( year != "" ) { ch_name += year ; }
          std::map<string, ChannelInfo_t>::iterator ichan_e = sumDataProc.channels.find( ch_name ) ;
          if ( ichan_e == sumDataProc.channels.end() ) { printf("\n\n *** can't find data %s channel in sumAllInfo.  I quit.\n\n", ch_name.c_str() ) ; gSystem -> Exit(-1); }
          ChannelInfo_t chan_e = ichan_e->second ;

          ch_name = "mu_B_SR_4b" ;
          if ( year != "" ) { ch_name += year ; }
          std::map<string, ChannelInfo_t>::iterator ichan_mu = sumDataProc.channels.find( ch_name ) ;
          if ( ichan_mu == sumDataProc.channels.end() ) { printf("\n\n *** can't find data mu_B_SR_4b channel in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ChannelInfo_t chan_mu = ichan_mu->second ;

          std::map<string, ShapeData_t>::iterator ishape_e = chan_e.shapes.find( mainHisto.Data() ) ;
          if ( ishape_e == chan_e.shapes.end() ) { printf("\n\n *** can't find data e_B_SR_4b channel shape in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ShapeData_t shape_e = ishape_e->second ;
          TH1* h_B_e = shape_e.histo() ;
          if ( h_B_e == 0x0 ) { printf("\n\n *** can't find data e_B_SR_4b channel shape histogram in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

          std::map<string, ShapeData_t>::iterator ishape_mu = chan_mu.shapes.find( mainHisto.Data() ) ;
          if ( ishape_mu == chan_mu.shapes.end() ) { printf("\n\n *** can't find data mu_B_SR_4b channel shape in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ShapeData_t shape_mu = ishape_mu->second ;
          TH1* h_B_mu = shape_mu.histo() ;
          if ( h_B_mu == 0x0 ) { printf("\n\n *** can't find data mu_B_SR_4b channel shape histogram in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

          h_B = (TH1*) h_B_e -> Clone("h_B_from_sumAllInfo") ;
          h_B -> Add( h_B_mu ) ;

          if ( verbose ) {
             printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hist from sumAllInfo for  e_B_SR_4b :  integral = %9.1f\n", h_B_e->Integral() ) ;
             printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hist from sumAllInfo for mu_B_SR_4b :  integral = %9.1f\n", h_B_mu->Integral() ) ;
          }
       }


      //--- C
       
       TH1* h_C(0x0) ;
       {
          string ch_name ;
          ch_name = "e_C_SR_4b" ;
          if ( year != "" ) { ch_name += year ; }
          std::map<string, ChannelInfo_t>::iterator ichan_e = sumDataProc.channels.find( ch_name ) ;
          if ( ichan_e == sumDataProc.channels.end() ) { printf("\n\n *** can't find data e_C_SR_4b channel in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ChannelInfo_t chan_e = ichan_e->second ;

          ch_name = "mu_C_SR_4b" ;
          if ( year != "" ) { ch_name += year ; }
          std::map<string, ChannelInfo_t>::iterator ichan_mu = sumDataProc.channels.find( ch_name ) ;
          if ( ichan_mu == sumDataProc.channels.end() ) { printf("\n\n *** can't find data mu_C_SR_4b channel in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ChannelInfo_t chan_mu = ichan_mu->second ;

          std::map<string, ShapeData_t>::iterator ishape_e = chan_e.shapes.find( mainHisto.Data() ) ;
          if ( ishape_e == chan_e.shapes.end() ) { printf("\n\n *** can't find data e_C_SR_4b channel shape in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ShapeData_t shape_e = ishape_e->second ;
          TH1* h_C_e = shape_e.histo() ;
          if ( h_C_e == 0x0 ) { printf("\n\n *** can't find data e_C_SR_4b channel shape histogram in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

          std::map<string, ShapeData_t>::iterator ishape_mu = chan_mu.shapes.find( mainHisto.Data() ) ;
          if ( ishape_mu == chan_mu.shapes.end() ) { printf("\n\n *** can't find data mu_C_SR_4b channel shape in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ShapeData_t shape_mu = ishape_mu->second ;
          TH1* h_C_mu = shape_mu.histo() ;
          if ( h_C_mu == 0x0 ) { printf("\n\n *** can't find data mu_C_SR_4b channel shape histogram in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

          h_C = (TH1*) h_C_e -> Clone("h_C_from_sumAllInfo") ;
          h_C -> Add( h_C_mu ) ;

          if ( verbose ) {
             printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hist from sumAllInfo for  e_C_SR_4b :  integral = %9.1f\n", h_C_e->Integral() ) ;
             printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hist from sumAllInfo for mu_C_SR_4b :  integral = %9.1f\n", h_C_mu->Integral() ) ;
          }
       }


      //--- D
       
       TH1* h_D(0x0) ;
       {
          string ch_name ;
          ch_name = "e_D_SR_4b" ;
          if ( year != "" ) { ch_name += year ; }
          std::map<string, ChannelInfo_t>::iterator ichan_e = sumDataProc.channels.find( ch_name ) ;
          if ( ichan_e == sumDataProc.channels.end() ) { printf("\n\n *** can't find data e_D_SR_4b channel in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ChannelInfo_t chan_e = ichan_e->second ;

          ch_name = "mu_D_SR_4b" ;
          if ( year != "" ) { ch_name += year ; }
          std::map<string, ChannelInfo_t>::iterator ichan_mu = sumDataProc.channels.find( ch_name ) ;
          if ( ichan_mu == sumDataProc.channels.end() ) { printf("\n\n *** can't find data mu_D_SR_4b channel in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ChannelInfo_t chan_mu = ichan_mu->second ;

          std::map<string, ShapeData_t>::iterator ishape_e = chan_e.shapes.find( mainHisto.Data() ) ;
          if ( ishape_e == chan_e.shapes.end() ) { printf("\n\n *** can't find data e_D_SR_4b channel shape in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ShapeData_t shape_e = ishape_e->second ;
          TH1* h_D_e = shape_e.histo() ;
          if ( h_D_e == 0x0 ) { printf("\n\n *** can't find data e_D_SR_4b channel shape histogram in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

          std::map<string, ShapeData_t>::iterator ishape_mu = chan_mu.shapes.find( mainHisto.Data() ) ;
          if ( ishape_mu == chan_mu.shapes.end() ) { printf("\n\n *** can't find data mu_D_SR_4b channel shape in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }
          ShapeData_t shape_mu = ishape_mu->second ;
          TH1* h_D_mu = shape_mu.histo() ;
          if ( h_D_mu == 0x0 ) { printf("\n\n *** can't find data mu_D_SR_4b channel shape histogram in sumAllInfo.  I quit.\n\n") ; gSystem -> Exit(-1); }

          h_D = (TH1*) h_D_e -> Clone("h_D_from_sumAllInfo") ;
          h_D -> Add( h_D_mu ) ;

          if ( verbose ) {
             printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hist from sumAllInfo for  e_D_SR_4b :  integral = %9.1f\n", h_D_e->Integral() ) ;
             printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  hist from sumAllInfo for mu_D_SR_4b :  integral = %9.1f\n", h_D_mu->Integral() ) ;
          }
       }

       resetNegativeBinsAndErrors( h_C, 1. ) ;
       resetNegativeBinsAndErrors( h_D, 1. ) ;
       resetNegativeBinsAndErrors( h_B, 0. ) ;

       double sum_C_val, sum_C_err ;
       double sum_D_val, sum_D_err ;

       sum_C_val = h_C -> IntegralAndError( 1, h_C->GetNbinsX(), sum_C_err ) ;
       sum_D_val = h_D -> IntegralAndError( 1, h_D->GetNbinsX(), sum_D_err ) ;

       double sum_C_over_D_val(0.) ;
       double sum_C_over_D_err(0.) ;
       if ( sum_D_val > 0 && sum_C_val > 0 ) {
          sum_C_over_D_val = sum_C_val / sum_D_val ;
          sum_C_over_D_err = sum_C_over_D_val * sqrt( pow( sum_C_err/sum_C_val, 2. ) + pow( sum_D_err/sum_D_val, 2. ) ) ;
       }
       if ( verbose ) {
          printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  from sumAllInfo, C/D = (%9.1f +/- %6.1f)/(%9.1f +/- %6.1f) = %6.3f +/- %6.3f\n",
             sum_C_val, sum_C_err, sum_D_val, sum_D_err, sum_C_over_D_val, sum_C_over_D_err ) ;
       }

      //-- use h_B from sum for the shape, but normalize B with the non-sum version

       hDD -> Add( h_B ) ;
       double integral_nonsum_B = hCtrl_SI -> Integral() ;
       double integral_sum_B = h_B -> Integral() ;
       if ( integral_sum_B > 0 ) {
          hDD -> Scale(  integral_nonsum_B / integral_sum_B ) ;
       } else {
          printf("\n\n *** integral of h_B from sumAllInfo is negative!!!  I quit.\n\n") ; gSystem -> Exit(-1) ;
       }

       hDD -> Scale( sum_C_over_D_val ) ;



    } else {

       hDD->Add(hCtrl_SI , 1.0);

       if(hCtrl_SB->Integral()<=0 || hCtrl_SI->Integral()<0 || hChan_SB->Integral()<0) alpha = 0;

       if ( verbose ) { printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  %s  final C/D = %6.3f\n", chData->first.c_str(), alpha ) ; }
       hDD->Scale(alpha);

    }


    hDD->SetTitle(DDProcName.Data());
    if ( verbose ) {
       printf(" --- verbose :  AllInfo_t::doBackgroundSubtraction :  %s  final prediction hist:  integral = %9.1f | ", chData->first.c_str(), hDD->Integral() ) ;
       for ( int i=1; i<hDD->GetNbinsX(); i++ ) {
          printf(" bin %2d = %9.1f +/- %6.1f | ", i, hDD->GetBinContent( i ) , hDD->GetBinError( i ) ) ;
       }
       printf("\n") ;
    }

    //save values for printout
    double valDD, valDD_err;
    //valDD = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-6){valDD=0.0; valDD_err=0.0;}
    valDD = hDD->IntegralAndError(1,hDD->GetXaxis()->GetNbins()+1,valDD_err); if(valDD<1E-3){valDD=0.0; valDD_err=0.0;}
    
    if(chDD->second.shapes[mainHisto.Data()].histo() == NULL){
      hDD->SetFillColor(634); hDD->SetLineColor(1); hDD->SetMarkerColor(634);
      hDD->SetFillStyle(1001);  hDD->SetLineWidth(1); hDD->SetMarkerStyle(20); hDD->SetLineStyle(1);
      chDD->second.shapes[mainHisto.Data()].uncShape[""] = hDD;
    }
    //remove all syst uncertainty
    chDD->second.shapes[mainHisto.Data()].clearSyst();
    //add syst uncertainty
    //chDD->second.shapes[mainHisto.Data()].uncScale[string("CMS_haa4b_sys_ddqcd_") + binName.Data() +"_"+chData->second.bin.c_str() + systpostfix.Data()] = valDD_err; //:1.8*valDD;
    chDD->second.shapes[mainHisto.Data()].uncScale[string("CMS_haa4b_sys_ddqcd_") + binName.Data() +"_"+chData->second.bin.c_str() + year.Data() + systpostfix.Data()] = valDD*datadriven_qcd_Syst; //:1.8*valDD;
    
    //    chDD->second.shapes[mainHisto.Data()].uncScale[string("CMS_haa4b_sys_ddqcd_") + binName.Data() + systpostfix.Data()] = ratioMC<0.5?valDD*datadriven_qcd_Syst:fabs(1.-ratioMC)*valDD;    

    //printout
    sprintf(Lcol    , "%s%s"  ,Lcol,    "|c");
    sprintf(Lchan   , "%s%25s",Lchan,   (string(" & $ ") + chData->second.channel+string(" - ")+chData->second.bin + string(" $ ")).c_str());
    sprintf(Lalph1  , "%s%25s",Lalph1,  (string(" &") + utils::toLatexRounded(alpha,alpha_err)).c_str());
    sprintf(Lyield  , "%s%25s",Lyield,  (string(" &") + utils::toLatexRounded(valDD,valDD_err,valDD*datadriven_qcd_Syst)).c_str());
    sprintf(LyieldMC, "%s%25s",LyieldMC,(string(" &") + utils::toLatexRounded(valMC,valMC_err)).c_str());
    sprintf(LalphMC  , "%s%25s",LalphMC,  (string(" &") + utils::toLatexRounded(alphaMC,alphaMC_err)).c_str()); 
    sprintf(LyieldPred  , "%s%25s",LyieldPred,  (string(" &") + utils::toLatexRounded(valDD_MC,valDD_MC_err)).c_str());     
    sprintf(RatioMC, "%s%25s",RatioMC, (string(" &") + utils::toLatexRounded(ratioMC,ratioMC_err)).c_str());
  } // end data channels




  procs["NonQCD"] = ProcessInfo_t(); //reset  

  //recompute total background  
  if ( verbose ) { printf("\n --- verbose : AllInfo_t::doBackgroundSubtraction : calling computeTotalBackground().\n") ; fflush(stdout) ; }
  computeTotalBackground(); 

  if(pFile){
    if (postfit){
      fprintf(pFile,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\tiny\n\\begin{center}\n\\caption{Data-driven QCD background estimation (using initial W and Top scale factors).}\n\\label{tab:table}\n");} else {
      fprintf(pFile,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\tiny\n\\begin{center}\n\\caption{Data-driven QCD background estimation.}\n\\label{tab:table}\n");
    }
    fprintf(pFile,"\\begin{tabular}{%s|}\\hline\n", Lcol);
    fprintf(pFile,"channel               %s\\\\\\hline\n", Lchan);
    fprintf(pFile,"$\\texttt{SF}_{qcd}$ measured    %s\\\\\n", Lalph1);
    fprintf(pFile,"QCD yield predicted (data)            %s\\\\\n", Lyield);
    fprintf(pFile,"QCD yield observed (MC)              %s\\\\\n", LyieldMC);
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"$\\texttt{SF}_{qcd}$ (MC)    %s\\\\\n", LalphMC);
    fprintf(pFile,"QCD yield predicted (MC)   %s\\\\\n", LyieldPred);
    fprintf(pFile,"QCD yield observed (MC)  %s\\\\\n", LyieldMC); 
    fprintf(pFile,"ratio MC (syst)        %s\\\\\n", RatioMC); 
    fprintf(pFile,"\\hline\n");   
    fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n\\end{document}\n");
  }



  if ( verbose ) { printf("\n\n --- verbose :  AllInfo_t::doBackgroundSubtraction : end \n\n") ; }

} // doBackgroundSubtraction



//
// Sum up all background processes and add this as a total process
//
void AllInfo_t::computeTotalBackground(){
  if ( verbose ) { printf("\n  --- verbose: computeTotalBackground : begin.\n" ) ; fflush(stdout) ; }
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
    if ( verbose ) { printf("    --- verbose: computeTotalBackground :  calling addProc( procInfo_Bckgs, it->second, false)  it->first = %s\n", (it->first).c_str() ) ; fflush(stdout) ; }
    addProc(procInfo_Bckgs, it->second, false);
  }
  //Compute total background systematics
  for(std::map<string, ProcessInfo_t>::iterator it=procs.begin(); it!=procs.end();it++){
    if(it->first=="total" || it->second.isBckg!=true)continue;
    if ( verbose ) { 
      printf("    --- verbose: computeTotalBackground :  calling addProc( procInfo_Bckgs, it->second, true)  it->first = %s\n", (it->first).c_str() ) ; fflush(stdout) ; }
    addProc(procInfo_Bckgs, it->second, true);
  }
  if ( verbose ) { printf(" ---  verbose: computeTotalBackground : end.\n\n" ) ; fflush(stdout) ; }

}


//
// Replace the Data process by TotalBackground
//
void AllInfo_t::blind() {
   if ( verbose ) { printf("\n  --- verbose : AllInfo_t::blind : begin.\n") ; fflush(stdout) ; }
  if(procs.find("total")==procs.end())computeTotalBackground();

  if(true){ //always replace data
    //if(procs.find("data")==procs.end()){ //true only if there is no "data" samples in the json file
    sorted_procs.push_back("data");           
    if ( verbose ) { printf("  verbose :  AllInfo_t::blind() :  before resetting data\n" ) ; procs["data"].printProcess() ; }
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
    if ( verbose ) { printf("  verbose :  AllInfo_t::blind() :  after resetting data\n" ) ; procs["data"].printProcess() ; }
  }
   if ( verbose ) { printf(" --- verbose : AllInfo_t::blind : end.\n\n") ; fflush(stdout) ; }
}



//---------------------------------------------------------------

void AllInfo_t::replaceHighSensitivityBinsWithBG() {

  if ( verbose ) { printf("\n  --- verbose : AllInfo_t::replaceHighSensitivityBinsWithBG : begin.\n") ; fflush(stdout) ; }
  if(procs.find("total")==procs.end())computeTotalBackground();

  std::map<string, ProcessInfo_t>::iterator itbg=procs.find("total");
  if ( itbg==procs.end() ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG : no total process???\n\n") ; return ; }
  std::map<string, ProcessInfo_t>::iterator idata=procs.find("data");
  if ( idata==procs.end() ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG : no data process???\n\n") ; return ; }

  ProcessInfo_t& total_proc = itbg -> second ;
  ProcessInfo_t& data_proc = idata -> second ;

  if ( verbose ) { printf(" AllInfo_t::replaceHighSensitivityBinsWithBG : before replacement.\n") ; total_proc.printProcess() ; data_proc.printProcess() ; }

  for ( std::map<string, ChannelInfo_t>::iterator ic = data_proc.channels.begin(); ic!= data_proc.channels.end(); ic++ ) {

     string chan_key = ic -> first ;
     ChannelInfo_t& data_chan = ic -> second ;

     if ( chan_key.find("SR")!=string::npos ) {

        std::map<string, ChannelInfo_t>::iterator itbgc = total_proc.channels.find( chan_key ) ;
        if ( itbgc == total_proc.channels.end() ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG : can't find channel %s in total BG!  bailing out.\n\n", chan_key.c_str() ) ; return ; }
        ChannelInfo_t& total_chan = itbgc -> second ;

        for ( std::map<string, ShapeData_t>::iterator is = data_chan.shapes.begin(); is!= data_chan.shapes.end(); is++ ) {

           string shape_key = is -> first ;
           ShapeData_t& data_shape = is -> second ;

           std::map<string, ShapeData_t>::iterator itbgs = total_chan.shapes.find( shape_key ) ;
           if ( itbgs == total_chan.shapes.end() ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG : can't find shape %s for channel %s in in total BG!  bailing out.\n\n", shape_key.c_str(), chan_key.c_str() ) ; return ; }
           ShapeData_t& total_shape = itbgs -> second ;

           TH1* data_hist = data_shape.histo() ;
           TH1* total_hist = total_shape.histo() ;

           if ( data_hist == 0x0 ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG : data hist is null pointer!!! bailing out.\n\n") ; return ; }
           if ( total_hist == 0x0 ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG : total BG hist is null pointer!!! bailing out.\n\n") ; return ; }

           if ( data_hist -> GetNbinsX() != 5 ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG :  was expecting 5 bins.  found %d in data hist.  bailing out.\n\n", data_hist -> GetNbinsX() ) ; }
           if ( total_hist -> GetNbinsX() != 5 ) { printf("\n\n *** AllInfo_t::replaceHighSensitivityBinsWithBG :  was expecting 5 bins.  found %d in total BG hist.  bailing out.\n\n", total_hist -> GetNbinsX() ) ; }

           data_hist -> SetBinContent( 4, total_hist -> GetBinContent( 4 ) ) ;
           data_hist -> SetBinContent( 5, total_hist -> GetBinContent( 5 ) ) ;


        } // is

     }


  } // ic

  if ( verbose ) { printf(" AllInfo_t::replaceHighSensitivityBinsWithBG : after replacement.\n") ; data_proc.printProcess() ; }


  if ( verbose ) { printf("\n  --- verbose : AllInfo_t::replaceHighSensitivityBinsWithBG : end.\n") ; fflush(stdout) ; }
} // AllInfo_t::replaceHighSensitivityBinsWithBG


//---------------------------------------------------------------



void AllInfo_t::getYieldsFromShape(FILE* pFile, std::vector<TString>& selCh, string histoName, FILE* pFileInc){
  if(!pFileInc)pFileInc=pFile;

  std::vector<string> VectorProc;
  std::map<string, bool> MapChannel;
  std::map<string, std::map<string, string> > MapProcChYields;         
  std::map<string, bool> MapChannelBin;
  std::map<string, std::map<string, string> > MapProcChYieldsBin;         
  std::map<string, std::map<string, std::vector<double> > > MapProcChBinYields;         
  std::map<string, std::map<string, std::vector<double> > > MapProcChBinErrors;         
  std::map<string, std::vector<int> > MapProcChBin;         
  std::vector<double> MapSignChBinYields;
  std::vector<double> MapBckgChBinYields;

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
      fflush(stdout) ;

      TH1* h = ch->second.shapes[histoName].histo();
      double valerr = 0.;
      double val  = 0.;
      if (h!=NULL) {
        int start_bin = (docut ? h->FindBin(sstyCut) : 1);
        val = h->IntegralAndError(start_bin,h->GetXaxis()->GetNbins(),valerr);
      }

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

      TString LabelText = TString("$") + ch->second.channel+ "\\ " +ch->second.bin + TString("$");
      LabelText.ReplaceAll("eq"," ="); LabelText.ReplaceAll("g =","\\geq"); LabelText.ReplaceAll("l =","\\leq"); LabelText.ReplaceAll("mu","\\mu"); LabelText.ReplaceAll("_","\\_");
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


    //All Channels
  fprintf(pFile,"\\documentclass{article}\n\\usepackage{graphicx}\n\\usepackage{geometry}\n\\geometry{\n\tleft=10mm,\n\tright=10mm,\n\ttop=10mm,\n\tbottom=10mm\n}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n\\resizebox{\\textwidth}{!}{\n ");
  fprintf(pFile, "\\begin{tabular}{|c|"); for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ fprintf(pFile, "c|"); } fprintf(pFile, "}\\hline\n");
  fprintf(pFile, "channel");   for(auto ch = MapChannel.begin(); ch!=MapChannel.end();ch++){ fprintf(pFile, " & %s", ch->first.c_str()); } fprintf(pFile, "\\\\\\hline\n");
  for(auto proc = VectorProc.begin();proc!=VectorProc.end(); proc++){
    if(*proc=="total")fprintf(pFile, "\\hline\n");
    auto ChannelYields = MapProcChYields.find(*proc);
    if(ChannelYields == MapProcChYields.end())continue;
    TString procName = (*proc).c_str(); if(procName.Contains("#")){procName = "$" + procName + "$";} procName.ReplaceAll("#","\\");
    fprintf(pFile, "%s ", procName.Data()); 
//    fprintf(pFile, "%s ", proc->c_str()); 
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
  fprintf(pFile,"\\end{tabular}\n}\n\\end{center}\n\\end{sidewaystable}\n\\end{document}\n");
  
    //All Bins
  fprintf(pFileInc,"\\documentclass{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");
  fprintf(pFileInc, "\\begin{tabular}{|c|"); for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ fprintf(pFileInc, "c|"); } fprintf(pFileInc, "}\\hline\n");
  fprintf(pFileInc, "channel");   for(auto ch = MapChannelBin.begin(); ch!=MapChannelBin.end();ch++){ fprintf(pFileInc, " & %s", ch->first.c_str()); } fprintf(pFileInc, "\\\\\\hline\n");
  for(auto proc = VectorProc.begin();proc!=VectorProc.end(); proc++){
    if(*proc=="total")fprintf(pFileInc, "\\hline\n");
    auto ChannelYields = MapProcChYieldsBin.find(*proc);
    if(ChannelYields == MapProcChYieldsBin.end())continue;
    TString procName = (*proc).c_str(); procName.ReplaceAll("#","\\"); procName = "$" + procName + "$";
    fprintf(pFile, "%s ", procName.Data()); 
//    fprintf(pFileInc, "%s ", proc->c_str()); 
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
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end()){
	  printf(" --- dropCtrlChannels::: process = %s, channel dropped : %s\n",procName.c_str(),ch->first.c_str());
	  it->second.channels.erase(ch); ch=it->second.channels.begin();}
      }
    }
  }



  //
  // drop background process that have a negligible yield
  //
  void AllInfo_t::dropSmallBckgProc(std::vector<TString>& selCh, string histoName, double threshold)
  {
   
    auto total = procs.find("total");
    if(total==procs.end()) {printf("dropSmallBckgProc: Error, cannot find process: total\n");return;}
    std::map<string, double> total_yields;
    std::map<string, std::map<string, double> > map_yields;
    for(std::map<string, ChannelInfo_t>::iterator ch = total->second.channels.begin(); ch!=total->second.channels.end(); ch++){
      if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
      TH1 *h=ch->second.shapes[histoName].histo();
      total_yields[ch->first] = 0.;
      if(h!=NULL) {total_yields[ch->first] = h->Integral();printf("dropsmallBckgProc, total in channel %s %f\n",ch->first.c_str(), total_yields[ch->first]);}
    }

    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];
      if(procName.compare("total") == 0) continue;
      if ( procName.compare("ddqcd") == 0 ) {
         printf("  AllInfo_t::dropSmallBckgProc:  excluding ddqcd proc from consideration.  Always keep it.\n") ; fflush(stdout) ;
         continue ;
      }
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      if(!it->second.isBckg)continue;
//      printf("dropSmallBckgProc, Process: %s\n",  procName.c_str());
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end();ch++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        
        map_yields[it->first][ch->first] = 0.;  
        TH1 *h=ch->second.shapes[histoName].histo();
        if (h!=NULL) map_yields[it->first][ch->first] = h->Integral(); 
      }
    }

    for(std::map<string, std::map<string, double> >::iterator p = map_yields.begin();p!=map_yields.end();p++){
      for(std::map<string, double>::iterator ch = p->second.begin();ch!=p->second.end();ch++){
        double tot = total_yields[ch->first];
        double yield = map_yields[p->first][ch->first];
        if(tot>0 && yield/tot<threshold){
          printf("Drop %s from the list of backgrounds in the channel %s because of negligible rate (%f of total bckq)\n", p->first.c_str(), ch->first.c_str(), yield/tot);
          procs.find(p->first)->second.channels.erase(procs.find(p->first)->second.channels.find(ch->first));
        }
      }
    }

//    double total = map_yields["total"];
//    for(std::map<string, double>::iterator Y=map_yields.begin();Y!=map_yields.end();Y++){
      //      if(Y->first.find("ddqcd")<std::string::npos)continue;//never drop this background
      //      if(Y->first.find("VV")<std::string::npos)continue;//never drop this background
      // if(Y->first.find("Vh")<std::string::npos)continue;//never drop this background 
      //      if(Y->first.find("t#bar{t}+#gammaZW")<std::string::npos)continue;//never drop this background    
//      if(Y->second/total<threshold){
//        printf("Drop %s from the list of backgrounds because of negligible rate (%f%% of total bckq)\n", Y->first.c_str(), Y->second/total);
//        for(std::vector<string>::iterator p=sorted_procs.begin(); p!=sorted_procs.end();p++){if((*p)==Y->first ){sorted_procs.erase(p);break;}}
//        procs.erase(procs.find(Y->first ));
//      }
//    }
  }



  //------------------------------------------------------------------------------------------------------


  //
  // Make a summary plot
  //
  void AllInfo_t::showShape(std::vector<TString>& selCh , TString histoName, TString SaveName)
  {

    if ( verbose ) {
       printf(" --- verbose : AllInfo_t::showShape :  begin \n") ;
       printf(" --- verbose : AllInfo_t::showShape :  channels : ") ;
       for ( int ci=0; ci<selCh.size(); ci++ ) { printf(" %s , ", selCh[ci].Data() ) ; }
       printf("\n") ;
       printf(" --- verbose : AllInfo_t::showShape :  histoName = %s , SaveName = %s\n", histoName.Data(), SaveName.Data() ) ;
       fflush(stdout) ;
    }
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

      if ( verbose ) { printf(" --- verbose : AllInfo_t::showShape :  proc = %s\n", process.Data() ) ; fflush(stdout) ; }

      if(it==procs.end())continue;
      //loop on channels for each process
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        if(std::find(selCh.begin(), selCh.end(), ch->second.channel)==selCh.end())continue;

        if(ch->second.shapes.find(histoName.Data())==(ch->second.shapes).end())continue;
        if(modeDD && ch->first.find("_A_")==std::string::npos) continue; //  only consider region A

        TH1* h = ch->second.shapes[histoName.Data()].histo();
        if (!h) continue;
        if ( verbose ) { printf(" --- verbose : AllInfo_t::showShape :  proc = %s , chan = %s, hist = %s\n", process.Data(), ch->first.c_str(), h->GetName() ) ; fflush(stdout) ; }
        //if(process.Contains("SOnly_S") ) h->Scale(10);  
        if(it->first=="total"){

          double syst_scale = std::max(0.0, ch->second.shapes[histoName.Data()].getScaleUncertainty());
	  //          double syst_shape = std::max(0.0, ch->second.shapes[histoName.Data()].getBinShapeUncertainty((it->first+ch->first).c_str()));
	  //          double syst = sqrt(pow(syst_scale,2)+pow(syst_shape,2)); 
	  //          double Uncertainty_scale = syst / h->Integral();
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

	  //	  if ( verbose ) { printf(" --- verbose : AllInfo_t::showShape :  proc = %s\n", process.Data() ) ; fflush(stdout) ; }

          continue;//otherwise it will fill the legend
        }else if(it->second.isBckg){                 
          if(map_stack.find(ch->first)==map_stack.end()){
            map_stack[ch->first] = new THStack((ch->first+"stack").c_str(),(ch->first+"stack").c_str());
            map_mc   [ch->first] = (TH1D*)h->Clone((ch->first+"mc").c_str()); //utils::root::checkSumw2((ch->first+"mc").c_str());
          }else{
            map_mc [ch->first]->Add(h); 
          }
          map_stack   [ch->first]->Add(h,"HIST");
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

      legA->Draw("same");    legA->SetTextFont(42);

      double iLumi=35.9;
      double iEcm=13;
      if(lumi > 0) iLumi = lumi;
   
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
      c[I]->SaveAs(LabelText+"_Shape"+year+".png");
      c[I]->SaveAs(LabelText+"_Shape"+year+".pdf");
      c[I]->SaveAs(LabelText+"_Shape"+year+".C");
      delete c[I];
      
      I++;
    }
    if ( verbose ) {
       printf(" --- verbose : AllInfo_t::showShape :  end \n") ;
       fflush(stdout) ;
    }
  } // showShape


  //------------------------------------------------------------------------------------------------------

  //
  // Make a summary plot
  //
  void AllInfo_t::showUncertainty(std::vector<TString>& selCh , TString histoName, TString SaveName)
  {
    string UncertaintyOnYield="";  char txtBuffer[4096];
    TFile *unc_f = TFile::Open("unc.root", "recreate");

    sprintf(txtBuffer,"\\documentclass{article}\n\\usepackage{graphicx}\n\\usepackage{geometry}\n\\geometry{\n\tleft=10mm,\n\tright=10mm,\n\ttop=10mm,\n\tbottom=10mm\n}\n\\usepackage[utf8]{inputenc}\n\\usepackage{rotating}\n\\begin{document}\n"); UncertaintyOnYield += txtBuffer;    

    //loop on sorted proc
    for(unsigned int p=0;p<sorted_procs.size();p++){
      int NLegEntry = 0;
      std::map<string, int               > map_legend;
      std::vector<TH1*>                    toDelete;             
      TLegend* legA  = new TLegend(0.03,0.89,0.97,0.95, "");
      legA->SetTextSize(0.015);

      string procName = sorted_procs[p];
      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      //if(it->first=="total" || it->first=="data")continue;  //only do samples which have systematics
      if(it->first=="data")continue;  //only do samples which have systematics  

      std::map<string, bool> mapUncType;
      std::map<string, std::map< string, double> > mapYieldPerBin;
      std::map<string, std::pair< double, double> > mapYieldInc;


      int NBins = it->second.channels.size()/selCh.size() > 1 ? it->second.channels.size()/selCh.size() : 2;
//      if((it->second.channels.size()-NBins*selCh.size())>0) NBins++;
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


        //-- Look for the total BG process.  If found, find the corresponding total background histogram.
        //     Will pass it to makeStatUnc so that it can decide whether to include BinByBin stat uncertainty based on err_i / sqrt( N(total)_i ) for each bin i.
	
        TH1* h_total(0x0) ;
        std::map<string, ProcessInfo_t>::iterator tpi = procs.find("total") ;
        if ( tpi != procs.end() ) {
           ProcessInfo_t total_proc = tpi->second ;
           std::map<string, ChannelInfo_t>::iterator tci = total_proc.channels.find( ch->first ) ;
           if ( tci != total_proc.channels.end() ) {
              ChannelInfo_t total_chan = tci -> second ;
              std::map<string, ShapeData_t>::iterator tsi = total_chan.shapes.find( histoName.Data() ) ;
              if ( tsi != total_chan.shapes.end() ) {
                 ShapeData_t total_shape = tsi -> second ;
                 h_total = total_shape.histo() ;
              } // tsi
           } // tci
        } // tpi
	

        //add the stat uncertainty is there;
        //ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );//add stat uncertainty to the uncertainty map;

        if((it->second.shortName).find("wh")!=std::string::npos)ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_wh")).Data(),systpostfix.Data(), h_total );// attention

        //      else if((it->second.shortName).find("qqH")!=std::string::npos)ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_qqH")).Data(),systpostfix.Data(), false );

        else ch->second.shapes[histoName.Data()].makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false, h_total );

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
          //systName.ReplaceAll("up","");
          //systName.ReplaceAll("down","");
          //if(systName.Index("jes")<0 && systName.Index("umet")<0 && systName.Index("resj")<0) continue;

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
              mapYieldPerBin[systName.Data()][ch->first] = ( 1 - (varYield/yield)); // fabs( 1 - (varYield/yield));
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
          hvar->Write(vh_tag + "_" + systName+"_Uncertainty_"+ch->second.channel+"_"+ch->second.bin+"_"+it->second.shortName);
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
      ////////////c1->SaveAs(SaveName+vh_tag+"_Uncertainty_"+it->second.shortName+".png");  //-- owen: temporarily disable this (file is huge)
      c1->SaveAs(SaveName+vh_tag+"_Uncertainty_"+it->second.shortName+".pdf");  //-- owen: temporarily disable this (file is huge)
      ////////////c1->SaveAs(SaveName+vh_tag+"_Uncertainty_"+it->second.shortName+".C");  //-- owen: temporarily disable this (file is huge)
      delete c1;             

      for(unsigned int i=0;i<toDelete.size();i++){delete toDelete[i];} //clear the objects


      //add inclusive uncertainty as a channel
      //
      for(auto systIt=mapYieldInc.begin(); systIt!=mapYieldInc.end(); systIt++){ mapYieldPerBin[systIt->first][" Inc"] = systIt->second.first/systIt->second.second;  }
      //print uncertainty on yield            
      sprintf(txtBuffer,"\\begin{table}[htp]\n\\tiny\n\\begin{center}\n\\caption{Uncertainty on the yield for the process: \\bf{%s}}\n\\label{tab:tablesys}\n\\resizebox{\\textwidth}{!}{\n ",it->first.c_str()); UncertaintyOnYield += txtBuffer; 
      sprintf(txtBuffer, "\\begin{tabular}{ccccccc} \n"); UncertaintyOnYield+= txtBuffer; 
      //      sprintf(txtBuffer, "\\begin{tabular}{ccccccc} \n {\\bf{%s}} & & & & & & \\\\ \\hline \n", it->first.c_str());  UncertaintyOnYield+= txtBuffer; 
      //      sprintf(txtBuffer, "\\multicolumn{%i}{'c'}{\\bf{%s}}\\\\ \n", I+1, it->first.c_str());  UncertaintyOnYield+= txtBuffer;
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
      sprintf(txtBuffer,"\\end{tabular}\n}\n\\end{center}\n\\end{table}\n");UncertaintyOnYield += txtBuffer; 
    }
    sprintf(txtBuffer,"\\end{document}\n");UncertaintyOnYield += txtBuffer;

    FILE* pFile = fopen(SaveName+vh_tag+"_Uncertainty.txt", "w");
    if(pFile){ fprintf(pFile, "%s\n", UncertaintyOnYield.c_str()); fclose(pFile);}
    unc_f->Close();
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

    if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  begin with histoName = %s, file = %s\n", histoName.c_str(), fout->GetName() ) ; fflush(stdout) ; }
    //order the proc first
    sortProc();

    //Loop on processes and channels
    for(unsigned int p=0;p<sorted_procs.size();p++){
      string procName = sorted_procs[p];

      //////if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  proc_name = %s\n", procName.c_str() ) ; fflush(stdout) ; }

      std::map<string, ProcessInfo_t>::iterator it=procs.find(procName);
      if(it==procs.end())continue;
      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;

        //////if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  chbin = %s , directory for next histograms.\n", chbin.Data() ) ; fflush(stdout) ; }

        if(!fout->GetDirectory(chbin)){fout->mkdir(chbin);}fout->cd(chbin);

        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        TH1* h = shapeInfo.histo();
        if(h==NULL)continue;




        //-- Look for the total BG process.  If found, find the corresponding total background histogram.
        //     Will pass it to makeStatUnc so that it can decide whether to include BinByBin stat uncertainty based on err_i / sqrt( N(total)_i ) for each bin i.

        TH1* h_total(0x0) ;
        std::map<string, ProcessInfo_t>::iterator tpi = procs.find("total") ;
        if ( tpi != procs.end() ) {
           ProcessInfo_t total_proc = tpi->second ;
           std::map<string, ChannelInfo_t>::iterator tci = total_proc.channels.find( ch->first ) ;
           if ( tci != total_proc.channels.end() ) {
              ChannelInfo_t total_chan = tci -> second ;
              std::map<string, ShapeData_t>::iterator tsi = total_chan.shapes.find( histoName ) ;
              if ( tsi != total_chan.shapes.end() ) {
                 ShapeData_t total_shape = tsi -> second ;
                 h_total = total_shape.histo() ;
              } // tsi
           } // tci
        } // tpi


        //shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), it->second.isSign );//add stat uncertainty to the uncertainty map;
        //shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false );//add stat uncertainty to the uncertainty map;
        //////if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  TH1 name : %s\n", h -> GetName() ) ; fflush(stdout) ; }

        //Fix
        if((it->second.shortName).find("ggH")!=std::string::npos)shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_ggH")).Data(),systpostfix.Data(), false, h_total );// attention
        else if((it->second.shortName).find("qqH")!=std::string::npos)shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+TString("_qqH")).Data(),systpostfix.Data(), false, h_total );
        else shapeInfo.makeStatUnc("_CMS_haa4b_", (TString("_")+ch->first+"_"+it->second.shortName).Data(),systpostfix.Data(), false, h_total );
        fout->cd(chbin);

        TString proc = it->second.shortName.c_str();
        //      printf("\n\n");
        //      printf("Process: %s and channel : %s , with N shapes= %d\n",proc.Data(),ch->first.c_str(),(int)shapeInfo.uncShape.size());

        for(std::map<string, TH1*  >::iterator unc=shapeInfo.uncShape.begin();unc!=shapeInfo.uncShape.end();unc++){
          TString syst   = unc->first.c_str();
          TH1*    hshape = unc->second;
          hshape->SetDirectory(0);

          if(syst==""){
            //central shape (for data call it data_obs)
            hshape->SetName(proc); 
            if(it->first=="data"){
              hshape->Write("data_obs");
            }else{
              hshape->Write(proc+postfix);
            }
          }else if(runSystematics && proc!="data" && (syst.EndsWith("Up") || syst.EndsWith("Down"))){
            //if empty histogram --> no variation is applied except for stat
            
            if(!syst.Contains("stat") && (hshape->Integral()<h->Integral()*0.01 || isnan((float)hshape->Integral()))){hshape->Reset(); hshape->Add(h,1); }

//          std::cout << "*****************" << proc << postfix << syst << "******************" << std::endl;
            if(hshape->Integral()<=0){
              for(int ibin=1; ibin<=hshape->GetXaxis()->GetNbins(); ibin++) //hmirrorshape->SetBinContent(ibin, 1E-10);     
                hshape->SetBinContent(ibin, 1E-10);
            }
            if(syst.Contains("CMS_haa4b_AbsoluteStat_jes") || syst.Contains("CMS_haa4b_RelativeJEREC1_jes") || syst.Contains("CMS_haa4b_RelativeJEREC2_jes") || syst.Contains("CMS_haa4b_RelativePtEC1_jes") || syst.Contains("CMS_haa4b_RelativePtEC2_jes") || syst.Contains("CMS_haa4b_RelativeSample_jes") || syst.Contains("CMS_haa4b_RelativeStatEC_jes") || syst.Contains("CMS_haa4b_RelativeStatFSR_jes") || syst.Contains("CMS_haa4b_RelativeStatHF_jes") || syst.Contains("CMS_haa4b_TimePtEta_jes") || syst.Contains("CMS_haa4b_umet") || syst.Contains("CMS_res_j") ){
                if(syst.Contains("Up")) syst = syst.ReplaceAll("Up", year+"Up");
                if(syst.Contains("Down")) syst = syst.ReplaceAll("Down", year+"Down");;
            }
            //write variation to file
            hshape->SetName(proc+syst);

            ///if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  hshape name just before Write = %s, Write argument = %s\n", hshape->GetName(), (proc+postfix+syst).Data() ) ; fflush(stdout) ; }

            hshape->Write(proc+postfix+syst);
          }else if(runSystematics && proc!="data"){
            //for one sided systematics the down variation mirrors the difference bin by bin
            hshape->SetName(proc+syst);

            //////if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  hshape name just before Write = %s, Write argument = %s\n", hshape->GetName(), (proc+postfix+syst+"Up").Data() ) ; fflush(stdout) ; }

            hshape->Write(proc+postfix+syst+"Up");
            TH1 *hmirrorshape=(TH1 *)hshape->Clone(proc+syst+"Down");
            for(int ibin=1; ibin<=hmirrorshape->GetXaxis()->GetNbins(); ibin++){
              double bin = 2*h->GetBinContent(ibin)-hmirrorshape->GetBinContent(ibin);
              if(bin<0)bin=0;
              hmirrorshape->SetBinContent(ibin,bin);
            }
            if(hmirrorshape->Integral()<=0) hmirrorshape->SetBinContent(1, 1E-10);

            //////if ( verbose ) { printf(" ---  verbose : AllInfo_t::saveHistoForLimit :  hshape name just before Write = %s, Write argument = %s\n", hmirrorshape->GetName(), (proc+postfix+syst+"Down").Data() ) ; fflush(stdout) ; }

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
        //      if(it->second.shortName.find("wlnu")!=string::npos){shapeInfo.uncScale["norm_wlnu"] = integral*0.10;}  
        //      if(it->second.shortName.find("ttbarjet")!=string::npos){shapeInfo.uncScale["norm_tt"] = integral*0.10;}  

        //if(it->second.shortName.find("zvv")!=string::npos){shapeInfo.uncScale["norm_zvv"] = integral*0.50;}    
        //if(it->second.shortName.find("vhbb")!=string::npos){shapeInfo.uncScale["norm_vhbb"] = integral*0.10;}     
        //if(it->second.shortName.find("singleto")!=string::npos){shapeInfo.uncScale["norm_singletop"] = integral*0.05;}
        //if(it->second.shortName.find("ttbargam")!=string::npos){shapeInfo.uncScale["norm_topgzw"] = integral*0.15;} 
        if(it->second.shortName.find("otherbkg")!=string::npos){shapeInfo.uncScale["norm_otherbkgds"] = integral*0.27;} 
        if(runZh) {
          if(it->second.shortName.find("wlnu")!=string::npos){shapeInfo.uncScale["norm_wjet"] = integral*0.02;}     
        } else {
          if(it->second.shortName.find("zll")!=string::npos){shapeInfo.uncScale["norm_zll"] = integral*0.02;}
        }

        if(it->second.shortName.find("ttbarlig")!=string::npos){shapeInfo.uncScale["norm_toplight"] = integral*0.06;} 
        //      if(it->second.shortName.find("ttbarcba")!=string::npos){shapeInfo.uncScale["norm_topcc"] = integral*0.50;} 
        if (!subFake){
          if(it->second.shortName.find("qcd")!=string::npos){shapeInfo.uncScale["norm_qcd"] = integral*0.50;} 
        }
        //uncertainties to be applied only in higgs analyses
        if(mass>0){

          //Introduce theory uncertaintly in signal x-section between 2016 and 2017/2018 samples due to change in PYTHIA tune:
          // https://hypernews.cern.ch/HyperNews/CMS/get/generators/4546/1.html
          if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["thxsec_wh"] = integral*0.011;}

          //bin migration at theory level
          // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#WH_Process
          if (runZh) {
            if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["QCDscale_wh"]  = integral*0.038;} //QCD scale
            if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["PDFscale_wh"]  = integral*0.016;} //PDF+as scale 
          } else {
            if(it->second.shortName.find("wh")!=string::npos ){shapeInfo.uncScale["QCDscale_wh"]  = integral*0.007;} //QCD scale
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
    if ( verbose ) { printf(" --- verbose : AllInfo_t::buildDataCards : histoName = %s, url = %s\n", histoName.c_str(), url.Data() ) ; fflush(stdout) ; }
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
        printf("pushing into sign_procs: %s\n", it->first.c_str());
      }
      if(it->second.isBckg) {clean_procs.push_back(procName); }

      for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
        TString chbin = ch->first;
        if(ch->second.shapes.find(histoName)==(ch->second.shapes).end())continue;
        allChannels[ch->first] = true;
        ShapeData_t& shapeInfo = ch->second.shapes[histoName];      
        for(std::map<string, double>::iterator unc=shapeInfo.uncScale.begin();unc!=shapeInfo.uncScale.end();unc++){
          if(unc->first=="")continue;
          if ( verbose ) { printf(" --- verbose : AllInfo_t::buildDataCards :  proc = %20s , chan = %20s , syst = %s\n", procName.c_str(), chbin.Data(), unc->first.c_str() ) ; }
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

      //-- owen: August 10, 2020:  remove CR5j completely.
      //      if (C->first.find("CR5j")!=string::npos) continue;

      if(C->first.find("emu")==string::npos) combinedcard += (C->first+"=").c_str()+dcName+" ";
      if(runZh) {
        if(C->first.find("ee"  )!=string::npos)eecard   += (C->first+"=").c_str()+dcName+" ";
        if(C->first.find("mumu")!=string::npos)mumucard += (C->first+"=").c_str()+dcName+" ";

        if(C->first.find("ee_A_CR"  )!=string::npos)ecrcard   += (C->first+"=").c_str()+dcName+" ";      
        if(C->first.find("mumu_A_CR")!=string::npos)mucrcard += (C->first+"=").c_str()+dcName+" ";  
        
        if(C->first.find("emu_A_CR"  )!=string::npos || C->first.find("emu_A_SR_3b")!=string::npos){ // emu_A_SR_4b is dropped in 2016 and 2017 zh
          eecard   += (C->first+"=").c_str()+dcName+" ";mumucard += (C->first+"=").c_str()+dcName+" ";
        }
        if(C->first.find("emu_A_SR_4b"  )!=string::npos && inFileUrl.Contains("2018")){
          eecard   += (C->first+"=").c_str()+dcName+" ";mumucard += (C->first+"=").c_str()+dcName+" ";
        }


      } else {
        if(C->first.find("e"  )!=string::npos)eecard   += (C->first+"=").c_str()+dcName+" ";    
        if(C->first.find("mu")!=string::npos)mumucard += (C->first+"=").c_str()+dcName+" ";  

        if(C->first.find("e_A_CR"  )!=string::npos)ecrcard   += (C->first+"=").c_str()+dcName+" ";  
        if(C->first.find("mu_A_CR")!=string::npos)mucrcard += (C->first+"=").c_str()+dcName+" "; 
      }

      std::vector<string> valid_procs;
      for(unsigned int j=0; j<clean_procs.size(); j++){
        if(clean_procs[j].find("Wh")!=string::npos) valid_procs.push_back(clean_procs[j]); //always include signal process in
        else if (procs[clean_procs[j]].channels[C->first].shapes[histoName].histo()!=NULL) valid_procs.push_back(clean_procs[j]);
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
      fprintf(pFile,"%55s  ", "bin");     for(unsigned int j=0; j<valid_procs.size(); j++){ 
        fprintf(pFile,"%8s ", "bin1")                     ;}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<valid_procs.size(); j++){ 
        fprintf(pFile,"%8s ", procs[valid_procs[j]].shortName.c_str());
      }  
      fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "process"); for(unsigned int j=0; j<valid_procs.size(); j++){ 
        fprintf(pFile,"%8i ", ((int)j)-(nsign-1)    );}  fprintf(pFile,"\n");
      fprintf(pFile,"%55s  ", "rate");    for(unsigned int j=0; j<valid_procs.size(); j++){ 
        double fval=0;      
        if (procs[valid_procs[j]].channels[C->first].shapes[histoName].histo()!=NULL) { 
          fval=procs[valid_procs[j]].channels[C->first].shapes[histoName].histo()->Integral();}    
        fprintf(pFile,"%8f ", fval);     
        //      fprintf(pFile,"%8f ", procs[valid_procs[j]].channels[C->first].shapes[histoName].histo()->Integral() );
      }
      fprintf(pFile,"\n");
      fprintf(pFile, "-------------------------------\n");

      for(std::map<string, bool>::iterator U=allSysts.begin(); U!=allSysts.end();U++){
        if ( verbose ) { printf(" --- verbose : AllInfo_t::buildDataCards : datacard %s :  allSysts : %s  %s \n", dcName.Data(), U->first.c_str(), U->second?"shape":"lnN" ) ; }
        if ( noCorrelatedStatUnc && U->first.find("CMS_haa4b_stat_") !=string::npos ) {
           if ( verbose ) { printf(" --- verbose : AllInfo_t::buildDataCards :   noCorrelatedStatUnc set to true.  Skipping this one.\n") ; }
           continue ;
        }
        char line[2048];
        sprintf(line,"%-45s %-10s ", U->first.c_str(), U->second?"shape":"lnN");
        bool isNonNull = false;
        for(unsigned int j=0; j<valid_procs.size(); j++){
          ShapeData_t& shapeInfo = procs[valid_procs[j]].channels[C->first].shapes[histoName];
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
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + b#bar{b}")!=valid_procs.end())  fprintf(pFile,"tt_norm_e%s rateParam bin1 ttbarbba 1\n",year.Data()); 
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + c#bar{c}")!=valid_procs.end())  fprintf(pFile,"tt_norm_e%s rateParam bin1 ttbarcba 1\n",year.Data());    
          if(std::find(valid_procs.begin(), valid_procs.end(), "Z#rightarrow ll")!=valid_procs.end()){  
            if (C->first.find("3b")!=string::npos) fprintf(pFile,"z_norm_3b_e%s rateParam bin1 zll 1 \n",year.Data());
            if (C->first.find("4b")!=string::npos) fprintf(pFile,"z_norm_4b_e%s rateParam bin1 zll 1 \n",year.Data());
          }
        } else if (C->first.find("mumu" )!=string::npos) {  
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + b#bar{b}")!=valid_procs.end())  fprintf(pFile,"tt_norm_mu%s rateParam bin1 ttbarbba 1\n",year.Data()); 
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + c#bar{c}")!=valid_procs.end())  fprintf(pFile,"tt_norm_mu%s rateParam bin1 ttbarcba 1\n",year.Data());    
          if(std::find(valid_procs.begin(), valid_procs.end(), "Z#rightarrow ll")!=valid_procs.end()){
            if (C->first.find("3b")!=string::npos) fprintf(pFile,"z_norm_3b_mu%s rateParam bin1 zll 1 \n",year.Data());
            if (C->first.find("4b")!=string::npos) fprintf(pFile,"z_norm_4b_mu%s rateParam bin1 zll 1 \n",year.Data());
          }
        } else if  (C->first.find("emu" )!=string::npos) {     
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + b#bar{b}")!=valid_procs.end())  {fprintf(pFile,"tt_norm_e%s rateParam bin1 ttbarbba 1\n",year.Data());fprintf(pFile,"tt_norm_mu%s rateParam bin1 ttbarbba 1\n",year.Data());}
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + c#bar{c}")!=valid_procs.end())  {fprintf(pFile,"tt_norm_e%s rateParam bin1 ttbarcba 1\n",year.Data());fprintf(pFile,"tt_norm_mu%s rateParam bin1 ttbarcba 1\n",year.Data());} 
        } 
      } else {
        if(C->first.find("e" )!=string::npos) {             
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + b#bar{b}")!=valid_procs.end())  fprintf(pFile,"tt_norm_e%s rateParam bin1 ttbarbba 1\n",year.Data());      
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + c#bar{c}")!=valid_procs.end())  fprintf(pFile,"tt_norm_e%s rateParam bin1 ttbarcba 1\n",year.Data());   
          if(std::find(valid_procs.begin(), valid_procs.end(), "W#rightarrow l#nu")!=valid_procs.end())  fprintf(pFile,"w_norm_e%s rateParam bin1 wlnu 1 \n",year.Data());
          //      fprintf(pFile,"w_norm_e rateParam bin1 wlnu 1 \n");     
        } else if (C->first.find("mu" )!=string::npos) {    
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + b#bar{b}")!=valid_procs.end())  fprintf(pFile,"tt_norm_mu%s rateParam bin1 ttbarbba 1\n",year.Data());            
          if(std::find(valid_procs.begin(), valid_procs.end(), "t#bar{t} + c#bar{c}")!=valid_procs.end())  fprintf(pFile,"tt_norm_mu%s rateParam bin1 ttbarcba 1\n",year.Data()); 
          if(std::find(valid_procs.begin(), valid_procs.end(), "W#rightarrow l#nu")!=valid_procs.end())  fprintf(pFile,"w_norm_mu%s rateParam bin1 wlnu 1 \n",year.Data()); 
          //      fprintf(pFile,"w_norm_mu rateParam bin1 wlnu 1 \n"); 
        }                  
      }
      
      fprintf(pFile, "-------------------------------\n");  
      if ( autoMCStats ) {
         fprintf( pFile, "* autoMCStats 5\n") ;
         fprintf(pFile, "-------------------------------\n");  
      }
      fclose(pFile);

      

    }

    FILE* pFile = fopen("combineCards"+vh_tag+".sh","w");   
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcard + " > " + "card_combined"+vh_tag+".dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + eecard       + " > " + "card_e"+vh_tag+".dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + mumucard     + " > " + "card_mu"+vh_tag+".dat").Data());
    
    fprintf(pFile,"%s;\n",(TString("sed -i '/tt_norm_mu/d' card_e"+vh_tag+".dat")).Data());   
    fprintf(pFile,"%s;\n",(TString("sed -i '/tt_norm_e/d' card_mu"+vh_tag+".dat")).Data());         

    fclose(pFile);         

    if ( verbose ) { printf(" --- verbose : AllInfo_t::buildDataCards : done.\n" ) ; fflush(stdout) ; }

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

      //      printf("\n\nProcess (get shape from file) = %s\n",proc.Data());
      //      printf("1:channelsANdShapes size= %d\n",(int)channelsAndShapes.size());    

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
        std::string xsec_str = procInfo.jsonObj["data"].daughters()[0].getString("xsec", "1.0");
        std::string delimiter = "*";
        size_t pos = 0;
        if((pos = xsec_str.find(delimiter)) != std::string::npos){
          double m1 = stod(xsec_str.substr(0,pos));
          double m2 = stod(xsec_str.erase(0,pos+delimiter.length()));
          procInfo.xsec = m1 * m2;
        }else{
          procInfo.xsec = stod(xsec_str);
        }
//        procInfo.xsec = procInfo.jsonObj["data"].daughters()[0].getDouble("xsec", 1);
//      std::cout << proc.Data() << ", xsec: " << procInfo.xsec << std::endl;
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
        TString ch_postfix= (year == "") ? ch : ch + year;

                //printf("channel= %s, bin= %s, shape name= %s, ch name= %s\n",chName.Data(), binName.Data(), shapeName.Data(), ch.Data());

        ChannelInfo_t& channelInfo = procInfo.channels[ch_postfix.Data()];
        channelInfo.bin        = binName.Data();
        channelInfo.channel    = chName.Data();
        ShapeData_t& shapeInfo = channelInfo.shapes[shapeName.Data()];

	//	printf("%s SYST SIZE=%i\n", (ch+"_"+shapeName).Data(), syst->GetNbinsX());

        for(int ivar = 1; ivar<=syst->GetNbinsX();ivar++){                
          TH1D* hshape   = NULL;

	  TString varName   = syst->GetXaxis()->GetBinLabel(ivar);
	  TString histoName = ch+"_"+shapeName+(isSignal?signalSufix:"")+varName ; 
          //if(isSignal && ivar==1)printf("Syst %i = %s\n", ivar, varName.Data()); 
	  //	  if(ivar>syst->GetNbinsX()) 
	  //printf("Histo %s  for syst:%s has Integral:%f\n", histoName.Data(), varName.Data(),);          

          TH2* hshape2D = (TH2*)pdir->Get(histoName ); 
          //      else {hshape = (TH1D*)pdir->Get(histoName ); }
          if(!hshape2D){
	    //	    printf("Histo %s is NULL for syst:%s\n", histoName.Data(), varName.Data());   
            
	    if(shapeName==histo && histoVBF!="" && ch.Contains("vbf")){   hshape2D = (TH2*)pdir->Get(TString("all_")+histoVBF+(isSignal?signalSufix:"")+varName);
            }else{                                                        hshape2D = (TH2*)pdir->Get(TString("all_")+shapeName+varName);
            }
          
            if(hshape2D){
              hshape2D->Reset();
            }else{  //if still no histo, skip this proc...
	      //	      printf("Histo %s does not exist for syst:%s\n", histoName.Data(), varName.Data());
              continue;
            }
          }

          //special treatment for side mass points
          int cutBinUsed = cutBin;
          if(shapeName == histo && !ch.Contains("vbf") && procMass==massL)cutBinUsed = indexcutML[channelInfo.bin];
          if(shapeName == histo && !ch.Contains("vbf") && procMass==massR)cutBinUsed = indexcutMR[channelInfo.bin];

          histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
          hshape   = hshape2D->ProjectionY(histoName,cutBinUsed,cutBinUsed);

	  //	  if(ivar>syst->GetNbinsX()) 
	  //	    printf("proc = %s . Histo %s  for syst:%s has Integral:%f\n", proc.Data(), histoName.Data(), varName.Data(),hshape->Integral()); 

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
	  /*
          if ( verbose ) {
             printf(" --- verbose : AllInfo_t::getShapeFromFile : " ) ;
             printf(" proc = %s ", proc.Data() ) ;
             printf(", shortName = %s ", shortName.Data() ) ;
             printf(", chName = %s ", chName.Data() ) ;
             printf(", binName = %s ", binName.Data() ) ;
             printf(", ch = %s ", ch.Data() ) ;
             printf(", ch_postfix = %s ", ch_postfix.Data() ) ;
             printf(", year = %s ", year.Data() ) ;
             printf(", histoName = %s ", histoName.Data() ) ;
             printf("\n") ;
             fflush(stdout) ;
          }
	  */
          if (postfit) {
            ///if(!(ch.Contains("_A_"))) {
              ///printf("W/Top NORMALIZATIONs: Process = %s and channel = %s\n\n",proc.Data(),ch.Data());
              if ( proc.Contains("t#bar{t} + b#bar{b}") ||
                   proc.Contains("t#bar{t} + c#bar{c}") ){   
                 double scale_factor(1.) ;
                 if ( chName.Contains("mu_") ) {
                   scale_factor = rfr_tt_norm_mu ;
                 } else {
                   scale_factor = rfr_tt_norm_e ;
                 }
                 TString hname( hshape->GetName() ) ;
                 if ( ! hname.Contains("shapes_") ) {
                    printf("     postfit:  %15s %20s %10s : %-40s : tt_norm = %6.3f\n", shortName.Data(), chName.Data(), year.Data(), hshape->GetName(), scale_factor ) ;
                    fflush(stdout) ;
                 }
                 hshape -> Scale( scale_factor ) ;
              }
              
              if ( proc.Contains("W#rightarrow l#nu") && !runZh ) {
                 double scale_factor(1.) ;
                 if ( chName.Contains("mu_") ) {
                   scale_factor = rfr_w_norm_mu ;
                 } else {
                   scale_factor = rfr_w_norm_e ;
                 }
                 TString hname( hshape->GetName() ) ;
                 if ( ! hname.Contains("shapes_") ) {
                    printf("     postfit:  %15s %20s %10s : %-40s : w_norm = %6.3f\n", shortName.Data(), chName.Data(), year.Data(), hshape->GetName(), scale_factor ) ;
                    fflush(stdout) ;
                 }
                 hshape -> Scale( scale_factor ) ;
              }

              if ( runZh && proc.Contains("Z#rightarrow ll") ) {
                 double scale_factor(1.) ;
                 if ( chName.Contains("mu_") ) {
                    if(ch.Contains("3b")) { scale_factor = rfr_z_norm_3b_mu ; }
                    if(ch.Contains("4b")) { scale_factor = rfr_z_norm_4b_mu ; }
                 } else {
                    if(ch.Contains("3b")) { scale_factor = rfr_z_norm_3b_e ; }
                    if(ch.Contains("4b")) { scale_factor = rfr_z_norm_4b_e ; }
                 }
                 TString hname( hshape->GetName() ) ;
                 if ( ! hname.Contains("shapes_") ) {
                    printf("     postfit:  %15s %20s %10s : %-40s : z_norm = %6.3f\n", shortName.Data(), chName.Data(), year.Data(), hshape->GetName(), scale_factor ) ;
                    fflush(stdout) ;
                 }
                 hshape -> Scale( scale_factor ) ;
              }

            ///} // _A_?

          } // postfit?

          if(isSignal)hshape->Scale(SignalRescale);


          //Do Renaming and cleaning
          varName.ReplaceAll("down","Down");
          varName.ReplaceAll("up","Up");

          if(varName==""){//does nothing
            //    }else if(varName.EndsWith("_jes")){varName.ReplaceAll("_jes","_CMS_scale_j");
            //    }else if(varName.BeginsWith("_umet")) { continue; //skip this one for now
          }else if(varName.BeginsWith("_jer")){varName.ReplaceAll("_jer","_CMS_res_j"); // continue;//skip res for now
          }else if(varName.BeginsWith("_les")){
            continue; // skip this one for now
            //      if(ch.Contains("e"  ))varName.ReplaceAll("_les","_CMS_scale_e");
            //      if(ch.Contains("mu"))varName.ReplaceAll("_les","_CMS_scale_m");
          }else if(varName.BeginsWith("_btag"  )){varName.ReplaceAll("_btag","_CMS_eff_b");
          }else if(varName.BeginsWith("_ctag"  )){varName.ReplaceAll("_ctag","_CMS_eff_c");
          }else if(varName.BeginsWith("_ltag"  )){varName.ReplaceAll("_ltag","_CMS_eff_mistag");
            //          }else if(varName.BeginsWith("_pu"    )){varName.ReplaceAll("_pu", "_CMS_haa4b_pu");
            //    }else if(varName.BeginsWith("_pdf" )){
            //      if (proc.Contains("wh")!=std::string::npos) {continue; }
          }else if(varName.BeginsWith("_bnorm"  )){continue; //skip this one
          }else{ varName="_CMS_haa4b"+varName;}

          hshape->SetTitle(proc+varName);
          if(shapeInfo.uncShape.find(varName.Data())==shapeInfo.uncShape.end()){
            shapeInfo.uncShape[varName.Data()] = hshape;
          }else{
            shapeInfo.uncShape[varName.Data()]->Add(hshape);
          }
        } // end ivar
      }
    }

  } // getShapeFromFile



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

            //      printf("Now going to rebinMainHisto of: %s\n",jetBin.Data());

            if(jetBin.Contains("3b")){

             //-----------------
              //double xbins[] = {-0.31, -0.09, 0.09, 0.19, 0.21, 0.35};   // last bins from Yuan
              //if(runZh)  {xbins[0]=-0.31;xbins[1]=-0.09;xbins[2]=0.09;xbins[3]=0.17; xbins[4]=0.19;xbins[5]=0.35;}  // last bins from Yuan
             //-----------------
              ///double xbins[] = {-0.31, -0.25, -0.07, 0.12, 0.19, 0.35};   // july 21, new trial bins
              ///if(runZh)  {xbins[0]=-0.31;xbins[1]=-0.09;xbins[2]=0.00;xbins[3]=0.09; xbins[4]=0.17;xbins[5]=0.35;}  // july 21, new trial bins
             //-----------------
              double xbins[] = {-0.31, -0.12,  0.00, 0.12, 0.19, 0.35};   // july 22, new trial bins
              if(runZh)  {xbins[0]=-0.31;xbins[1]=-0.09;xbins[2]=0.00;xbins[3]=0.09; xbins[4]=0.17;xbins[5]=0.35;}  // july 22, (same as 21), new trial bins
             //-----------------

              int nbins=sizeof(xbins)/sizeof(double);    
              unc->second = histo->Rebin(nbins-1, histo->GetName(), (double*)xbins);  
              utils::root::fixExtremities(unc->second, false, true); 
            }else if(jetBin.Contains("4b")){ 

             //-----------------
              ///double xbins[] = {-0.31, -0.09, 0.07, 0.11, 0.17, 0.35};     // last bins from Yuan
              ///if(runZh) {xbins[0]=-0.31;xbins[1]=-0.09;xbins[2]=0.09;xbins[3]=0.15; xbins[4]=0.21;xbins[5]=0.35;}     // last bins from Yuan
             //-----------------
              ///double xbins[] = {-0.31, -0.20, -0.10, -0.02, 0.10, 0.35};  // july 21, new trial bins
              ///if(runZh) {xbins[0]=-0.31;xbins[1]=-0.15;xbins[2]=-0.07;xbins[3]=0.00; xbins[4]=0.09;xbins[5]=0.35;}  // july 21, new trial bins
             //-----------------
              double xbins[] = {-0.31, -0.14, -0.06,  0.02, 0.10, 0.35};  // july 22, new trial bins
              if(runZh) {xbins[0]=-0.31;xbins[1]=-0.15;xbins[2]=-0.07;xbins[3]=0.00; xbins[4]=0.09;xbins[5]=0.35;}  // july 22, (same as 21), new trial bins
             //-----------------
              
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
            addChannel(ch->second, ch2->second, false); // add nominal 
            addChannel(ch->second, ch2->second, true, false); // add systematics
            it->second.channels.erase(ch2);  
            ch2=ch;
          }
          ch->second.bin = NewName;
        }

        //also update the map keys
        std::map<string, ChannelInfo_t> newMap;
        for(std::map<string, ChannelInfo_t>::iterator ch = it->second.channels.begin(); ch!=it->second.channels.end(); ch++){
          //newMap[ch->second.channel+ch->second.bin] = ch->second;
          newMap[ch->second.channel+"_"+ch->second.bin] = ch->second;
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
	if(it->second.isData)continue; //only do this for MC
	//        if(it->second.isData && procName!="Instr. MET")continue; //only do this for MC or InstrMET (which also contains MC and negative bins)
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
              if(unc->second->GetBinContent(binx)<=0){unc->second->SetBinContent(binx, 1E-6); };
            }
          }
        }
      }

    //recompute total background
    computeTotalBackground();
  }


  void AllInfo_t::printInventory() {

     printf("\n\n\n =============== AllInfo_t::printInventory :  begin\n\n") ;

     for ( std::map<string, ProcessInfo_t>::iterator ip = procs.begin(); ip!= procs.end(); ip++ ) {

        string proc_key = ip -> first ;
        ProcessInfo_t proc = ip -> second ;

        proc.printProcess() ;

     } // ip

     printf("\n\n\n =============== AllInfo_t::printInventory :  end\n\n") ;

  } // AllInfo_t::printInventory






