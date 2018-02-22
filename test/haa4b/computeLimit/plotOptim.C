#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"

#include "UserCode/bsmhiggs_fwk/interface/tdrstyle.h"
#include "UserCode/bsmhiggs_fwk/src/tdrstyle.C"
#include "UserCode/bsmhiggs_fwk/interface/RootUtils.h"

TGraph* multiplyGraph(TGraph* A, TGraph* B){   double x, y;   for(int i=0;i<A->GetN();i++){A->GetPoint(i, x, y); A->SetPoint(i, x, y*B->Eval(x));}    return A;   }
TGraph* divideGraph (TGraph* A, TGraph* B){ double x, y; for(int i=0;i<A->GetN();i++){A->GetPoint(i, x, y); A->SetPoint(i, x, y/B->Eval(x));} return A; } 


TGraph* makeGraphFromTxt(std::string dataFile, int cutIndex){

  FILE* pFile = fopen(dataFile.c_str(), "r");
  if(!pFile){ printf("Couldn't open file %s to read limit values\n", dataFile.c_str()); return NULL; }

  TGraph* toReturn = new TGraph(999); int N=0;
  char line [4096];

  stringstream ss;
  ss << cutIndex;
  std::string cut = ss.str();

  while(fgets(line, 4096, pFile)){
    if(std::string(line).find("------")==0)continue; //skip line starting by //

    char str0[10];
    char cmass[5];char icmass[2];
    char climit[10];char iclimit[10];
    char cindex[5];char icindex[5];

    int icut;
    int mA; double rlimit;

    //Read formatted data from string (line)
    sscanf(line,"%s %s %s %s %s %s %s",&cmass[0],&icmass[1],&str0[0],&climit[0],&iclimit[0],&cindex[0],&icindex[0]);

    sscanf(&icmass[1],"%d",&mA);
    sscanf(&iclimit[0],"%lf",&rlimit);
    sscanf(&icindex[0],"%d",&icut);

    if (icut!=cutIndex) continue; // skip lines of Index != cutIndex
    //    printf ("Mass %d -> limit: %lf , Index: %d\n",mA,rlimit,icut);
    toReturn->SetPoint(N, (double)mA, rlimit); N++; 

  }fclose(pFile);
  toReturn->Set(N);
  return toReturn;

}

void scaleGraph(TGraph* Limit, double scale){ 
   for(int i=0;i<Limit->GetN();i++){ 
      Limit->SetPoint(i, Limit->GetX()[i], Limit->GetY()[i]*scale);
   }   
}

void printLimits(FILE* pFile, TGraph* graph, double Mmin=12, double Mmax=60){
   double previous = graph->Eval(Mmin);
   fprintf(pFile, "Exclude ");
   bool opened = false;
   for(double M=Mmin;M<=Mmax;M+=1.0){
      double NEW = graph->Eval(M);
      //printf("%04f - %6.3f\n", M, NEW);
      if(previous>1 && NEW<=1){fprintf(pFile,"[%f,",M); opened=true;}
      if(previous<1 && NEW>=1){fprintf(pFile,"%f]",M-1); opened=false;}
      if(M==Mmax    && opened){fprintf(pFile,">%f]",Mmax); opened=false;}
      previous = NEW;
   }fprintf(pFile,"\n");
}



void plotOptim(string outputDir="./", string inputs="", TString inputXSec="", bool strengthLimit=true, bool blind=false, double energy=7, double luminosity=5.035, TString legendName="Wh channels")
{
   setTDRStyle();  
   gStyle->SetPadTopMargin   (0.05);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);
   
   TCanvas *c = new TCanvas("c", "c",1000,700);  ;
   c->SetGridx();
   c->SetGridy();

   TH1F* framework = new TH1F(inputs.c_str(),"Graph",1,10,60); //3000);
   framework->SetStats(false);
   framework->SetTitle("");
   framework->GetXaxis()->SetTitle("M_{a} [GeV]");
   framework->GetYaxis()->SetTitleOffset(1.70);

   framework->GetYaxis()->SetTitle("#mu = #sigma_{95%} / #sigma_{th}");
   framework->GetYaxis()->SetRangeUser(1E-2,1E3);
   c->SetLogy(true);

   framework->GetXaxis()->SetLabelOffset(0.007);
   framework->GetXaxis()->SetLabelSize(0.03);
   framework->GetXaxis()->SetTitleOffset(1.0);
   framework->GetXaxis()->SetTitleFont(42);
   framework->GetXaxis()->SetTitleSize(0.035);
   framework->GetYaxis()->SetLabelFont(42);
   framework->GetYaxis()->SetLabelOffset(0.007);
   framework->GetYaxis()->SetLabelSize(0.03);
   framework->GetYaxis()->SetTitleOffset(1.3);
   framework->GetYaxis()->SetTitleFont(42);
   framework->GetYaxis()->SetTitleSize(0.035);
   framework->Draw();

   TLine* SMLine = new TLine(framework->GetXaxis()->GetXmin(),1.0,framework->GetXaxis()->GetXmax(),1.0);
   SMLine->SetLineWidth(2); SMLine->SetLineStyle(1); SMLine->SetLineColor(4);      
   SMLine->Draw("same C");

   TLegend* LEG = new TLegend(0.4,0.55,0.85,0.95);
   LEG->SetHeader("");
   LEG->SetFillColor(0);
   LEG->SetFillStyle(0);
   LEG->SetTextFont(42);
   LEG->SetBorderSize(0);

   TString ival[12]={"-0.4","-0.05","-0.025","0.0","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2"};

   TGraph* gLimit[12];
 
   stringstream ss; 
  
   int i(0);
   //scan through the BDT cut indices
   for (int icut=1; icut<12; icut++) {
     
     gLimit[icut] = makeGraphFromTxt(inputs, icut);
     gLimit[icut]->SetMarkerColor(icut);
     gLimit[icut]->SetLineColor(icut);
     
     gLimit[icut]->Draw("same C");
     
     ss << icut;
     std::string cut = ss.str();  
     LEG->AddEntry(gLimit[icut], "median expected; BDT cut= "+ival[i+1],"L");
     i++;
   }
   
   LEG->SetTextSize(0.025);

   LEG->Draw();
   c->RedrawAxis();
   
   inputs.erase(inputs.length()-4);
   c->SaveAs((inputs+"_scan.pdf").c_str()); 
   c->SaveAs((inputs+"_scan.png").c_str()); 
}




