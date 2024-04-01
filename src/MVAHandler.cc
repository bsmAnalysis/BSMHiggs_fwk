#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"

//
MVAHandler::MVAHandler()
{
}

//
MVAEvtContainer &MVAHandler::getEvent()
{
  return evSummary_;
}

//
void MVAHandler::resetStruct()
{
  //catagory type
  evSummary_.isEven = false;
  //  evSummary_.isOdd = false;
  evSummary_.is3b = false;
  evSummary_.is4b = false; 
  //W boson related only related var
  evSummary_.WpT = -1.0;
  //Higgs boson related only var
  evSummary_.Hmass = -1.0; 
  evSummary_.HpT = -1.0; 
  evSummary_.bbdRAve = -1.0; 
  evSummary_.bbdMMin = -1.0;
  evSummary_.HHt = -1.0;
  //dr W and Higgs 
  evSummary_.WHdR = -1.0;
  evSummary_.lepPt = -1.0;
  evSummary_.pfMET = -1.0;
  evSummary_.MTw = -1.0;
  evSummary_.ljDR = -1.0;
  //weight
  evSummary_.weight = 0.0;
  evSummary_.xsecWeight = 1.0;
  evSummary_.lheNJets = -1;
}

//
void MVAHandler::getEntry( 
			  bool isEven, 
                          bool is3b, bool is4b, 
                          float Wpt, //W only
                          float Hmass, float HpT, float bbdRAve, float bbdMMin, float HHt, //Higgs only
                          float WHdR, //W and H
			  float lepPt, float pfMET, float MTw, 
			  float ljDR,
                          float weight, float xsecWeight,
                          int lheNJets
                         )
{
  resetStruct();
  evSummary_.isEven = isEven;
  //  evSummary_.isOdd = isOdd;
  //catagory type
  evSummary_.is3b = is3b;
  evSummary_.is4b = is4b;
  //W boson related only related var
  evSummary_.WpT = Wpt;
  //Higgs boson related only var
  evSummary_.Hmass = Hmass;
  evSummary_.HpT = HpT;
  evSummary_.bbdRAve = bbdRAve;
  evSummary_.bbdMMin = bbdMMin;
  evSummary_.HHt = HHt;
  //dr W and Higgs 
  evSummary_.WHdR = WHdR;
  // few additions
  evSummary_.lepPt = lepPt;
  evSummary_.pfMET = pfMET;
  evSummary_.MTw = MTw;
  evSummary_.ljDR = ljDR;
  //weight
  evSummary_.weight = weight;
  evSummary_.xsecWeight = xsecWeight;  
  //AUX
  evSummary_.lheNJets = lheNJets;
  return ;
}

//
bool MVAHandler::initTree(TString mvaout)
{
  //write mode, to mva tree
  MVAofile = TFile::Open( mvaout, "recreate");
  to3b_e = new TTree("TribMVA_e","TribMVA");
  to4b_e = new TTree("QuabMVA_e","QuabMVA");
  toSignal_e = new TTree("H4bMVA_e","H4bMVA");

  to3b_o = new TTree("TribMVA_o","TribMVA");   
  to4b_o = new TTree("QuabMVA_o","QuabMVA"); 
  toSignal_o = new TTree("H4bMVA_o","H4bMVA");  

  toSignal_e->Branch("WpT",  &evSummary_.WpT,  "WpT/F"); 
  toSignal_e->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");   
  toSignal_e->Branch("HpT",  &evSummary_.HpT,  "HpT/F");      
  toSignal_e->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F"); 
  toSignal_e->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");  
  toSignal_e->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  toSignal_e->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F"); 
  toSignal_e->Branch("lepPt", &evSummary_.lepPt, "lepPt/F");
  toSignal_e->Branch("pfMET", &evSummary_.pfMET, "pfMET/F"); 
  toSignal_e->Branch("MTw", &evSummary_.MTw, "MTw/F");
  toSignal_e->Branch("ljDR", &evSummary_.ljDR, "ljDR/F"); 
  toSignal_e->Branch("weight",  &evSummary_.weight,  "weight/F");
  toSignal_e->Branch("xsecWeight",  &evSummary_.xsecWeight,  "xsecWeight/F"); 
  toSignal_e->Branch("lheNJets",  &evSummary_.lheNJets,  "lheNJets/I"); 

  toSignal_o->Branch("WpT",  &evSummary_.WpT,  "WpT/F"); 
  toSignal_o->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");   
  toSignal_o->Branch("HpT",  &evSummary_.HpT,  "HpT/F");      
  toSignal_o->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F"); 
  toSignal_o->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");  
  toSignal_o->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  toSignal_o->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F"); 
  toSignal_o->Branch("lepPt", &evSummary_.lepPt, "lepPt/F");
  toSignal_o->Branch("pfMET", &evSummary_.pfMET, "pfMET/F"); 
  toSignal_o->Branch("MTw", &evSummary_.MTw, "MTw/F");
  toSignal_o->Branch("ljDR", &evSummary_.ljDR, "ljDR/F"); 
  toSignal_o->Branch("weight",  &evSummary_.weight,  "weight/F");
  toSignal_o->Branch("xsecWeight",  &evSummary_.xsecWeight,  "xsecWeight/F"); 
  toSignal_o->Branch("lheNJets",  &evSummary_.lheNJets,  "lheNJets/I"); 
  
  to3b_e->Branch("WpT",  &evSummary_.WpT,  "WpT/F");
  to3b_e->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");
  to3b_e->Branch("HpT",  &evSummary_.HpT,  "HpT/F");
  to3b_e->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F");
  //to3b_e->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");
  to3b_e->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  to3b_e->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F");
  to3b_e->Branch("lepPt", &evSummary_.lepPt, "lepPt/F");
  to3b_e->Branch("pfMET", &evSummary_.pfMET, "pfMET/F");
  to3b_e->Branch("MTw", &evSummary_.MTw, "MTw/F");  
  to3b_e->Branch("ljDR", &evSummary_.ljDR, "ljDR/F");  
  to3b_e->Branch("weight",  &evSummary_.weight,  "weight/F");
  to3b_e->Branch("xsecWeight",  &evSummary_.xsecWeight,  "xsecWeight/F");    
  to3b_e->Branch("lheNJets",  &evSummary_.lheNJets,  "lheNJets/I");

  to3b_o->Branch("WpT",  &evSummary_.WpT,  "WpT/F");
  to3b_o->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");
  to3b_o->Branch("HpT",  &evSummary_.HpT,  "HpT/F");
  to3b_o->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F");
  //to3b_o->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");
  to3b_o->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  to3b_o->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F");
  to3b_o->Branch("lepPt", &evSummary_.lepPt, "lepPt/F");
  to3b_o->Branch("pfMET", &evSummary_.pfMET, "pfMET/F");
  to3b_o->Branch("MTw", &evSummary_.MTw, "MTw/F");  
  to3b_o->Branch("ljDR", &evSummary_.ljDR, "ljDR/F");  
  to3b_o->Branch("weight",  &evSummary_.weight,  "weight/F");
  to3b_o->Branch("xsecWeight",  &evSummary_.xsecWeight,  "xsecWeight/F");    
  to3b_o->Branch("lheNJets",  &evSummary_.lheNJets,  "lheNJets/I");

  to4b_e->Branch("WpT",  &evSummary_.WpT,  "WpT/F");
  to4b_e->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");
  to4b_e->Branch("HpT",  &evSummary_.HpT,  "HpT/F");
  to4b_e->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F");
  to4b_e->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");
  to4b_e->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  to4b_e->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F");
  to4b_e->Branch("lepPt", &evSummary_.lepPt, "lepPt/F");   
  to4b_e->Branch("pfMET", &evSummary_.pfMET, "pfMET/F");      
  to4b_e->Branch("MTw", &evSummary_.MTw, "MTw/F"); 
  to4b_e->Branch("ljDR", &evSummary_.ljDR, "ljDR/F");  
  to4b_e->Branch("weight",  &evSummary_.weight,  "weight/F");
  to4b_e->Branch("xsecWeight",  &evSummary_.xsecWeight,  "xsecWeight/F"); 
  to4b_e->Branch("lheNJets",  &evSummary_.lheNJets,  "lheNJets/I");

  to4b_o->Branch("WpT",  &evSummary_.WpT,  "WpT/F");
  to4b_o->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");
  to4b_o->Branch("HpT",  &evSummary_.HpT,  "HpT/F");
  to4b_o->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F");
  to4b_o->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");
  to4b_o->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  to4b_o->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F");
  to4b_o->Branch("lepPt", &evSummary_.lepPt, "lepPt/F");   
  to4b_o->Branch("pfMET", &evSummary_.pfMET, "pfMET/F");      
  to4b_o->Branch("MTw", &evSummary_.MTw, "MTw/F"); 
  to4b_o->Branch("ljDR", &evSummary_.ljDR, "ljDR/F");  
  to4b_o->Branch("weight",  &evSummary_.weight,  "weight/F");
  to4b_o->Branch("xsecWeight",  &evSummary_.xsecWeight,  "xsecWeight/F"); 
  to4b_o->Branch("lheNJets",  &evSummary_.lheNJets,  "lheNJets/I");
  
  return true;
}

//
void MVAHandler::fillTree()
{

  if ( evSummary_.is3b || evSummary_.is4b )
    {
      if (evSummary_.isEven) {
	if ( toSignal_e ) toSignal_e->Fill();
      } else {
	//      if (evSummary_.isOdd) {
	if ( toSignal_o ) toSignal_o->Fill();
      }
    }

  // Fill MVA separately for 3b and 4b
  if ( evSummary_.is3b && evSummary_.is4b )
  {
    std::cout << "One event can not be both in 3 and 4 b cat! Please check!" << std::endl;
  }
  else if ( evSummary_.is3b && !evSummary_.is4b ) 
  {
    // if ( to3b_e ) to3b_e->Fill();
     if (evSummary_.isEven) {
	if ( to3b_e ) to3b_e->Fill();
     } else {
       //      if (evSummary_.isOdd) {
	if ( to3b_o ) to3b_o->Fill();
     }
    return ;
  }
  else if ( !evSummary_.is3b && evSummary_.is4b )
  {
    //  if ( to4b_e ) to4b_e->Fill();
     if (evSummary_.isEven) {
	if ( to4b_e ) to4b_e->Fill();
     } else {
       //      if (evSummary_.isOdd) {
       if ( to4b_o ) to4b_o->Fill();
     }
    return ;
  }
  else return ;
  return ;
}

//
void MVAHandler::writeTree()
{
  //TFile *MVAofile=TFile::Open( outURL, "recreate");  
  toSignal_e->Write(); toSignal_o->Write();
  to3b_e->Write(); to3b_o->Write();
  to4b_e->Write(); to4b_o->Write();
  MVAofile->Close();
  return ;
}
//
MVAHandler::~MVAHandler()
{
}
