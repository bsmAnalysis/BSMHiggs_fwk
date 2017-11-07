#include "UserCode/bsmhiggs_fwk/interface/MVAHandler.h"

//
MVAHandler::MVAHandler()
{
}

//
void MVAHandler::resetStruct()
{
  //catagory type
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
}

//
bool MVAHandler::initTree(TTree *t3b, TTree *t4b)
{
  if ( t3b == 0 ) return false;
  to3b = t3b;

  to3b->Branch("WpT",  &evSummary_.WpT,  "WpT/F");
  to3b->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");
  to3b->Branch("HpT",  &evSummary_.HpT,  "HpT/F");
  to3b->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F");
  to3b->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");
  to3b->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  to3b->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F");

  if ( t4b == 0 ) return false;
  to4b = t4b;

  to4b->Branch("WpT",  &evSummary_.WpT,  "WpT/F");
  to4b->Branch("Hmass",  &evSummary_.Hmass,  "Hmass/F");
  to4b->Branch("HpT",  &evSummary_.HpT,  "HpT/F");
  to4b->Branch("bbdRAve",  &evSummary_.bbdRAve,  "bbdRAve/F");
  to4b->Branch("bbdMMin",  &evSummary_.bbdMMin,  "bbdMMin/F");
  to4b->Branch("HHt",  &evSummary_.HHt,  "HHt/F");
  to4b->Branch("WHdR",  &evSummary_.WHdR,  "WHdR/F");

  return true;
}

//
void MVAHandler::fillTree()
{
  if ( evSummary_.is3b && !evSummary_.is4b )
  {
    std::cout << "One event can not be both in 3 and 4 b cat! Please check!" << std::endl;
  }
  else if ( evSummary_.is3b && !evSummary_.is4b )
  {
    if ( to3b ) to3b->Fill();
    return ;
  }
  else if ( !evSummary_.is3b && evSummary_.is4b )
  {
    if ( to4b ) to4b->Fill();
    return ;
  }
  else return ;
  return ;
}

//
MVAHandler::~MVAHandler()
{
}
