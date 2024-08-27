#include "UserCode/bsmhiggs_fwk/interface/MVAHandler_wh.h"
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
  evSummary_.Njets=0;
  evSummary_.HTjets=0;
  evSummary_.H_pt=0;
  evSummary_.W_pt=0;
  evSummary_.lepton_pt=0;
  evSummary_.bjet1_pt=0;
  evSummary_.MET=0;
  evSummary_.H_M=0;
  evSummary_.MTW=0;
  evSummary_.dphiWH=0;
  evSummary_.dphiMetLepton=0;
  evSummary_.dphiMetJetMin=0;
  evSummary_.btag1=0;
  evSummary_.btag2=0;
  evSummary_.btag3=0;
  evSummary_.weight=0;

}


void MVAHandler::getEntry(
			  float Njets,
			  float HTjets,
			  float H_pt,
			  float W_pt,
			  float lepton_pt,
			  float bjet1_pt,
			  float MET,
			  float H_M,
			  float MTW,
			  float dphiWH,
			  float dphiMetLepton,
			  float dphiMetJetMin,
			  float btag1,
			  float btag2,
			  float btag3,
			  float weight
			  ) 
 
{
  
  resetStruct();
  
  evSummary_.Njets=Njets;
  evSummary_.HTjets=HTjets;
  evSummary_.H_pt=H_pt;
  evSummary_.W_pt=W_pt;
  evSummary_.lepton_pt=lepton_pt;
  evSummary_.bjet1_pt=bjet1_pt;
  evSummary_.MET=MET;
  evSummary_.H_M=H_M;
  evSummary_.MTW=MTW;
  evSummary_.dphiWH=dphiWH;
  evSummary_.dphiMetLepton=dphiMetLepton;
  evSummary_.dphiMetJetMin=dphiMetJetMin;
  evSummary_.btag1=btag1;
  evSummary_.btag2=btag2;
  evSummary_.btag3=btag3;
  evSummary_.weight=weight;


  return ;
}

//
bool MVAHandler::initTree()
{
  
  t1 = new TTree("t1","trMVA");
  t1->SetDirectory(0);
  std::cout << "check init tree"<<std::endl;
  
  t1->Branch("Njets",&evSummary_.Njets,"Njets/F");
  t1->Branch("HTjets",&evSummary_.HTjets,"HTjets/F");
  t1->Branch("H_pt",&evSummary_.H_pt,"H_pt/F");
  t1->Branch("W_pt",&evSummary_.W_pt,"W_pt/F");
  t1->Branch("lepton_pt",&evSummary_.lepton_pt,"lepton_pt/F");
  t1->Branch("bjet1_pt",&evSummary_.bjet1_pt,"bjet1_pt/F");
  t1->Branch("MET",&evSummary_.MET,"MET/F");  
  t1->Branch("H_M",&evSummary_.H_M,"H_M/F");
  t1->Branch("MTW",&evSummary_.MTW,"MTW/F");
  t1->Branch("dphiWH",&evSummary_.dphiWH,"dphiWH/F");
  t1->Branch("dphiMetLepton",&evSummary_.dphiMetLepton,"dphiMetLepton/F");
  t1->Branch("dphiMetJetMin",&evSummary_.dphiMetJetMin,"dphiMetJetMin/F");
  t1->Branch("btag1",&evSummary_.btag1,"btag1/F");
  t1->Branch("btag2",&evSummary_.btag2,"btag2/F");
  t1->Branch("btag3",&evSummary_.btag3,"btag3/F");
  t1->Branch("weight",&evSummary_.weight,"weight/F");
                                                                                       
  std::cout << "branches check!" << std::endl;
  return true;
}

//
void MVAHandler::fillTree()
{
  t1->Fill();
  
  std::cout << "tree fill check!" << std::endl;

  return;   
}

void MVAHandler::writeTree(TString mvaout)
{
  TFile *MVAofile=TFile::Open( mvaout, "recreate");  
  t1->SetDirectory(0);
  t1->Write();
 
  MVAofile->Close();
  std::cout << "write tree check!" << std::endl;
 
  return;
}
//
MVAHandler::~MVAHandler()
{
}
