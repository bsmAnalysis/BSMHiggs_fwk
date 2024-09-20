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
  evSummary_.isSR1=false;
  evSummary_.isSR2=false;
  evSummary_.Njets=0;
  evSummary_.Nlepton=0;
  evSummary_.HT=0;
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
  evSummary_.DRbbav=0;
  evSummary_.DeltaMassbbMin=0;
  evSummary_.mbbuntag=0;
  evSummary_.btag1=0;
  evSummary_.btag2=0;
  evSummary_.btag3=0;
  evSummary_.weight=0;

}


void MVAHandler::getEntry(
			  bool isSR1,
			  bool isSR2,
			  float Njets,
			  float Nlepton,
			  float HT,
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
			  float DRbbav,
			  float DeltaMassbbMin,
			  float mbbuntag,
			  float btag1,
			  float btag2,
			  float btag3,
			  float weight
			  ) 
 
{
  
  resetStruct();

  evSummary_.isSR1=isSR1;
  evSummary_.isSR2=isSR2;
  evSummary_.Njets=Njets;
  evSummary_.Nlepton=Nlepton;
  evSummary_.HT=HT;
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
  evSummary_.DRbbav=DRbbav;
  evSummary_.DeltaMassbbMin=DeltaMassbbMin;
  evSummary_.mbbuntag=mbbuntag;
  evSummary_.btag1=btag1;
  evSummary_.btag2=btag2;
  evSummary_.btag3=btag3;
  evSummary_.weight=weight;


  return ;
}

//
bool MVAHandler::initTree()
{
  tSR0 = new TTree("tSR0","trMVA");
  tSR0->SetDirectory(0);
  std::cout << "Check init tree" << std::endl;
  
  tSR0->Branch("Njets",&evSummary_.Njets,"Njets/F");
  tSR0->Branch("Nlepton",&evSummary_.Nlepton,"Nlepton/F");
  tSR0->Branch("HT",&evSummary_.HT,"HT/F");
  tSR0->Branch("H_pt",&evSummary_.H_pt,"H_pt/F");
  tSR0->Branch("W_pt",&evSummary_.W_pt,"W_pt/F");
  tSR0->Branch("lepton_pt",&evSummary_.lepton_pt,"lepton_pt/F");
  tSR0->Branch("bjet1_pt",&evSummary_.bjet1_pt,"bjet1_pt/F");
  tSR0->Branch("MET",&evSummary_.MET,"MET/F");  
  tSR0->Branch("H_M",&evSummary_.H_M,"H_M/F");
  tSR0->Branch("MTW",&evSummary_.MTW,"MTW/F");
  tSR0->Branch("dphiWH",&evSummary_.dphiWH,"dphiWH/F");
  tSR0->Branch("dphiMetLepton",&evSummary_.dphiMetLepton,"dphiMetLepton/F");
  tSR0->Branch("dphiMetJetMin",&evSummary_.dphiMetJetMin,"dphiMetJetMin/F");
  tSR0->Branch("DRbbav",&evSummary_.DRbbav,"DRbbav/F");
  tSR0->Branch("DeltaMassbbMin",&evSummary_.DeltaMassbbMin,"DeltaMassbbMin/F"); 
  tSR0->Branch("btag1",&evSummary_.btag1,"btag1/F");
  tSR0->Branch("btag2",&evSummary_.btag2,"btag2/F");
  tSR0->Branch("btag3",&evSummary_.btag3,"btag3/F");
  tSR0->Branch("weight",&evSummary_.weight,"weight/F");
  

  tSR1 = new TTree("tSR1","trMVA1");
  tSR1->SetDirectory(0);
  std::cout << "Check init tree for SR1" << std::endl;
  
  tSR1->Branch("Njets",&evSummary_.Njets,"Njets/F");
  tSR1->Branch("Nlepton",&evSummary_.Nlepton,"Nlepton/F");
  tSR1->Branch("HT",&evSummary_.HT,"HT/F");
  tSR1->Branch("H_pt",&evSummary_.H_pt,"H_pt/F");
  tSR1->Branch("W_pt",&evSummary_.W_pt,"W_pt/F");
  tSR1->Branch("lepton_pt",&evSummary_.lepton_pt,"lepton_pt/F");
  tSR1->Branch("bjet1_pt",&evSummary_.bjet1_pt,"bjet1_pt/F");
  tSR1->Branch("MET",&evSummary_.MET,"MET/F");  
  tSR1->Branch("H_M",&evSummary_.H_M,"H_M/F");
  tSR1->Branch("MTW",&evSummary_.MTW,"MTW/F");
  tSR1->Branch("dphiWH",&evSummary_.dphiWH,"dphiWH/F");
  tSR1->Branch("dphiMetLepton",&evSummary_.dphiMetLepton,"dphiMetLepton/F");
  tSR1->Branch("dphiMetJetMin",&evSummary_.dphiMetJetMin,"dphiMetJetMin/F");
  tSR1->Branch("DRbbav",&evSummary_.DRbbav,"DRbbav/F");
  tSR1->Branch("DeltaMassbbMin",&evSummary_.DeltaMassbbMin,"DeltaMassbbMin/F");
  tSR1->Branch("btag1",&evSummary_.btag1,"btag1/F");
  tSR1->Branch("btag2",&evSummary_.btag2,"btag2/F");
  tSR1->Branch("btag3",&evSummary_.btag3,"btag3/F");
  tSR1->Branch("weight",&evSummary_.weight,"weight/F");
                                                                                      
  
  tSR2 = new TTree("tSR2","trMVA2");
  tSR2->SetDirectory(0);
  std::cout << "check init tree for SR2" << std::endl;
  
  tSR2->Branch("Njets",&evSummary_.Njets,"Njets/F");
  tSR2->Branch("Nlepton",&evSummary_.Nlepton,"Nlepton/F");
  tSR2->Branch("HT",&evSummary_.HT,"HT/F");
  tSR2->Branch("H_pt",&evSummary_.H_pt,"H_pt/F");
  tSR2->Branch("W_pt",&evSummary_.W_pt,"W_pt/F");
  tSR2->Branch("lepton_pt",&evSummary_.lepton_pt,"lepton_pt/F");
  tSR2->Branch("bjet1_pt",&evSummary_.bjet1_pt,"bjet1_pt/F");
  tSR2->Branch("MET",&evSummary_.MET,"MET/F");  
  tSR2->Branch("H_M",&evSummary_.H_M,"H_M/F");
  tSR2->Branch("MTW",&evSummary_.MTW,"MTW/F");
  tSR2->Branch("dphiWH",&evSummary_.dphiWH,"dphiWH/F");
  tSR2->Branch("dphiMetLepton",&evSummary_.dphiMetLepton,"dphiMetLepton/F");
  tSR2->Branch("dphiMetJetMin",&evSummary_.dphiMetJetMin,"dphiMetJetMin/F");
  tSR2->Branch("DRbbav",&evSummary_.DRbbav,"DRbbav/F");
  tSR2->Branch("DeltaMassbbMin",&evSummary_.DeltaMassbbMin,"DeltaMassbbMin/F");
  tSR2->Branch("btag1",&evSummary_.btag1,"btag1/F");
  tSR2->Branch("btag2",&evSummary_.btag2,"btag2/F");
  tSR2->Branch("btag3",&evSummary_.btag3,"btag3/F");
  tSR2->Branch("weight",&evSummary_.weight,"weight/F");
                                                                                       
  std::cout << "Branches check!" << std::endl;
  return true;
}

//
void MVAHandler::fillTree()
{
  tSR0->Fill();
  
  if(evSummary_.isSR1)
    {
      tSR1->Fill();
    }

  else if(evSummary_.isSR2)
    {
      tSR2->Fill();
    }
  
  std::cout << "Tree fill check!" << std::endl;
  return;   
}

void MVAHandler::writeTree(TString mvaout)
{
  TFile *MVAofile=TFile::Open( mvaout, "recreate");

  tSR0->SetDirectory(0);
  tSR0->Write();
  
  tSR1->SetDirectory(0);
  tSR1->Write();

  tSR2->SetDirectory(0);
  tSR2->Write();
 
  MVAofile->Close();
  std::cout << "Write tree check!" << std::endl;
 
  return;
}
//
MVAHandler::~MVAHandler()
{
}
