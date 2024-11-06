#include "UserCode/bsmhiggs_fwk/interface/MVAHandler_vbfh.h"
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
  evSummary_.higgs_m=0;
  evSummary_.higgs_pt=0;
  evSummary_.higgs_eta=0;
  evSummary_.costheta0=0;
  evSummary_.dR_higgs_q1=0;
  evSummary_.dR_higgs_q2=0;
  evSummary_.phi_qq_higgs=0;
  evSummary_.avDeltaR=0;
  evSummary_.minDeltaM=0;
  evSummary_.MassOfbbj=0;
  evSummary_.DeltaEta=0;
  evSummary_.qq_m=0;
  evSummary_.ProductEta=0;
  evSummary_.qq_dR=0;
  evSummary_.qq_dphi=0;
  evSummary_.alpha_qq=0;
  evSummary_.N_jet=0;
  evSummary_.N_jet_eta_cut=0;
  evSummary_.jet1_pt=0;
  evSummary_.jet2_pt=0;
  evSummary_.jet3_pt=0;
  evSummary_.jet4_pt=0;
  evSummary_.jet5_pt=0;
  evSummary_.jet6_pt=0;
  evSummary_.jet1_eta=0;
  evSummary_.jet2_eta=0;
  evSummary_.jet3_eta=0;
  evSummary_.jet4_eta=0;
  evSummary_.jet5_eta=0;
  evSummary_.jet6_eta=0;
  evSummary_.H_T=0;
  evSummary_.H_z=0;
  evSummary_.H_Tvec=0;
  evSummary_.pt_rest=0;
  evSummary_.E_rest=0;
  evSummary_.N_bjet=0;
  evSummary_.bjet1_pt=0;
  evSummary_.bjet2_pt=0;
  evSummary_.bjet3_pt=0;
  evSummary_.bjet4_pt=0;
  evSummary_.bjet1_eta=0;
  evSummary_.bjet2_eta=0;
  evSummary_.bjet3_eta=0;
  evSummary_.bjet4_eta=0;
  evSummary_.bjet1_btag=0;
  evSummary_.bjet2_btag=0;
  evSummary_.bjet3_btag=0;
  evSummary_.bjet4_btag=0;
  evSummary_.N_qjet=0;
  evSummary_.qjet1_pt=0;
  evSummary_.qjet2_pt=0;
  evSummary_.qjet1_eta=0;
  evSummary_.qjet2_eta=0;
  evSummary_.N_untaggedjet=0;
  evSummary_.MET_pt=0;
  evSummary_.MET_phi=0;
  evSummary_.weight=0;

}


void MVAHandler::getEntry(float higgs_m,
			  float higgs_pt,
			  float higgs_eta,
			  float costheta0,
			  float dR_higgs_q1,
			  float dR_higgs_q2,
			  float phi_qq_higgs,
			  float avDeltaR,
			  float minDeltaM,
			  float MassOfbbj,
			  float DeltaEta,
			  float qq_m,
			  float ProductEta,
			  float qq_dR,
			  float qq_dphi,
			  float alpha_qq,
			  float N_jet,
			  float N_jet_eta_cut,
			  float jet1_pt,
			  float jet2_pt,
			  float jet3_pt,
			  float jet4_pt,
			  float jet5_pt,
			  float jet6_pt,
			  float jet1_eta,
			  float jet2_eta,
			  float jet3_eta,
			  float jet4_eta,
			  float jet5_eta,
			  float jet6_eta,
			  float H_T,
			  float H_z,
			  float H_Tvec,
			  float pt_rest,
			  float E_rest,
			  float N_bjet,
			  float bjet1_pt,
			  float bjet2_pt,
			  float bjet3_pt,
			  float bjet4_pt,
			  float bjet1_eta,
			  float bjet2_eta,
			  float bjet3_eta,
			  float bjet4_eta,
			  float bjet1_btag,
			  float bjet2_btag,
			  float bjet3_btag,
			  float bjet4_btag,
			  float N_qjet,
			  float qjet1_pt,
			  float qjet2_pt,
			  float qjet1_eta,
			  float qjet2_eta,
			  float N_untaggedjet,
			  float MET_pt,
			  float MET_phi,
			  float weight
			  ) 
 
{
  
  resetStruct();

  evSummary_.higgs_m=higgs_m;
  evSummary_.higgs_pt=higgs_pt;
  evSummary_.higgs_eta=higgs_eta;
  evSummary_.costheta0=costheta0;
  evSummary_.dR_higgs_q1=dR_higgs_q1;
  evSummary_.dR_higgs_q2=dR_higgs_q2;
  evSummary_.phi_qq_higgs=phi_qq_higgs;
  evSummary_.avDeltaR=avDeltaR;
  evSummary_.minDeltaM=minDeltaM;
  evSummary_.MassOfbbj=MassOfbbj;
  evSummary_.DeltaEta=DeltaEta;
  evSummary_.qq_m=qq_m;
  evSummary_.ProductEta=ProductEta;
  evSummary_.qq_dR=qq_dR;
  evSummary_.qq_dphi=qq_dphi;
  evSummary_.alpha_qq=alpha_qq;
  evSummary_.N_jet=N_jet;
  evSummary_.N_jet_eta_cut=N_jet_eta_cut;
  evSummary_.jet1_pt=jet1_pt;
  evSummary_.jet2_pt=jet2_pt;
  evSummary_.jet3_pt=jet3_pt;
  evSummary_.jet4_pt=jet4_pt;
  evSummary_.jet5_pt=jet5_pt;
  evSummary_.jet6_pt=jet6_pt;
  evSummary_.jet1_eta=jet1_eta;
  evSummary_.jet2_eta=jet2_eta;
  evSummary_.jet3_eta=jet3_eta;
  evSummary_.jet4_eta=jet4_eta;
  evSummary_.jet5_eta=jet5_eta;
  evSummary_.jet6_eta=jet6_eta;
  evSummary_.H_T=H_T;
  evSummary_.H_z=H_z;
  evSummary_.H_Tvec=H_Tvec;
  evSummary_.pt_rest=pt_rest;
  evSummary_.E_rest=E_rest;
  evSummary_.N_bjet=N_bjet;
  evSummary_.bjet1_pt=bjet1_pt;
  evSummary_.bjet2_pt=bjet2_pt;
  evSummary_.bjet3_pt=bjet3_pt;
  evSummary_.bjet4_pt=bjet4_pt;
  evSummary_.bjet1_eta=bjet1_eta;
  evSummary_.bjet2_eta=bjet2_eta;
  evSummary_.bjet3_eta=bjet3_eta;
  evSummary_.bjet4_eta=bjet4_eta;
  evSummary_.bjet1_btag=bjet1_btag;
  evSummary_.bjet2_btag=bjet2_btag;
  evSummary_.bjet3_btag=bjet3_btag;
  evSummary_.bjet4_btag=bjet4_btag;
  evSummary_.N_qjet=N_qjet;
  evSummary_.qjet1_pt=qjet1_pt;
  evSummary_.qjet2_pt=qjet2_pt;
  evSummary_.qjet1_eta=qjet1_eta;
  evSummary_.qjet2_eta=qjet2_eta;
  evSummary_.N_untaggedjet=N_untaggedjet;
  evSummary_.MET_pt=MET_pt;
  evSummary_.MET_phi=MET_phi;
  evSummary_.weight=weight;

  return ;
}

//
bool MVAHandler::initTree()
{
  
  t1 = new TTree("t1","trMVA");
  t1->SetDirectory(0);
  std::cout << "check init tree"<<std::endl;

  t1->Branch("higgs_m",        &evSummary_.higgs_m,           "higgs_m/F");
  t1->Branch("higgs_pt",       &evSummary_.higgs_pt,          "higgs_pt/F");
  t1->Branch("higgs_eta",      &evSummary_.higgs_eta,         "higgs_eta/F");
  t1->Branch("costheta0",      &evSummary_.costheta0,         "costheta0/F");
  t1->Branch("dR_higgs_q1",    &evSummary_.dR_higgs_q1,       "dR_higgs_q1/F");
  t1->Branch("dR_higgs_q2",    &evSummary_.dR_higgs_q2,       "dR_higgs_q2/F");
  t1->Branch("phi_qq_higgs",   &evSummary_.phi_qq_higgs,      "phi_qq_higgs/F");
  t1->Branch("avDeltaR",       &evSummary_.avDeltaR,          "avDeltaR/F");
  t1->Branch("minDeltaM",      &evSummary_.minDeltaM,         "minDeltaM/F");
  t1->Branch("MassOfbbj",      &evSummary_.MassOfbbj,         "MassOfbbj/F");
  t1->Branch("DeltaEta",       &evSummary_.DeltaEta,          "DeltaEta/F");
  t1->Branch("qq_m",           &evSummary_.qq_m,              "qq_m/F");
  t1->Branch("ProductEta",     &evSummary_.ProductEta,        "ProductEta/F");
  t1->Branch("qq_dR",          &evSummary_.qq_dR,             "qq_dR/F");
  t1->Branch("qq_dphi",        &evSummary_.qq_dphi,           "qq_dphi/F");
  t1->Branch("alpha_qq",       &evSummary_.alpha_qq,          "alpha_qq/F");
  t1->Branch("N_jet",          &evSummary_.N_jet,             "N_jet/F");
  t1->Branch("N_jet_eta_cut",  &evSummary_.N_jet_eta_cut,     "N_jet_eta_cut/F");
  t1->Branch("jet1_pt",        &evSummary_.jet1_pt,           "jet1_pt/F");
  t1->Branch("jet2_pt",        &evSummary_.jet2_pt,           "jet2_pt/F");
  t1->Branch("jet3_pt",        &evSummary_.jet3_pt,           "jet3_pt/F");
  t1->Branch("jet4_pt",        &evSummary_.jet4_pt,           "jet4_pt/F");
  t1->Branch("jet5_pt",        &evSummary_.jet5_pt,           "jet5_pt/F");
  t1->Branch("jet6_pt",        &evSummary_.jet6_pt,           "jet6_pt/F");
  t1->Branch("jet1_eta",       &evSummary_.jet1_eta,          "jet1_eta/F");
  t1->Branch("jet2_eta",       &evSummary_.jet2_eta,          "jet2_eta/F");
  t1->Branch("jet3_eta",       &evSummary_.jet3_eta,          "jet3_eta/F");
  t1->Branch("jet4_eta",       &evSummary_.jet4_eta,          "jet4_eta/F");
  t1->Branch("jet5_eta",       &evSummary_.jet5_eta,          "jet5_eta/F");
  t1->Branch("jet6_eta",       &evSummary_.jet6_eta,          "jet6_eta/F");
  t1->Branch("H_T",            &evSummary_.H_T,               "H_T/F");
  t1->Branch("H_z",            &evSummary_.H_z,               "H_z/F");
  t1->Branch("H_Tvec",         &evSummary_.H_Tvec,            "H_Tvec/F");
  t1->Branch("pt_rest",        &evSummary_.pt_rest,           "pt_rest/F");
  t1->Branch("E_rest",         &evSummary_.E_rest,            "E_rest/F");
  t1->Branch("N_bjet",         &evSummary_.N_bjet,            "N_bjet/F");
  t1->Branch("bjet1_pt",       &evSummary_.bjet1_pt,          "bjet1_pt/F");
  t1->Branch("bjet2_pt",       &evSummary_.bjet2_pt,          "bjet2_pt/F");
  t1->Branch("bjet3_pt",       &evSummary_.bjet3_pt,          "bjet3_pt/F");
  t1->Branch("bjet4_pt",       &evSummary_.bjet4_pt,          "bjet4_pt/F");
  t1->Branch("bjet1_eta",      &evSummary_.bjet1_eta,         "bjet1_eta/F");
  t1->Branch("bjet2_eta",      &evSummary_.bjet2_eta,         "bjet2_eta/F");
  t1->Branch("bjet3_eta",      &evSummary_.bjet3_eta,         "bjet3_eta/F");
  t1->Branch("bjet4_eta",      &evSummary_.bjet4_eta,         "bjet4_eta/F");
  t1->Branch("bjet1_btag",     &evSummary_.bjet1_btag,        "bjet1_btag/F");
  t1->Branch("bjet2_btag",     &evSummary_.bjet2_btag,        "bjet2_btag/F");
  t1->Branch("bjet3_btag",     &evSummary_.bjet3_btag,        "bjet3_btag/F");
  t1->Branch("bjet4_btag",     &evSummary_.bjet4_btag,        "bjet4_btag/F");
  t1->Branch("N_qjet",         &evSummary_.N_qjet,            "N_qjet/F");
  t1->Branch("qjet1_pt",       &evSummary_.qjet1_pt,          "qjet1_pt/F");
  t1->Branch("qjet2_pt",       &evSummary_.qjet2_pt,          "qjet2_pt/F");
  t1->Branch("qjet1_eta",      &evSummary_.qjet1_eta,         "qjet1_eta/F");
  t1->Branch("qjet2_eta",      &evSummary_.qjet2_eta,         "qjet2_eta/F");
  t1->Branch("N_untaggedjet",  &evSummary_.N_untaggedjet,     "N_untaggedjet/F");
  t1->Branch("MET_pt",         &evSummary_.MET_pt,            "MET_pt/F");
  t1->Branch("MET_phi",        &evSummary_.MET_phi,           "MET_phi/F");
  t1->Branch("weight",         &evSummary_.weight,            "weight/F");

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
