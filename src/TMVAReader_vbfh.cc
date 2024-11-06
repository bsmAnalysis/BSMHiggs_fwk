#include "UserCode/bsmhiggs_fwk/interface/TMVAReader_vbfh.h"

void TMVAReader::InitTMVAReader()
{
  myreader = new TMVA::Reader( "!Color:!Silent" );
  return ;
}


void TMVAReader::SetupMVAReader( std::string methodName, std::string modelPath )
{
  if(methodName.find("VBFhAnalysisClass") != std::string::npos)
    {
      if(modelPath.find("VBFHaa4bMVA") != std::string::npos)
	{
	  myreader->AddVariable("higgs_m",        &higgs_m);
	  myreader->AddVariable("higgs_pt",       &higgs_pt);
	  myreader->AddVariable("higgs_eta",      &higgs_eta);
	  myreader->AddVariable("costheta0",      &costheta0);
	  myreader->AddVariable("dR_higgs_q1",    &dR_higgs_q1);
	  myreader->AddVariable("dR_higgs_q2",    &dR_higgs_q2);
	  myreader->AddVariable("phi_qq_higgs",   &phi_qq_higgs);
	  myreader->AddVariable("avDeltaR",       &avDeltaR);
	  myreader->AddVariable("minDeltaM",      &minDeltaM);
	  myreader->AddVariable("MassOfbbj",      &MassOfbbj);
	  myreader->AddVariable("DeltaEta",       &DeltaEta);
	  myreader->AddVariable("qq_m",           &qq_m);              
	  myreader->AddVariable("ProductEta",     &ProductEta);
	  myreader->AddVariable("qq_dphi",        &qq_dphi);
	  myreader->AddVariable("qq_dR",          &qq_dR);
	  myreader->AddVariable("alpha_qq",       &alpha_qq);
	  myreader->AddVariable("jet1_eta",       &jet1_eta);
	  myreader->AddVariable("jet1_pt",        &jet1_pt);
	  myreader->AddVariable("jet2_pt",        &jet2_pt);
	  myreader->AddVariable("jet3_pt",        &jet3_pt);
	  myreader->AddVariable("jet4_pt",        &jet4_pt);
	  myreader->AddVariable("jet5_pt",        &jet5_pt);
	  myreader->AddVariable("N_jet_eta_cut",  &N_jet_eta_cut);
	  myreader->AddVariable("H_T",            &H_T);
	  myreader->AddVariable("H_z",            &H_z);
	  myreader->AddVariable("H_Tvec",         &H_Tvec);
	  myreader->AddVariable("bjet1_pt",       &bjet1_pt);
	  myreader->AddVariable("bjet2_pt",       &bjet2_pt);
	  myreader->AddVariable("bjet3_pt",       &bjet3_pt);
	  myreader->AddVariable("bjet1_eta",      &bjet1_eta);
	  myreader->AddVariable("bjet2_eta",      &bjet2_eta);
	  myreader->AddVariable("bjet3_eta",      &bjet3_eta);
	  myreader->AddVariable("bjet1_btag",     &bjet1_btag);
	  myreader->AddVariable("bjet2_btag",     &bjet2_btag);
	  myreader->AddVariable("bjet3_btag",     &bjet3_btag);
	  myreader->AddVariable("N_bjet",         &N_bjet);	  
	  myreader->AddVariable("qjet1_pt",       &qjet1_pt);
	  myreader->AddVariable("qjet2_pt",       &qjet2_pt);
	  myreader->AddVariable("qjet1_eta",      &qjet1_eta);
	  myreader->AddVariable("qjet2_eta",      &qjet2_eta);
	  myreader->AddVariable("pt_rest",        &pt_rest);
	  myreader->AddVariable("E_rest",         &E_rest);
	  myreader->AddVariable("N_untaggedjet",  &N_untaggedjet);	  
	  myreader->AddVariable("MET_pt",         &MET_pt);
	  myreader->AddVariable("MET_phi",        &MET_phi);
	}
    }

  myreader->BookMVA(methodName.c_str(), modelPath.c_str());
  return ;
}


float TMVAReader::GenReMVAReader(float thishiggs_m, float thishiggs_pt, float thishiggs_eta, float thiscostheta0, float thisdR_higgs_q1, float thisdR_higgs_q2, float thisphi_qq_higgs, float thisavDeltaR, float thisminDeltaM, float thisMassOfbbj, float thisDeltaEta, float thisqq_m, float thisProductEta, float thisqq_dphi, float thisqq_dR, float thisalpha_qq, float thisjet1_eta, float thisjet1_pt, float thisjet2_pt, float thisjet3_pt, float thisjet4_pt, float thisjet5_pt, float thisN_jet_eta_cut, float thisH_T, float thisH_z, float thisH_Tvec, float thisbjet1_pt, float thisbjet2_pt, float thisbjet3_pt, float thisbjet1_eta, float thisbjet2_eta, float thisbjet3_eta, float thisbjet1_btag, float thisbjet2_btag, float thisbjet3_btag, float thisN_bjet, float thisqjet1_pt, float thisqjet2_pt, float thisqjet1_eta, float thisqjet2_eta, float thispt_rest, float thisE_rest, float thisN_untaggedjet, float thisMET_pt, float thisMET_phi, std::string methodName)

{
  higgs_m=thishiggs_m; higgs_pt=thishiggs_pt; higgs_eta=thishiggs_eta; costheta0=thiscostheta0; dR_higgs_q1=thisdR_higgs_q1; dR_higgs_q2=thisdR_higgs_q2; phi_qq_higgs=thisphi_qq_higgs; avDeltaR=thisavDeltaR; minDeltaM=thisminDeltaM; MassOfbbj=thisMassOfbbj; DeltaEta=thisDeltaEta; qq_m=thisqq_m; ProductEta=thisProductEta; qq_dphi=thisqq_dphi; qq_dR=thisqq_dR; alpha_qq=thisalpha_qq; jet1_eta=thisjet1_eta; jet1_pt=thisjet1_pt; jet2_pt=thisjet2_pt; jet3_pt=thisjet3_pt; jet4_pt=thisjet4_pt; jet5_pt=thisjet5_pt; N_jet_eta_cut=thisN_jet_eta_cut; H_T=thisH_T; H_z=thisH_z; H_Tvec=thisH_Tvec; bjet1_pt=thisbjet1_pt; bjet2_pt=thisbjet2_pt; bjet3_pt=thisbjet3_pt; bjet1_eta=thisbjet1_eta; bjet2_eta=thisbjet2_eta; bjet3_eta=thisbjet3_eta; bjet1_btag=thisbjet1_btag; bjet2_btag=thisbjet2_btag; bjet3_btag=thisbjet3_btag; N_bjet=thisN_bjet; qjet1_pt=thisqjet1_pt; qjet2_pt=thisqjet2_pt; qjet1_eta=thisqjet1_eta; qjet2_eta=thisqjet2_eta; pt_rest=thispt_rest; E_rest=thisE_rest; N_untaggedjet=thisN_untaggedjet; MET_pt=thisMET_pt; MET_phi=thisMET_phi;
  
  float mvaOut = myreader->EvaluateMVA( methodName.c_str() );
  return mvaOut;

}


void TMVAReader::CloseMVAReader()
{
  delete myreader;
  return ;
}
