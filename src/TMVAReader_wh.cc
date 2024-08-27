#include "UserCode/bsmhiggs_fwk/interface/TMVAReader_wh.h"

void TMVAReader::InitTMVAReader()
{
  myreader = new TMVA::Reader( "!Color:!Silent" );
  return ;
}


void TMVAReader::SetupMVAReader( std::string methodName, std::string modelPath )
{
  if(methodName.find("WhAnalysisClass") != std::string::npos)
    {
      if(modelPath.find("WhAnalysis") != std::string::npos)
	{
	  myreader->AddVariable("Njets"     , &Njets);
	  myreader->AddVariable("HTjets"    , &HTjets);
	  myreader->AddVariable("H_pt"      , &H_pt);
	  myreader->AddVariable("W_pt"      , &W_pt);
	  myreader->AddVariable("lepton_pt" , &lepton_pt);  
	  myreader->AddVariable("bjet1_pt"  , &bjet1_pt);
	  myreader->AddVariable("MET"       , &MET);              
	  myreader->AddVariable("H_M"       , &H_M);
	  myreader->AddVariable("MTW"       , &MTW);
	  myreader->AddVariable("dphiWH"    , &dphiWH);
	  myreader->AddVariable("dphiMetJetMin", &dphiMetJetMin);
	  myreader->AddVariable("btag1"     , &btag1);
	  myreader->AddVariable("btag2"     , &btag2);
	  myreader->AddVariable("btag3"     , &btag3);
	}
    }

  myreader->BookMVA(methodName.c_str(), modelPath.c_str());
  return ;
}


float TMVAReader::GenReMVAReader(float thisNjets, float thisHTjets, float thisH_pt, float thisW_pt, float thislepton_pt, float thisbjet1_pt, float thisMET, float thisH_M, float thisMTW, float thisdphiWH, float thisdphiMetJetMin, float thisbtag1, float thisbtag2, float thisbtag3, std::string methodName)

{
  Njets=thisNjets; HTjets=thisHTjets; H_pt=thisH_pt; W_pt=thisW_pt; lepton_pt=thislepton_pt; bjet1_pt=thisbjet1_pt; MET=thisMET; H_M=thisH_M; MTW=thisMTW; dphiWH=thisdphiWH; dphiMetJetMin=thisdphiMetJetMin; btag1=thisbtag1; btag2=thisbtag2; btag3=thisbtag3;

  float mvaOut = myreader->EvaluateMVA( methodName.c_str() );
  return mvaOut;

}


void TMVAReader::CloseMVAReader()
{
  delete myreader;
  return ;
}
