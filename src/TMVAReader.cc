#include "UserCode/bsmhiggs_fwk/interface/TMVAReader.h"

void TMVAReader::InitTMVAReader()
{
  myreader = new TMVA::Reader( "!Color:!Silent" );
  return ;
}


void TMVAReader::SetupMVAReader( std::string methodName, std::string modelPath )
{
  if ( methodName.find("TribMVA") != std::string::npos )
  {
    myreader->AddVariable( "WpT"    , &WpT );
    myreader->AddVariable( "Hmass"  , &Hmass );
    myreader->AddVariable( "HpT"    , &HpT );
    myreader->AddVariable( "bbdRAve", &bbdRAve );
    myreader->AddVariable( "HHt"    , &HHt );
    myreader->AddVariable( "WHdR"   , &WHdR );
  }
  else if ( methodName.find("QuabMVA") != std::string::npos )
  {
    myreader->AddVariable( "WpT"    , &WpT );
    myreader->AddVariable( "Hmass"  , &Hmass );
    myreader->AddVariable( "HpT"    , &HpT );
    myreader->AddVariable( "bbdRAve", &bbdRAve );
    myreader->AddVariable( "bbdMMin", &bbdMMin );
    myreader->AddVariable( "HHt"    , &HHt );
    myreader->AddVariable( "WHdR"   , &WHdR );
  }
  else
  {
    return ;
  }

  myreader->BookMVA( methodName.c_str(), modelPath.c_str() );
  return ;
}

float TMVAReader::GenReMVAReader(
                                 float thisWpT,
                                 float thisHmass, float thisHpT, float thisbbdRAve, float thisbbdMMin, float thisHHt,
                                 float thisWHdR,
                                 std::string methodName
                                )
{
  WpT = thisWpT; 
  Hmass= thisHmass; HpT = thisHpT; bbdRAve = thisbbdRAve; bbdMMin = thisbbdMMin; HHt = thisHHt;
  WHdR = thisWHdR;

  //float mvaOut = myreader->GetProba( methodName.c_str(), 0.2);
  float mvaOut = myreader->EvaluateMVA( methodName.c_str() );
  return mvaOut;
}
  
void TMVAReader::CloseMVAReader()
{
  delete myreader;
  return ;
}
