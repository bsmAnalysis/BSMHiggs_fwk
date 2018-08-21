#include "UserCode/bsmhiggs_fwk/interface/EventCategory.h"


using namespace std;

EventCategory::EventCategory(int mode_)
{
    mode = mode_;
    if(mode==1) {
        NStates = 9;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "CR_nonTT_3b"; // Nb=0
	EvtCategoryLabel[1] = "CR_nonTT_4b"; // Nb=0
        EvtCategoryLabel[2] = "CR_3b"; // Nb=2 -> W CR
	EvtCategoryLabel[3] = "CR_4b"; // Nb=2 -> tt+light CR
        EvtCategoryLabel[4] = "SR_3b"; // Nb=3,4
	EvtCategoryLabel[5] = "SR_4b"; // Nb=3,4
	EvtCategoryLabel[6] = "CR_3b5j"; // Nb=3,4 -> tt+light CR
	EvtCategoryLabel[7] = "CR_4b5j"; // Nb=3,4 -> tt+HF CR
	EvtCategoryLabel[8] = "CR_5b"; // Nb=3,4
    } else if(mode==2) {
        NStates = 18;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "0b_3j";
        EvtCategoryLabel[1] = "0b_4j";
	EvtCategoryLabel[2] = "0b_geq5j";
	EvtCategoryLabel[3] = "1b_3j";
        EvtCategoryLabel[4] = "1b_4j";
	EvtCategoryLabel[5] = "1b_geq5j";
	EvtCategoryLabel[6] = "2b_3j";
        EvtCategoryLabel[7] = "2b_4j";
	EvtCategoryLabel[8] = "2b_geq5j";
	EvtCategoryLabel[9] = "3b_3j";
        EvtCategoryLabel[10] = "3b_4j";
	EvtCategoryLabel[11] = "3b_geq5j";
	EvtCategoryLabel[12] = "4b_3j";
        EvtCategoryLabel[13] = "4b_4j";
	EvtCategoryLabel[14] = "4b_geq5j";
	EvtCategoryLabel[15] = "5b_3j";
        EvtCategoryLabel[16] = "5b_4j";
	EvtCategoryLabel[17] = "5b_geq5j";
    } else {
        mode = 0;
        NStates = 1;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "";
    }
}

EventCategory::~EventCategory() {}

//
int EventCategory::Get(const PhysicsEvent_t& phys, std::vector<LorentzVector>* bJetsP4, std::vector<LorentzVector>* JetsP4)
{
    if(!JetsP4 || !bJetsP4) {
        printf("NoJetCollection given to EventCategory::Get\n");
        exit(0);
    }
    
    //Nj jets
    //NJets = 0;
    std::vector<LorentzVector> jets = *JetsP4;
    NJets=jets.size();
    
    //Nb jets
    // NbJets = 0;
    std::vector<LorentzVector> bjets = *bJetsP4;
    NbJets=bjets.size();

    switch(mode) {
    case 1: {
      if (NbJets==5) return 8;
      if(NbJets==4) {
	if (NJets==3 || NJets==4) return 5;
	if (NJets>=5) return 7;
      }
      if(NbJets==3) {
	if (NJets==3 || NJets==4) return 4;
	if (NJets>=5) return 6;
      }
      if(NbJets==2) { // Top CR
	if (NJets>=4) { return 3;}
	else if (NJets==3) { return 2; }
      }
      if(NbJets==0) { // W CR
	if (NJets>=4) { return 1;}
	else if (NJets==3) {return 0;}
      }
      return -1;
    }
    break;
    case 2: {
      //  if(NJets>=2) return 1;
      if (NbJets==0) {
	if (NJets==3) return 0;
	if (NJets==4) return 1;
	if (NJets>=5) return 2;
      }
      if (NbJets==1) {
	if (NJets==3) return 3;
	if (NJets==4) return 4;
	if (NJets>=5) return 5;
      }
      if (NbJets==2) {
	if (NJets==3) return 6;
	if (NJets==4) return 7;
	if (NJets>=5) return 8;
      }
      if (NbJets==3) {
	if (NJets==3) return 9;
	if (NJets==4) return 10;
	if (NJets>=5) return 11;
      }
      if (NbJets==4) {
	if (NJets==3) return 12;
	if (NJets==4) return 13;
	if (NJets>=5) return 14;
      }
      if (NbJets==5) {
	if (NJets==3) return 15;
	if (NJets==4) return 16;
	if (NJets>=5) return 17;
      }
      return -1;
    }
    break;
    case 0:
    default: {
        return 0;
    }
    break;
    }
}

//
TString EventCategory::GetLabel(int CategoryType)
{
    if(mode==0) {
        if(CategoryType<=2) return EvtCategoryLabel[0];
        return EvtCategoryLabel[1];
    } else {
        return EvtCategoryLabel[CategoryType];
    }
}

//
TString EventCategory::GetLabel(const PhysicsEvent_t& phys)
{
    return GetLabel(Get(phys));
}
