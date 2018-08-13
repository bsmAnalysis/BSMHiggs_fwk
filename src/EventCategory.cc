#include "UserCode/bsmhiggs_fwk/interface/EventCategory.h"


using namespace std;

EventCategory::EventCategory(int mode_)
{
    mode = mode_;
    if(mode==1) {
        NStates = 6;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "CR_nonTT_3b"; // Nb=0
	EvtCategoryLabel[1] = "CR_nonTT_4b"; // Nb=0
        EvtCategoryLabel[2] = "CR_3b"; // Nb=1,2
	EvtCategoryLabel[3] = "CR_4b"; // Nb=1,2
        EvtCategoryLabel[4] = "SR_3b"; // Nb=3,4
	EvtCategoryLabel[5] = "SR_4b"; // Nb=3,4
    } else if(mode==2) {
        NStates = 2;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "lesq1jets";
        EvtCategoryLabel[1] = "geq2jets";
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
    if(!JetsP4) {
        printf("NoJetCollection given to EventCategory::Get\n");
        exit(0);
    }

    //Nb jets
    // NbJets = 0;
    std::vector<LorentzVector> bjets = *bJetsP4;
    NbJets=bjets.size();
    
    // for(size_t ijet=0; ijet<bjets.size(); ijet++) {
    //     if(bjets[ijet].pt()<=20)continue;
    //     NbJets++;
    // }

    //Nj jets
    //NJets = 0;
    std::vector<LorentzVector> jets = *JetsP4;
    NJets=jets.size();
    
    // for(size_t ijet=0; ijet<jets.size(); ijet++) {
    //     if(jets[ijet].pt()<=20)continue;
    //     NJets++;
    // }
    
    switch(mode) {
    case 1: {
      if(NbJets==4) {
	if (NJets==4) return 5;
      }
      if(NbJets==3) {
	if (NJets==3 || NJets==4) return 4;
      }
      if(NbJets==2) { // Top CR
	if (NJets==4) { return 3;}
	else if (NJets==3) { return 2; }
      }
      if(NbJets==0) { // W CR
	if (NJets==4) { return 1;}
	else if (NJets==3) {return 0;}
      }
      return -1;
    }
    break;
    case 2: {
        if(NJets>=2) return 1;
        return 0;
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
