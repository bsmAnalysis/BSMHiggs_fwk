#ifndef ____EventCategory__
#define ____EventCategory__

#include <vector>
#include <iostream>
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"
#include "UserCode/bsmhiggs_fwk/interface/DataEvtSummaryHandler.h"

class EventCategory {
public:
    int mode;
    int NStates;
    TString* EvtCategoryLabel;


    /// Constructor
    EventCategory(int mode=0);

    /// Destructor
    virtual ~EventCategory();

    int GetLabelSize() {
        return NStates;
    }

    int Get(const PhysicsEvent_t& phys, std::vector<LorentzVector>* bJetsP4=NULL, std::vector<LorentzVector>* JetsP4=NULL);
    TString GetLabel(int CategoryType);
    TString GetLabel(const PhysicsEvent_t& phys);


private:
    int NJets;
    int NbJets;

};


#endif /* defined(____EventCategory__) */
