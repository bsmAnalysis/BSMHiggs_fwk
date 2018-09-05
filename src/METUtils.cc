//
//  METUtils.cpp
//

#include "UserCode/bsmhiggs_fwk/interface/METUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;


namespace METUtils {

double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass)
{
    if(assumeSameMass) {
        LorentzVector sum=visible+invisible;
        double tMass = TMath::Power(TMath::Sqrt(TMath::Power(visible.pt(),2)+pow(visible.mass(),2))+TMath::Sqrt(TMath::Power(invisible.pt(),2)+pow(visible.mass(),2)),2);
        tMass-=TMath::Power(sum.pt(),2);
        return TMath::Sqrt(tMass);
    } else {
        double dphi=fabs(deltaPhi(invisible.phi(),visible.phi()));
        return TMath::Sqrt(2*invisible.pt()*visible.pt()*(1-TMath::Cos(dphi)));
    }
    return -1;
}

double response(LorentzVector &Z, LorentzVector &MET)
{
        TVector2 z(Z.px(),Z.py());
        TVector2 met(MET.px(),MET.py());
        TVector2 sum = z+met;
        sum *= -1;
        double cos_ = (sum*z)/(sum.Mod()*z.Mod());
        return cos_*sum.Mod()/z.Mod();
}


//Jet energy resoltuion, 13TeV scale factors, updated on 30/08/2018
PhysicsObject_Jet smearedJet(const PhysicsObject_Jet &origJet, double genJetPt, int mode)
{
    if(genJetPt<=0) return origJet;

    //smearing factors are described in https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    double eta=fabs(origJet.eta());

    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
    // Spring16_25nsV10 (80X, 2016, BCD+GH PromtReco) DATA/MC SFs 
    double ptSF(1.0), ptSF_up(1.0), ptSF_down(1.0);
    if(eta<0.5)			   {
        ptSF=1.109;
        ptSF_up=ptSF+0.008;
        ptSF_down=ptSF-0.008;
    } else if(eta>=0.5 && eta<0.8) {
        ptSF=1.138;
        ptSF_up=ptSF+0.013;
        ptSF_down=ptSF-0.013;
    } else if(eta>=0.8 && eta<1.1) {
        ptSF=1.114;
        ptSF_up=ptSF+0.013;
        ptSF_down=ptSF-0.013;
    } else if(eta>=1.1 && eta<1.3) {
        ptSF=1.123;
        ptSF_up=ptSF+0.024;
        ptSF_down=ptSF-0.024;
    } else if(eta>=1.3 && eta<1.7) {
        ptSF=1.084;
        ptSF_up=ptSF+0.011;
        ptSF_down=ptSF-0.011;
    } else if(eta>=1.7 && eta<1.9) {
        ptSF=1.082;
        ptSF_up=ptSF+0.035;
        ptSF_down=ptSF-0.035;
    } else if(eta>=1.9 && eta<2.1) {
        ptSF=1.140;
        ptSF_up=ptSF+0.047;
        ptSF_down=ptSF-0.047;
    } else if(eta>=2.1 && eta<2.3) {
        ptSF=1.067;
        ptSF_up=ptSF+0.053;
        ptSF_down=ptSF-0.053;
    } else if(eta>=2.3 && eta<2.5) {
        ptSF=1.177;
        ptSF_up=ptSF+0.041;
        ptSF_down=ptSF-0.041;
    } else if(eta>=2.5 && eta<2.8) {
        ptSF=1.364;
        ptSF_up=ptSF+0.039;
        ptSF_down=ptSF-0.039;
    } else if(eta>=2.8 && eta<3.0) {
        ptSF=1.857;
        ptSF_up=ptSF+0.071;
        ptSF_down=ptSF-0.071;
    } else if(eta>=3.0 && eta<3.2) {
        ptSF=1.328;
        ptSF_up=ptSF+0.022;
        ptSF_down=ptSF-0.022;
    } else if(eta>=3.2 && eta<5.0) {
        ptSF=1.16;
        ptSF_up=ptSF+0.029;
        ptSF_down=ptSF-0.029;
    }

    if(mode==1) ptSF = ptSF_up;
    if(mode==2) ptSF = ptSF_down;


    ptSF=max(0.,(genJetPt+ptSF*(origJet.pt()-genJetPt)))/origJet.pt();                      //deterministic version
    //ptSF=max(0.,(genJetPt+gRandom->Gaus(ptSF,ptSF_err)*(origJet.pt()-genJetPt)))/origJet.pt();  //deterministic version
    if( ptSF<=0 /*|| isnan(ptSF)*/ ) return origJet;

    double px(origJet.px()*ptSF), py(origJet.py()*ptSF), pz(origJet.pz()*ptSF), mass(origJet.mass()*ptSF);
    double en = sqrt(mass*mass+px*px+py*py+pz*pz);

    PhysicsObject_Jet toReturn = origJet;
    toReturn.SetCoordinates(px, py, pz, en);
    //cout << "eta: " << eta << "\t" << toReturn.eta() << endl;
    return toReturn;
}


//
void computeVariation(PhysicsObjectJetCollection& jets,
                      PhysicsObjectLeptonCollection& leptons,
                      LorentzVector& met,
                      std::vector<PhysicsObjectJetCollection>& jetsVar,
                      LorentzVectorCollection& metsVar,
                      JetCorrectionUncertainty *jecUnc)
{
    jetsVar.clear();
    metsVar.clear();

    int vars[]= {JER, JER_UP,JER_DOWN, JES_UP,JES_DOWN, UMET_UP,UMET_DOWN, LES_UP,LES_DOWN};
    for(size_t ivar=0; ivar<sizeof(vars)/sizeof(int); ivar++) {
        PhysicsObjectJetCollection newJets;
        LorentzVector newMET(met), jetDiff(0,0,0,0), lepDiff(0,0,0,0), unclustDiff(0,0,0,0), clusteredFlux(0,0,0,0);
        for(size_t ijet=0; ijet<jets.size(); ijet++) {
            if(ivar==JER || ivar==JER_UP || ivar==JER_DOWN) {
                PhysicsObject_Jet iSmearJet=METUtils::smearedJet(jets[ijet],jets[ijet].genPt,ivar);
                jetDiff += (iSmearJet-jets[ijet]);
                newJets.push_back( iSmearJet );
            } else if(ivar==JES_UP || ivar==JES_DOWN) {
                double varSign=(ivar==JES_UP ? 1.0 : -1.0 );
                double jetScale(1.0);
                try {
                    jecUnc->setJetEta(jets[ijet].eta());
                    jecUnc->setJetPt(jets[ijet].pt());
                    jetScale = 1.0 + varSign*fabs(jecUnc->getUncertainty(true));
                } catch(std::exception &e) {
                    cout << "[METUtils::computeVariation]" << e.what() << ijet << " " << jets[ijet].pt() << endl;
                }
                PhysicsObject_Jet iScaleJet(jets[ijet]);
                iScaleJet *= jetScale;
                jetDiff += (iScaleJet-jets[ijet]);
                newJets.push_back(iScaleJet);
            } else if(ivar==UMET_UP || ivar==UMET_DOWN)  clusteredFlux += jets[ijet];
        }

        if(ivar==UMET_UP || ivar==UMET_DOWN || ivar==LES_UP || ivar==LES_DOWN) {
            for(size_t ilep=0; ilep<leptons.size(); ilep++) {
                if(ivar==UMET_UP || ivar==UMET_DOWN)  clusteredFlux +=leptons[ilep];
                else if(ivar==LES_UP  || ivar==LES_DOWN) {
                    LorentzVector iScaleLepton=leptons[ilep];
                    double varSign=(ivar==LES_UP ? 1.0 : -1.0);
                    if(fabs(leptons[ilep].id)==13)          iScaleLepton *= (1.0+varSign*0.002);
                    else if(fabs(leptons[ilep].id)==11) {
                        if(fabs(leptons[ilep].eta())<1.442) iScaleLepton *= (1.0+varSign*0.006);
                        else                                iScaleLepton *= (1.0+varSign*0.015);
                    }
                    lepDiff += (iScaleLepton-leptons[ilep]);
                }
            }
        }

        //vary unclustered component
        if(ivar==UMET_UP || ivar==UMET_DOWN) {
            unclustDiff=(met+clusteredFlux);
            double varSign=(ivar==UMET_UP ? 1.0 : -1.0);
            unclustDiff *= (varSign*0.10); //10% variation of residule recoil
        }

        //add new met
        newMET -= jetDiff;
        newMET -= lepDiff;
        newMET -= unclustDiff;
        metsVar.push_back(newMET);

        //add new jets (if some change has occured)
        jetsVar.push_back(newJets);


    }
}


//
LorentzVector applyMETXYCorr(LorentzVector met, bool isMC, int nvtx)
{
    double corX = 0.;
    double corY = 0.;


if(isMC){
  corX = -0.5091044 + (-0.0600439*nvtx);
  corY = 0.2873434 + (0.0183584*nvtx);
}
else{
  corX = -1.7325953 + (-0.1804222*nvtx);
  corY = 0.9452573 + (0.1277629*nvtx);
}

    double px = met.px()-corX;
    double py = met.py()-corY;
    return LorentzVector(px,py,0,sqrt(px*px+py*py));
}



}
