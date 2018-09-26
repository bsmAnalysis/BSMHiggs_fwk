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
  if (mode==0) return origJet;

    if(genJetPt<=0) return origJet;

    //smearing factors are described in https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    double eta=fabs(origJet.eta());

    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
    // Spring16_25nsV10 (80X, 2016, BCD+GH PromtReco) DATA/MC SFs 
    double ptSF(1.0), ptSF_up(1.0), ptSF_down(1.0);
    if(eta<0.522)			   {
      //      ptSF=1.1595;
      ptSF_up=ptSF+0.0645;
      ptSF_down=ptSF-0.0645;
    } else if(eta>=0.522 && eta<0.783) {
      //    ptSF=1.1948;
      ptSF_up=ptSF+0.0652;
      ptSF_down=ptSF-0.0652;
    } else if(eta>=0.783 && eta<1.131) {
      //      ptSF=1.1464;
      ptSF_up=ptSF+0.0632;
      ptSF_down=ptSF-0.0632;
    } else if(eta>=1.131 && eta<1.305) {
      //      ptSF=1.1609;
      ptSF_up=ptSF+0.1025;
      ptSF_down=ptSF-0.1025;
    } else if(eta>=1.305 && eta<1.740) {
      //      ptSF=1.1278;
      ptSF_up=ptSF+0.0986;
      ptSF_down=ptSF-0.0986;
    } else if(eta>=1.740 && eta<1.930) {
      //      ptSF=1.1000;
      ptSF_up=ptSF+0.1079;
      ptSF_down=ptSF-0.1079;
    } else if(eta>=1.930 && eta<2.043) {
      //      ptSF=1.1426;
      ptSF_up=ptSF+0.1214;
      ptSF_down=ptSF-0.1214;
    } else if(eta>=2.043 && eta<2.322) {
      //      ptSF=1.1512;
      ptSF_up=ptSF+0.1140;
      ptSF_down=ptSF-0.1140;
    } else if(eta>=2.322 && eta<2.5) {
      //      ptSF=1.2963;
      ptSF_up=ptSF+0.2371;
      ptSF_down=ptSF-0.2371;
    } else if(eta>=2.5 && eta<2.853) {
      //      ptSF=1.3418;
      ptSF_up=ptSF+0.2091;
      ptSF_down=ptSF-0.2091;
    } else if(eta>=2.853 && eta<2.964) {
      //      ptSF=1.7788;
      ptSF_up=ptSF+0.2008;
      ptSF_down=ptSF-0.2008;
    } else if(eta>=2.964 && eta<3.139) {
      //      ptSF=1.1869;
      ptSF_up=ptSF+0.1243;
      ptSF_down=ptSF-0.1243;
    } else if(eta>=3.139 && eta<5.191) {
      //      ptSF=1.1922;
      ptSF_up=ptSF+0.1488;
      ptSF_down=ptSF-0.1488;
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
		      std::vector<JetCorrectionUncertainty*> &jecUnc)
		      //                      JetCorrectionUncertainty *jecUnc)
{
    jetsVar.clear();
    metsVar.clear();

    int vars[]= {JER, JER_UP,JER_DOWN, UMET_UP,UMET_DOWN, LES_UP,LES_DOWN};
    //  int vars[]= {JER, JER_UP,JER_DOWN, JES_UP,JES_DOWN, UMET_UP,UMET_DOWN, LES_UP,LES_DOWN};
    for(size_t ivar=0; ivar<sizeof(vars)/sizeof(int); ivar++) {
        PhysicsObjectJetCollection newJets;
        LorentzVector newMET(met), jetDiff(0,0,0,0), lepDiff(0,0,0,0), unclustDiff(0,0,0,0), clusteredFlux(0,0,0,0);
        for(size_t ijet=0; ijet<jets.size(); ijet++) {
	  if(ivar==JER || ivar==JER_UP || ivar==JER_DOWN) {  
                PhysicsObject_Jet iSmearJet=METUtils::smearedJet(jets[ijet],jets[ijet].genPt,ivar);
                jetDiff += (iSmearJet-jets[ijet]);
                newJets.push_back( iSmearJet );
	  // } else if(ivar==JES_UP || ivar==JES_DOWN) {

	  //   for (int isrc = 0; isrc < nsrc; isrc++) {
	  //   // sub-loop to decorrelate JES

	  //     JetCorrectionUncertainty *unc = jetUnc[isrc];

	  //     double varSign=(ivar==JES_UP ? 1.0 : -1.0 );
	  //     double jetScale(1.0);
	  //     try {
	  // 	unc->setJetEta(jets[ijet].eta());
	  // 	unc->setJetPt(jets[ijet].pt());
	  // 	jetScale = 1.0 + varSign*fabs(unc->getUncertainty(true));
	  //     } catch(std::exception &e) {
	  // 	cout << "[METUtils::computeVariation]" << e.what() << ijet << " " << jets[ijet].pt() << endl;
	  //     }
	  //     PhysicsObject_Jet iScaleJet(jets[ijet]);
	  //     iScaleJet *= jetScale;
	  //     jetDiff += (iScaleJet-jets[ijet]);
	  //     newJets.push_back(iScaleJet);

	  //   }// end subloop on JES
	  } else if(ivar==UMET_UP || ivar==UMET_DOWN)  clusteredFlux += jets[ijet];
	  
	  
	  if(ivar==UMET_UP || ivar==UMET_DOWN || ivar==LES_UP || ivar==LES_DOWN) {    
	    //  if(ivar==UMET_UP || ivar==UMET_DOWN)  clusteredFlux += jets[ijet];  
	    // Reset jets  
	    PhysicsObject_Jet iScaleJet(jets[ijet]); 
	    newJets.push_back(iScaleJet);      
	  }

        } // end loop on ijet

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

    } // end ivariation

    // Now add JES as recommended from here: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources
    // Instantiate JES uncertainty sources
    const int nsrc = 27;
    
    for (int isrc = 0; isrc < nsrc; isrc++) {

      for(size_t ivar=0; ivar<2; ivar++) { //ivar<sizeof(jesvars)/sizeof(int); ivar++) {
	PhysicsObjectJetCollection newJets;
	LorentzVector newMET(met), jetDiff(0,0,0,0);
	
	for(size_t ijet=0; ijet<jets.size(); ijet++) {
	  
	  // sub-loop to decorrelate JES
	  JetCorrectionUncertainty *unc = jecUnc.at(isrc);
	  
	  //	  bool varSign=(ivar==0 ? true : false );    
	  double varSign=(ivar==0 ? 1.0 : -1.0 );
	  double jetScale(1.0);
	  try {
	    unc->setJetEta(jets[ijet].eta());
	    unc->setJetPt(jets[ijet].pt());
	    //jetScale = 1.0 + fabs(unc->getUncertainty(varSign)); 
	    jetScale = 1.0 + varSign*fabs(unc->getUncertainty(true));
	  } catch(std::exception &e) {
	    cout << "[METUtils::computeVariation]" << e.what() << ijet << " " << jets[ijet].pt() << endl;
	  }
	  
	  PhysicsObject_Jet iScaleJet(jets[ijet]);
	  iScaleJet *= jetScale;
	  jetDiff += (iScaleJet-jets[ijet]);
	  newJets.push_back(iScaleJet);
	  
	  // up OR down end
	} // ijet loop 

	//add new met
	newMET -= jetDiff;
	metsVar.push_back(newMET);   
	//add new jets (if some change has occured)
	jetsVar.push_back(newJets);
      
      } // loop on up and down
    } //JES nsrc sources end
    
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
