//
//  BSMPhysicsEvent.cc
//
//#define YEAR_2017
#include "UserCode/bsmhiggs_fwk/interface/BSMPhysicsEvent.h"

using namespace std;


//
PhysicsEvent_t getPhysicsEventFrom(DataEvtSummary_t &ev)
{
    PhysicsEvent_t phys;

    phys.run=ev.run;
    phys.event=ev.event;
    phys.lumi=ev.lumi;

    phys.nvtx = ev.nvtx;

    // Leptons
    size_t nlep(0);
    for(Int_t i=0; i<ev.mn; i++) {
        LorentzVector P4(ev.mn_px[i],ev.mn_py[i],ev.mn_pz[i],ev.mn_en[i]);
        if(P4.pt()>0) {
            phys.leptons.push_back(PhysicsObject_Lepton(P4, ev.mn_id[i]));
            phys.leptons[nlep].setLeptonIDInfo(ev.mn_IsLoose[i], ev.mn_IsTight[i], ev.mn_IsMedium[i], ev.mn_IsSoft[i],ev.mn_IsHighPt[i],
					       ev.mn_passId[i], ev.mn_passIdLoose[i], ev.mn_passSoftMuon[i],
                                               ev.en_passVeto[i], ev.en_passLoose[i], ev.en_passMedium[i], ev.en_passTight[i],
					       ev.en_passId[i], ev.en_passIdLoose[i],
                                               ev.ta_dm[i]
					       );
	    phys.leptons[nlep].setLeptonVar(ev.en_cor_en[i],ev.en_EtaSC[i],ev.en_R9[i],ev.en_gainSeed[i],
					    ev.mn_validMuonHits[i],ev.mn_trkLayersWithMeasurement[i],ev.mn_pixelLayersWithMeasurement[i]
					    );
            phys.leptons[nlep].setLeptonIsoInfo(ev.mn_pileupIsoR03[i],ev.mn_chargedIsoR03[i],ev.mn_photonIsoR03[i],ev.mn_neutralHadIsoR03[i],ev.mn_passIso[i], ev.mn_relIso[i], ev.mn_trkrelIso[i],
                                                ev.en_pileupIso[i],ev.en_chargedIso[i],ev.en_photonIso[i],ev.en_neutralHadIso[i],ev.en_relIsoWithEA[i],ev.en_passIso[i],
                                                ev.ta_IsLooseIso[i], ev.ta_IsMediumIso[i], ev.ta_IsTightIso[i], ev.en_relIso[i]
                                               );

            nlep++;
        }
    }

    for(Int_t i=0; i<ev.en; i++) {
        LorentzVector P4(ev.en_px[i],ev.en_py[i],ev.en_pz[i],ev.en_en[i]);
        if(P4.pt()>0) {
            phys.leptons.push_back(PhysicsObject_Lepton(P4, ev.en_id[i]));
            phys.leptons[nlep].setLeptonIDInfo(ev.mn_IsLoose[i], ev.mn_IsTight[i], ev.mn_IsMedium[i], ev.mn_IsSoft[i],ev.mn_IsHighPt[i],
					       ev.mn_passId[i], ev.mn_passIdLoose[i], ev.mn_passSoftMuon[i],
                                               ev.en_passVeto[i], ev.en_passLoose[i], ev.en_passMedium[i], ev.en_passTight[i],
					       ev.en_passId[i], ev.en_passIdLoose[i],
                                               ev.ta_dm[i]
                                              );
	    phys.leptons[nlep].setLeptonVar(ev.en_cor_en[i],ev.en_EtaSC[i],ev.en_R9[i],ev.en_gainSeed[i],   
					    ev.mn_validMuonHits[i],ev.mn_trkLayersWithMeasurement[i],ev.mn_pixelLayersWithMeasurement[i] 
					    );
            phys.leptons[nlep].setLeptonIsoInfo(ev.mn_pileupIsoR03[i],ev.mn_chargedIsoR03[i],ev.mn_photonIsoR03[i],ev.mn_neutralHadIsoR03[i],ev.mn_passIso[i], ev.mn_relIso[i], ev.mn_trkrelIso[i],
                                                ev.en_pileupIso[i],ev.en_chargedIso[i],ev.en_photonIso[i],ev.en_neutralHadIso[i],ev.en_relIsoWithEA[i],ev.en_passIso[i],
                                                ev.ta_IsLooseIso[i], ev.ta_IsMediumIso[i], ev.ta_IsTightIso[i], ev.en_relIso[i]
                                               );
#ifdef YEAR_2017
            phys.leptons[nlep].setLeptonScaleFac(ev.en_enSmearNrSigma[i],ev.en_enScaleValue[i],
                                                 ev.en_enScaleStatUp[i],ev.en_enScaleStatDown[i],ev.en_enScaleSystUp[i],ev.en_enScaleSystDown[i],ev.en_enScaleGainUp[i],ev.en_enScaleGainDown[i],ev.en_enSigmaRhoUp[i],ev.en_enSigmaRhoDown[i],ev.en_enSigmaPhiDown[i]);
#endif

            nlep++;
        }
    }


    for(Int_t i=0; i<ev.ta; i++) {
        LorentzVector P4(ev.ta_px[i],ev.ta_py[i],ev.ta_pz[i],ev.ta_en[i]);
        if(P4.pt()>0) {
            phys.leptons.push_back(PhysicsObject_Lepton(P4, ev.ta_id[i]));
            phys.leptons[nlep].setLeptonIDInfo(ev.mn_IsLoose[i], ev.mn_IsTight[i], ev.mn_IsMedium[i], ev.mn_IsSoft[i],ev.mn_IsHighPt[i],
					       ev.mn_passId[i], ev.mn_passIdLoose[i], ev.mn_passSoftMuon[i], 
                                               ev.en_passVeto[i], ev.en_passLoose[i], ev.en_passMedium[i], ev.en_passTight[i],
					       ev.en_passId[i], ev.en_passIdLoose[i],
                                               ev.ta_dm[i]
                                              );
            phys.leptons[nlep].setLeptonIsoInfo(ev.mn_pileupIsoR03[i],ev.mn_chargedIsoR03[i],ev.mn_photonIsoR03[i],ev.mn_neutralHadIsoR03[i],ev.mn_passIso[i], ev.mn_relIso[i], ev.mn_trkrelIso[i],
                                                ev.en_pileupIso[i],ev.en_chargedIso[i],ev.en_photonIso[i],ev.en_neutralHadIso[i],ev.en_relIsoWithEA[i],ev.en_passIso[i],
                                                ev.ta_IsLooseIso[i], ev.ta_IsMediumIso[i], ev.ta_IsTightIso[i], ev.en_relIso[i]
                                               );

            nlep++;
        }
    }





    // MET
    phys.met 	 = LorentzVector( ev.met_pt*cos(ev.met_phi), ev.met_pt*sin(ev.met_phi), 0, ev.met_pt );
    phys.metNoHF = LorentzVector( ev.metNoHF_pt*cos(ev.metNoHF_phi), ev.metNoHF_pt*sin(ev.metNoHF_phi), 0, ev.metNoHF_pt );

    for(Int_t i=0; i<11; i++) {
      phys.imet[i] = LorentzVector( ev.imet_pt[i]*cos(ev.imet_phi[i]), ev.imet_pt[i]*sin(ev.imet_phi[i]), 0, ev.imet_pt[i] ); 
      phys.variedMet.push_back( phys.imet[i] );
    }



    // Jet
    size_t njet(0);
    for(Int_t i=0; i<ev.jet; i++) {
        LorentzVector P4( ev.jet_px[i],ev.jet_py[i],ev.jet_pz[i],ev.jet_en[i] );
        if(P4.pt()>0) {
	  phys.jets.push_back( PhysicsObject_Jet( P4, ev.jet_puId[i],ev.jet_PFLoose[i],ev.jet_PFTight[i] ) );
	  phys.jets[i].setBtagInfo(ev.jet_btag0[i], ev.jet_btag1[i]); //,ev.jet_btag1[i],ev.jet_btag2[i],ev.jet_btag3[i],ev.jet_btag4[i],ev.jet_btag5[i],ev.jet_btag6[i],ev.jet_btag7[i]);
	  phys.jets[i].setGenInfo(ev.jet_partonFlavour[i], ev.jet_hadronFlavour[i], ev.jet_mother_id[i], ev.jet_parton_px[i], ev.jet_parton_py[i], ev.jet_parton_pz[i], ev.jet_parton_en[i], ev.jet_genpt[i]);
	  
	  njet++;
        }
    }

    // fat Jet
    /*
    size_t nfatjet(0);
    for(Int_t i=0; i<ev.fjet; i++) {
      LorentzVector P4( ev.fjet_px[i],ev.fjet_py[i],ev.fjet_pz[i],ev.fjet_en[i] );
      if (P4.pt()>0) {  
	phys.fatjets.push_back( PhysicsObject_FatJet(P4) );
	phys.fatjets[i].setBtagInfo(ev.fjet_btag0[i]);
	phys.fatjets[i].setSubjetInfo(ev.fjet_prunedM[i], ev.fjet_softdropM[i], ev.fjet_tau1[i], ev.fjet_tau2[i], ev.fjet_tau3[i]);
	phys.fatjets[i].setSubjets(ev.fjet_subjet_count[i], ev.fjet_subjets_px[i], ev.fjet_subjets_py[i], ev.fjet_subjets_pz[i], ev.fjet_subjets_en[i]);
	
	phys.fatjets[i].setGenInfo(ev.fjet_partonFlavour[i], ev.fjet_hadronFlavour[i], ev.fjet_mother_id[i], ev.fjet_parton_px[i], ev.fjet_parton_py[i], ev.fjet_parton_pz[i], ev.fjet_parton_en[i]);

	nfatjet++;
      }
    }
    */

    // secondary vertices 
    size_t nsv(0);
    for (Int_t i=0; i<ev.sv; i++) {
      LorentzVector P4( ev.sv_px[i], ev.sv_py[i], ev.sv_pz[i], ev.sv_en[i] );
      if (P4.pt()>0) {
	phys.svs.push_back( PhysicsObject_SV(P4, ev.sv_chi2[i], ev.sv_ndof[i]) );
	phys.svs[i].setSVinfo(ev.sv_ntrk[i], ev.sv_dxy[i], ev.sv_dxyz[i], ev.sv_dxyz_signif[i], ev.sv_cos_dxyz_p[i]);
	phys.svs[i].setGenInfo(ev.sv_mc_nbh_moms[i], ev.sv_mc_nbh_daus[i], ev.sv_mc_mcbh_ind[i]);

	nsv++;
      }
    }


    //generator level particles
    for(Int_t ipart=0; ipart<ev.nmcparticles; ipart++) {
        LorentzVector p4(ev.mc_px[ipart],ev.mc_py[ipart],ev.mc_pz[ipart],ev.mc_en[ipart]);
	phys.genparticles.push_back(PhysicsObject(p4,ev.mc_id[ipart],ev.mc_mom[ipart],ev.mc_momidx[ipart],ev.mc_status[ipart]) );  
        if(ev.mc_status[ipart]==2 && abs(ev.mc_id[ipart])==15) phys.genleptons.push_back( PhysicsObject(p4,ev.mc_id[ipart],ev.mc_mom[ipart],ev.mc_momidx[ipart],ev.mc_status[ipart]) );
	//        if(ev.mc_status[ipart]!=1) continue;
        switch( ev.mc_id[ipart] ) {
        case 12:
        case -12:
        case 14:
        case -14:
        case 16:
        case -16: {
	  phys.genneutrinos.push_back( PhysicsObject(p4,ev.mc_id[ipart],ev.mc_mom[ipart],ev.mc_momidx[ipart],ev.mc_status[ipart]) );
        }
        break;
	case 35:
        case 36: {
	  phys.genHiggs.push_back( PhysicsObject(p4,ev.mc_id[ipart],ev.mc_mom[ipart],ev.mc_momidx[ipart],ev.mc_status[ipart]) );
        }
        break;
	case 4:
	case -4:
	case 5:
	case -5:
	case 6:
	case -6: {
	  phys.genpartons.push_back (PhysicsObject(p4,ev.mc_id[ipart],ev.mc_mom[ipart],ev.mc_momidx[ipart],ev.mc_status[ipart]) );
	}
	break;  
        case 11:
        case -11:
        case 13:
        case -13:
	case 15:
	case -15:
        {
	  phys.genleptons.push_back( PhysicsObject(p4,ev.mc_id[ipart],ev.mc_mom[ipart],ev.mc_momidx[ipart],ev.mc_status[ipart]) );
        }
        break;
        }
    }

    //add generator jets separately
    for(Int_t ipart=0; ipart<ev.nmcjparticles; ipart++) {
        if(ev.mcj_id[ipart]!=1 || ev.mcj_status[ipart]!=0) continue; //genjet id set
        LorentzVector p4(ev.mcj_px[ipart],ev.mcj_py[ipart],ev.mcj_pz[ipart],ev.mcj_en[ipart]);

        bool overlap(false);
        //check overlap with previous jets
        for(size_t igj=0; igj<phys.genjets.size(); igj++) {
            if(deltaR(p4,phys.genjets[igj])<0.1) {
                overlap=true;
                break;
            }
        }
        if(overlap) continue;

        //check overlap with any other gen particle
        for(Int_t jpart=0; jpart<ev.nmcjparticles; jpart++) {
            if(ev.mcj_id[jpart]==1 && ev.mcj_status[jpart]==0) continue; // not genjet
            if(ev.mcj_status[jpart]!=1) continue; // no intermediator
            LorentzVector jp4(ev.mcj_px[jpart],ev.mcj_py[jpart],ev.mcj_pz[jpart],ev.mcj_en[jpart]);
            if(deltaR(p4,jp4)<0.1) {
                overlap=true;
                break;
            }
        }
        if(overlap) continue;

        phys.genjets.push_back( PhysicsObject(p4,ev.mcj_id[ipart],ev.mcj_mom[ipart],-1,-1) );
    }


    return phys;
}


//
int getLeptonId(int id)
{
  if (id==MUON) return MU;
  if (id==ELECTRON) return E;
  return UNKNOWN;

}

int getDileptonId(int id1, int id2)
{
    if( id1==MUON     && id2==MUON)     return MUMU;
    if( id1==ELECTRON && id2==ELECTRON) return EE;
    if( ( id1==MUON && id2==ELECTRON) || ( id1==ELECTRON && id2==MUON)  ) return EMU;
    return UNKNOWN;
}


bool isDYToLL(int id1, int id2)
{
    if(id1==11 || id2==-11) return true;
    else if(id1==-11 || id2==11) return true;
    else if(id1==13 || id2==-13) return true;
    else if(id1==-13 || id2==13) return true;
    else return false;
}

bool isDYToTauTau(int id1, int id2)
{
    if(id1==15 || id2==-15) return true;
    else if(id1==-15 || id2==15) return true;
    else return false;
}

float getSFfrom1DHist(double xval, TH1F* h_)
{

  if(h_==NULL) {
    cout << "[getSFfrom1DHist]: empty hist! " << endl;
    return 1;
  }
  int xbins = h_->GetXaxis()->GetNbins();
  if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
  if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);


  int binx = h_->GetXaxis()->FindBin(xval);
  float sf_ = h_->GetBinContent(binx);

  if(sf_==0.) return 1.;
  else return sf_;
}

float getSFfrom2DHist(double xval, double yval, TH2F* h_)
{

  if(h_==NULL) {
    cout << "[getSFfrom2DHist]: empty hist! " << endl;
    return 1;
  }
  int xbins = h_->GetXaxis()->GetNbins();
  if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
  if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

  int ybins = h_->GetYaxis()->GetNbins();
  if(yval > h_->GetYaxis()->GetBinUpEdge(ybins)    ) yval = h_->GetYaxis()->GetBinUpEdge(ybins);
  if(yval < h_->GetYaxis()->GetBinLowEdge(1)       ) yval = h_->GetYaxis()->GetBinLowEdge(1);

  int binx = h_->GetXaxis()->FindBin(xval);
  int biny = h_->GetYaxis()->FindBin(yval);
  float sf_ = h_->GetBinContent(binx,biny);

  if(sf_==0.) return 1.;
  else return sf_;
}


float getNLOEWKZZWeight(double trailing_pt)
{
    double ewk_w = 1.0;

    if(     trailing_pt<60.)  ewk_w = 1.-( 4.0/100.);
    else if(trailing_pt<80.)  ewk_w = 1.-( 5.0/100.);
    else if(trailing_pt<100.) ewk_w = 1.-( 6.3/100.);
    else if(trailing_pt<120.) ewk_w = 1.-( 7.6/100.);
    else if(trailing_pt<140.) ewk_w = 1.-( 9.2/100.);
    else if(trailing_pt<160.) ewk_w = 1.-(10.0/100.);
    else if(trailing_pt<180.) ewk_w = 1.-(11.4/100.);
    else if(trailing_pt<200.) ewk_w = 1.-(12.8/100.);
    else if(trailing_pt<220.) ewk_w = 1.-(14.2/100.);
    else if(trailing_pt<240.) ewk_w = 1.-(15.6/100.);
    else if(trailing_pt<260.) ewk_w = 1.-(17.0/100.);
    else if(trailing_pt<280.) ewk_w = 1.-(18.4/100.);
    else if(trailing_pt<300.) ewk_w = 1.-(20.0/100.);
    else if(trailing_pt<320.) ewk_w = 1.-(21.2/100.);
    else if(trailing_pt<340.) ewk_w = 1.-(22.4/100.);
    else if(trailing_pt<360.) ewk_w = 1.-(23.6/100.);
    else if(trailing_pt<380.) ewk_w = 1.-(24.8/100.);
    else                      ewk_w = 1.-(26.0/100.);

    return ewk_w;
}


float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ)
{

    float k = 0.0;
    k += 1.513834489150 * (abs(GENdPhiZZ)>0.0 && abs(GENdPhiZZ)<=0.1);
    k += 1.541738780180 * (abs(GENdPhiZZ)>0.1 && abs(GENdPhiZZ)<=0.2);
    k += 1.497829632510 * (abs(GENdPhiZZ)>0.2 && abs(GENdPhiZZ)<=0.3);
    k += 1.534956782920 * (abs(GENdPhiZZ)>0.3 && abs(GENdPhiZZ)<=0.4);
    k += 1.478217033060 * (abs(GENdPhiZZ)>0.4 && abs(GENdPhiZZ)<=0.5);
    k += 1.504330859290 * (abs(GENdPhiZZ)>0.5 && abs(GENdPhiZZ)<=0.6);
    k += 1.520626246850 * (abs(GENdPhiZZ)>0.6 && abs(GENdPhiZZ)<=0.7);
    k += 1.507013090030 * (abs(GENdPhiZZ)>0.7 && abs(GENdPhiZZ)<=0.8);
    k += 1.494243156250 * (abs(GENdPhiZZ)>0.8 && abs(GENdPhiZZ)<=0.9);
    k += 1.450536096150 * (abs(GENdPhiZZ)>0.9 && abs(GENdPhiZZ)<=1.0);
    k += 1.460812521660 * (abs(GENdPhiZZ)>1.0 && abs(GENdPhiZZ)<=1.1);
    k += 1.471603622200 * (abs(GENdPhiZZ)>1.1 && abs(GENdPhiZZ)<=1.2);
    k += 1.467700038200 * (abs(GENdPhiZZ)>1.2 && abs(GENdPhiZZ)<=1.3);
    k += 1.422408690640 * (abs(GENdPhiZZ)>1.3 && abs(GENdPhiZZ)<=1.4);
    k += 1.397184022730 * (abs(GENdPhiZZ)>1.4 && abs(GENdPhiZZ)<=1.5);
    k += 1.375593447520 * (abs(GENdPhiZZ)>1.5 && abs(GENdPhiZZ)<=1.6);
    k += 1.391901318370 * (abs(GENdPhiZZ)>1.6 && abs(GENdPhiZZ)<=1.7);
    k += 1.368564350560 * (abs(GENdPhiZZ)>1.7 && abs(GENdPhiZZ)<=1.8);
    k += 1.317884804290 * (abs(GENdPhiZZ)>1.8 && abs(GENdPhiZZ)<=1.9);
    k += 1.314019950800 * (abs(GENdPhiZZ)>1.9 && abs(GENdPhiZZ)<=2.0);
    k += 1.274641749910 * (abs(GENdPhiZZ)>2.0 && abs(GENdPhiZZ)<=2.1);
    k += 1.242346606820 * (abs(GENdPhiZZ)>2.1 && abs(GENdPhiZZ)<=2.2);
    k += 1.244727403840 * (abs(GENdPhiZZ)>2.2 && abs(GENdPhiZZ)<=2.3);
    k += 1.146259351670 * (abs(GENdPhiZZ)>2.3 && abs(GENdPhiZZ)<=2.4);
    k += 1.107804993520 * (abs(GENdPhiZZ)>2.4 && abs(GENdPhiZZ)<=2.5);
    k += 1.042053646740 * (abs(GENdPhiZZ)>2.5 && abs(GENdPhiZZ)<=2.6);
    k += 0.973608545141 * (abs(GENdPhiZZ)>2.6 && abs(GENdPhiZZ)<=2.7);
    k += 0.872169942668 * (abs(GENdPhiZZ)>2.7 && abs(GENdPhiZZ)<=2.8);
    k += 0.734505279177 * (abs(GENdPhiZZ)>2.8 && abs(GENdPhiZZ)<=2.9);
    k += 1.163152837230 * (abs(GENdPhiZZ)>2.9 && abs(GENdPhiZZ)<=3.1416);

    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;
}
