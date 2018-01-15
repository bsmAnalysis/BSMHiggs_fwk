#include <string>
#include <map>

const std::map<std::string, int> mStat = 
  {
    //  Map between dtag and total number of events per process (Npos - Nneg)
    // including Negative weight events

    //h->aa->4b Signal
    { "MC13TeV_Wh_amass12", 377600 }, //3.62818e-06 },
    { "MC13TeV_Wh_amass15", 379120 }, //3.61363e-06 },
    { "MC13TeV_Wh_amass20", 383440 },
    { "MC13TeV_Wh_amass25", 389360 },
    { "MC13TeV_Wh_amass30", 389560 },
    { "MC13TeV_Wh_amass40", 398000 },
    { "MC13TeV_Wh_amass50", 399000 },
    { "MC13TeV_Wh_amass60", 388720 },

    //Single top
    { "MC13TeV_SingleT_s_2016", 622990 }, //0.193558 },
    { "MC13TeV_SingleT_atW_2016", 6933094 }, //0.184169 },
    { "MC13TeV_SingleT_tW_2016", 6952830 }, //0.183646 },
    { "MC13TeV_SingleT_at_2016", 37904304 }, //67240808 }, //0.0765989 },
    { "MC13TeV_SingleT_t_2016", 67240808 }, //37904304 }, //0.0725545 },
    //TTJets
    { "MC13TeV_TTJets_slt_2016", 61973977 }, //0.105736*0.1929 },
    { "MC13TeV_TTJets_slt_ext1_2016", 61973977 }, //0.105736*0.8071 },
    { "MC13TeV_TTJets_slat_2016", 60210394 }, // 0.108833*0.1984 },
    { "MC13TeV_TTJets_slat_ext1_2016", 60210394 }, //0.109933*0.8016 },
    { "MC13TeV_TTJets_dl_2016", 30444678 }, //0.104012*0.2002 },
    { "MC13TeV_TTJets_dl_ext1_2016", 30444678 }, // 0.104012*0.7998 },
    //WJets
    { "MC13TeV_WJets_2016", 86731806 }, //30.7868*0.3425 },
    { "MC13TeV_WJets_ext2_2016", 86731806 }, // 30.7868*0.6575 },
    //DY
    { "MC13TeV_DYJetsToLL_10to50_2016", 51433009 }, //12.9777*0.434 },
    { "MC13TeV_DYJetsToLL_10to50_ext1_2016", 51433009 }, // 12.9777*0.566 },
    { "MC13TeV_DYJetsToLL_50toInf_2016", 81781064 }, //2.52854 },
    //TG, TTG, TTW, TTZ
    { "MC13TeV_TGJets_2016", 368562 }, //0.288736*0.25 },
    { "MC13TeV_TGJets_ext1_2016", 368562 }, //0.288736*0.75 },
    { "MC13TeV_TTGJets_2016", 1577833 }, // 0.0840393 },
    { "MC13TeV_TTWJetslnu_2016", 1603527 }, //0.00456969 },
    { "MC13TeV_TTZJets2l2nu_2016", 927976 }, //0.00977476 },
    //Di-boson
    { "MC13TeV_WW2l2nu_2016", 1999000 }, // 0.218503 },
    { "MC13TeV_WWlnu2q_2016", 8997800 }, //0.199298*0.5 },
    { "MC13TeV_WWlnu2q_ext1_2016", 8997800 }, //0.199298*0.5 },
    { "MC13TeV_WZ_2016", 3871065 }, //0.436678*0.25 },
    { "MC13TeV_WZ_ext1_2016", 3871065 }, //0.436678*0.75 },
    { "MC13TeV_ZZ_2016", 1988098 }, //0.298089*0.5 },
    { "MC13TeV_ZZ_ext1_2016", 1988098 }, //0.298089*0.5 },
    //Tri-boson
    { "MC13TeV_ZZZ_2016", 213197 }, //0.0355368 },
    { "MC13TeV_WZZ_2016", 216366 }, //0.0267381 },
    { "MC13TeV_WWZ_2016", 221468 }, //0.00922509 },
    { "MC13TeV_WWW_4F_2016", 210538 }, //0.00235191 },
    //QCD
    { "MC13TeV_QCD_Pt50To80_EMEnr_2016", 23474171 },
    { "MC13TeV_QCD_Pt80To120_EMEnr_2016", 35841783 }, // 350.246 },
    { "MC13TeV_QCD_Pt120To170_EMEnr_2016", 35817281 }, //63.0513 },
    { "MC13TeV_QCD_Pt170To300_EMEnr_2016", 11540163 }, //58.4617 },
    { "MC13TeV_QCD_Pt300ToInf_EMEnr_2016", 7373633 }, //6.56669 },
    { "MC13TeV_QCD_Pt20ToInf_MuEnr_2016", 22094081 } //491.35 }

  };
