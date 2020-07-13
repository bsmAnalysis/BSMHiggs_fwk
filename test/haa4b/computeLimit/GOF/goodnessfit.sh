# comment out Wh or Zh part
###################### ZH #####################
cd /afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/cards_SB13TeV_SM_Zh_2016_noSoftb/0060/

# run distribution for ee channel
combine -M GoodnessOfFit -d workspace_e.root --algo=saturated --setParametersForFit mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --setParametersForEval mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --freezeParameters r --setParameters r=0
combine -M GoodnessOfFit -d workspace_e.root --algo=saturated --setParametersForFit mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --setParametersForEval mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --freezeParameters r --setParameters r=0,mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 -t 500 --toysFrequentist

# run distributioni for mumu channel
#combine -M GoodnessOfFit -d workspace_mu.root --algo=saturated --setParametersForFit mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --setParametersForEval mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --freezeParameters r --setParameters r=0
#combine -M GoodnessOfFit -d workspace_mu.root --algo=saturated --setParametersForFit mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --setParametersForEval mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --freezeParameters r --setParameters r=0,mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 -t 500 --toysFrequentist


##################### WH ######################
#cd /afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/cards_SB13TeV_SM_Wh_2018_noSoftb/0060/

#combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --setParametersForEval mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --freezeParameters r --setParameters r=0
#combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --setParametersForEval mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --freezeParameters r --setParameters r=0,mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 -t 500 --toysFrequentist

