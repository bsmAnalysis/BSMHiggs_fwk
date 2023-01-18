# comment out Wh or Zh part
###################### ZH #####################
#cd ../cards_SB13TeV_SM_$1$2_noSoftb/0060/    
cd ../cards_SB13TeV_SM_$1$2_noSoftb/
cards_dir=`pwd`

cd 0060/

if [[ $cards_dir == *"Wh_all"* ]]; then 
    printf "\n Found Wh_all\n\n"
    rangeString="tt_norm_e_2016=0.1,4.0:tt_norm_mu_2016=0.1,4.0:w_norm_e_2016=0.1,4:w_norm_mu_2016=0.1,4:tt_norm_e_2017=0.1,4.0:tt_norm_mu_2017=0.1,4.0:w_norm_e_2017=0.1,4:w_norm_mu_2017=0.1,4:tt_norm_e_2018=0.1,4.0:tt_norm_mu_2018=0.1,4.0:w_norm_e_2018=0.1,4:w_norm_mu_2018=0.1,4"
elif [[ $cards_dir == *"Zh_all"* ]]; then
    printf "\n Found Zh_all\n\n"
    rangeString="z_norm_3b_e_2016=0.1,4.0:z_norm_3b_mu_2016=0.1,4.0:z_norm_4b_e_2016=0.1,4.0:z_norm_4b_mu_2016=0.1,4.0:tt_norm_e_2016=0.1,4.0:tt_norm_mu_2016=0.1,4.0:z_norm_3b_e_2017=0.1,4.0:z_norm_3b_mu_2017=0.1,4.0:z_norm_4b_e_2017=0.1,4.0:z_norm_4b_mu_2017=0.1,4.0:tt_norm_e_2017=0.1,4.0:tt_norm_mu_2017=0.1,4.0:z_norm_3b_e_2018=0.1,4.0:z_norm_3b_mu_2018=0.1,4.0:z_norm_4b_e_2018=0.1,4.0:z_norm_4b_mu_2018=0.1,4.0:tt_norm_e_2018=0.1,4.0:tt_norm_mu_2018=0.1,4.0"
elif [[ $cards_dir == *"Wh_"* ]]; then
    printf "\n Found Wh_\n\n"
    rangeString="tt_norm_e=0.1,4.0:tt_norm_mu=0.1,4.0:w_norm_e=0.1,4:w_norm_mu=0.1,4"
elif [[ $cards_dir == *"Zh_"* ]]; then
    printf "\n Found Zh_\n\n"
    rangeString="z_norm_3b_e=0.1,4.0:z_norm_3b_mu=0.1,4.0:z_norm_4b_e=0.1,4.0:z_norm_4b_mu=0.1,4.0:tt_norm_e=0.1,4.0:tt_norm_mu=0.1,4.0"
elif [[ $cards_dir == *"SM_2"* ]]; then
    printf "\n Found SM_2\n\n"
    rangeString="tt_norm_e=0.1,4.0:tt_norm_mu=0.1,4.0:w_norm_e=0.1,4:w_norm_mu=0.1,4:z_norm_3b_e=0.1,4.0:z_norm_3b_mu=0.1,4.0:z_norm_4b_e=0.1,4.0:z_norm_4b_mu=0.1,4.0"
elif [[ $cards_dir == *"SM_all"* ]]; then
    printf "\n Found SM_all\n\n"
    rangeString="tt_norm_e_2016=0.1,4.0:tt_norm_mu_2016=0.1,4.0:w_norm_e_2016=0.1,4:w_norm_mu_2016=0.1,4:z_norm_3b_e_2016=0.1,4.0:z_norm_3b_mu_2016=0.1,4.0:z_norm_4b_e_2016=0.1,4.0:z_norm_4b_mu_2016=0.1,4.0:tt_norm_e_2017=0.1,4.0:tt_norm_mu_2017=0.1,4.0:w_norm_e_2017=0.1,4:w_norm_mu_2017=0.1,4:z_norm_3b_e_2017=0.1,4.0:z_norm_3b_mu_2017=0.1,4.0:z_norm_4b_e_2017=0.1,4.0:z_norm_4b_mu_2017=0.1,4.0:tt_norm_e_2018=0.1,4.0:tt_norm_mu_2018=0.1,4.0:w_norm_e_2018=0.1,4:w_norm_mu_2018=0.1,4:z_norm_3b_e_2018=0.1,4.0:z_norm_3b_mu_2018=0.1,4.0:z_norm_4b_e_2018=0.1,4.0:z_norm_4b_mu_2018=0.1,4.0"
else
    printf "\n\n *** can't figure out what kind of cards directory this is:  %s\n\n" $cards_dir
    exit -1
fi

printf "\n\n rangeString = %s\n\n" $rangeString

# run distribution for ee channel
combine -M GoodnessOfFit -d workspace.root --algo=saturated -n _result_sb --rMin -5.0 --rMax 5.0 --setParameterRanges $rangeString --cminDefaultMinimizerStrategy 0
#--freezeParameters r --setParameters r=0 
#--setParametersForFit mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --setParametersForEval mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --freezeParameters r --setParameters r=0
combine -M GoodnessOfFit -d workspace.root --algo=saturated -n result_toy_sb -t 200 --toysFrequentist   --rMin -5.0 --rMax 5.0 --setParameterRanges $rangeString --cminDefaultMinimizerStrategy 0
#--freezeParameters r --setParameters r=0 
#--setParametersForFit mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --setParametersForEval mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 --freezeParameters r --setParameters r=0,mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 -t 500 --toysFrequentist

# run distributioni for mumu channel
#combine -M GoodnessOfFit -d workspace_mu.root --algo=saturated --setParametersForFit mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --setParametersForEval mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --freezeParameters r --setParameters r=0
#combine -M GoodnessOfFit -d workspace_mu.root --algo=saturated --setParametersForFit mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --setParametersForEval mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --freezeParameters r --setParameters r=0,mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1 -t 500 --toysFrequentist


##################### WH ######################
#cd /afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/cards_SB13TeV_SM_Wh_2018_noSoftb/0060/

#combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --setParametersForEval mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --freezeParameters r --setParameters r=0
#combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --setParametersForEval mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 --freezeParameters r --setParameters r=0,mask_e_A_SR_3b=1,mask_mu_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_4b=1 -t 500 --toysFrequentist

