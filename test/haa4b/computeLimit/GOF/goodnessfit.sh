# cd to your card directory
cd /afs/cern.ch/work/y/yuanc/Analysis/H2a4b/CMSSW_10_2_13/src/UserCode/bsmhiggs_fwk/test/haa4b/computeLimit/JOBS/SB13TeV_SM_Zh_backup/test/out_60_new


######## ZH #########
tt_e=`cat simfit_m60_e.txt | grep 'tt_norm_e' | awk '{print $4;}'`;
tt_mu=`cat simfit_m60_mu.txt | grep 'tt_norm_mu' | awk '{print $4;}'`;
v_3b_e=`cat simfit_m60_e.txt | grep 'z_norm_3b_e' | awk '{print $4;}'`;
v_4b_e=`cat simfit_m60_e.txt | grep 'z_norm_4b_e' | awk '{print $4;}'`;
v_3b_mu=`cat simfit_m60_mu.txt | grep 'z_norm_3b_mu' | awk '{print $4;}'`;
v_4b_mu=`cat simfit_m60_mu.txt | grep 'z_norm_4b_mu' | awk '{print $4;}'`;

combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1,mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --setParametersForEval mask_ee_A_SR_3b=0,mask_ee_A_SR_4b=0,mask_mumu_A_SR_3b=0,mask_mumu_A_SR_4b=0 --freezeParameters r,tt_norm_e,z_norm_3b_e,z_norm_4b_e,tt_norm_mu,z_norm_3b_mu,z_norm_4b_mu --setParameters r=0,tt_norm_e=$tt_e,z_norm_3b_e=$v_3b_e,z_norm_4b_e=$v_4b_e,tt_norm_mu=$tt_mu,z_norm_3b_mu=$v_3b_mu,z_norm_4b_mu=$v_4b_mu
combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1,mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1 --setParametersForEval mask_ee_A_SR_3b=0,mask_ee_A_SR_4b=0,mask_mumu_A_SR_3b=0,mask_mumu_A_SR_4b=0 --freezeParameters r,tt_norm_e,z_norm_3b_e,z_norm_4b_e,tt_norm_mu,z_norm_3b_mu,z_norm_4b_mu --setParameters r=0,mask_ee_A_SR_3b=1,mask_ee_A_SR_4b=1,mask_mumu_A_SR_3b=1,mask_mumu_A_SR_4b=1,tt_norm_e=$tt_e,z_norm_3b_e=$v_3b_e,z_norm_4b_e=$v_4b_e,tt_norm_mu=$tt_mu,z_norm_3b_mu=$v_3b_mu,z_norm_4b_mu=$v_4b_mu -t 100 --toysFrequentist # -s -1

######## WH #########
#tt_e=`cat simfit_m60_e.txt | grep 'tt_norm_e' | awk '{print $4;}'`;
#tt_mu=`cat simfit_m60_mu.txt | grep 'tt_norm_mu' | awk '{print $4;}'`;
#v_e=`cat simfit_m60_e.txt | grep 'w_norm_e' | awk '{print $4;}'`;
#v_mu=`cat simfit_m60_mu.txt | grep 'w_norm_mu' | awk '{print $4;}'`;

#combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_e_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_3b=1,mask_mu_A_SR_4b=1 --setParametersForEval mask_e_A_SR_3b=0,mask_e_A_SR_4b=0,mask_mu_A_SR_3b=0,mask_mu_A_SR_4b=0 --freezeParameters r,tt_norm_e,w_norm_e,tt_norm_mu,w_norm_mu --setParameters r=0,tt_norm_e=$tt_e,w_norm_e=$v_e,tt_norm_mu=$tt_mu,w_norm_mu=$v_mu
#combine -M GoodnessOfFit -d workspace.root --algo=saturated --setParametersForFit mask_e_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_3b=1,mask_mu_A_SR_4b=1 --setParametersForEval mask_e_A_SR_3b=0,mask_e_A_SR_4b=0,mask_mu_A_SR_3b=0,mask_mu_A_SR_4b=0 --freezeParameters r,tt_norm_e,w_norm_e,tt_norm_mu,w_norm_mu --setParameters r=0,mask_e_A_SR_3b=1,mask_e_A_SR_4b=1,mask_mu_A_SR_3b=1,mask_mu_A_SR_4b=1,tt_norm_e=$tt_e,w_norm_e=$v_e,tt_norm_mu=$tt_mu,w_norm_mu=$v_mu -t 100 --toysFrequentist #-s -1


