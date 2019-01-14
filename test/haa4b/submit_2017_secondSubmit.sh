#!/usr/bin/env bash

#--------------------------------------------------
# Global Code 
#--------------------------------------------------

if [[ $# -eq 0 ]]; then 
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
## Run Analysis over samples
    printf "\n\t%-5s  %-40s\n"  "0"  "completely clean up the directory" 
    printf "\n\t%-5s  %-40s\n"  "1.0"  "run 'runNtuplizer' on all MINIAOD samples"
#    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runhaaAnalysis' on all Ntuples" 

## Merge Results
#    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples" 
#    printf "\n\t%-5s  %-40s\n"  "3.0"  "make plots and combine root files" 

## Make plots in mcbased(_blind), datadriven(_blind) cases
#    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for mcbased analysis"  
#    printf "\n\t%-5s  %-40s\n"  "3.01"  "make root file for input to the Limits"  
#    printf "\n\t%-5s  %-40s\n"  "3.02"  "make plots for QCD analysis"
#    printf "\n\t%-5s  %-40s\n"  "3.2"  "make plots with data-driven bkgs"
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
#arguments='crab3'; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 

#--------------------------------------------------
# Global Variables
#--------------------------------------------------

SUFFIX=_2018_12_06 #_2018_03_03

#SUFFIX=$(date +"_%Y_%m_%d") 
MAINDIR=$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b
#JSON=$MAINDIR/samples2017.json
JSON=$MAINDIR/samples2017_secondSubmit.json 
#GOLDENJSON=$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/data/json/

RESULTSDIR=$MAINDIR/results$SUFFIX

if [[ $arguments == *"crab3"* ]]; then STORAGEDIR='';
else STORAGEDIR=/eos/cms/store/user/georgia/results$SUFFIX ; fi

PLOTSDIR=$MAINDIR/plots${SUFFIX}
PLOTTER=$MAINDIR/plotter${SUFFIX}
 
####################### Settings for Ntuple Analysis ##################
NTPL_INPUT=results$SUFFIX

#NTPL_JSON=$MAINDIR/samples2017.json
NTPL_JSON=$MAINDIR/samples2017_secondSubmit.json 
NTPL_OUTDIR=$MAINDIR/results_Ntpl$SUFFIX
#NTPL_OUTDIR=/eos/cms/store/user/georgia/results_Ntpl$SUFFIX #only for Data
RUNLOG=$NTPL_OUTDIR/LOGFILES/runSelection.log

queue='8nh'   

#IF CRAB3 is provided in argument, use crab submission instead of condor/lsf 
if [[ $arguments == *"crab3"* ]]; then queue='crab3' ;fi  

################################################# STEPS between 0 and 1
if [[ $step == 0 ]]; then   
        #analysis cleanup
    echo "Really delete directory "$RESULTSDIR" ?" 
    echo "ALL DATA WILL BE LOST! [N/y]?"
    read answer
    if [[ $answer == "y" ]];
    then
	echo "CLEANING UP..."
	rm -rdf $RESULTSDIR $PLOTSDIR LSFJOB_* core.* *.sh.e* *.sh.o*
    fi
fi #end of step0
if [[ $step == 0.1 ]]; then
    echo "Really delete directory "$NTPL_OUTDIR" ?"
    echo "ALL DATA WILL BE LOST! [N/y]?"
    read answer
    if [[ $answer == "y" ]];
    then
	echo "CLEANING UP..."
	rm -rdf $NTPL_OUTDIR LSFJOB_* core.* *.sh.e* *.sh.o*
    fi
fi

###  ############################################## STEPS between 1 and 2
if [[ $step > 0.999 &&  $step < 2 ]]; then
   if [[ $step == 1.0 ]]; then
       echo "2017 JOB SUBMISSION for Ntuplization using full CMSSW fwk"
       echo -e "Input: " $JSON "\nOutput: " $RESULTSDIR
#       runAnalysisOverSamples.py -j $JSON -o $RESULTSDIR  -c $MAINDIR/../fullAnalysis_cfg_2017.py.templ -l results$SUFFIX -p "@verbose=False" --key haa_signal -s crab 
       runAnalysisOverSamples.py -j $JSON -o $RESULTSDIR  -c $MAINDIR/../fullAnalysis_cfg_2017.py.templ -l results$SUFFIX -p "@verbose=False" --key haa_mcbased -s crab 
   fi    

#   if [[ $step == 1.1 ]]; then  #submit jobs for h->aa->XXYY analysis
#       echo "JOB SUBMISSION for BSM h->aa Analysis"
#       echo -e "Input: " $NTPL_JSON "\n Output: " $NTPL_OUTDIR
#       ## if the output directory does not exist, create it:
#       if [ ! -d "$NTPL_OUTDIR" ]; then
#	   mkdir $NTPL_OUTDIR
#       fi
#       runLocalAnalysisOverSamples.py -e runhaaAnalysis -g $RUNLOG -j $NTPL_JSON -o $NTPL_OUTDIR -d $NTPL_INPUT -c $MAINDIR/../runNtplAnalysis_cfg.py.templ -p "@runSystematics=True @runMVA=False @reweightTopPt=False @usemetNoHF=False @verbose=False @useDeepCSV=True @runQCD=False" -s $queue #-r True
#   fi
fi

###  ############################################## STEPS between 2 and 3
#if [[ $step > 1.999 && $step < 3 ]]; then
#   if [[ $step == 2 ]]; then    #extract integrated luminosity of the processed lumi blocks
#	echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
#	mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
#	mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
#	mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
#	mergeJSON.py --output=$RESULTSDIR/json_muEG.json       $RESULTSDIR/Data*_MuEG*.json
#	mergeJSON.py --output=$RESULTSDIR/json_in.json  $GOLDENJSON/Cert_*.txt
#	echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
#	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json 
#	echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
#	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json 
#	echo "MISSING LUMI BLOCKS IN MUON EGAMMA DATASET"
#	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_muEG.json 
#
#	echo "COMPUTE INTEGRATED LUMINOSITY"
#	export LD_LIBRARY_PATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root/lib
#        export PYTHONPATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root/lib 
#        export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH 
#        export ROOTSYS=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root 
#        export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH 
##	export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
#
#	pip uninstall brilws -y 
#	pip install --upgrade --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time
#
#	if [[ $JSON =~ "2016" ]]; then
#	    brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
#	else
#	    brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
#	fi
#	tail -n 3 $RESULTSDIR/LUMI.txt  
#     fi
#  fi     
#
####  ############################################## STEPS between 3 and 4
#if [[ $step > 2.999 && $step < 4 ]]; then
#    if [ -f $NTPL_OUTDIR/LUMI.txt ]; then
#      INTLUMI=`tail -n 3 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
#    else
#	if [[ $JSON =~ "2016" ]]; then  
#	    INTLUMI=35866.932
#            echo "Please run step==2 above to calculate int. luminosity for 2016 data!" 
#        else
#	    echo "Please run step==2 above to calculate int. luminosity!"
#        fi                                                                                                                   
#	echo "WARNING: $RESULTSDIR/LUMI.txt file is missing so use fixed integrated luminosity value, this might be different than the dataset you ran on"
#    fi
#    
#    if [[ $step == 3 || $step == 3.0 ]]; then  # make plots and combined root files
#	echo "MAKE SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
#	echo "Input DIR = "$NTPL_OUTDIR
#        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $NTPL_OUTDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased $arguments        
#   fi        
#
#    if [[ $step == 3 || $step == 3.01 ]]; then  # make plots and combined root files for limits only
#	echo "MAKE SUMMARY ROOT FILE FOR LIMITS, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI" 
#	echo "Input DIR = "$NTPL_OUTDIR  
#	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}_v1/ --outFile ${PLOTTER}_forLimits_v1.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased --only "(all_optim_systs|all_optim_cut|(e|mu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_(bdt)(|_shapes)(|_umetup|_umetdown|_jerup|_jerdown|_jesup|_jesdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_pdfup|_pdfdown|_btagup|_btagdown))" $arguments 
#    fi
#
#    if [[ $step == 3.02 ]]; then # make plots for data-driven QCD bkg
#	echo "MAKE SUMMARY ROOT FILE, for data-driven QCD estimate"
#	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outFile ${PLOTTER}_qcd.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased --only "(all_optim_systs|all_optim_cut|(e|mu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_(bdt)(|_shapes)(|_umetup|_umetdown|_jerup|_jerdown|_jesup|_jesdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_pdfup|_pdfdown|_btagup|_btagdown))" $arguments
#	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased_qcd/ --outFile ${PLOTTER}_qcd.root  --json $JSON --plotExt .pdf --key haa_mcbased --fileOption READ --noLog --only "(e|mu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_bdt" $arguments
#    fi
#
##(ht|pfmet|ptw|mtw|higgsPt|higgsMass|dRave|dmmin|dphijmet|dphiWh|
#    if [[ $step == 3 || $step == 3.1 ]]; then  # make plots and combine root files for mcbased study    
#	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased --only "((e|mu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_(pfmet|ht|mtw|ptw|dphiWh|dRave|dmmin|higgsMass|higgsPt|dphijmet|dphijmet1|dphijmet12|bdt))|(all_(eventflow|nvtx_raw|nvtxwgt_raw)|((presel1|presel2)_(nvtx_raw|nvtxwgt_raw))|(e|mu)(|_presf|_nj2)_leadlep_(pt|eta)_raw|raw_dphijmet|raw_j1_dphijmet|raw_minj1j2_dphijmet|nj_njets_raw|(jet_(b1|b2|b3|b4)_jet_(pt|eta)_raw)|nb_nbjets_raw|(DeepCSV_(b1|b2|b3|b4)_jet_(pt|eta)_raw)|nb_soft_nbjets_raw|(softb_(b1|b2|b3|b4)_softjet_pt_raw)|(softb_(b1|b2|b3|b4)_jet_eta_raw)|raw_pfmet|raw_mtw|raw_ptw|raw_softb_ntrk|raw_softb_dxy|raw_softb_dxyz_signif|raw_softb_cos|sv_b_dR_raw|nb_(LOOSE|MEDIUM|TIGHT)_nbjets_raw|(sel1|sel2|sel3)_evt_cat|DeepCSV_(b1|b2|b3|b4|high)_b_discrim)" $arguments 
#	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased/ --outFile ${PLOTTER}.root  --json $JSON --plotExt .png --plotExt .pdf --key haa_mcbased --fileOption READ --noLog --signalScale 100 $arguments 
#        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased_blind/ --outFile ${PLOTTER}.root  --json $JSON --plotExt .png --plotExt .pdf --key haa_mcbased --fileOption READ --noLog --signalScale 100 --blind 0.1 --only "(e|mu)_A_SR_(3b|4b)_bdt" $arguments
#    fi
#
#    if [[ $step == 3 || $step == 3.2 ]]; then # make plots and combine root files for data-driven backgrounds 
#	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key haa_datadriven $arguments
#        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven/  --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key haa_datadriven --fileOption READ $arguments
#    fi	
#fi

