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
    printf "\n\t%-5s  %-40s\n"  "1.1"  "run 'runhaaAnalysis' on all Ntuples" 

## Merge Results
    printf "\n\t%-5s  %-40s\n"  "2"  "compute integrated luminosity from processed samples" 
    printf "\n\t%-5s  %-40s\n"  "3.0"  "make plots and combine root files" 

## Make plots in mcbased(_blind), datadriven(_blind) cases
    printf "\n\t%-5s  %-40s\n"  "3.1"  "make plots for mcbased analysis"  
    printf "\n\t%-5s  %-40s\n"  "3.01"  "make root file for input to the Limits"  
    printf "\n\t%-5s  %-40s\n"  "3.02"  "make plots for QCD analysis"
    printf "\n\t%-5s  %-40s\n"  "3.2"  "make plots with data-driven bkgs"
fi

step=$1   #variable that store the analysis step to run

#Additional arguments to take into account
arguments=''; for var in "${@:2}"; do arguments=$arguments" "$var; done
#arguments='crab3'; for var in "${@:2}"; do arguments=$arguments" "$var; done
if [[ $# -ge 4 ]]; then echo "Additional arguments will be considered: "$arguments ;fi 

#--------------------------------------------------
# Global Variables
#--------------------------------------------------

YEAR=2016
CHANNEL=WH

do_syst=True # Always run with Systematics, unless its QCD mode

if [[ $CHANNEL == "WH_QCD" ]]; then 
    doQCD=True ; do_syst=False
else doQCD=False ; fi

if [[ $CHANNEL == "ZH" ]]; then doZH=True ; 
else doZH=False ; fi

if [[ $YEAR == "2016" ]]; then SUFFIX=_2023_11_15 
else SUFFIX=_2020_02_05 ; fi

MAINDIR=$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b

# Json and python template for all years
JSON=$MAINDIR/samples$YEAR.json
NTPL_JSON=$MAINDIR/samples$YEAR.json
FULLANALYSISCFG=$MAINDIR/../fullAnalysis_cfg_$YEAR.py.templ
RUNNTPLANALYSISCFG=$MAINDIR/../runNtplAnalysis_cfg_$YEAR.py.templ
    
#SUFFIX=$(date +"_%Y_%m_%d") 
GOLDENJSON=$CMSSW_BASE/src/UserCode/bsmhiggs_fwk/data/json/

RESULTSDIR=$MAINDIR/results_$YEAR$SUFFIX 

if [[ $arguments == *"crab3"* ]]; then STORAGEDIR='';
else STORAGEDIR=/eos/cms/store/user/e/esiamark/results$SUFFIX ; fi

PLOTSDIR=$MAINDIR/plots_${CHANNEL}_${YEAR}${SUFFIX}
PLOTTER=$MAINDIR/plotter_${CHANNEL}_${YEAR}${SUFFIX}
 
####################### Settings for Ntuple Analysis ##################
if [[ $YEAR == "2016" ]]; then
    NTPL_INPUT=/eos/user/e/esiamark/results_2016$SUFFIX
elif [[ $YEAR == "2017" ]]; then
    NTPL_INPUT=/eos/cms/store/user/e/esiamark/results_2017$SUFFIX 
elif [[ $YEAR == "2018" ]]; then
    NTPL_INPUT=/eos/cms/store/user/e/esiamark/results_2018$SUFFIX   
fi

ZPtSF_OUT=$MAINDIR/VPtSF_${CHANNEL}_$YEAR$SUFFIX
BTAG_NTPL_OUTDIR=$MAINDIR/btag_SFs/$YEAR/btag_Ntpl$SUFFIX
NTPL_OUTDIR=$MAINDIR/results_Ntpl_${CHANNEL}_$YEAR${SUFFIX}
#NTPL_OUTDIR=$EOSDIR/results_Ntpl_${CHANNEL}_$YEAR$SUFFIX

RUNLOG=$NTPL_OUTDIR/LOGFILES/runSelection.log
queue='workday'   

TopSF_INPUT=$MAINDIR/PrefitPlots_${YEAR}${CHANNEL}_noSoftb/TEST_ht
TopSF_OUT=$MAINDIR/TopPtSF_$YEAR$CHANNEL
if [[ $CHANNEL == "WH" ]] ; then
    vh_tag="wh" # wh or zh channel when computing Top Pt weights
else vh_tag="zh" ; fi

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
       echo "JOB SUBMISSION for Ntuplization using full CMSSW fwk"
       echo -e "Input: " $JSON "\nOutput: " $RESULTSDIR
       #runAnalysisOverSamples.py -j $JSON -o $RESULTSDIR  -c $MAINDIR/../fullAnalysis_cfg.py.templ -l results$SUFFIX -p "@verbose=False" --key haa_mcbased -s crab 
       runAnalysisOverSamples.py -j $JSON -o $RESULTSDIR  -c $FULLANALYSISCFG -l results_$YEAR$SUFFIX -p "@verbose=False" --key haa_signal -s crab
       # Ntuplize 2016 signal samples under 94X:
       #runAnalysisOverSamples.py -j $MAINDIR/samples2016.json -o $RESULTSDIR  -c $MAINDIR/../fullAnalysis_cfg_2016Signal.py.templ -l results$SUFFIX -p "@verbose=False" --key haa_signal  -s crab
   fi    

   if [[ $step == 1.01 ]]; then  # compute BTagging efficiency
       echo "Merge and Calculate BTag Efficiency in MC"
       echo -e "Input: " $NTPL_INPUT "\n Output: " $BTAG_NTPL_OUTDIR
       ## if the output directory does not exist, create it:
       if [ ! -d "$BTAG_NTPL_OUTDIR" ]; then
	   mkdir $BTAG_NTPL_OUTDIR
       fi
       computeBTagRatio.py -j $NTPL_JSON -d $NTPL_INPUT -o $BTAG_NTPL_OUTDIR -t MC13TeV_DY # -t MC13TeV_TTTo
   fi
   
   if [[ $step == 1.02 ]]; then  # compute Z pt SFs from LO and NLO DY samples
#       NTPL_OUTDIR=results_Ntpl_ZH_2016_2020_06_19_DY_nozpt
       echo "Merge and Calculate Z pt SFs:"
       echo -e "Input: " $NTPL_OUTDIR "\n Output: " $ZPtSF_OUT
       ## if the output directory does not exist, create it:
       if [ ! -d "$ZPtSF_OUT" ]; then
	   mkdir $ZPtSF_OUT
       fi
       computeDYZPtSF.py -j $NTPL_JSON -d $NTPL_OUTDIR -o $ZPtSF_OUT -t MC13TeV_DY -f True  
   fi
   
   if [[ $step == 1.03 ]]; then  # compute top Pt reweighting
       echo "Merge and Calculate Top pt SFs:"
       echo -e "Input: " $TopSF_INPUT "\n Output: " $TopSF_OUT
       ## if the output directory does not exist, create it:
       if [ ! -d "$TopSF_OUT" ]; then
	   mkdir $TopSF_OUT
       fi
        computeTopSF.py -d $TopSF_INPUT -o $TopSF_OUT -c $vh_tag -l 160 -u 600 -a False -f True
   fi
   
   if [[ $step == 1.1 ]]; then  #submit jobs for h->aa->XXYY analysis
       echo "JOB SUBMISSION for BSM h->aa Analysis"
       echo -e "Input: " $NTPL_JSON "\n Output: " $NTPL_OUTDIR
       ## if the output directory does not exist, create it:
       if [ ! -d "$NTPL_OUTDIR" ]; then
	   mkdir $NTPL_OUTDIR
       fi
	# @btagSFMethod=1: Jet-by-jet updating of b-tagging status
	# @btagSFMethod=2: Event reweighting using discriminator-dependent scale factors
       runLocalAnalysisOverSamples.py -e runhaaAnalysis -b $BTAG_NTPL_OUTDIR -g $RUNLOG -j $NTPL_JSON -o $NTPL_OUTDIR -d $NTPL_INPUT -c $RUNNTPLANALYSISCFG -p "@runSystematics=$do_syst @runMVA=False @reweightDYZPt=True @reweightDYdR16=True @reweightTopPt=True @usemetNoHF=False @verbose=False @useDeepCSV=True @useWNJet=False @runQCD=$doQCD @runZH=$doZH @btagSFMethod=1" -s $queue #-t MC13TeV_TTTo #-r true
   fi
fi

###  ############################################## STEPS between 2 and 3
if [[ $step > 1.999 && $step < 3 ]]; then
   if [[ $step == 2 ]]; then    #extract integrated luminosity of the processed lumi blocks
	echo "MISSING LUMI WILL APPEAR AS DIFFERENCE LUMI ONLY IN in.json"
	mergeJSON.py --output=$RESULTSDIR/json_all.json        $RESULTSDIR/Data*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleMu.json   $RESULTSDIR/Data*_DoubleMu*.json
	mergeJSON.py --output=$RESULTSDIR/json_doubleEl.json   $RESULTSDIR/Data*_DoubleElectron*.json
	mergeJSON.py --output=$RESULTSDIR/json_muEG.json       $RESULTSDIR/Data*_MuEG*.json
	mergeJSON.py --output=$RESULTSDIR/json_in.json  $GOLDENJSON/Cert_*.txt
	echo "MISSING LUMI BLOCKS IN DOUBLE MU DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleMu.json 
	echo "MISSING LUMI BLOCKS IN DOUBLE ELECTRON DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_doubleEl.json 
	echo "MISSING LUMI BLOCKS IN MUON EGAMMA DATASET"
	compareJSON.py --diff $RESULTSDIR/json_in.json $RESULTSDIR/json_muEG.json 

	echo "COMPUTE INTEGRATED LUMINOSITY"
	export LD_LIBRARY_PATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root/lib
        export PYTHONPATH=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root/lib 
        export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH 
        export ROOTSYS=/afs/cern.ch/cms/lumi/brilconda-1.1.7/root 
        export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH 
#	export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH

	pip uninstall brilws -y 
	pip install --upgrade --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time

	if [[ $JSON =~ "2016" ]]; then
	    brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/Normtags/normtag_DATACERT.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
	else
	    brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json -i $RESULTSDIR/json_all.json -u /pb -o $RESULTSDIR/LUMI.txt 
	fi
	tail -n 3 $RESULTSDIR/LUMI.txt  
     fi
  fi     

###  ############################################## STEPS between 3 and 4
if [[ $step > 2.999 && $step < 4 ]]; then
    if [ -f $NTPL_OUTDIR/LUMI.txt ]; then
      INTLUMI=`tail -n 3 $RESULTSDIR/LUMI.txt | cut -d ',' -f 6`
    else
	if [[ $JSON =~ "2016" ]]; then  
	    INTLUMI=36330.00 #35866.932
            echo "Please run step==2 above to calculate int. luminosity for 2016 data!" 
	else
            if [[ $JSON =~ "2017" ]]; then
		INTLUMI=41529.152
	        echo "Please run step==2 above to calculate int. luminosity for 2017 data!"
            else
		if [[ $JSON =~ "2018" ]]; then
		    INTLUMI=59740.565
		    echo "Please run step==2 above to calculate int. luminosity for 2018 data!"
		else
	            echo "Please run step==2 above to calculate int. luminosity!"
		fi
            fi                                                                                                                   
	fi
	echo "WARNING: $RESULTSDIR/LUMI.txt file is missing so use fixed integrated luminosity value, this might be different than the dataset you ran on"
    fi
    
    if [[ $step == 3 || $step == 3.0 ]]; then  # make plots and combined root files
	echo "MAKE SUMMARY ROOT FILE, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI"
	echo "Input DIR = "$NTPL_OUTDIR
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $NTPL_OUTDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased $arguments        
   fi        

    if [[ $step == 3 || $step == 3.01 ]]; then  # make plots and combined root files for limits only
	echo "MAKE SUMMARY ROOT FILE FOR LIMITS, BASED ON AN INTEGRATED LUMINOSITY OF $INTLUMI" 
	echo "Input DIR = "$NTPL_OUTDIR  
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outFile ${PLOTTER}_forLimits.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased --only "all_optim_systs|all_optim_cut|(e|mu|ee|mumu|emu)_(A|B|C|D)_(SR|CR)_(3b|4b)_(bdt|higgsMass|higgsPt|ht|pfmet|ptw|mtw|dphiWh|dRave|dmmin|dphijmet|dphilepmet|lep_pt_raw)(|_shapes)(|_umetup|_umetdown|_jerup|_jerdown|_(SubTotalPileUp|SubTotalRelative|SubTotalPt|SubTotalScale|SubTotalAbsolute|SubTotalMC|AbsoluteStat|AbsoluteScale|AbsoluteMPFBias|Fragmentation|SinglePionECAL|SinglePionHCAL|FlavorQCD|TimePtEta|RelativeJEREC1|RelativeJEREC2|RelativeJERHF|RelativePtBB|RelativePtEC1|RelativePtEC2|RelativePtHF|RelativeBal|RelativeSample|RelativeFSR|RelativeStatFSR|RelativeStatEC|RelativeStatHF|PileUpDataMC|PileUpPtRef|PileUpPtBB|PileUpPtEC1|PileUpPtEC2|PileUpPtHF)_jes(up|down)|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_pdfup|_pdfdown|_nloEWK_up|_nloEWK_down|_lesup|_lesdown|_softbup|_softbdown|_(b|c|l)tag(up|down))" $arguments
    fi

    if [[ $step == 3.02 ]]; then # make plots for data-driven QCD bkg
	echo "MAKE SUMMARY ROOT FILE, for data-driven QCD estimate"
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outFile ${PLOTTER}_qcd.root  --json $JSON --noPlot --fileOption RECREATE --key haa_mcbased --only "(all_optim_systs|all_optim_cut|(e|mu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_(bdt|higgsMass|higgsPt|ht|pfmet|ptw|mtw|dphiWh|dRave|dmmin|dphijmet|lep_pt_raw)(|_shapes)(|_umetup|_umetdown|_jerup|_jerdown|_jesup|_jesdown|_scale_mup|_scale_mdown|_stat_eup|_stat_edown|_stat_eup|_stat_edown|_sys_eup|_sys_edown|_GS_eup|_GS_edown|_resRho_eup|_resRho_edown|_puup|_pudown|_pdfup|_pdfdown|_btagup|_btagdown))" $arguments
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased_qcd/ --outFile ${PLOTTER}_qcd.root  --json $JSON --plotExt .pdf --key haa_mcbased --fileOption READ --noLog --only "(e|mu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_bdt" $arguments
    fi

#(ht|pfmet|ptw|mtw|higgsPt|higgsMass|dRave|dmmin|dphijmet|dphiWh|
    if [[ $step == 3 || $step == 3.1 ]]; then  # make plots and combine root files for mcbased study    
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption RECREATE  --key haa_mcbased --only "((e|mu|ee|mumu|emu)_(A|B|C|D)_(SR|CR|CR5j)_(3b|4b)_(pfmet|pfmetphi|ht|mtw|ptw|dphiWh|dRave|dmmin|higgsMass|higgsPt|nbjets_raw|nbtags_raw|jet_(b1|b2|b3|b4)_jet_(pt|eta|phi)_raw|DeepCSV_(b1|b2|b3|b4)_jet_(pt|eta|phi)_raw|zmass_raw|dphijmet|dphijmet1|dphijmet12|(lead|)lep_(pt|eta)_raw|lep_reliso|bdt|(nvtx_raw|nvtxwgt_raw)))|(e|mu|ee|mumu|emu)_eventflow|all_(nvtx_raw|nvtxwgt_raw)|(e|mu|ee|mumu|emu|e_noniso|mu_noniso|mutrk|mutrk_noniso)_lep_reliso|(e|mu|ee|mumu|emu)(|_metmt|_nj2)_lep_(pt|eta|phi)_raw|raw_dphijmet|raw_j1_dphijmet|raw_minj1j2_dphijmet|(e|mu|ee|mumu|emu)_nj_njets_raw|(e|mu|ee|mumu|emu)_jet_(b1|b2|b3|b4)_jet_(pt|eta|phi)_raw|(e|mu|ee|mumu|emu)_nb_nbjets_raw|(e|mu|ee|mumu|emu)_DeepCSV_(b1|b2|b3|b4)_jet_(pt|eta|phi)_raw|(e|mu|ee|mumu|emu)_alljets_jet_(pt|eta|phi)_raw|(e|mu|ee|mumu|emu)_dRlj_raw|raw(e|mu|ee|mumu|emu)_(pfmet|mtw|ptw)|(e|mu|ee|mumu|emu)_nb_(LOOSE|MEDIUM|TIGHT)_nbjets_raw|(e|mu|ee|mumu|emu)_(sel1|sel2)_evt_cat|(e|mu|ee|mumu|emu)_DeepCSV_(b1|b2|b3|b4)_b_discrim" $arguments 
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased/ --outFile ${PLOTTER}.root  --json $JSON --plotExt .png --plotExt .pdf --key haa_mcbased --fileOption READ --noLog  --signalScale 1 $arguments 
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased_log/ --outFile ${PLOTTER}.root  --json $JSON --plotExt .png --plotExt .pdf --key haa_mcbased  --fileOption READ
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir ${NTPL_OUTDIR}/ --outDir $PLOTSDIR/mcbased_blind/ --outFile ${PLOTTER}.root  --json $JSON --plotExt .png --plotExt .pdf --key haa_mcbased --fileOption READ --noLog --signalScale 10 --blind 0.1 --only "(e|mu|ee|mumu)_A_SR_(3b|4b)_bdt" $arguments
    fi

    if [[ $step == 3 || $step == 3.2 ]]; then # make plots and combine root files for data-driven backgrounds 
	runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outFile ${PLOTTER}.root  --json $JSON --noPlot --fileOption UPDATE   --key haa_datadriven $arguments
        runPlotter --iEcm 13 --iLumi $INTLUMI --inDir $RESULTSDIR/ --outDir $PLOTSDIR/datadriven/  --outFile ${PLOTTER}.root  --json $JSON --no2D --plotExt .png --plotExt .pdf  --key haa_datadriven --fileOption READ $arguments
    fi	
fi

