#!/bin/sh -f

#----------------------------------------------------------- $1 ZH, WH ---------------------------------------------------
#---------------------------- all variables except bdt -------------------------------
#shapes=("bdt")
shapes=("higgsMass" "higgsPt" "ht" "pfmet" "ptw" "mtw" "dphiWh" "dRave" "dmmin" "dphilepmet" "lep_pt_raw") # "bdt")
blind="" #--SRblind"
flds=("PrefitPlots_$1ZH_noSoftb" "PrefitPlots_$1WH_noSoftb") 

PLOTTER_ZH=plotter_ZH_$1_2020_02_05_forLimits.root
PLOTTER_WH=plotter_WH_$1_2020_02_05_forLimits.root

SAMPLES=samples$1.json

if [[ $1 =~ "2016" ]]; then  
    INTLUMI=36330.0
    SAMPLES=samples2016_legacy.json
    PLOTTER_ZH=plotter_ZH_$1_2020_06_19_forLimits_v1.root #_v2?  
    PLOTTER_WH=plotter_WH_$1_2020_06_19_forLimits.root 
elif [[ $1 =~ "2017" ]]; then
    INTLUMI=41529.152
elif [[ $1 =~ "2018" ]]; then
    INTLUMI=59740.565
else
    echo "Please specify the year to calculate int. luminosity!"
fi


## ZH 
if [[ $2 =~ "ZH" ]]; then

    fld="${flds[0]}"
    mkdir -p "${fld}" && cd "${fld}"
    count=1
    for shape in ${shapes[@]}
    do
	echo "${count}. **********************************${shape}***********************************"
	dir="TEST_${shape}/"
	rm -rf "${dir}"
	mkdir -p "${dir}";
	cd "${dir}";
	
	computeLimit "${blind}" --noLogy --runZh --lumi $INTLUMI --signalScale 10 --m 60 --histo "${shape}_shapes" --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$PLOTTER_ZH --syst --simfit --shape --index 1 --shapeMin -9999 --shapeMax 9999  --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$SAMPLES --key haa_mcbased --dropBckgBelow 0.0001  --systpostfix _13TeV ;
	
	cd ..
	count=${count}+1
    done
    cd ..

    exit

fi

#exit

## WH
fld="${flds[1]}"
mkdir -p "${fld}" && cd "${fld}"
count=1
for shape in ${shapes[@]}
do
    echo "${count}. **********************************${shape}***********************************"
    dir="TEST_${shape}/"
    rm -rf "${dir}"
    mkdir -p "${dir}";
    cd "${dir}";

    computeLimit "${blind}" --noLogy --lumi $INTLUMI --verbose --signalScale 10  --m 60 --histo "${shape}_shapes" --modeDD --shape --subFake --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$PLOTTER_WH  --syst --simfit --index 1  --shapeMin -9999 --shapeMax 9999  --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$SAMPLES --sumInputFile $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/all_plotter_forLimits.root  --key haa_mcbased  --systpostfix _13TeV --dropBckgBelow 0.0015  ;
    
    cd ..
    count=${count}+1
done
cd ..

exit
#----------------------------  bdt -------------------------------

shapes=("bdt")
blind="--SRblind"
flds=("PrefitPlotsBlind_$1ZH_DYSF_noSoftb" "PrefitPlotsBlind_$1WH_noSoftb")

## ZH 
fld="${flds[0]}"
mkdir -p "${fld}" && cd "${fld}"
count=1
for shape in ${shapes[@]}
do
    echo "${count}. **********************************${shape}***********************************"
    dir="TEST_${shape}/"
    rm -rf "${dir}"
    mkdir -p "${dir}";
    cd "${dir}";

    computeLimit "${blind}" --noLogy --lumi $INTLUMI --runZh --signalScale 10 --m 60 --histo "${shape}_shapes" --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$PLOTTER_ZH --syst --simfit --shape --index 1 --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$SAMPLES --key haa_mcbased --dropBckgBelow 0.01  --systpostfix _13TeV ;

    cd ..
    count=${count}+1
done
cd ..

## WH 
fld="${flds[1]}"
mkdir -p "${fld}" && cd "${fld}"
count=1
for shape in ${shapes[@]}
do
    echo "${count}. **********************************${shape}***********************************"
    dir="TEST_${shape}/"
    rm -rf "${dir}"
    mkdir -p "${dir}";
    cd "${dir}";

    computeLimit "${blind}" --noLogy --lumi $INTLUMI --signalScale 10 --m 60 --histo "${shape}_shapes" --modeDD --shape --subFake --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$PLOTTER_WH --syst --simfit --index 1 --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/$SAMPLES --key haa_mcbased  --systpostfix _13TeV  --dropBckgBelow 0.015  ;

    cd ..
    count=${count}+1
done
cd ..
