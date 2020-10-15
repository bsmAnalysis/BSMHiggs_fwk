#!/bin/sh -f

#----------------------------------------------------------- $1 ZH, WH ---------------------------------------------------
#---------------------------- all variables except bdt -------------------------------
shapes=("higgsMass" "higgsPt" "ht" "pfmet" "ptw" "mtw" "dphiWh" "dRave" "dmmin" "dphijmet" "lep_pt_raw")
blind=""
flds=("PrefitPlots_$1ZH_noSoftb" "PrefitPlots_$1WH_noSoftb") 

if [[ $1 =~ "2016" ]]; then  
    INTLUMI=35866.932
elif [[ $1 =~ "2017" ]]; then
    INTLUMI=41529.152
elif [[ $1 =~ "2018" ]]; then
    INTLUMI=59740.565
else
    echo "Please specify the year to calculate int. luminosity!"
fi


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

    computeLimit "${blind}" --noLogy --runZh --lumi $INTLUMI --signalScale 10 --m 60 --histo "${shape}_shapes" --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_$1_ZH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --simfit --shape --index 1 --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples$1.json --key haa_mcbased --dropBckgBelow 0.01  --systpostfix _13TeV ;

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

    computeLimit "${blind}" --noLogy --lumi $INTLUMI --signalScale 10 --m 60 --histo "${shape}_shapes" --modeDD --shape --subFake --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_$1_WH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --simfit --index 1 --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples$1.json --key haa_mcbased  --systpostfix _13TeV --dropBckgBelow 0.015  ;

    cd ..
    count=${count}+1
done
cd ..

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

    computeLimit "${blind}" --noLogy --lumi $INTLUMI --runZh --signalScale 10 --m 60 --histo "${shape}_shapes" --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_$1_ZH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --simfit --shape --index 1 --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples$1.json --key haa_mcbased --dropBckgBelow 0.01  --systpostfix _13TeV ;

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

    computeLimit "${blind}" --noLogy --lumi $INTLUMI --signalScale 10 --m 60 --histo "${shape}_shapes" --modeDD --shape --subFake --in $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/plotter_$1_WH_Sys_noSoftb_forLimits-replaced-signal-mc.root --syst --simfit --index 1 --bins 3b,4b --json $CMSSW_BASE/src/UserCode/bsmhiggs_fwk/test/haa4b/samples$1.json --key haa_mcbased  --systpostfix _13TeV  --dropBckgBelow 0.015  ;

    cd ..
    count=${count}+1
done
cd ..
