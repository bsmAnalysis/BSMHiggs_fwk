rm -rf TEST
mkdir -p TEST;
cd TEST;

computeLimit --m 60 --histo bdt_shapes --in $PWD/../../plotter_2018_03_03_forLimits.root --syst --controlR --shape --index 1 --json $PWD/../../samples2016.json  --shapeMin -9999 --shapeMax 9999 --bins 3b,4b --systpostfix _13TeV --dropBckgBelow 0.0 ;
sh combineCards.sh;

## Only Top:
#combineCards.py E_CR_4b=haa4b_60_13TeV_E_CR_4b.dat MU_CR_4b=haa4b_60_13TeV_MU_CR_4b.dat E_SR_4b=haa4b_60_13TeV_E_SR_4b.dat MU_SR_4b=haa4b_60_13TeV_MU_SR_4b.dat > card_combined.dat
#combineCards.py E_CR_3b=haa4b_60_13TeV_E_CR_3b.dat MU_CR_3b=haa4b_60_13TeV_MU_CR_3b.dat E_SR_3b=haa4b_60_13TeV_E_SR_3b.dat MU_SR_3b=haa4b_60_13TeV_MU_SR_3b.dat > card_combined.dat  

## Only NonTT:
#combineCards.py E_CR_nonTT_4b=haa4b_60_13TeV_E_CR_nonTT_4b.dat MU_CR_nonTT_4b=haa4b_60_13TeV_MU_CR_nonTT_4b.dat  E_SR_4b=haa4b_60_13TeV_E_SR_4b.dat MU_SR_4b=haa4b_60_13TeV_MU_SR_4b.dat > card_combined.dat 
#combineCards.py E_CR_nonTT_3b=haa4b_60_13TeV_E_CR_nonTT_3b.dat MU_CR_nonTT_3b=haa4b_60_13TeV_MU_CR_nonTT_3b.dat  E_SR_3b=haa4b_60_13TeV_E_SR_3b.dat MU_SR_3b=haa4b_60_13TeV_MU_SR_3b.dat > card_combined.dat

## All 4b:
combineCards.py E_CR_nonTT_4b=haa4b_60_13TeV_E_CR_nonTT_4b.dat MU_CR_nonTT_4b=haa4b_60_13TeV_MU_CR_nonTT_4b.dat E_CR_4b=haa4b_60_13TeV_E_CR_4b.dat MU_CR_4b=haa4b_60_13TeV_MU_CR_4b.dat E_SR_4b=haa4b_60_13TeV_E_SR_4b.dat MU_SR_4b=haa4b_60_13TeV_MU_SR_4b.dat > card_combined.dat  

## All 3b:
#combineCards.py E_CR_nonTT_3b=haa4b_60_13TeV_E_CR_nonTT_3b.dat MU_CR_nonTT_3b=haa4b_60_13TeV_MU_CR_nonTT_3b.dat E_CR_3b=haa4b_60_13TeV_E_CR_3b.dat MU_CR_3b=haa4b_60_13TeV_MU_CR_3b.dat E_SR_3b=haa4b_60_13TeV_E_SR_3b.dat MU_SR_3b=haa4b_60_13TeV_MU_SR_3b.dat > card_combined.dat 


#text2workspace.py haa4b_60_13TeV_E_CR_4b.dat -o workspace.root --PO verbose 
text2workspace.py card_combined.dat -o workspace.root --PO verbose --channel-masks
#combine -M FitDiagnostics workspace.root --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --minos all --rMin=-2 --rMax=2 --robustFit 1 
combine -M FitDiagnostics workspace.root --saveShapes --saveWithUncertainties --plots --saveNormalizations --setParameters mask_E_SR_3b=1,mask_MU_SR_3b=1,mask_E_SR_4b=1,mask_MU_SR_4b=1

python ../print.py -u fitDiagnostics.root > log.txt
python ../diffNuisances.py -A -a fitDiagnostics.root -g outputfile.root >> log.txt

#combine -M Asymptotic -m 60 --run expected card_combined.dat -v 3 > COMB.log;

cd ..
