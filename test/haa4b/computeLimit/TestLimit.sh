rm -rf TEST
mkdir -p TEST;
cd TEST;
computeLimit --m 60 --histo bdt_shapes --in $PWD/../../plotter_2018_03_03_forLimits.root --controlR --shape --index 1 --json $PWD/../../samples2016.json  --shapeMin -9999 --shapeMax 9999 --bins 4b --systpostfix _13TeV --dropBckgBelow 0.0 ;
sh combineCards.sh;

#combineCards.py E_CR_4b=haa4b_60_13TeV_E_CR_4b.dat MU_CR_4b=haa4b_60_13TeV_MU_CR_4b.dat > card_CR_4b_combined.dat
#combineCards.py E_CR_4b=haa4b_60_13TeV_E_CR_nonTT_4b.dat MU_CR_4b=haa4b_60_13TeV_MU_CR_nonTT_4b.dat > card_CR_nonTT_4b_combined.dat  

text2workspace.py card_combined.dat -o workspace.root --PO verbose  
combine -M FitDiagnostics workspace.root --plots --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL --rMin=-2 --rMax=2 --robustFit 1 

python ../print.py -u fitDiagnostics.root

#combine -M Asymptotic -m 60 --run expected card_combined.dat > COMB.log;

cd ..
